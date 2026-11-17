# Prepare the Beijing air-quality data for the RI and RS experiment.
source("real_data_experiment.R")

public_dataset <- function(name, data, features, schemes, units, original_rows, audit) {
  metadata <- setdiff(names(data), unique(c(features$Full, features$Excl)))
  stopifnot(!anyNA(data), all(vapply(data[features$Full], is.numeric, logical(1))),
            all(is.finite(data$response)), !anyNA(units), length(units) == nrow(data))
  list(name = name, data = data, units = units, response = "response",
       original_rows = original_rows, schemes = schemes,
       exclusions = lapply(features, function(x) setdiff(names(data), x)),
       features = features, metadata = metadata, audit = audit)
}

load_beijing <- function(directory) {
  files <- sort(list.files(directory, pattern = "^PRSA_Data_.*\\.csv$", full.names = TRUE))
  stopifnot(length(files) == 12L)
  raw <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
  raw$date <- as.Date(sprintf("%04d-%02d-%02d", raw$year, raw$month, raw$day))
  directions <- c("N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
                  "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW")
  angle <- (match(raw$wd, directions) - 1) * pi / 8
  raw$wind_sin <- sin(angle); raw$wind_cos <- cos(angle)
  columns <- c("PM2.5", "TEMP", "PRES", "DEWP", "RAIN", "WSPM", "wind_sin", "wind_cos")
  ids <- split(seq_len(nrow(raw)), paste(raw$station, raw$date, sep = "|"))
  # Each daily target/covariate needs at least 18 observed hours. No response or
  # covariate interpolation across days or partitions is performed.
  daily <- lapply(ids, function(i) {
    z <- raw[i, , drop = FALSE]
    count <- colSums(!is.na(z[columns]))
    if (any(count < 18L)) return(NULL)
    data.frame(date = z$date[1], station = z$station[1],
      as.list(vapply(z[columns], mean, numeric(1), na.rm = TRUE)), check.names = FALSE)
  })
  daily <- do.call(rbind, daily)
  daily <- daily[order(daily$date, daily$station), , drop = FALSE]
  stopifnot(all(daily$PM2.5 >= 0), all(is.finite(as.matrix(daily[columns]))))
  day <- as.integer(daily$date - as.Date("2013-03-01"))
  block <- day %/% 28L
  # Remove two boundary days on either side of each block, for every station.
  # Test periods are interleaved held-out blocks, not future forecasting.
  keep <- day %% 28L >= 2L & day %% 28L <= 25L
  d <- daily[keep, , drop = FALSE]; block <- block[keep]
  month <- as.integer(format(d$date, "%m"))
  yday <- as.integer(format(d$date, "%j"))
  data <- data.frame(response = log1p(d$PM2.5), TEMP = d$TEMP, PRES = d$PRES,
    DEWP = d$DEWP, RAIN = d$RAIN, WSPM = d$WSPM,
    wind_sin = d$wind_sin, wind_cos = d$wind_cos,
    annual_sin = sin(2 * pi * yday / 365.25), annual_cos = cos(2 * pi * yday / 365.25),
    year_index = as.integer(format(d$date, "%Y")) - 2013L,
    station = as.character(d$station),
    season = c("Winter", "Spring", "Summer", "Autumn")[(month %% 12L) %/% 3L + 1L],
    date = as.character(d$date), block = as.character(block), stringsAsFactors = FALSE)
  weather <- c("TEMP", "PRES", "DEWP", "RAIN", "WSPM", "wind_sin", "wind_cos")
  public_dataset("beijing", data,
    features = list(Full = c(weather, "annual_sin", "annual_cos", "year_index"), Excl = weather),
    schemes = list(Station = "station", Season = "season", Station_Season = c("station", "season")),
    units = block, original_rows = rownames(d),
    audit = list(original_hourly_N = nrow(raw), n_stations = length(unique(raw$station)),
      original_station_days = length(ids), sufficient_coverage_days = nrow(daily),
      coverage_removed = length(ids) - nrow(daily), boundary_days_removed = sum(!keep),
      final_N = nrow(data), n_split_units = length(unique(block)),
      response = "log1p(daily mean PM2.5); at least 18 observed hours per target/covariate",
      split = "50% train / 20% reference / 30% test by common 28-day blocks across all sites; two-day boundary purge",
      estimand = "Prediction for held-out calendar blocks, conditional on contemporaneous weather; not future forecasting",
      days_per_block = table(block), first_date = min(d$date), last_date = max(d$date)))
}
