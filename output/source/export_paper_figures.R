#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args <- commandArgs(trailingOnly = TRUE)
if (!nzchar(Sys.getenv("XDG_CACHE_HOME"))) Sys.setenv(XDG_CACHE_HOME = file.path(tempdir(), "font-cache"))
dir.create(Sys.getenv("XDG_CACHE_HOME"), recursive = TRUE, showWarnings = FALSE)
input <- if (length(args)) args[1] else "result/results_main/main_results.rds"
out <- if (length(args) > 1L) args[2] else "output"
scope <- if (length(args) > 2L) args[3] else "all"
stopifnot(scope %in% c("all", "cases", "additional"))
study <- readRDS(input)
stopifnot(isTRUE(study$complete), study$task_count == 64L, study$config$nloop == 50L)
case_map <- c(11L, 12L, 13L, 14L, 15L, 3L, 16L)
selected_cases <- if (scope == "additional") setdiff(1:16, case_map) else case_map
expected_points <- length(selected_cases) * 72L
methods <- c("CIRG", "KM", "BLM", "CPF", "OBS", "LM")
raw <- subset(study$raw, case %in% selected_cases & method %in% methods)
raw$paper_case <- match(raw$case, case_map)
stopifnot(nrow(raw) == expected_points * 50L, all(is.finite(raw$MSPE)))
key <- function(z) interaction(z$model, z$case, z$R_setting, z$Var_a_setting, z$method, drop = TRUE)
points <- do.call(rbind, lapply(split(raw, key(raw)), function(z) {
  stopifnot(nrow(z) == 50L, identical(sort(as.integer(z$replication)), 1:50))
  value <- z$MSPE
  avg <- mean(value); mcse <- sd(value) / sqrt(length(value))
  data.frame(model = z$model[1], paper_case = z$paper_case[1], code_case = z$case[1],
    R = z$R_setting[1], Var_a = z$Var_a_setting[1], method = z$method[1],
    n = length(value), mean = avg, MCSE = mcse,
    mean_lower = avg - 1.96 * mcse, mean_upper = avg + 1.96 * mcse,
    median = median(value), q05 = unname(quantile(value, .05, type = 7)),
    q95 = unname(quantile(value, .95, type = 7)),
    n_converged = sum(z$converged), n_nonconverged = sum(!z$converged))
}))
rownames(points) <- NULL
stopifnot(nrow(points) == expected_points, all(points$n == 50L),
          all(points$mean_lower <= points$mean & points$mean <= points$mean_upper),
          all(points$q05 <= points$median & points$median <= points$q95))
existing <- study$summary[study$summary$case %in% selected_cases & study$summary$method %in% methods, ]
old_key <- paste(existing$model, existing$case, existing$R_setting, existing$Var_a_setting, existing$method)
new_key <- paste(points$model, points$code_case, points$R, points$Var_a, points$method)
expected <- existing[match(new_key, old_key), ]
stopifnot(!anyNA(expected$MSPE_mean),
          isTRUE(all.equal(points$mean, expected$MSPE_mean, check.attributes = FALSE)),
          isTRUE(all.equal(points$MCSE, expected$MSPE_MCSE, check.attributes = FALSE)))

colors <- c(CIRG = "#D55E00", KM = "#0072B2", BLM = "#009E73", CPF = "#56B4E9", OBS = "#CC79A7", LM = "#252525")
shapes <- c(CIRG = 16, KM = 17, BLM = 15, CPF = 3, OBS = 8, LM = 18)
line_types <- c(CIRG = "solid", KM = "22", BLM = "dotdash", CPF = "longdash", OBS = "dashed", LM = "dotted")
offsets <- setNames(seq(-.15, .15, length.out = length(methods)), methods)
points$method <- factor(points$method, levels = methods)
points$case_label <- if (scope == "additional") {
  factor(paste("Code CASE", points$code_case), levels = paste("Code CASE", selected_cases))
} else factor(paste("CASE", points$paper_case), levels = paste("CASE", 1:7))
points$var_label <- factor(sprintf("Var(a) = %.2f", points$Var_a), levels = c("Var(a) = 0.50", "Var(a) = 2.25"))
points$x <- match(points$R, c(10, 20, 50)) + offsets[as.character(points$method)]
points$panel <- factor(paste0(points$model, " | ", points$var_label,
  ifelse(points$model == "RS", " | Var(b) = 0.10", "")),
  levels = c("RI | Var(a) = 0.50", "RI | Var(a) = 2.25",
             "RS | Var(a) = 0.50 | Var(b) = 0.10", "RS | Var(a) = 2.25 | Var(b) = 0.10"))
for (folder in c("pdf", "figures", "source")) dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
statistics_name <- if (scope == "additional") "additional_figure_statistics.rds" else "paper_figure_statistics.rds"
saveRDS(list(points = points, case_map = case_map, selected_cases = selected_cases, source = normalizePath(input),
             source_md5 = tools::md5sum(input), config = study$config),
        file.path(out, "source", statistics_name))

make_plot <- function(data, interval = "mean_ci", compact = FALSE) {
  z <- data
  if (interval == "mean_ci") {
    z$center <- z$mean; z$lower <- z$mean_lower; z$upper <- z$mean_upper
    y_label <- "Mean test MSPE"
  } else {
    z$center <- z$median; z$lower <- z$q05; z$upper <- z$q95
    y_label <- "Median test MSPE"
  }
  ggplot(z, aes(x, center, color = method, shape = method, linetype = method, group = method)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = .035, linewidth = .26, alpha = .60) +
    geom_line(linewidth = .37) +
    geom_point(size = if (compact) 1.25 else 1.6, stroke = .45) +
    scale_color_manual(values = colors, breaks = methods) +
    scale_shape_manual(values = shapes, breaks = methods) +
    scale_linetype_manual(values = line_types, breaks = methods) +
    scale_x_continuous(breaks = 1:3, labels = c("10", "20", "50"), limits = c(.68, 3.32), expand = expansion(mult = 0)) +
    scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(.06, .08))) +
    labs(x = "Number of data-generating groups (R)", y = y_label, color = NULL, shape = NULL, linetype = NULL) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = .6, size = 1.8)),
           shape = guide_legend(nrow = 1, byrow = TRUE), linetype = guide_legend(nrow = 1, byrow = TRUE)) +
    theme_bw(base_size = if (compact) 8.5 else 9, base_family = "Arial") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#E7E7E7", linewidth = .25),
      panel.border = element_rect(color = "#777777", linewidth = .32),
      strip.background = element_rect(fill = "#F3F3F3", color = "#777777", linewidth = .3),
      strip.text = element_text(size = 8.5, color = "#111111"),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      strip.placement = "outside", strip.switch.pad.grid = grid::unit(2, "pt"),
      axis.text = element_text(color = "#222222", size = 8),
      axis.title = element_text(size = if (compact) 8.5 else 9),
      axis.ticks = element_line(linewidth = .25), axis.ticks.length = grid::unit(2, "pt"),
      panel.spacing = grid::unit(if (compact) .11 else .18, "in"),
      legend.position = "bottom", legend.key.width = grid::unit(.33, "in"),
      legend.key.height = grid::unit(.16, "in"), legend.spacing.x = grid::unit(3, "pt"),
      legend.text = element_text(size = 8.5), legend.margin = margin(3, 0, 0, 0),
      plot.title = element_text(size = 10, face = "bold", hjust = .5, margin = margin(b = 7)),
      plot.margin = margin(7, 8, 5, 7))
}

make_case_plot <- function(data) {
  panel_labels <- function(labels) {
    ids <- match(labels$panel, levels(data$panel))
    list(
      c("Random-intercept model (RI)", "Random-intercept model (RI)",
        "Random-slope model (RS)", "Random-slope model (RS)")[ids],
      c("Var(a) = 0.50", "Var(a) = 2.25", "Var(a) = 0.50", "Var(a) = 2.25")[ids]
    )
  }
  make_plot(data) +
    facet_wrap(vars(panel), ncol = 2, scales = "free_y", labeller = panel_labels) +
    labs(y = "Test mean squared prediction error (MSPE)") +
    theme(strip.text = element_text(size = 8.5, face = "bold", color = "#111111",
                                   margin = margin(3.5, 2, 3.5, 2)),
          strip.background = element_rect(fill = "#EFEFEF", color = "#777777", linewidth = .3),
          panel.spacing.x = grid::unit(.18, "in"),
          panel.spacing.y = grid::unit(.06, "in"))
}

manifest_path <- file.path(out, "source", "figure_manifest.rds")
exports <- if (scope == "cases" && file.exists(manifest_path)) readRDS(manifest_path) else list()
export <- function(p, name, width, height, interval, expected_points) {
  pdf_path <- file.path(out, "pdf", paste0(name, ".pdf"))
  png_path <- file.path(out, "figures", paste0(name, ".png"))
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "Arial", onefile = TRUE)
  print(p); invisible(dev.off())
  ragg::agg_png(png_path, width = width, height = height, units = "in", res = 360, background = "white")
  print(p); invisible(dev.off())
  built <- ggplot_build(p)
  stopifnot(nrow(built$data[[1]]) == expected_points,
            all(is.finite(built$data[[1]]$ymin)), all(is.finite(built$data[[1]]$ymax)))
  exports[[name]] <<- list(pdf = pdf_path, png = png_path, interval = interval, points = expected_points)
  cat("Saved", name, ":", expected_points, "points, each with lower and upper bounds.\n")
}
if (scope == "all") for (model in c("RI", "RS")) {
  z <- points[points$model == model, ]
  title <- if (model == "RI") "Random-intercept model (RI), Var(b) = 0" else "Random-slope model (RS), Var(b) = 0.10"
  for (kind in c("mean_ci", "empirical_range")) {
    p <- make_plot(z, kind, compact = TRUE) +
      facet_grid(case_label ~ var_label, scales = "free_y", switch = "y", axes = "all_x", axis.labels = "all_x") +
      labs(title = title)
    filename <- if (kind == "mean_ci") paste0("main_CASE1-7_", model, "_mean_95MC") else paste0("robustness_CASE1-7_", model, "_median_Q05-Q95")
    export(p, filename, 7.1, 8.8, kind, 252L)
  }
}
for (code_case in selected_cases) {
  z <- points[points$code_case == code_case, ]
  p <- make_case_plot(z)
  filename <- if (scope == "additional") sprintf("CODE_CASE%02d_mean_95MC", code_case) else
    sprintf("CASE%02d_mean_95MC", match(code_case, case_map))
  export(p, filename, 7.1, 5.0, "mean_ci", 72L)
}
saveRDS(exports, manifest_path)
cat("Validated:", expected_points, "distinct parameter/method points, 50 repetitions each; no failed/nonconverged rows were dropped.\n")
