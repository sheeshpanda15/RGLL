# Source ROG_method_comparison.R before fitting. Data preparation is independent
# of the simulation generator and never supplies test responses to an estimator.
real_data_schemes <- function(dataset) {
  if (dataset == "student") {
    columns <- c("studytime", "reason", "goout", "parent_edu")
    pairs <- combn(columns, 2L, simplify = FALSE)
    schemes <- c(as.list(columns), pairs)
  } else if (dataset == "wage") {
    columns <- c("R8496500", "R7898500", "R8497000")
    schemes <- c(as.list(columns), combn(columns, 2L, simplify = FALSE), list(columns))
  } else stop("Unknown dataset: ", dataset)
  setNames(schemes, vapply(schemes, paste, character(1), collapse = "+"))
}

student_identity_columns <- c("school", "sex", "age", "address", "famsize", "Pstatus",
  "Medu", "Fedu", "Mjob", "Fjob", "reason", "nursery", "internet")

row_signature <- function(data, columns) {
  do.call(paste, c(data[columns], sep = "|"))
}

load_real_dataset <- function(dataset, directory = "real_data/data") {
  if (dataset == "student") {
    por <- read.csv(file.path(directory, "student-por.csv"), sep = ";", stringsAsFactors = FALSE)
    mat <- read.csv(file.path(directory, "student-mat.csv"), sep = ";", stringsAsFactors = FALSE)
    por$subject <- 1L; mat$subject <- 2L
    data <- rbind(por, mat)
    units <- match(row_signature(data, student_identity_columns),
                   unique(row_signature(data, student_identity_columns)))
    data$parent_edu <- pmax(2L, as.integer(data$Medu + data$Fedu))
    response <- "G3"
    exclusions <- list(
      Full = c("G3", "G1", "G2", "school", "sex"),
      Excl = c("G3", "G1", "G2", "subject", "parent_edu", "sex", "school", "address",
               "Pstatus", "studytime", "reason", "goout", "Medu", "Fedu")
    )
  } else if (dataset == "wage") {
    data <- read.csv(file.path(directory, "2004.csv"), stringsAsFactors = FALSE)
    data[data < 0] <- NA
    data$R8378703 <- NULL
    original_rows <- which(complete.cases(data))
    data <- data[original_rows, , drop = FALSE]
    units <- data$R0000100
    response <- "R8495200"
    exclusions <- list(
      Full = c("R8495200", "R0000100", "R0173600"),
      Excl = c("R8495200", "R0000100", "R0173600", "R0214800", "R8496500",
               "R7898500", "R8497000")
    )
  } else stop("Unknown dataset: ", dataset)
  stopifnot(!anyNA(data), all(is.finite(data[[response]])), !anyNA(units))
  list(name = dataset, data = data, units = units, response = response,
       original_rows = if (dataset == "student") seq_len(nrow(data)) else original_rows,
       exclusions = exclusions, schemes = real_data_schemes(dataset))
}

real_data_split <- function(units, replication, seed = 42L) {
  stopifnot(replication >= 1L)
  split_seed <- as.integer(seed + (replication - 1L) * 104729L)
  set.seed(split_seed)
  unique_units <- unique(units)
  n <- length(unique_units)
  development <- sort(sample.int(n, round(.7 * n)))
  train <- sort(sample(development, round(.5 * n)))
  reference <- setdiff(development, train)
  test <- setdiff(seq_len(n), development)
  id <- match(units, unique_units)
  indices <- list(train = which(id %in% train), reference = which(id %in% reference),
                  test = which(id %in% test))
  stopifnot(sum(lengths(indices)) == length(units), all(lengths(indices) > 0L),
            !length(intersect(units[indices$train], units[indices$reference])),
            !length(intersect(units[indices$train], units[indices$test])),
            !length(intersect(units[indices$reference], units[indices$test])))
  c(indices, list(seed = split_seed))
}

prepare_real_split <- function(dataset, variant, split) {
  columns <- setdiff(names(dataset$data), dataset$exclusions[[variant]])
  stopifnot(variant %in% names(dataset$exclusions), length(columns) > 0L)
  matrices <- lapply(split[c("train", "reference", "test")], function(idx)
    matrix(0, length(idx), length(columns), dimnames = list(NULL, columns)))
  preprocessing <- vector("list", length(columns)); names(preprocessing) <- columns
  for (column in columns) {
    value <- dataset$data[[column]]
    levels <- NULL
    if (is.character(value) || is.factor(value)) {
      levels <- sort(unique(as.character(value[split$train])))
      value <- match(as.character(value), levels)
      value[is.na(value)] <- median(value[split$train])
    }
    value <- as.numeric(value)
    lo <- min(value[split$train]); hi <- max(value[split$train])
    for (name in names(matrices)) {
      matrices[[name]][, column] <- if (hi > lo) 2 * (value[split[[name]]] - lo) / (hi - lo) - 1 else 0
    }
    preprocessing[[column]] <- list(min = lo, max = hi, levels = levels)
  }
  keep <- vapply(preprocessing, function(z) z$max > z$min, logical(1))
  matrices <- lapply(matrices, function(z) z[, keep, drop = FALSE])
  y <- as.numeric(dataset$data[[dataset$response]])
  center <- mean(y[split$train]); scale <- sd(y[split$train])
  stopifnot(is.finite(scale), scale > 0, all(vapply(matrices, function(z) all(is.finite(z)), logical(1))))
  list(X_train = matrices$train, X_reference = matrices$reference, X_test = matrices$test,
       y_train = (y[split$train] - center) / scale,
       y_test = (y[split$test] - center) / scale,
       y_test_original = y[split$test], y_center = center, y_scale = scale,
       columns = columns[keep], constant_columns_removed = columns[!keep],
       preprocessing = preprocessing, split = split)
}

predict_administrative_group <- function(X_train, y_train, X_test, labels_train,
                                         labels_test, model, config) {
  levels <- sort(unique(labels_train))
  train <- match(labels_train, levels)
  test <- match(labels_test, levels)
  fit <- fit_lmm_from_labels(X_train, y_train, train, model, config$em_tol, config$em_max)
  # Entirely new administrative groups have no estimable BLUP; use the marginal
  # prediction for those rows rather than dropping them from the test set.
  pred <- drop(cbind(1, X_test) %*% fit$beta)
  known <- !is.na(test)
  if (any(known)) pred[known] <- predict_lmm_from_labels(fit, X_test[known, , drop = FALSE], test[known], model)
  list(fit = fit, pred = pred, K = length(levels), new_group_fraction = mean(!known))
}

fit_real_method <- function(method, data, model, seed, config, observed, labels = NULL) {
  X <- data$X_train; y <- data$y_train; Xt <- data$X_test
  p <- ncol(X); tau <- 5L * (p + 1L)
  labeled <- function(train, test) {
    fit <- fit_predict_labeled(X, y, Xt, train, test, model, config$em_tol, config$em_max)
    list(fit = fit$fit, pred = fit$pred, K = fit$K)
  }
  if (method == "LM") {
    fit <- fit_pooled_lm(X, y)
    return(list(fit = fit, pred = predict_pooled_lm(fit, Xt), K = 1L))
  }
  if (method == "OBS") return(labeled(observed$train, observed$test))
  if (method == "KM") {
    grid <- thin_integer_grid(seq.int(2L, max(2L, floor(nrow(X) / tau))), config$k_grid_points)
    km <- fit_kmeans_bic(X, grid, seed + 401L, tau = tau,
                        num_init = config$kmeans_num_init, max_iters = config$kmeans_max_iters)
    return(labeled(km$cluster$labels, assign_kmeans(Xt, km$cluster$centroids)))
  }
  if (method == "CPF") {
    cpf <- fit_cpf_bic(X, y, observed$train, lambda_points = config$cpf_lambda_points,
                       theta = 1, gamma = 3, bic_c = 10, tol = 2e-3, max_iter = config$cpf_max_iter)
    answer <- labeled(cpf$group_map[observed$train], cpf$group_map[observed$test])
    answer$first_stage_converged <- cpf$converged
    return(answer)
  }
  if (method == "BLM") {
    blm <- fit_blm_discretization(X, y, observed$train, model, seed = seed + 607L,
                                  ridge = 1, nstart = config$blm_nstart, iter.max = 50L)
    return(labeled(blm$group_map[observed$train], blm$group_map[observed$test]))
  }
  if (method == "CIRG") {
    search <- cirg_search(X, y, data$X_reference, initial_Cn = 2L, lambda = config$lambda,
      tau = tau, model_type = model, max_iter = config$sa_max, seed = seed + 701L,
      em_tol = config$em_tol, em_max_iter = config$em_max, cirg_criterion = "IMSPE",
      kmeans_num_init = config$kmeans_num_init, kmeans_max_iters = config$kmeans_max_iters)
    best <- search$best
    pred <- if (model == "RI") predict_ri_soft(best$imspe_fit, Xt, best$cluster$params) else
      predict_rs_soft(best$imspe_fit, Xt, best$cluster$params)
    return(list(fit = best$imspe_fit, pred = pred, K = best$K, objective = best$objective,
                IMSPE_selected = best$IMSPE))
  }
  if (startsWith(method, "ADMIN:")) {
    return(predict_administrative_group(X, y, Xt, labels$train, labels$test, model, config))
  }
  stop("Unsupported method: ", method)
}

score_real_method <- function(method, data, model, seed, config, observed, labels = NULL) {
  warnings <- character(0)
  t0 <- proc.time()[[3L]]
  answer <- tryCatch(withCallingHandlers(
    fit_real_method(method, data, model, seed, config, observed, labels),
    warning = function(w) {
      warnings <<- unique(c(warnings, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  row <- data.frame(method_id = method, method = if (method == "CIRG") "PGR" else method,
    model = model, p = ncol(data$X_train), tau = 5L * (ncol(data$X_train) + 1L),
    N_train = nrow(data$X_train), N_reference = nrow(data$X_reference), N_test = nrow(data$X_test),
    MSPE = NA_real_, MSPE_original = NA_real_, K = NA_integer_, converged = NA,
    iterations = NA_integer_, first_stage_converged = NA, new_group_fraction = NA_real_,
    objective = NA_real_, IMSPE_selected = NA_real_, runtime = proc.time()[[3L]] - t0,
    status = "failed", detail = NA_character_, warnings = paste(warnings, collapse = " | "))
  if (inherits(answer, "error")) {
    row$detail <- conditionMessage(answer)
    return(row)
  }
  if (length(answer$pred) != length(data$y_test) || any(!is.finite(answer$pred))) {
    row$detail <- "Nonfinite or incorrectly sized test predictions"
    return(row)
  }
  row$MSPE <- mean((data$y_test - answer$pred)^2)
  row$MSPE_original <- mean((data$y_test_original - (data$y_center + data$y_scale * answer$pred))^2)
  row$K <- answer$K
  row$converged <- isTRUE(answer$fit$converged)
  row$iterations <- as.integer(answer$fit$iterations)
  for (name in c("first_stage_converged", "new_group_fraction", "objective", "IMSPE_selected"))
    if (!is.null(answer[[name]])) row[[name]] <- answer[[name]]
  row$status <- "ok"
  row
}

run_real_replication <- function(dataset, replication, config) {
  split <- real_data_split(dataset$units, replication, config$seed)
  rows <- list(); preprocessing <- list()
  for (variant in config$variants) {
    data <- prepare_real_split(dataset, variant, split)
    preprocessing[[variant]] <- data[c("columns", "constant_columns_removed", "preprocessing", "y_center", "y_scale")]
    observed <- make_pseudo_groups(data$X_train, data$X_test, K = config$pseudo_groups, seed = split$seed + 313L)
    for (model in config$models) {
      for (method in config$methods) {
        row <- score_real_method(method, data, model, split$seed, config, observed)
        row$dataset <- dataset$name; row$variant <- variant; row$replication <- replication
        row$split_seed <- split$seed
        rows[[length(rows) + 1L]] <- row
      }
      if (config$include_admin) for (name in names(dataset$schemes)) {
        labels <- row_signature(dataset$data, dataset$schemes[[name]])
        row <- score_real_method(paste0("ADMIN:", name), data, model, split$seed, config, observed,
          labels = list(train = labels[split$train], test = labels[split$test]))
        row$dataset <- dataset$name; row$variant <- variant; row$replication <- replication
        row$split_seed <- split$seed
        rows[[length(rows) + 1L]] <- row
      }
    }
  }
  raw <- do.call(rbind, rows)
  cat(sprintf("[REAL] %s replication=%d/%d rows=%d failures=%d\n",
              dataset$name, replication, config$nloop, nrow(raw), sum(raw$status != "ok")))
  list(raw = raw, split = split, preprocessing = preprocessing)
}

summarize_real_results <- function(raw, nloop) {
  groups <- split(raw, interaction(raw$dataset, raw$variant, raw$model, raw$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    valid <- z$status == "ok" & is.finite(z$MSPE)
    n <- sum(valid)
    mean <- if (n) mean(z$MSPE[valid]) else NA_real_
    mcse <- if (n > 1L) sd(z$MSPE[valid]) / sqrt(n) else NA_real_
    data.frame(dataset = z$dataset[1], variant = z$variant[1], model = z$model[1],
      method = z$method[1], n_expected = nloop, n_rows = nrow(z), n_valid = n,
      mean_MSPE = mean, MCSE = mcse, lower_95MC = mean - 1.96 * mcse, upper_95MC = mean + 1.96 * mcse,
      median_MSPE = if (n) median(z$MSPE[valid]) else NA_real_,
      mean_MSPE_original = if (n) mean(z$MSPE_original[valid]) else NA_real_,
      n_converged = sum(z$converged %in% TRUE), n_nonconverged = sum(z$converged %in% FALSE),
      n_failures = sum(!valid), mean_K = if (n) mean(z$K[valid]) else NA_real_,
      mean_runtime_seconds = mean(z$runtime), mean_N_train = mean(z$N_train),
      mean_N_reference = mean(z$N_reference), mean_N_test = mean(z$N_test),
      mean_new_group_fraction = if (all(is.na(z$new_group_fraction))) NA_real_ else mean(z$new_group_fraction, na.rm = TRUE))
  }))
}

paired_real_comparisons <- function(raw) {
  key <- c("dataset", "variant", "model", "replication")
  pgr <- raw[raw$method == "PGR", c(key, "MSPE", "status")]
  comparison <- merge(raw[raw$method != "PGR", ], pgr, by = key, suffixes = c("", "_PGR"))
  comparison$difference <- comparison$MSPE_PGR - comparison$MSPE
  groups <- split(comparison, interaction(comparison$dataset, comparison$variant,
                                          comparison$model, comparison$method, drop = TRUE))
  do.call(rbind, lapply(groups, function(z) {
    valid <- z$status == "ok" & z$status_PGR == "ok" & is.finite(z$difference)
    n <- sum(valid); delta <- z$difference[valid]
    avg <- if (n) mean(delta) else NA_real_
    se <- if (n > 1) sd(delta) / sqrt(n) else NA_real_
    data.frame(dataset = z$dataset[1], variant = z$variant[1], model = z$model[1],
      comparator = z$method[1], n_pairs = n, mean_PGR_minus_comparator = avg,
      MCSE_difference = se, lower_95MC = avg - 1.96 * se, upper_95MC = avg + 1.96 * se,
      PGR_win_fraction = if (n) mean(delta < 0) else NA_real_)
  }))
}

real_method_comparison <- function(raw, nloop,
                                   methods = c("PGR", "KM", "BLM", "CPF", "OBS", "LM")) {
  main <- raw[raw$method %in% methods, , drop = FALSE]
  index <- which(raw$method %in% methods)
  key <- interaction(main$dataset, main$variant, main$model, main$replication, drop = TRUE)
  groups <- split(seq_len(nrow(main)), key)
  keep <- integer(0); excluded <- list()
  for (ids in groups) {
    block <- main[ids, , drop = FALSE]
    stopifnot(!anyDuplicated(block$method))
    bad <- block$method[block$status != "ok" | !is.finite(block$MSPE)]
    missing <- setdiff(methods, block$method)
    if (!length(bad) && !length(missing)) {
      keep <- c(keep, ids)
    } else {
      excluded[[length(excluded) + 1L]] <- data.frame(dataset = block$dataset[1],
        variant = block$variant[1], model = block$model[1], replication = block$replication[1],
        failed_or_missing_methods = paste(c(bad, missing), collapse = ","))
    }
  }
  keep <- sort(keep)
  matched <- main[keep, , drop = FALSE]
  if (!nrow(matched)) stop("No common successful replications for the requested method comparison.")
  list(row_indices = index[keep], summary = summarize_real_results(matched, nloop),
       paired = paired_real_comparisons(matched),
       excluded_splits = if (length(excluded)) do.call(rbind, excluded) else
         data.frame(dataset = character(), variant = character(), model = character(),
                    replication = integer(), failed_or_missing_methods = character()))
}
