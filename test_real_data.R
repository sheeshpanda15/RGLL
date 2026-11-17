#!/usr/bin/env Rscript
source("real_data_experiment.R")
for (name in c("student", "wage")) {
  dataset <- load_real_dataset(name)
  stopifnot(nrow(dataset$data) == if (name == "student") 1044L else 5247L)
  if (name == "student") stopifnot(length(unique(dataset$units)) == 662L)
  for (replication in 1:50) {
    split <- real_data_split(dataset$units, replication)
    stopifnot(identical(split, real_data_split(dataset$units, replication)))
    for (variant in c("Full", "Excl")) {
      prepared <- prepare_real_split(dataset, variant, split)
      tau <- 5L * (ncol(prepared$X_train) + 1L)
      stopifnot(nrow(prepared$X_train) >= 2L * tau,
                all(prepared$X_train >= -1 - 1e-10), all(prepared$X_train <= 1 + 1e-10),
                abs(mean(prepared$y_train)) < 1e-10, abs(sd(prepared$y_train) - 1) < 1e-10)
    }
  }
  split <- real_data_split(dataset$units, 1L)
  original <- prepare_real_split(dataset, "Full", split)
  changed <- dataset
  heldout <- c(split$reference, split$test)
  changed$data[[dataset$response]][heldout] <- 1e9
  updated <- prepare_real_split(changed, "Full", split)
  stopifnot(identical(original$X_train, updated$X_train),
            identical(original$y_train, updated$y_train),
            identical(original$preprocessing, updated$preprocessing),
            identical(original$y_center, updated$y_center), identical(original$y_scale, updated$y_scale))
  numeric_column <- if (name == "student") "age" else "R8497200"
  changed <- dataset
  changed$data[[numeric_column]][heldout] <- 1e9
  updated <- prepare_real_split(changed, "Full", split)
  stopifnot(identical(original$X_train, updated$X_train),
            identical(original$preprocessing, updated$preprocessing))
  cat("PASS:", name, "50 unit-disjoint splits; N_train >= 2*tau; training-only scaling and category encoding.\n")
}
stopifnot(!any(grepl("y_test", deparse(body(fit_real_method)), fixed = TRUE)))
cat("PASS: fitting dispatcher does not access held-out responses.\n")
if (file.exists("result/real_data_latest/real_data_results.rds")) {
  original <- readRDS("result/real_data_latest/real_data_results.rds")$raw
  fixture <- original[original$dataset == "student" & original$variant == "Excl" &
    original$model == "RI" & original$replication %in% 1:2 &
    original$method %in% c("PGR", "KM", "BLM", "CPF", "OBS", "LM"), ]
  fixture$status[fixture$method == "KM" & fixture$replication == 2L] <- "failed"
  fixture$MSPE[fixture$method == "KM" & fixture$replication == 2L] <- NA_real_
  fixture$converged[fixture$method == "OBS" & fixture$replication == 1L] <- FALSE
  fixture$MSPE[fixture$method == "LM" & fixture$replication == 1L] <- 1e9
  matched <- real_method_comparison(fixture, 2L)
  stopifnot(length(matched$row_indices) == 6L,
    all(fixture$replication[matched$row_indices] == 1L),
    nrow(matched$excluded_splits) == 1L, matched$excluded_splits$replication == 2L,
    all(matched$summary$n_valid == 1L),
    matched$summary$n_nonconverged[matched$summary$method == "OBS"] == 1L,
    matched$summary$mean_MSPE[matched$summary$method == "LM"] == 1e9)
  cat("PASS: common-split comparison removes entire failed blocks, retains nonconverged and high-loss finite fits.\n")
}
