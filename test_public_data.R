#!/usr/bin/env Rscript
# Check the standalone Beijing experiment from the repository root.
source("beijing_experiment/code/public_data_experiment.R", chdir = TRUE)
dataset <- load_beijing("beijing_experiment/data")
for (replication in 1:50) {
  split <- real_data_split(dataset$units, replication)
  stopifnot(identical(split, real_data_split(dataset$units, replication)))
  for (variant in c("Full", "Excl")) {
    prepared <- prepare_real_split(dataset, variant, split)
    stopifnot(nrow(prepared$X_train) >= 10 * (ncol(prepared$X_train) + 1),
      all(is.finite(prepared$X_train)), all(is.finite(prepared$X_test)),
      !any(colnames(prepared$X_train) %in% c("response", "station", "state", "city", "source", "date", "block")))
  }
  dates <- lapply(split[c("train", "reference", "test")], function(i)
    unique(as.numeric(as.Date(dataset$data$date[i]))))
  for (pair in combn(names(dates), 2L, simplify = FALSE))
    stopifnot(min(abs(outer(dates[[pair[1]]], dates[[pair[2]]], "-"))) >= 5)
}
split <- real_data_split(dataset$units, 1L)
original <- prepare_real_split(dataset, "Full", split)
changed <- dataset
heldout <- c(split$reference, split$test)
changed$data$response[heldout] <- 1e9
column <- dataset$features$Full[1]
changed$data[[column]][heldout] <- 1e9
updated <- prepare_real_split(changed, "Full", split)
stopifnot(identical(original$X_train, updated$X_train),
  identical(original$y_train, updated$y_train),
  identical(original$preprocessing, updated$preprocessing),
  identical(original$y_center, updated$y_center), identical(original$y_scale, updated$y_scale))
cat("PASS: Beijing", "50 unit-disjoint splits, feature masks and training-only transforms.\n")
stopifnot(!any(grepl("y_test", deparse(body(fit_real_method)), fixed = TRUE)))
cat("PASS: estimator dispatcher never reads held-out responses.\n")
