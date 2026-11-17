#!/usr/bin/env Rscript

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
stopifnot(length(script) == 1L)
root <- dirname(normalizePath(script, mustWork = TRUE))
data_dir <- file.path(root, "data")
out <- Sys.getenv("RESULT_DIR", file.path(root, "results"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)

packages <- c("Rcpp", "RcppArmadillo", "ClusterR", "MASS")
available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
if (any(!available)) stop("Missing R packages: ", paste(packages[!available], collapse = ", "))

config <- list(
  nloop = as.integer(Sys.getenv("NLOOP", "50")),
  n_cores = as.integer(Sys.getenv("N_CORES", "4")),
  seed = 42L, variants = c("Full", "Excl"), models = c("RI", "RS"),
  methods = c("LM", "KM", "OBS", "CPF", "BLM", "CIRG"),
  include_admin = FALSE, lambda = 0, criterion = "IMSPE", sa_max = 25L,
  em_max = 300L, em_tol = 1e-5, k_grid_points = 8L, kmeans_num_init = 1L,
  kmeans_max_iters = 25L, cpf_lambda_points = 6L, cpf_max_iter = 300L,
  blm_nstart = 5L, pseudo_groups = 20L
)
stopifnot(config$nloop >= 1L, config$n_cores >= 1L)
if (.Platform$OS.type == "windows") config$n_cores <- 1L
group_config <- config
group_config$variants <- "Full"
group_config$methods <- character(0)
group_config$include_admin <- TRUE

setwd(file.path(root, "code"))
source("public_data_experiment.R")
source("ROG_method_comparison.R")
dataset <- load_beijing(data_dir)
cat("Prepared", nrow(dataset$data), "station-day observations.\n")

fits <- parallel::mclapply(seq_len(config$nloop), function(replication) {
  methods <- run_real_replication(dataset, replication, config)
  groups <- run_real_replication(dataset, replication, group_config)
  stopifnot(identical(methods$split, groups$split),
            identical(methods$preprocessing$Full, groups$preprocessing$Full))
  methods$raw <- rbind(methods$raw, groups$raw)
  methods
}, mc.cores = config$n_cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
errors <- vapply(fits, function(fit) !is.list(fit) || is.null(fit$raw), logical(1))
if (any(errors)) stop("Worker errors in repetitions: ", paste(which(errors), collapse = ", "))
raw <- do.call(rbind, lapply(fits, function(fit) fit$raw))

# Table 2 uses common successful splits across the methods and both models
# within each covariate setting. Keep finite predictions from capped fits.
valid <- raw$status == "ok" & is.finite(raw$MSPE)
split_id <- interaction(raw$variant, raw$replication, drop = TRUE)
use <- ave(valid, split_id, FUN = all)
if (!any(use)) stop("No common successful repetitions.")
summary <- summarize_real_results(raw[use, ], config$nloop)
summary <- summary[, c("dataset", "variant", "model", "method", "n_valid",
                       "mean_MSPE", "MCSE", "n_nonconverged")]
excluded <- unique(raw[!use, c("variant", "replication")])
write.csv(raw, file.path(out, "raw.csv"), row.names = FALSE)
write.csv(summary, file.path(out, "table2_summary.csv"), row.names = FALSE)
write.csv(excluded, file.path(out, "excluded_splits.csv"), row.names = FALSE)
saveRDS(list(
  raw = raw, summary = summary, config = config, group_config = group_config,
  splits = lapply(fits, function(fit) fit$split),
  preprocessing = lapply(fits, function(fit) fit$preprocessing),
  audit = dataset$audit, source_md5 = tools::md5sum(list.files()),
  data_md5 = tools::md5sum(list.files(data_dir, full.names = TRUE)),
  session = sessionInfo()
), file.path(out, "beijing_results.rds"))
cat("Saved", nrow(raw), "fits and", nrow(summary), "table entries to", out, "\n")
