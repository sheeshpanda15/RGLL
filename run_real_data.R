#!/usr/bin/env Rscript
options(width = 180)
setwd(Sys.getenv("PROJECT_DIR", getwd()))
source("prepare_real_data.R")
source("real_data_experiment.R")
get_integer <- function(name, default) as.integer(Sys.getenv(name, as.character(default)))
get_list <- function(name, default) strsplit(Sys.getenv(name, default), ",", fixed = TRUE)[[1L]]
config <- list(nloop = get_integer("NLOOP", 50L), n_cores = get_integer("N_CORES", 4L),
  seed = 42L, datasets = get_list("DATASETS", "student,wage"),
  variants = get_list("VARIANTS", "Full,Excl"), models = get_list("MODELS", "RI"),
  methods = c("LM", "KM", "OBS", "CPF", "BLM", "CIRG"),
  include_admin = identical(tolower(Sys.getenv("INCLUDE_ADMIN", "false")), "true"),
  lambda = 0, criterion = "IMSPE", sa_max = 25L, em_max = 300L, em_tol = 1e-5,
  k_grid_points = 8L, kmeans_num_init = 1L, kmeans_max_iters = 25L,
  cpf_lambda_points = 6L, cpf_max_iter = 300L, blm_nstart = 5L, pseudo_groups = 20L,
  split = "unit-level 50% labeled train, 20% independent unlabeled reference, 30% test",
  implementation = "real_data_latest_v1")
stopifnot(config$nloop >= 1L, config$n_cores >= 1L, length(config$datasets) > 0L,
          all(config$datasets %in% c("student", "wage")), !anyDuplicated(config$datasets),
          all(config$variants %in% c("Full", "Excl")), !anyDuplicated(config$variants),
          all(config$models %in% c("RI", "RS")), !anyDuplicated(config$models))
inputs <- prepare_real_data_inputs(datasets = config$datasets)
directory <- Sys.getenv("RESULT_DIR", "result/real_data_RI")
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
source_files <- c("ROG_aligned_fixed.R", "ROG_method_comparison.R",
                  "lmm_fast_aligned.cpp", "prepare_real_data.R", "real_data_experiment.R", "run_real_data.R")
config$source_md5 <- tools::md5sum(source_files)
cat("[REAL] Loading the unchanged main-study R/C++ estimators.\n")
source("ROG_method_comparison.R")
cat("[REAL]", config$nloop, "repeated splits;", config$n_cores,
    "workers; no checkpoints; criterion=IMSPE lambda=0.\n")
results <- list()
for (name in config$datasets) {
  dataset <- load_real_dataset(name)
  cat("[REAL] Dataset=", name, " N=", nrow(dataset$data), " units=", length(unique(dataset$units)), "\n", sep = "")
  runs <- parallel::mclapply(seq_len(config$nloop), function(replication)
    run_real_replication(dataset, replication, config),
    mc.cores = config$n_cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
  errors <- vapply(runs, inherits, logical(1), "try-error")
  if (any(errors)) stop("Replication worker failed: ", paste(runs[errors], collapse = "\n"))
  raw <- do.call(rbind, lapply(runs, `[[`, "raw"))
  summary <- summarize_real_results(raw, config$nloop)
  complete <- all(summary$n_valid == config$nloop & summary$n_rows == config$nloop)
  result <- list(raw = raw, summary = summary, paired = paired_real_comparisons(raw),
    comparison = real_method_comparison(raw, config$nloop),
    complete = complete, config = config, inputs = inputs,
    dataset = list(name = name, N = nrow(dataset$data), n_units = length(unique(dataset$units)),
                   original_rows = dataset$original_rows, schemes = dataset$schemes, exclusions = dataset$exclusions),
    splits = lapply(runs, `[[`, "split"), preprocessing = lapply(runs, `[[`, "preprocessing"))
  saveRDS(result, file.path(directory, paste0(name, "_results.rds")))
  write.csv(result$summary, file.path(directory, paste0(name, "_summary.csv")), row.names = FALSE)
  results[[name]] <- result
  cat("[REAL] Saved final dataset result:", name, " complete=", complete, "\n")
}
result <- list(raw = do.call(rbind, lapply(results, `[[`, "raw")),
  summary = do.call(rbind, lapply(results, `[[`, "summary")),
  paired = do.call(rbind, lapply(results, `[[`, "paired")),
  complete = all(vapply(results, `[[`, logical(1), "complete")), config = config, inputs = inputs)
result$comparison <- real_method_comparison(result$raw, config$nloop)
saveRDS(result, file.path(directory, "real_data_results.rds"))
write.csv(result$raw, file.path(directory, "real_data_raw.csv"), row.names = FALSE)
write.csv(result$summary, file.path(directory, "real_data_summary.csv"), row.names = FALSE)
write.csv(result$comparison$summary, file.path(directory, "method_comparison_summary.csv"), row.names = FALSE)
write.csv(result$comparison$paired, file.path(directory, "paired_comparisons.csv"), row.names = FALSE)
write.csv(result$comparison$excluded_splits, file.path(directory, "comparison_excluded_splits.csv"), row.names = FALSE)
failures <- subset(result$raw, status != "ok")
write.csv(failures, file.path(directory, "failures.csv"), row.names = FALSE)
cat("[REAL] FINISHED; rows=", nrow(result$raw), " failures=", nrow(failures),
    " complete=", result$complete, " output=", normalizePath(directory), "\n", sep = "")
if (!result$complete) stop("Some methods failed. Final outputs preserve all failures; inspect failures.csv.")
