#!/usr/bin/env Rscript
# Collect completed main-study outputs without loading or compiling solvers.
args <- commandArgs(trailingOnly = TRUE)
component <- if (length(args) >= 2L) args[2] else "base"
if (length(args) > 2L || component != "base") stop("Usage: Rscript collect_results.R [RESULT_DIR] [base]")
default_dir <- "result/results_main"
directory <- if (length(args)) args[1] else Sys.getenv("RESULT_DIR", default_dir)
source("experiment_plan.R")
# Load only the pure R summarizer; collection does not compile or run solvers.
for (expr in parse("ROG_method_comparison.R")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("summarize_comparison"))) eval(expr)
}
files <- list.files(file.path(directory, "tasks"), full.names = TRUE)
files <- files[endsWith(files, ".rds") & startsWith(basename(files), paste0(toupper(component), "_"))]
if (!length(files)) stop("No final task results in: ", directory)
first <- readRDS(files[1])
if (!isTRUE(first$complete) || is.null(first$config)) stop("Invalid final result: ", files[1])
study <- first$config
plan <- experiment_tasks(study$nloop)
names <- mapply(task_result_name, plan$component, plan$model, plan$case, plan$R_list,
                plan$Var_a, plan$rep_first, plan$rep_last, USE.NAMES = FALSE)
missing <- setdiff(names, basename(files))
extra <- setdiff(basename(files), names)
if (length(missing) || length(extra)) {
  cat("INCOMPLETE: missing", length(missing), "task results; unexpected files:", length(extra), "\n")
  cat(paste(head(c(missing, extra), 12L), collapse = "\n"), "\n")
  quit(status = 2L)
}
results <- lapply(file.path(directory, "tasks", names), readRDS)
for (i in seq_along(results)) {
  ans <- results[[i]]; task <- plan[i, ]
  if (!isTRUE(ans$complete) || !identical(ans$config, study))
    stop("Incomplete or different experiment configuration: ", names[i])
  for (field in colnames(plan)) {
    if (!identical(as.character(ans$task[[field]]), as.character(task[[field]])))
      stop("Wrong task metadata in: ", names[i])
  }
  methods <- study$base_methods
  rs <- as.integer(strsplit(task$R_list, ",", fixed = TRUE)[[1]])
  var_a <- as.numeric(strsplit(as.character(task$Var_a), ",", fixed = TRUE)[[1]])
  expected <- expand.grid(method = methods, replication = seq.int(task$rep_first, task$rep_last),
                          R_setting = rs, Var_a_setting = var_a, stringsAsFactors = FALSE)
  raw <- ans$raw
  key <- function(x) paste(x$method, x$replication, x$R_setting, x$Var_a_setting, sep = "/")
  if (nrow(raw) != nrow(expected) || anyDuplicated(key(raw)) ||
      !setequal(key(raw), key(expected)) ||
      any(raw$model != task$model) || any(raw$case != task$case))
    stop("Missing, duplicate or incorrectly labeled observations in: ", names[i])
}
raw <- do.call(rbind, lapply(results, `[[`, "raw"))
groups <- split(raw, interaction(raw$model, raw$case, raw$R_setting, raw$Var_a_setting, drop = TRUE))
summary <- do.call(rbind, lapply(groups, function(d) {
  ans <- summarize_comparison(d)
  for (field in c("model", "case", "R_setting", "Var_a_setting")) ans[[field]] <- d[[field]][1]
  ans
}))
combined <- list(raw = raw, summary = summary, complete = TRUE, config = study,
                 component = component, task_count = nrow(plan))
outfile <- file.path(directory, "main_results.rds")
tmp <- tempfile(".result-", tmpdir = directory)
saveRDS(combined, tmp)
if (!file.rename(tmp, outfile)) { unlink(tmp); stop("Cannot replace final result: ", outfile) }
cat("COMPLETE:", nrow(plan), "tasks,", nrow(raw), "rows. Saved:", outfile, "\n")
