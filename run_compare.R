#!/usr/bin/env Rscript
options(width = 180)
setwd(Sys.getenv("PROJECT_DIR", getwd()))
cat("[RGLL] Loading main-study R/C++ sources.\n")
source("ROG_method_comparison.R")
source("experiment_plan.R")

get_i <- function(name, default) as.integer(Sys.getenv(name, as.character(default)))
get_d <- function(name, default) as.numeric(Sys.getenv(name, as.character(default)))
get_s <- function(name, default) Sys.getenv(name, default)
numbers <- function(name, default) {
  text <- Sys.getenv(name, "")
  if (!nzchar(text)) return(default)
  as.numeric(strsplit(text, ",", fixed = TRUE)[[1]])
}
model <- toupper(get_s("MODEL", "RI"))
case_id <- get_i("CASE", 1L)
var_a_values <- numbers("VAR_A_LIST", .5)
r_values <- numbers("R_LIST", c(10L, 20L, 50L))
nloop <- get_i("NLOOP", 50L)
first <- get_i("REP_FIRST", 1L); last <- get_i("REP_LAST", nloop)
mode <- get_s("RUN_MODE", "base")
if (!model %in% c("RI", "RS") || !case_id %in% 1:16 ||
    mode != "base" || length(var_a_values) != 1L || !is.finite(var_a_values) || var_a_values < 0)
  stop("One main-study task must specify one model, CASE and nonnegative Var.a value.")
if (is.na(nloop) || is.na(first) || is.na(last) || nloop < 1L ||
    first < 1L || first > last || last > nloop) stop("Invalid replication range.")
if (!length(r_values) || any(!is.finite(r_values)) || anyDuplicated(r_values))
  stop("Invalid R_LIST.")
replications <- seq.int(first, last)
n_cores <- get_i("N_CORES", 1L)
base_methods <- c("LM", "KM", "OBS", "CPF", "BLMX", "BLM", "CIRG")
methods <- base_methods

study <- list(
  N = get_i("N", 2500L), p = get_i("P", 50L), nloop = nloop,
  var_e = 9, rs_var_b = get_d("VAR_B", .1),
  label_policy = get_s("LABEL_POLICY", "unknown"), mis_type = get_s("MIS_TYPE", "none"),
  rho = get_d("RHO", 0), tau = get_i("TAU", 5L * (get_i("P", 50L) + 1L)),
  lambda = get_d("LAMBDA", 0), criterion = get_s("CIRG_CRITERION", "IMSPE"),
  sa_max = get_i("SA_MAX", 25L), em_max = get_i("EM_MAX", 300L),
  k_grid_points = get_i("K_GRID_POINTS", 8L),
  kmeans_num_init = get_i("KMEANS_NUM_INIT", 1L), kmeans_max_iters = get_i("KMEANS_MAX_ITERS", 25L),
  cpf_lambda_points = get_i("CPF_LAMBDA_POINTS", 6L), cpf_max_iter = get_i("CPF_MAX_ITER", 300L),
  blm_nstart = get_i("BLM_NSTART", 5L),
  base_methods = base_methods,
  sources = tools::md5sum(c("ROG_aligned_fixed.R", "ROG_method_comparison.R", "literature_baselines.R",
                           "lmm_fast_aligned.cpp", "run_compare.R", "experiment_plan.R"))
)
cat("[RGLL] Component=", mode, " model=", model, " CASE=", case_id,
    " R=", paste(r_values, collapse = ","), " Var.a=", paste(var_a_values, collapse = ","),
    " repetitions=", first, "-", last, " N_CORES=", n_cores, "\n", sep = "")
runs <- list()
for (r in r_values) for (var_a in var_a_values) {
  ans <- run_method_comparison(
    N_all = study$N, p = study$p, R = r, Var.e = study$var_e, Var.a = var_a,
    Var.b = if (model == "RS") study$rs_var_b else 0,
    nloop = nloop, dist_x = paste0("case", case_id), groupsize = "large",
    mis_type = study$mis_type, rho = study$rho, tau = study$tau, lambda = study$lambda,
    cirg_criterion = study$criterion, label_policy = study$label_policy,
    sa_max_iter = study$sa_max, em_tol = 1e-5, em_max_iter = study$em_max,
    k_grid_points = study$k_grid_points, kmeans_num_init = study$kmeans_num_init,
    kmeans_max_iters = study$kmeans_max_iters, cpf_lambda_points = study$cpf_lambda_points,
    cpf_max_iter = study$cpf_max_iter, blm_nstart = study$blm_nstart,
    n_cores = n_cores, methods = methods, seed = 12345L,
    replication_ids = replications, keep_extra = FALSE, verbose = TRUE)
  if (!ans$complete) stop("Task incomplete; no final result will be saved.")
  ans$raw$model <- model; ans$raw$case <- case_id
  ans$raw$R_setting <- r; ans$raw$Var_a_setting <- var_a
  ans$summary$model <- model; ans$summary$case <- case_id
  ans$summary$R_setting <- r; ans$summary$Var_a_setting <- var_a
  runs[[length(runs) + 1L]] <- ans
}
result <- list(
  raw = do.call(rbind, lapply(runs, `[[`, "raw")),
  summary = do.call(rbind, lapply(runs, `[[`, "summary")),
  diagnostics = lapply(runs, `[[`, "diagnostics"),
  complete = TRUE, config = study,
  task = list(component = mode, model = model, case = case_id, Var_a = paste(var_a_values, collapse = ","),
              R_list = paste(r_values, collapse = ","), rep_first = first, rep_last = last,
              backend = "cpu"),
  implementation = "main_study_v1"
)
directory <- get_s("RESULT_DIR", "result/results_main")
outfile <- file.path(directory, "tasks",
  task_result_name(mode, model, case_id, paste(r_values, collapse = ","), paste(var_a_values, collapse = ","), first, last))
atomic_save_rds(result, outfile)
cat("[RGLL] COMPLETE. Saved one final result:", outfile, "\n")
