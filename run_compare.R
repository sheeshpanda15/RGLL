#!/usr/bin/env Rscript
options(width = 220)
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)
source("ROG_method_comparison.R")

get_i <- function(name, default) as.integer(Sys.getenv(name, unset = as.character(default)))
get_d <- function(name, default) as.numeric(Sys.getenv(name, unset = as.character(default)))
get_s <- function(name, default) Sys.getenv(name, unset = default)
get_num_list <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  as.numeric(trimws(strsplit(raw, ",", fixed = TRUE)[[1]]))
}
get_int_list <- function(name, default) as.integer(get_num_list(name, default))

MODEL <- toupper(get_s("MODEL", "RI"))
CASE <- get_i("CASE", 1)
CASE_LIST <- get_int_list("CASE_LIST", CASE)
LABEL_POLICY <- tolower(get_s("LABEL_POLICY", "unknown"))
if (!LABEL_POLICY %in% c("unknown", "observed")) {
  stop("LABEL_POLICY must be unknown or observed.")
}
MIS_TYPE <- tolower(get_s("MIS_TYPE", if (LABEL_POLICY == "observed") "contam" else "none"))
RHO <- get_d("RHO", 0.25)
NLOOP <- get_i("NLOOP", 20)
N <- get_i("N", 2500)
P <- get_i("P", 50)
R <- get_i("R", 20)
R_LIST_DEFAULT <- if (nzchar(Sys.getenv("R", unset = ""))) R else c(10L, 20L, 50L)
VAR_A_SINGLE <- get_d("VAR_A", 2.25)
VAR_A_LIST_DEFAULT <- if (nzchar(Sys.getenv("VAR_A", unset = ""))) VAR_A_SINGLE else c(0.5, 2.25)
R_LIST <- get_int_list("R_LIST", R_LIST_DEFAULT)
VAR_A_LIST <- get_num_list("VAR_A_LIST", VAR_A_LIST_DEFAULT)
TAU <- get_i("TAU", 5 * (P + 1L))
LAMBDA <- get_d("LAMBDA", 0)
CIRG_CRITERION <- toupper(get_s("CIRG_CRITERION", "IMSPE"))
if (CIRG_CRITERION %in% c("IOPT", "I_OPT", "I-OPT", "I_OPTIMAL", "I-OPTIMAL")) {
  CIRG_CRITERION <- "I"
}
if (!CIRG_CRITERION %in% c("IMSPE", "I")) {
  stop("CIRG_CRITERION must be IMSPE or I.")
}
EM_MAX <- get_i("EM_MAX", 300)
SA_MAX <- get_i("SA_MAX", 25)
K_GRID_POINTS <- get_i("K_GRID_POINTS", 8)
KMEANS_NUM_INIT <- get_i("KMEANS_NUM_INIT", 1)
KMEANS_MAX_ITERS <- get_i("KMEANS_MAX_ITERS", 25)
CPF_LAMBDA_POINTS <- get_i("CPF_LAMBDA_POINTS", 6)
CPF_MAX_ITER <- get_i("CPF_MAX_ITER", 300)
BLM_NSTART <- get_i("BLM_NSTART", 5)
N_CORES_DEFAULT <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (!is.finite(N_CORES_DEFAULT) || N_CORES_DEFAULT < 1L) N_CORES_DEFAULT <- 1L
N_CORES <- get_i("N_CORES", N_CORES_DEFAULT)
VAR_B <- if (MODEL == "RS") get_d("VAR_B", 0.1) else 0
default_methods <- if (LABEL_POLICY == "observed") {
  "ORACLE,OBS,LM,KM,CPF,BLM,CIRG"
} else {
  "LM,KM,OBS,CPF,BLM,CIRG"
}
methods <- strsplit(get_s("METHODS", default_methods), ",", fixed = TRUE)[[1]]
methods <- trimws(methods)
methods <- toupper(methods[nzchar(methods)])
methods <- setdiff(methods, "GMM")
if (!length(methods)) {
  stop("METHODS has no runnable entries after removing GMM.")
}

dir.create("results_comparison", showWarnings = FALSE, recursive = TRUE)
tag_rho <- if (MIS_TYPE == "contam") sprintf("_rho%02d", round(100 * RHO)) else ""
tag_vb <- if (MODEL == "RS") sprintf("_varb%s", format(VAR_B, trim = TRUE, scientific = FALSE)) else ""
tag_lambda <- sprintf("_lambda%s", format(LAMBDA, trim = TRUE, scientific = FALSE))
tag_criterion <- if (CIRG_CRITERION == "IMSPE") "" else "_iopt"
tag_rgrid <- paste(R_LIST, collapse = "-")
tag_vagrid <- paste(format(VAR_A_LIST, trim = TRUE, scientific = FALSE), collapse = "-")

cat("============================================================\n")
cat("Method comparison\n")
cat("MODEL=", MODEL, " cases=", paste(CASE_LIST, collapse=","),
    " label_policy=", LABEL_POLICY, " misspec=", MIS_TYPE, " rho=", RHO, "\n", sep="")
cat("N=", N, " p=", P, " R=", paste(R_LIST, collapse=","), " Var.a=", paste(VAR_A_LIST, collapse=","), "\n", sep="")
cat("nloop=", NLOOP, " tau=", TAU, " lambda=", LAMBDA,
    " CIRG_CRITERION=", CIRG_CRITERION, "\n", sep="")
cat("Budgets: SA_MAX=", SA_MAX, " EM_MAX=", EM_MAX,
    " K_GRID_POINTS=", K_GRID_POINTS,
    " KMEANS_NUM_INIT=", KMEANS_NUM_INIT,
    " KMEANS_MAX_ITERS=", KMEANS_MAX_ITERS,
    " CPF_LAMBDA_POINTS=", CPF_LAMBDA_POINTS,
    " CPF_MAX_ITER=", CPF_MAX_ITER,
    " BLM_NSTART=", BLM_NSTART,
    " N_CORES=", N_CORES, "\n", sep="")
cat("Methods: ", paste(methods, collapse=", "), "\n", sep="")
cat("============================================================\n")

t0 <- proc.time()[3]
all_runs <- list()
pos <- 0L
for (case_id in CASE_LIST) {
  for (R_value in R_LIST) {
    for (Var_a in VAR_A_LIST) {
      pos <- pos + 1L
      outfile <- sprintf(
        "results_comparison/%s_case%d_R%d_vara%s%s%s%s_%s%s.rds",
        MODEL, case_id, R_value,
        format(Var_a, trim = TRUE, scientific = FALSE),
        tag_vb, tag_lambda, tag_criterion, MIS_TYPE, tag_rho
      )
      cat("\nRunning case=", case_id, " R=", R_value, " Var.a=", Var_a, "\n", sep="")
      ans <- run_method_comparison(
        N_all = N, p = P, R = R_value, Var.e = 9, Var.a = Var_a, Var.b = VAR_B,
        nloop = NLOOP, dist_x = paste0("case", case_id), groupsize = "large",
        mis_type = MIS_TYPE, rho = RHO, tau = TAU, lambda = LAMBDA,
        cirg_criterion = CIRG_CRITERION,
        label_policy = LABEL_POLICY,
        sa_max_iter = SA_MAX, em_tol = 1e-5, em_max_iter = EM_MAX,
        k_grid_points = K_GRID_POINTS,
        kmeans_num_init = KMEANS_NUM_INIT,
        kmeans_max_iters = KMEANS_MAX_ITERS,
        cpf_lambda_points = CPF_LAMBDA_POINTS,
        cpf_max_iter = CPF_MAX_ITER,
        blm_nstart = BLM_NSTART,
        n_cores = N_CORES,
        methods = methods, seed = 12345L
      )
      ans$raw$case <- case_id
      ans$raw$R_setting <- R_value
      ans$raw$Var_a_setting <- Var_a
      ans$summary$case <- case_id
      ans$summary$R_setting <- R_value
      ans$summary$Var_a_setting <- Var_a
      saveRDS(ans, outfile)
      cat("Saved:", outfile, "\n")
      all_runs[[pos]] <- ans
    }
  }
}
elapsed <- (proc.time()[3] - t0) / 60
combined <- list(
  raw = do.call(rbind, lapply(all_runs, `[[`, "raw")),
  summary = do.call(rbind, lapply(all_runs, `[[`, "summary"))
)
combined_file <- sprintf("results_comparison/%s_cases%s_R%s_vara%s%s_grid%s_%s%s.rds",
                         MODEL, paste(CASE_LIST, collapse="-"),
                         tag_rgrid, tag_vagrid, tag_vb,
                         paste0(tag_lambda, tag_criterion), MIS_TYPE, tag_rho)
saveRDS(combined, combined_file)
cat("\nFinished. Elapsed minutes:", elapsed, "\n")
cat("Saved combined:", combined_file, "\n\n")
print(combined$summary)
