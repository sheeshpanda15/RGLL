# ROG_aligned_all_cases.R
# -----------------------------------------------------------------------------
# Full Case 1--9 runner for the paper-aligned CIRG/RASC implementation.
#
# This file intentionally contains no method implementation.  It sources the
# canonical functions from ROG_aligned_fixed.R so every entry point uses the same
# RASC clustering, SGA prediction, model-matched IMSPE, and tau defaults.
# -----------------------------------------------------------------------------

source("ROG_aligned_fixed.R")

getenv_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  as.numeric(value)
}

getenv_int <- function(name, default) as.integer(getenv_num(name, default))

dir.create("results", showWarnings = FALSE, recursive = TRUE)

P <- getenv_int("P", 50)
R <- getenv_int("R", 20)
N <- getenv_int("N", 2500)
NLOOP <- getenv_int("NLOOP", 50)
VAR_E <- getenv_num("VAR_E", 9)
VAR_A <- getenv_num("VAR_A", 2.25)
VAR_B <- getenv_num("VAR_B", 0.1)
TAU <- getenv_int("TAU", 5 * (P + 1L))
SA_MAX <- getenv_int("SA_MAX_ITER", 80)
EM_MAX <- getenv_int("EM_MAX_ITER", 500)
SEED <- getenv_int("SEED", 12345)

all_results_RI <- vector("list", 9L)
all_results_RS <- vector("list", 9L)
names(all_results_RI) <- paste0("case", 1:9)
names(all_results_RS) <- paste0("case", 1:9)

for (case_id in 1:9) {
  dx <- paste0("case", case_id)

  cat("\n============================================================\n")
  cat("Running ", dx, " : Random Intercept\n", sep = "")
  cat("============================================================\n")

  result_RI <- Comp(
    N_all = N,
    p = P,
    R = R,
    Var.e = VAR_E,
    nloop = NLOOP,
    dist_x = dx,
    dist_a = "N.ori",
    groupsize = "large",
    lambda = 1,
    tau = TAU,
    Var.a = VAR_A,
    em_tol = 1e-5,
    em_max_iter = EM_MAX,
    sa_max_iter = SA_MAX,
    seed = SEED
  )

  all_results_RI[[dx]] <- result_RI
  saveRDS(result_RI, file.path("results", paste0("RI_", dx, ".rds")))

  cat("\n============================================================\n")
  cat("Running ", dx, " : Random Slope\n", sep = "")
  cat("============================================================\n")

  result_RS <- Comp_RS(
    N_all = N,
    p = P,
    R = R,
    Var.e = VAR_E,
    nloop = NLOOP,
    dist_x = dx,
    dist_a = "N.ori",
    groupsize = "large",
    lambda = 1,
    tau = TAU,
    Var.a = VAR_A,
    Var.b = VAR_B,
    em_tol = 1e-5,
    em_max_iter = EM_MAX,
    sa_max_iter = SA_MAX,
    seed = SEED
  )

  all_results_RS[[dx]] <- result_RS
  saveRDS(result_RS, file.path("results", paste0("RS_", dx, ".rds")))

  saveRDS(all_results_RI, file.path("results", "all_RI_results.rds"))
  saveRDS(all_results_RS, file.path("results", "all_RS_results.rds"))

  cat("Finished ", dx, "\n", sep = "")
}

cat("\n============================================================\n")
cat("ALL CASE 1--9 EXPERIMENTS FINISHED\n")
cat("============================================================\n")
