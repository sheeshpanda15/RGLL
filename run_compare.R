#!/usr/bin/env Rscript
options(width = 220)
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)
source("ROG_method_comparison.R")

get_i <- function(name, default) as.integer(Sys.getenv(name, unset = as.character(default)))
get_d <- function(name, default) as.numeric(Sys.getenv(name, unset = as.character(default)))
get_s <- function(name, default) Sys.getenv(name, unset = default)

MODEL <- toupper(get_s("MODEL", "RI"))
CASE <- get_i("CASE", 1)
MIS_TYPE <- tolower(get_s("MIS_TYPE", "contam"))
RHO <- get_d("RHO", 0.25)
NLOOP <- get_i("NLOOP", 50)
N <- get_i("N", 2500)
P <- get_i("P", 50)
R <- get_i("R", 20)
TAU <- get_i("TAU", P + 1)
EM_MAX <- get_i("EM_MAX", 500)
SA_MAX <- get_i("SA_MAX", 80)
VAR_B <- if (MODEL == "RS") get_d("VAR_B", 0.1) else 0
methods <- strsplit(get_s("METHODS", "ORACLE,OBS,LM,KM,GMM,CPF,BLM,CIRG"), ",", fixed = TRUE)[[1]]
methods <- trimws(methods)

dir.create("results_comparison", showWarnings = FALSE, recursive = TRUE)
tag_rho <- if (MIS_TYPE == "contam") sprintf("_rho%02d", round(100 * RHO)) else ""
outfile <- sprintf("results_comparison/%s_case%d_%s%s.rds", MODEL, CASE, MIS_TYPE, tag_rho)

cat("============================================================\n")
cat("Method comparison\n")
cat("MODEL=", MODEL, " case=", CASE, " misspec=", MIS_TYPE, " rho=", RHO, "\n", sep="")
cat("N=", N, " p=", P, " R=", R, " nloop=", NLOOP, " tau=", TAU, "\n", sep="")
cat("Methods: ", paste(methods, collapse=", "), "\n", sep="")
cat("============================================================\n")

t0 <- proc.time()[3]
ans <- run_method_comparison(
  N_all = N, p = P, R = R, Var.e = 9, Var.a = 2.25, Var.b = VAR_B,
  nloop = NLOOP, dist_x = paste0("case", CASE), groupsize = "large",
  mis_type = MIS_TYPE, rho = RHO, tau = TAU,
  sa_max_iter = SA_MAX, em_tol = 1e-5, em_max_iter = EM_MAX,
  methods = methods, seed = 12345L
)
elapsed <- (proc.time()[3] - t0) / 60
saveRDS(ans, outfile)
cat("\nFinished. Elapsed minutes:", elapsed, "\n")
cat("Saved:", outfile, "\n\n")
print(ans$summary)
