source("ROG_cases10_11.R")
getenv_num <- function(name, default) {
  x <- Sys.getenv(name, unset = "")
  if (!nzchar(x)) return(default)
  as.numeric(x)
}
getenv_int <- function(name, default) as.integer(getenv_num(name, default))
methods <- strsplit(Sys.getenv("METHODS", unset = "ORACLE,UNIT,OBS,LM,KM,GMM,CPF,BLM,CIRG-U"), ",", fixed = TRUE)[[1]]
methods <- trimws(methods)
ans <- run_case11_comparison(
  nloop = getenv_int("NLOOP", 50),
  p = getenv_int("P", 50),
  R_unit = getenv_int("R_UNIT", 100),
  n_test_per_unit = getenv_int("N_TEST_PER_UNIT", 20),
  obs_bins = getenv_int("OBS_BINS", 4),
  proxy_sd = getenv_num("PROXY_SD", 0.35),
  Var.a = getenv_num("VAR_A", 2.25),
  Var.b = getenv_num("VAR_B", 0),
  Var.e = getenv_num("VAR_E", 9),
  tau = getenv_int("TAU", 5 * (getenv_int("P", 50) + 1L)),
  initial_Cn = getenv_int("INITIAL_K", 4),
  sa_max_iter = getenv_int("SA_MAX_ITER", 80),
  em_tol = getenv_num("EM_TOL", 1e-5),
  em_max_iter = getenv_int("EM_MAX_ITER", 500),
  seed = getenv_int("SEED", 12345),
  methods = methods,
  keep_extra = identical(Sys.getenv("KEEP_EXTRA", unset = "0"), "1")
)
dir.create("results_cases10_11", showWarnings = FALSE, recursive = TRUE)
model <- if (getenv_num("VAR_B", 0) > 0) "RS" else "RI"
outfile <- sprintf("results_cases10_11/%s_case11.rds", model)
saveRDS(ans, outfile)
print(ans$summary)
cat("Saved:", outfile, "\n")
