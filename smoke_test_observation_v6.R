setwd(Sys.getenv("PROJECT_DIR", unset = getwd()))
source("ROG_cases10_11_observation.R")

cat("\n--- observation-level Case 10 smoke ---\n")
a10 <- run_case10_observation(
  nloop=1, p=5, N_test=180, K_type=3, modes_per_type=2,
  Var.b=0, tau=6, K_grid=2:6, initial_Cn=2,
  sa_max_iter=4, em_tol=1e-4, em_max_iter=200,
  methods=c("ORACLE","OBS","LM","KM","CPF","MIXREG","CIRG")
)
print(a10$summary)

cat("\n--- observation-level Case 11 smoke ---\n")
a11 <- run_case11_observation(
  nloop=1, p=5, N_test=180, obs_bins=3,
  Var.b=0, tau=6, K_grid=2:6, initial_Cn=2,
  sa_max_iter=4, em_tol=1e-4, em_max_iter=200,
  methods=c("ORACLE","OBS","LM","KM","CPF","MIXREG","CIRG")
)
print(a11$summary)

cat("\nBLM is intentionally omitted: its defining classification object is a panel individual, not a single observation.\n")
