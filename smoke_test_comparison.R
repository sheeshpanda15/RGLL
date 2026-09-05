setwd(Sys.getenv("PROJECT_DIR", unset = getwd()))
source("ROG_method_comparison.R")
cat("\n--- comparison smoke test (RI, small) ---\n")
ans <- run_method_comparison(
  N_all = 200, p = 5, R = 4, Var.e = 9, Var.a = 2.25, Var.b = 0,
  nloop = 1, dist_x = "case2", groupsize = "small",
  mis_type = "contam", rho = 0.25, tau = 6,
  sa_max_iter = 3, em_tol = 1e-4, em_max_iter = 200,
  methods = c("ORACLE", "OBS", "LM", "KM", "CPF", "BLM", "CIRG")
)
print(ans$summary)
cat("\nNOTE: GMM is omitted from this smoke test so mclust is not required.\n")
