source("ROG_cases10_11.R")
cat("\n--- Case 10 RI smoke ---\n")
a10 <- run_case10_comparison(
  nloop = 1, p = 5, K_type = 3, modes_per_type = 2, units_per_mode = 2,
  n_test_per_unit = 8, Var.b = 0, tau = 6, K_grid = 2:6,
  initial_Cn = 2, sa_max_iter = 4, em_tol = 1e-4, em_max_iter = 200,
  methods = c("ORACLE","OBS","LM","KM","CPF","BLM","CIRG-ROW","CIRG-U")
)
print(a10$summary)

cat("\n--- Case 11 RI smoke ---\n")
a11 <- run_case11_comparison(
  nloop = 1, p = 5, R_unit = 12, n_test_per_unit = 8,
  obs_bins = 3, Var.b = 0, tau = 6, K_grid = 2:6,
  initial_Cn = 2, sa_max_iter = 4, em_tol = 1e-4, em_max_iter = 200,
  methods = c("ORACLE","UNIT","OBS","LM","KM","CPF","BLM","CIRG-ROW","CIRG-U")
)
print(a11$summary)
