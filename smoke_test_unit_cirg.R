source("ROG_cases10_11.R")

cat("\n=== Case 10: old row-level vs new unit-level CIRG ===\n")
a10 <- run_case10_comparison(
  nloop = 1, p = 5, K_type = 3, modes_per_type = 2, units_per_mode = 2,
  n_test_per_unit = 8, Var.b = 0, tau = 6, K_grid = 2:6,
  initial_Cn = 2, sa_max_iter = 12, em_tol = 1e-4, em_max_iter = 200,
  methods = c("ORACLE","OBS","LM","KM","BLM","CIRG-ROW","CIRG-U"),
  keep_extra = TRUE
)
print(a10$summary)

cat("\n=== Case 11: old row-level vs new unit-level CIRG ===\n")
a11 <- run_case11_comparison(
  nloop = 1, p = 5, R_unit = 12, n_test_per_unit = 8,
  obs_bins = 3, Var.b = 0, tau = 6, K_grid = 2:6,
  initial_Cn = 2, sa_max_iter = 12, em_tol = 1e-4, em_max_iter = 200,
  methods = c("ORACLE","UNIT","OBS","LM","KM","BLM","CIRG-ROW","CIRG-U"),
  keep_extra = TRUE
)
print(a11$summary)

cat("\nDiagnostic expectation: CIRG-U unit_consistency should equal 1.\n")
