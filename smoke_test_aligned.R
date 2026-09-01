# smoke_test_aligned.R
# Run from the directory containing ROG_aligned.R and lmm_fast_aligned.cpp.
source("ROG_aligned.R")

set.seed(1)
ans <- run_cirg_simulation(
  N_all = 200,
  p = 3,
  R = 5,
  Var.e = 9,
  Var.a = 2.25,
  Var.b = 0.1,
  nloop = 1,
  dist_x = "case1",
  groupsize = "small",
  lambda = 1,
  tau = 5 * (3 + 1),
  initial_Cn = 2,
  sa_max_iter = 5,
  em_max_iter = 50,
  keep_search = TRUE,
  verbose = TRUE
)

print(ans$raw)
print(ans$summary)
stopifnot(all(is.finite(ans$raw$MSPE)))
stopifnot(all(ans$raw$selected_K >= 2))
cat("Aligned smoke test completed.\n")
