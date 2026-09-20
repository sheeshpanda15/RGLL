setwd(Sys.getenv("PROJECT_DIR", unset = getwd()))
source("ROG_method_comparison.R")
cat("\n--- comparison smoke test (RI, small) ---\n")
ans <- run_method_comparison(
  N_all = 200, p = 5, R = 4, Var.e = 9, Var.a = 2.25, Var.b = 0,
  nloop = 1, dist_x = "case2", groupsize = "small",
  mis_type = "none", rho = 0, tau = 6, lambda = 0,
  sa_max_iter = 3, em_tol = 1e-4, em_max_iter = 200,
  label_policy = "unknown",
  methods = c("LM", "KM", "OBS", "CPF", "BLM", "CIRG")
)
print(ans$summary)
cat("\nNOTE: OBS/CPF/BLM use random pseudo labels in LABEL_POLICY=unknown.\n")

cat("\n--- prediction-focused DGP checks (cases 12--16) ---\n")
C_train <- rep(60L, 10L)
C_test <- rep(20L, 10L)
for (case_id in 12:16) {
  generator <- get(paste0("generate_case", case_id, "_comparison_data"))
  dat <- generator(C_train, C_test, p = 12, beta = rep(1, 12),
                   var_a = 2.25, var_b = 0.1, var_e = 9,
                   seed = 1000L + case_id)
  stopifnot(
    nrow(dat$X_train) == sum(C_train),
    nrow(dat$X_ref) == sum(C_test),
    nrow(dat$X_test) == sum(C_test),
    all(is.finite(dat$X_train)),
    all(is.finite(dat$resp$y_train)),
    isTRUE(all.equal(stats::var(dat$resp$a), 2.25, tolerance = 1e-8)),
    isTRUE(all.equal(mean(apply(dat$resp$B, 2, stats::var)), 0.1,
                     tolerance = 1e-8))
  )
  cat("case", case_id, dat$scenario, "OK\n")
}
