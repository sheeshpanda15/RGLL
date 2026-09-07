# lambda_diagnostic.R
# -----------------------------------------------------------------------------
# Paired quick diagnostic for CIRG residual augmentation.
# This is intentionally smaller than the production grid so it can run locally.
# -----------------------------------------------------------------------------

setwd(Sys.getenv("PROJECT_DIR", unset = getwd()))
source("ROG_method_comparison.R")

N <- as.integer(Sys.getenv("N", unset = "1000"))
P <- as.integer(Sys.getenv("P", unset = "20"))
R <- as.integer(Sys.getenv("R", unset = "20"))
NLOOP <- as.integer(Sys.getenv("NLOOP", unset = "5"))
VAR_A <- as.numeric(Sys.getenv("VAR_A", unset = "2.25"))
VAR_B_RS <- as.numeric(Sys.getenv("VAR_B", unset = "0.1"))
SA_MAX <- as.integer(Sys.getenv("SA_MAX_ITER", unset = "20"))
EM_MAX <- as.integer(Sys.getenv("EM_MAX_ITER", unset = "300"))
TAU <- as.integer(Sys.getenv("TAU", unset = as.character(5 * (P + 1L))))
CASE_LIST <- as.integer(strsplit(Sys.getenv("CASE_LIST", unset = "1,2,3,4,5,6,7,8,9"),
                                 ",", fixed = TRUE)[[1]])
LAMBDA_LIST <- as.numeric(strsplit(Sys.getenv("LAMBDA_LIST", unset = "0,1"),
                                   ",", fixed = TRUE)[[1]])

dir.create("results_comparison", showWarnings = FALSE, recursive = TRUE)

raw_all <- list()
summary_all <- list()
pos <- 0L

for (model in c("RI", "RS")) {
  var_b <- if (model == "RS") VAR_B_RS else 0
  for (lambda in LAMBDA_LIST) {
    for (case_id in CASE_LIST) {
      cat("\nRunning model=", model, " case=", case_id,
          " lambda=", lambda, "\n", sep = "")
      ans <- run_method_comparison(
        N_all = N,
        p = P,
        R = R,
        Var.e = 9,
        Var.a = VAR_A,
        Var.b = var_b,
        nloop = NLOOP,
        dist_x = paste0("case", case_id),
        groupsize = "large",
        mis_type = "contam",
        rho = 0.25,
        tau = TAU,
        lambda = lambda,
        sa_max_iter = SA_MAX,
        em_tol = 1e-5,
        em_max_iter = EM_MAX,
        methods = "CIRG",
        seed = 12345L
      )
      ans$raw$model <- model
      ans$raw$case <- case_id
      ans$raw$lambda <- lambda
      ans$summary$model <- model
      ans$summary$case <- case_id
      ans$summary$lambda <- lambda
      pos <- pos + 1L
      raw_all[[pos]] <- ans$raw
      summary_all[[pos]] <- ans$summary
    }
  }
}

raw <- do.call(rbind, raw_all)
summary <- do.call(rbind, summary_all)

wide <- reshape(summary[, c("model", "case", "lambda", "MSPE_mean",
                            "selected_K_mean", "group_ARI_mean",
                            "convergence_rate")],
                idvar = c("model", "case"),
                timevar = "lambda",
                direction = "wide")

if (all(c("MSPE_mean.0", "MSPE_mean.1") %in% names(wide))) {
  wide$MSPE_delta_lambda0_minus_1 <- wide$MSPE_mean.0 - wide$MSPE_mean.1
  wide$MSPE_pct_lambda0_vs_1 <- 100 * (wide$MSPE_mean.0 / wide$MSPE_mean.1 - 1)
}

saveRDS(list(raw = raw, summary = summary, wide = wide),
        "results_comparison/lambda_diagnostic.rds")
write.csv(raw, "results_comparison/lambda_diagnostic_raw.csv", row.names = FALSE)
write.csv(summary, "results_comparison/lambda_diagnostic_summary.csv", row.names = FALSE)
write.csv(wide, "results_comparison/lambda_diagnostic_wide.csv", row.names = FALSE)

cat("\nSummary:\n")
print(summary[, c("model", "case", "lambda", "MSPE_mean",
                  "selected_K_mean", "group_ARI_mean", "convergence_rate")])
cat("\nPaired lambda=0 vs lambda=1:\n")
print(wide)
