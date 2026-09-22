#!/usr/bin/env Rscript

options(width = 180)
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)
source("ROG_method_comparison.R")

output_dir <- file.path("result", "results_oracle")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- as.integer(Sys.getenv("N_CORES", unset = "8"))
if (!is.finite(n_cores) || n_cores < 1L) n_cores <- 1L

for (model in c("RI", "RS")) {
  var_b <- if (model == "RS") 0.1 else 0
  for (case_id in seq_len(16L)) {
    for (R_value in c(10L, 20L, 50L)) {
      for (var_a in c(0.5, 2.25)) {
        suffix <- if (model == "RS") "_varb0.1" else ""
        outfile <- file.path(
          output_dir,
          sprintf("%s_case%d_R%d_vara%s%s_oracle.rds",
                  model, case_id, R_value,
                  format(var_a, trim = TRUE, scientific = FALSE), suffix)
        )
        if (file.exists(outfile)) {
          existing <- readRDS(outfile)
          if (identical(unique(existing$raw$method), "ORACLE") && nrow(existing$raw) == 20L) {
            next
          }
        }

        message("Oracle: ", model, " case=", case_id,
                " R=", R_value, " Var(a)=", var_a)
        ans <- run_method_comparison(
          N_all = 2500, p = 50, R = R_value,
          Var.e = 9, Var.a = var_a, Var.b = var_b,
          nloop = 20, dist_x = paste0("case", case_id),
          groupsize = "large", mis_type = "none", rho = 0,
          tau = 5 * 51, lambda = 0, cirg_criterion = "IMSPE",
          label_policy = "unknown", n_cores = n_cores,
          methods = "ORACLE", seed = 12345L
        )
        ans$raw$case <- case_id
        ans$raw$R_setting <- R_value
        ans$raw$Var_a_setting <- var_a
        ans$summary$case <- case_id
        ans$summary$R_setting <- R_value
        ans$summary$Var_a_setting <- var_a
        saveRDS(ans, outfile)
      }
    }
  }
}

message("Saved Oracle-only benchmarks to: ", normalizePath(output_dir))
