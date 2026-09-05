#!/usr/bin/env Rscript
options(width = 220)
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)
source("ROG_method_comparison.R")

get_i <- function(name, default) as.integer(Sys.getenv(name, unset = as.character(default)))
get_d <- function(name, default) as.numeric(Sys.getenv(name, unset = as.character(default)))
get_s <- function(name, default) Sys.getenv(name, unset = default)
get_num_list <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  as.numeric(trimws(strsplit(raw, ",", fixed = TRUE)[[1]]))
}
get_int_list <- function(name, default) as.integer(get_num_list(name, default))

MODEL <- toupper(get_s("MODEL", "RI"))
CASE <- get_i("CASE", 1)
CASE_LIST <- get_int_list("CASE_LIST", CASE)
MIS_TYPE <- tolower(get_s("MIS_TYPE", "contam"))
RHO <- get_d("RHO", 0.25)
NLOOP <- get_i("NLOOP", 50)
N <- get_i("N", 2500)
P <- get_i("P", 50)
R <- get_i("R", 20)
R_LIST_DEFAULT <- if (nzchar(Sys.getenv("R", unset = ""))) R else c(10L, 20L, 50L)
VAR_A_SINGLE <- get_d("VAR_A", 2.25)
VAR_A_LIST_DEFAULT <- if (nzchar(Sys.getenv("VAR_A", unset = ""))) VAR_A_SINGLE else c(0.5, 2.25)
R_LIST <- get_int_list("R_LIST", R_LIST_DEFAULT)
VAR_A_LIST <- get_num_list("VAR_A_LIST", VAR_A_LIST_DEFAULT)
TAU <- get_i("TAU", 5 * (P + 1L))
EM_MAX <- get_i("EM_MAX", 500)
SA_MAX <- get_i("SA_MAX", 80)
VAR_B <- if (MODEL == "RS") get_d("VAR_B", 0.1) else 0
methods <- strsplit(get_s("METHODS", "ORACLE,OBS,LM,KM,GMM,CPF,BLM,CIRG"), ",", fixed = TRUE)[[1]]
methods <- trimws(methods)

dir.create("results_comparison", showWarnings = FALSE, recursive = TRUE)
tag_rho <- if (MIS_TYPE == "contam") sprintf("_rho%02d", round(100 * RHO)) else ""
tag_vb <- if (MODEL == "RS") sprintf("_varb%s", format(VAR_B, trim = TRUE, scientific = FALSE)) else ""

cat("============================================================\n")
cat("Method comparison\n")
cat("MODEL=", MODEL, " cases=", paste(CASE_LIST, collapse=","), " misspec=", MIS_TYPE, " rho=", RHO, "\n", sep="")
cat("N=", N, " p=", P, " R=", paste(R_LIST, collapse=","), " Var.a=", paste(VAR_A_LIST, collapse=","), "\n", sep="")
cat("nloop=", NLOOP, " tau=", TAU, "\n", sep="")
cat("Methods: ", paste(methods, collapse=", "), "\n", sep="")
cat("============================================================\n")

t0 <- proc.time()[3]
all_runs <- list()
pos <- 0L
for (case_id in CASE_LIST) {
  for (R_value in R_LIST) {
    for (Var_a in VAR_A_LIST) {
      pos <- pos + 1L
      outfile <- sprintf(
        "results_comparison/%s_case%d_R%d_vara%s%s_%s%s.rds",
        MODEL, case_id, R_value,
        format(Var_a, trim = TRUE, scientific = FALSE),
        tag_vb, MIS_TYPE, tag_rho
      )
      cat("\nRunning case=", case_id, " R=", R_value, " Var.a=", Var_a, "\n", sep="")
      ans <- run_method_comparison(
        N_all = N, p = P, R = R_value, Var.e = 9, Var.a = Var_a, Var.b = VAR_B,
        nloop = NLOOP, dist_x = paste0("case", case_id), groupsize = "large",
        mis_type = MIS_TYPE, rho = RHO, tau = TAU,
        sa_max_iter = SA_MAX, em_tol = 1e-5, em_max_iter = EM_MAX,
        methods = methods, seed = 12345L
      )
      ans$raw$case <- case_id
      ans$raw$R_setting <- R_value
      ans$raw$Var_a_setting <- Var_a
      ans$summary$case <- case_id
      ans$summary$R_setting <- R_value
      ans$summary$Var_a_setting <- Var_a
      saveRDS(ans, outfile)
      cat("Saved:", outfile, "\n")
      all_runs[[pos]] <- ans
    }
  }
}
elapsed <- (proc.time()[3] - t0) / 60
combined <- list(
  raw = do.call(rbind, lapply(all_runs, `[[`, "raw")),
  summary = do.call(rbind, lapply(all_runs, `[[`, "summary"))
)
combined_file <- sprintf("results_comparison/%s_cases%s_grid_%s%s.rds",
                         MODEL, paste(CASE_LIST, collapse="-"), MIS_TYPE, tag_rho)
saveRDS(combined, combined_file)
cat("\nFinished. Elapsed minutes:", elapsed, "\n")
cat("Saved combined:", combined_file, "\n\n")
print(combined$summary)
