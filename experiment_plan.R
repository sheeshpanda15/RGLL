# One independent job per model, CASE and Var.a setting.
experiment_tasks <- function(nloop = 50L) {
  if (length(nloop) != 1L || !is.finite(nloop) || nloop < 1L || nloop != as.integer(nloop))
    stop("NLOOP must be a positive integer.")
  base <- expand.grid(Var_a = c(.5, 2.25), case = 1:16, model = c("RI", "RS"),
                      stringsAsFactors = FALSE)
  base$component <- "base"; base$R_list <- "10,20,50"
  base$rep_first <- 1L; base$rep_last <- nloop
  columns <- c("component", "model", "case", "Var_a", "R_list", "rep_first", "rep_last")
  grid <- base[, columns]
  rownames(grid) <- NULL
  grid
}

task_result_name <- function(component, model, case, R_list, var_a, rep_first, rep_last) {
  sprintf("%s_%s_C%02d_R%s_A%s_rep%03d-%03d.rds", toupper(component), model, case,
          gsub(",", "-", R_list, fixed = TRUE),
          gsub(".", "p", gsub(",", "-", as.character(var_a), fixed = TRUE), fixed = TRUE),
          rep_first, rep_last)
}
