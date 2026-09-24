#!/usr/bin/env Rscript

# Export exactly two complete figure sets: with and without Oracle.
# The selected manuscript cases come first, followed by the remaining cases.
suppressPackageStartupMessages(library(ggplot2))

definitions <- c("case_titles", "method_colors", "method_shapes",
                 "method_linetypes", "make_case_plot")
loaded <- character()
for (expression in parse("plot_case_comparison.R")) {
  if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
      is.symbol(expression[[2]])) {
    name <- as.character(expression[[2]])
    if (name %in% definitions) {
      eval(expression)
      loaded <- c(loaded, name)
    }
  }
}
stopifnot(setequal(loaded, definitions))

original_cases <- c(11L, 12L, 13L, 14L, 15L, 3L, 16L,
                    1L, 2L, 4L, 5L, 6L, 7L, 8L, 9L, 10L)
stopifnot(setequal(original_cases, seq_len(16L)))
skip_existing <- identical(Sys.getenv("SKIP_EXISTING", "false"), "true")

for (include_oracle in c(FALSE, TRUE)) {
output_dir <- if (include_oracle) "result/figures_with_oracle" else "result/figures_without_oracle"
d <- read.csv(file.path(output_dir, "case_mspe_summary.csv"))
method_order <- c("CIRG", "KM", "BLM", "CPF", "OBS", "LM")
if (include_oracle) method_order <- c("ORACLE", method_order)
d$method <- factor(d$method, levels = method_order)
d$model <- factor(d$model, levels = c("Random-intercept model (RI)",
                                      "Random-slope model (RS)"))
d$R_factor <- factor(d$R, levels = c(10, 20, 50))
d$Var_a_label <- factor(d$Var_a_label,
                       levels = c("Var(a) = 0.50", "Var(a) = 2.25"))

for (new_case in seq_along(original_cases)) {
  old_case <- original_cases[new_case]
  stopifnot(nrow(d[d$case == old_case, ]) == 12L * length(method_order))
  stem <- file.path(output_dir, sprintf("case_%02d_mspe", new_case))
  if (skip_existing && all(file.exists(paste0(stem, c(".pdf", ".png"))))) next
  plot <- make_case_plot(old_case, d, include_oracle = include_oracle) +
    labs(title = sprintf("Case %d: %s", new_case, case_titles[old_case]))
  ggsave(paste0(stem, ".pdf"), plot, width = 10, height = 7,
         units = "in", device = cairo_pdf)
  ggsave(paste0(stem, ".png"), plot, width = 10, height = 7,
         units = "in", dpi = 300, bg = "white")
}

writeLines(c(
  if (include_oracle) "Complete figure set: all 16 cases, with infeasible Oracle benchmark." else
    "Complete figure set: all 16 cases, feasible methods, without Oracle.",
  "Cases 1-7: selected manuscript cases. Cases 8-16: remaining original cases.",
  "PDF: vector format for LaTeX. PNG: 300 dpi for image uploads.",
  "Only the case number in each title has changed; plot data and styling are preserved.",
  "",
  "New case | Original case | Description",
  sprintf("%d | %d | %s", seq_along(original_cases), original_cases,
          case_titles[original_cases])
), file.path(output_dir, "case_number_mapping.txt"))

cat("Complete set of 16 renumbered PDF/PNG pairs available in", output_dir, "\n")
}
