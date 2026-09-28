#!/usr/bin/env Rscript

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

make_paper_clean <- function(plot) {
  plot +
    labs(
      title = NULL,
      subtitle = NULL,
      caption = NULL
    ) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    )
}

make_paper_textfree <- function(plot) {
  make_paper_clean(plot) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

export_clean_set <- function(include_oracle, textfree = FALSE) {
  source_dir <- if (include_oracle) {
    "result/figures_with_oracle"
  } else {
    "result/figures_without_oracle"
  }
  output_dir <- if (include_oracle && textfree) {
    "result/figures_with_oracle_paper_textfree"
  } else if (include_oracle) {
    "result/figures_with_oracle_paper_clean"
  } else if (textfree) {
    "result/figures_without_oracle_paper_textfree"
  } else {
    "result/figures_without_oracle_paper_clean"
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  d <- read.csv(file.path(source_dir, "case_mspe_summary.csv"))
  method_order <- c("CIRG", "KM", "BLM", "CPF", "OBS", "LM")
  if (include_oracle) method_order <- c("ORACLE", method_order)
  d$method <- factor(d$method, levels = method_order)
  d$model <- factor(d$model, levels = c("Random-intercept model (RI)",
                                        "Random-slope model (RS)"))
  d$R_factor <- factor(d$R, levels = c(10, 20, 50))
  d$Var_a_label <- factor(d$Var_a_label,
                          levels = c("Var(a) = 0.50", "Var(a) = 2.25"))

  plots <- vector("list", length(original_cases))
  for (new_case in seq_along(original_cases)) {
    old_case <- original_cases[new_case]
    stopifnot(nrow(d[d$case == old_case, ]) == 12L * length(method_order))
    plot <- make_case_plot(old_case, d, include_oracle = include_oracle)
    plots[[new_case]] <- if (textfree) make_paper_textfree(plot) else make_paper_clean(plot)

    stem <- file.path(output_dir, sprintf("case_%02d_mspe", new_case))
    ggsave(paste0(stem, ".pdf"), plots[[new_case]], width = 10, height = 7,
           units = "in", device = cairo_pdf)
    ggsave(paste0(stem, ".png"), plots[[new_case]], width = 10, height = 7,
           units = "in", dpi = 300, bg = "white")
  }

  grDevices::cairo_pdf(
    file.path(output_dir, "all_cases_mspe.pdf"),
    width = 10, height = 7, onefile = TRUE
  )
  for (plot in plots) print(plot)
  grDevices::dev.off()

  writeLines(c(
    if (include_oracle && textfree) {
      "Paper text-free figure set: all 16 cases, with infeasible Oracle benchmark."
    } else if (include_oracle) {
      "Paper-clean figure set: all 16 cases, with infeasible Oracle benchmark."
    } else if (textfree) {
      "Paper text-free figure set: all 16 cases, feasible methods, without Oracle."
    } else {
      "Paper-clean figure set: all 16 cases, feasible methods, without Oracle."
    },
    if (textfree) {
      "All plot text has been removed, including axis tick labels."
    } else {
      paste(
        "Explanatory titles, subtitles, and captions have been removed.",
        "Method legends, model and Var(a) facet labels, axis titles, and tick labels are retained."
      )
    },
    "The case numbering follows export_paper_figures.R.",
    "",
    "New case | Original case | Description",
    sprintf("%d | %d | %s", seq_along(original_cases), original_cases,
            case_titles[original_cases])
  ), file.path(output_dir, "case_number_mapping.txt"))

  cat("Paper-clean PDF/PNG figures available in", output_dir, "\n")
}

export_clean_set(FALSE)
export_clean_set(TRUE)
export_clean_set(FALSE, textfree = TRUE)
export_clean_set(TRUE, textfree = TRUE)
