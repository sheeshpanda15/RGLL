#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
setwd(project_dir)

result_dir <- file.path("result", "results_comparison")
oracle_dir <- file.path("result", "results_oracle")
output_without_oracle <- file.path("result", "figures_without_oracle")
output_with_oracle <- file.path("result", "figures_with_oracle")

files <- list.files(
  result_dir,
  pattern = "^(RI|RS)_case[0-9]+_R[0-9]+_.*[.]rds$",
  full.names = TRUE
)
if (length(files) != 192L) {
  warning("Expected 192 individual result files; found ", length(files), ".")
}

base_raw <- do.call(rbind, lapply(files, function(path) {
  x <- readRDS(path)$raw
  x$model <- ifelse(x$var_b > 0, "RS", "RI")
  x
}))
oracle_files <- list.files(
  oracle_dir,
  pattern = "^(RI|RS)_case[0-9]+_R[0-9]+_.*[.]rds$",
  full.names = TRUE
)
if (length(oracle_files) != 192L) {
  stop("Expected 192 Oracle result files; found ", length(oracle_files), ". Run add_oracle_benchmark.R first.")
}
oracle_raw <- do.call(rbind, lapply(oracle_files, function(path) {
  x <- readRDS(path)$raw
  x$model <- ifelse(x$var_b > 0, "RS", "RI")
  x
}))
rownames(base_raw) <- NULL
rownames(oracle_raw) <- NULL

summarize_for_plot <- function(raw, method_order) {
  cell_id <- interaction(
    raw$model, raw$case, raw$R_setting, raw$Var_a_setting, raw$method,
    drop = TRUE
  )
  summary_rows <- lapply(split(raw, cell_id), function(d) {
    data.frame(
      model = d$model[1],
      case = d$case[1],
      R = d$R_setting[1],
      Var_a = d$Var_a_setting[1],
      method = d$method[1],
      n = sum(is.finite(d$MSPE)),
      mean_MSPE = mean(d$MSPE, na.rm = TRUE),
      MCSE = stats::sd(d$MSPE, na.rm = TRUE) / sqrt(sum(is.finite(d$MSPE)))
    )
  })
  d <- do.call(rbind, summary_rows)
  rownames(d) <- NULL
  d$lower <- d$mean_MSPE - 1.96 * d$MCSE
  d$upper <- d$mean_MSPE + 1.96 * d$MCSE
  d$method <- factor(d$method, levels = method_order)
  d$model <- factor(
    d$model,
    levels = c("RI", "RS"),
    labels = c("Random-intercept model (RI)", "Random-slope model (RS)")
  )
  d$R_factor <- factor(d$R, levels = c(10, 20, 50))
  d$Var_a_label <- factor(
    d$Var_a,
    levels = c(0.5, 2.25),
    labels = c("Var(a) = 0.50", "Var(a) = 2.25")
  )
  d
}

methods_without_oracle <- c("CIRG", "KM", "BLM", "CPF", "OBS", "LM")
methods_with_oracle <- c("ORACLE", methods_without_oracle)
plot_data_without_oracle <- summarize_for_plot(base_raw, methods_without_oracle)
plot_data_with_oracle <- summarize_for_plot(
  rbind(base_raw, oracle_raw), methods_with_oracle
)

case_titles <- c(
  "Independent uniform covariates",
  "Correlated Gaussian covariates",
  "Group-shifted uniform covariates",
  "Group-shifted Gaussian covariates",
  "Heavy-tailed t covariates",
  "Log-normal covariates",
  "Two-component Gaussian mixture",
  "Mean covariate shift",
  "Covariance shift",
  "Type structure with nuisance modes",
  "Structured latent predictive effects",
  "Continuous predictive heterogeneity",
  "Density-prediction conflict",
  "Overlapping predictive states",
  "Rare high-loss state",
  "Target-distribution shift"
)

method_colors <- c(
  ORACLE = "#7A5195",
  CIRG = "#D55E00",
  KM = "#0072B2",
  BLM = "#009E73",
  CPF = "#56B4E9",
  OBS = "#CC79A7",
  LM = "#000000"
)
method_shapes <- c(ORACLE = 1, CIRG = 16, KM = 17, BLM = 15, CPF = 3, OBS = 8, LM = 18)
method_linetypes <- c(
  ORACLE = "solid", CIRG = "solid", KM = "dashed", BLM = "dotdash",
  CPF = "longdash", OBS = "twodash", LM = "dotted"
)

make_case_plot <- function(case_id, plot_data, include_oracle) {
  d <- plot_data[plot_data$case == case_id, ]
  d$emphasis <- d$method == "CIRG"
  dodge <- position_dodge(width = 0.32)

  ggplot(d, aes(
    x = R_factor, y = mean_MSPE, group = method,
    color = method, shape = method, linetype = method
  )) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      width = 0.10, linewidth = 0.35, alpha = 0.75,
      position = dodge
    ) +
    geom_line(aes(linewidth = emphasis), position = dodge) +
    geom_point(aes(size = emphasis), position = dodge, stroke = 0.8) +
    facet_wrap(
      vars(model, Var_a_label),
      ncol = 2,
      scales = "free_y",
      labeller = label_value
    ) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    scale_shape_manual(values = method_shapes, drop = FALSE) +
    scale_linetype_manual(values = method_linetypes, drop = FALSE) +
    scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 1.05), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 2.0, `TRUE` = 2.8), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))) +
    labs(
      title = sprintf("Case %d: %s", case_id, case_titles[case_id]),
      subtitle = "Mean test MSPE with 95% Monte Carlo error bars; each panel has an independent y-scale",
      x = "Number of data-generating groups (R)",
      y = "Test mean squared prediction error (MSPE)",
      color = NULL, shape = NULL, linetype = NULL,
      caption = if (include_oracle) {
        "Oracle uses the true generating groups and is an infeasible benchmark; all other methods use no true-group labels."
      } else {
        NULL
      }
    ) +
    guides(
      color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(linewidth = 0.8, size = 2.4)),
      shape = guide_legend(nrow = 1, byrow = TRUE),
      linetype = guide_legend(nrow = 1, byrow = TRUE)
    ) +
    theme_bw(base_size = 11, base_family = "sans") +
    theme(
      plot.title = element_text(size = 15, face = "bold", margin = margin(b = 4)),
      plot.subtitle = element_text(size = 10, color = "grey30", margin = margin(b = 10)),
      plot.caption = element_text(size = 8.5, color = "grey30", hjust = 0, margin = margin(t = 6)),
      strip.text = element_text(size = 10.5, face = "bold"),
      strip.background = element_rect(fill = "grey94", color = "grey45", linewidth = 0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.3),
      panel.border = element_rect(color = "grey35", linewidth = 0.45),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9.5, color = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 5),
      legend.key.width = grid::unit(1.6, "lines"),
      plot.margin = margin(10, 12, 8, 10)
    )
}

save_version <- function(plot_data, output_dir, include_oracle) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    plot_data[order(plot_data$case, plot_data$model, plot_data$Var_a,
                    plot_data$R, plot_data$method), ],
    file.path(output_dir, "case_mspe_summary.csv"),
    row.names = FALSE
  )

  plots <- lapply(
    seq_len(16L), make_case_plot,
    plot_data = plot_data, include_oracle = include_oracle
  )
  for (case_id in seq_len(16L)) {
    ggsave(
      file.path(output_dir, sprintf("case_%02d_mspe.pdf", case_id)),
      plots[[case_id]], width = 10, height = 7, units = "in",
      device = cairo_pdf
    )
  }

  grDevices::cairo_pdf(
    file.path(output_dir, "all_cases_mspe.pdf"),
    width = 10, height = 7, onefile = TRUE
  )
  for (p in plots) print(p)
  grDevices::dev.off()

  cat("Saved 16 individual PDFs and combined PDF to:",
      normalizePath(output_dir), "\n")
}

save_version(plot_data_without_oracle, output_without_oracle, FALSE)
save_version(plot_data_with_oracle, output_with_oracle, TRUE)
