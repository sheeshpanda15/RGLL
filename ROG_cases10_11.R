# ROG_cases10_11.R
# -----------------------------------------------------------------------------
# Two additional simulation scenarios designed to study the actual target of
# CIRG: prediction-oriented regrouping when the observed grouping is either
# over-refined relative to a latent predictive partition (Case 10), or when
# heterogeneity is genuinely continuous and must be discretized (Case 11).
#
# This file extends ROG_method_comparison.R. It does NOT alter the existing
# cases 1--9.
#
# Case 10: 60 observed administrative units are nested inside 6 latent
# predictive types. X has 12 geometric modes (2 per type), while the random
# effects are shared within each latent type. Thus pure geometric clustering
# may over-split, whereas a prediction criterion may prefer pooling the paired
# modes. This is deliberately a regrouping problem rather than an oracle-recovery
# problem.
#
# Case 11: each repeated-measure unit has a continuous latent coordinate z.
# Covariate distributions and random effects vary smoothly with z. There is no
# true discrete K. Observed groups are a coarse noisy discretization. The goal
# is to choose a useful finite approximation for prediction.
# -----------------------------------------------------------------------------

if (!exists("run_method_comparison", mode = "function")) {
  if (!file.exists("ROG_method_comparison.R")) stop("ROG_method_comparison.R not found.")
  source("ROG_method_comparison.R")
}
if (!exists("unit_cirg_search", mode = "function")) {
  if (!file.exists("ROG_unit_cirg.R")) stop("ROG_unit_cirg.R not found.")
  source("ROG_unit_cirg.R")
}

# ----------------------------- general helpers -------------------------------

ar_cov <- function(p, rho = 0.25, variances = rep(1, p)) {
  idx <- seq_len(p)
  R <- outer(idx, idx, function(i, j) rho^abs(i - j))
  S <- diag(sqrt(variances), p) %*% R %*% diag(sqrt(variances), p)
  (S + t(S)) / 2
}

balanced_unit_sizes <- function(R_unit, n_test_per_unit, train_multiplier = 3L) {
  C_test <- rep(as.integer(n_test_per_unit), R_unit)
  C_train <- as.integer(train_multiplier) * C_test
  list(train = C_train, test = C_test, ref = C_test)
}

labels_from_units <- function(unit_type, C) {
  rep(as.integer(unit_type), as.integer(C))
}

standardize_direction <- function(v) {
  s <- sqrt(sum(v^2))
  if (!is.finite(s) || s < 1e-12) return(v)
  v / s
}

continuous_group_r2 <- function(z, labels) {
  if (is.null(z) || all(is.na(z)) || is.null(labels)) return(NA_real_)
  z <- as.numeric(z)
  labels <- as.integer(factor(labels))
  tot <- sum((z - mean(z))^2)
  if (!is.finite(tot) || tot < 1e-14) return(NA_real_)
  fitted <- ave(z, labels, FUN = mean)
  1 - sum((z - fitted)^2) / tot
}

unit_consistency <- function(unit_labels, estimated_labels) {
  if (is.null(unit_labels) || is.null(estimated_labels)) return(NA_real_)
  u <- as.integer(factor(unit_labels))
  g <- as.integer(factor(estimated_labels))
  vals <- vapply(split(g, u), function(x) max(tabulate(x)) / length(x), numeric(1))
  mean(vals)
}

pseudo_oracle_fit <- function(beta, var_a, var_b, var_e, model_type = c("RI", "RS")) {
  model_type <- match.arg(model_type)
  q <- length(beta) + 1L
  list(
    beta = c(1, beta),
    Var.a = var_a,
    Var.b = if (model_type == "RS") var_b else NA_real_,
    Var.e = var_e,
    G = if (model_type == "RS") diag(c(var_a, rep(var_b, length(beta)))) else NULL,
    u_hat = if (model_type == "RS") matrix(0, 0, q) else numeric(0),
    converged = TRUE,
    iterations = 0L
  )
}

# Map a unit-level classification back to observations.
map_unit_partition <- function(unit_group_map, unit_labels) {
  unit_group_map[as.integer(unit_labels)]
}

# Common evaluator for a pre-generated dataset. The dataset list must contain:
# X_train, X_ref, X_test, y_train, y_test, beta, model_type, var_a_metric,
# var_b_metric, var_e, observed_train/test, unit_train/test.
# Optional: oracle_train/test labels, oracle_pred, true_partition_test, z_test,
# G_true.
evaluate_generated_dataset <- function(dat,
                                       methods = c("ORACLE", "OBS", "LM", "KM", "GMM", "CPF", "BLM", "CIRG-U"),
                                       tau = 5 * (ncol(dat$X_train) + 1L),
                                       K_grid = NULL,
                                       initial_Cn = 2L,
                                       sa_max_iter = 80,
                                       em_tol = 1e-5,
                                       em_max_iter = 500,
                                       seed = 12345L,
                                       gmm_models = c("EII", "VII", "EEI", "VEI"),
                                       verbose = FALSE) {
  methods <- unique(toupper(methods))
  X_train <- dat$X_train; X_ref <- dat$X_ref; X_test <- dat$X_test
  y_train <- dat$y_train; y_test <- dat$y_test
  beta <- dat$beta; model_type <- dat$model_type
  var_a <- dat$var_a_metric; var_b <- dat$var_b_metric; var_e <- dat$var_e
  G_true <- if (!is.null(dat$G_true)) dat$G_true else NULL
  true_te <- if (!is.null(dat$true_partition_test)) dat$true_partition_test else NULL
  z_te <- if (!is.null(dat$z_test)) dat$z_test else NULL
  unit_te <- if (!is.null(dat$unit_test)) dat$unit_test else NULL

  if (is.null(K_grid)) {
    maxK <- min(20L, floor(nrow(X_train) / max(2L, as.integer(tau))))
    K_grid <- seq.int(2L, max(2L, maxK))
  }
  K_grid <- sort(unique(as.integer(K_grid)))
  K_grid <- K_grid[K_grid >= 2L]

  rows <- list(); extra <- list()
  add <- function(d, labels_test = NULL) {
    d$z_R2 <- continuous_group_r2(z_te, labels_test)
    d$unit_consistency <- unit_consistency(unit_te, labels_test)
    rows[[length(rows) + 1L]] <<- d
  }

  if ("ORACLE" %in% methods) {
    t0 <- proc.time()[3]
    if (!is.null(dat$oracle_pred)) {
      fit <- pseudo_oracle_fit(beta, var_a, var_b, var_e, model_type)
      pred <- dat$oracle_pred
      rt <- proc.time()[3] - t0
      d <- metric_row("ORACLE", fit, pred, y_test, beta, var_a, var_b, var_e,
                      G_true, if (!is.null(dat$oracle_K)) dat$oracle_K else NA_integer_,
                      rt, if (!is.null(true_te) && !is.null(dat$oracle_test)) adjusted_rand_index(dat$oracle_test, true_te) else NA_real_)
      add(d, if (!is.null(dat$oracle_test)) dat$oracle_test else NULL)
    } else {
      fp <- fit_predict_labeled(X_train, y_train, X_test,
                                dat$oracle_train, dat$oracle_test,
                                model_type, em_tol, em_max_iter)
      rt <- proc.time()[3] - t0
      ari <- if (!is.null(true_te)) adjusted_rand_index(fp$labels_test, true_te) else NA_real_
      d <- metric_row("ORACLE", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                      G_true, fp$K, rt, ari)
      add(d, fp$labels_test)
    }
  }

  if ("UNIT" %in% methods) {
    t0 <- proc.time()[3]
    fp <- fit_predict_labeled(X_train, y_train, X_test,
                              dat$unit_train, dat$unit_test,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(fp$labels_test, true_te) else NA_real_
    d <- metric_row("UNIT", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                    G_true, fp$K, rt, ari)
    add(d, fp$labels_test)
  }

  if ("OBS" %in% methods) {
    t0 <- proc.time()[3]
    fp <- fit_predict_labeled(X_train, y_train, X_test,
                              dat$observed_train, dat$observed_test,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(fp$labels_test, true_te) else NA_real_
    d <- metric_row("OBS", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                    G_true, fp$K, rt, ari)
    add(d, fp$labels_test)
  }

  if ("LM" %in% methods) {
    t0 <- proc.time()[3]
    fit <- fit_pooled_lm(X_train, y_train)
    pred <- predict_pooled_lm(fit, X_test)
    rt <- proc.time()[3] - t0
    d <- metric_row("LM", fit, pred, y_test, beta, var_a, var_b, var_e,
                    G_true, 1L, rt, NA_real_)
    add(d, rep(1L, nrow(X_test)))
  }

  if ("KM" %in% methods) {
    t0 <- proc.time()[3]
    km <- fit_kmeans_bic(X_train, K_grid, seed + 401L, tau = tau)
    trlab <- km$cluster$labels
    telab <- assign_kmeans(X_test, km$cluster$centroids)
    fp <- fit_predict_labeled(X_train, y_train, X_test, trlab, telab,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(telab, true_te) else NA_real_
    d <- metric_row("KM", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                    G_true, fp$K, rt, ari)
    add(d, telab); extra$KM <- km
  }

  if ("GMM" %in% methods) {
    if (requireNamespace("mclust", quietly = TRUE)) {
      t0 <- proc.time()[3]
      gm <- fit_gmm_bic(X_train, K_grid, seed + 503L, gmm_models)
      trlab <- gm$labels
      telab <- predict_gmm_labels(gm, X_test)
      enough <- model_type == "RI" || min(tabulate(trlab, nbins = max(trlab))) >= ncol(X_train) + 1L
      if (enough) {
        fp <- fit_predict_labeled(X_train, y_train, X_test, trlab, telab,
                                  model_type, em_tol, em_max_iter)
        rt <- proc.time()[3] - t0
        ari <- if (!is.null(true_te)) adjusted_rand_index(telab, true_te) else NA_real_
        d <- metric_row("GMM", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                        G_true, fp$K, rt, ari)
        add(d, telab)
      } else warning("GMM produced an RS group smaller than p+1; row skipped.")
      extra$GMM <- gm
    } else warning("Package mclust unavailable; GMM skipped.")
  }

  # CPF and BLM operate on repeated-measure units, not the possibly coarse OBS
  # labels. This mirrors their subgroup/discretization use of unit-level data.
  if ("CPF" %in% methods) {
    t0 <- proc.time()[3]
    cpf <- fit_cpf_bic(X_train, y_train, dat$unit_train,
                       theta = 1, gamma = 3, bic_c = 10,
                       tol = 1e-3, max_iter = 1000)
    cpf_tr <- map_unit_partition(cpf$group_map, dat$unit_train)
    cpf_te <- map_unit_partition(cpf$group_map, dat$unit_test)
    fp <- fit_predict_labeled(X_train, y_train, X_test, cpf_tr, cpf_te,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(cpf_te, true_te) else NA_real_
    d <- metric_row("CPF", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                    G_true, fp$K, rt, ari)
    add(d, cpf_te); extra$CPF <- cpf
  }

  if ("BLM" %in% methods) {
    t0 <- proc.time()[3]
    blm <- fit_blm_discretization(X_train, y_train, dat$unit_train,
                                  model_type, seed = seed + 607L, ridge = 1)
    blm_tr <- map_unit_partition(blm$group_map, dat$unit_train)
    blm_te <- map_unit_partition(blm$group_map, dat$unit_test)
    fp <- fit_predict_labeled(X_train, y_train, X_test, blm_tr, blm_te,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(blm_te, true_te) else NA_real_
    d <- metric_row("BLM", fp$fit, fp$pred, y_test, beta, var_a, var_b, var_e,
                    G_true, fp$K, rt, ari)
    add(d, blm_te); extra$BLM <- blm
  }

  # Improved CIRG-U: regroup whole repeated-measure units.  "CIRG" is kept
  # as an alias for backward compatibility.  CIRG-ROW retains the old
  # observation-level algorithm for ablation only.
  if (any(methods %in% c("CIRG", "CIRG-U", "CIRG_U"))) {
    if (is.null(dat$unit_ref)) stop("CIRG-U requires dat$unit_ref.")
    t0 <- proc.time()[3]
    uc <- uc_compact_units(dat$unit_train, dat$unit_ref, dat$unit_test)
    search <- unit_cirg_search(
      X_train, y_train, uc$train,
      X_ref, uc$ref,
      model_type = model_type,
      K_grid = K_grid,
      initial_Cn = initial_Cn,
      tau = tau,
      min_units = 2L,
      x_pc = min(5L, ncol(X_train)),
      effect_pc = min(5L, ncol(X_train)),
      ridge = 1,
      effect_weight = 1,
      max_iter = sa_max_iter,
      seed = seed + 701L,
      em_tol = em_tol,
      em_max_iter = em_max_iter,
      verbose = verbose
    )
    pp <- predict_unit_cirg(search, X_test, uc$test, model_type)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(pp$labels, true_te) else NA_real_
    d <- metric_row("CIRG-U", pp$fit, pp$pred, y_test, beta, var_a, var_b, var_e,
                    G_true, pp$K, rt, ari, search$best$objective)
    add(d, pp$labels); extra$CIRG_U <- search
  }

  if (any(methods %in% c("CIRG-ROW", "CIRG_ROW"))) {
    t0 <- proc.time()[3]
    search <- cirg_search(X_train, y_train, X_ref,
                          initial_Cn = initial_Cn, tau = tau,
                          model_type = model_type,
                          max_iter = sa_max_iter, seed = seed + 1701L,
                          em_tol = em_tol, em_max_iter = em_max_iter,
                          verbose = verbose)
    best <- search$best
    fit <- best$imspe_fit
    pred <- if (model_type == "RI") {
      predict_ri_soft(fit, X_test, best$cluster$params)
    } else {
      predict_rs_soft(fit, X_test, best$cluster$params)
    }
    telab <- predict_soft_labels(X_test, best$cluster$params)
    rt <- proc.time()[3] - t0
    ari <- if (!is.null(true_te)) adjusted_rand_index(telab, true_te) else NA_real_
    d <- metric_row("CIRG-ROW", fit, pred, y_test, beta, var_a, var_b, var_e,
                    G_true, best$K, rt, ari, best$objective)
    add(d, telab); extra$CIRG_ROW <- search
  }

  raw <- do.call(rbind, rows)
  raw$scenario <- dat$scenario
  raw$N_test <- nrow(X_test)
  raw$R_unit <- length(unique(dat$unit_train))
  list(metrics = raw, extra = extra)
}

summarize_structured_comparison <- function(raw) {
  num <- c("MSPE", "beta_MSE", "beta_SSE", "intercept_MSE",
           "Var_a_hat", "Var_b_hat", "Var_e_hat", "Var_a_MSE",
           "Var_b_MSE", "Var_e_MSE", "G_MSE", "selected_K", "runtime",
           "iterations", "group_ARI", "z_R2", "unit_consistency")
  pieces <- lapply(split(raw, interaction(raw$scenario, raw$method, drop = TRUE)), function(d) {
    out <- data.frame(scenario = d$scenario[1], method = d$method[1], n = nrow(d),
                      convergence_rate = mean(d$converged))
    for (nm in num) {
      x <- d[[nm]]
      out[[paste0(nm, "_mean")]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      out[[paste0(nm, "_MCSE")]] <- if (sum(is.finite(x)) <= 1L) NA_real_ else sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
    }
    out
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

# ------------------------------- CASE 10 -------------------------------------

# Case 10 deliberately separates geometric structure from predictive grouping.
# There are K_type latent predictive types. Each type has modes_per_type X-modes
# that differ along nuisance directions but share the SAME random effect. Each
# X-mode is further split into several observed administrative units. Thus:
#   observed units > geometric modes > predictive types.
# A clustering criterion may prefer the geometric modes, while prediction may
# benefit from pooling modes that share the same random effect.
generate_case10 <- function(p = 50,
                            K_type = 6,
                            modes_per_type = 2,
                            units_per_mode = 5,
                            n_test_per_unit = 30,
                            train_multiplier = 3L,
                            var_a = 2.25,
                            var_b = 0,
                            var_e = 9,
                            type_separation = 2.0,
                            nuisance_separation = 1.25,
                            unit_jitter_sd = 0.08,
                            seed = 12345L) {
  set.seed(seed)
  beta <- rep(1, p)
  model_type <- if (var_b > 0) "RS" else "RI"
  M <- K_type * modes_per_type
  R_unit <- M * units_per_mode
  sz <- balanced_unit_sizes(R_unit, n_test_per_unit, train_multiplier)

  mode_to_type <- rep(seq_len(K_type), each = modes_per_type)
  unit_to_mode <- rep(seq_len(M), each = units_per_mode)
  unit_to_type <- mode_to_type[unit_to_mode]

  # Predictive/type direction occupies many coordinates so it remains visible
  # in p=50 without making the groups perfectly separable.
  v_type <- numeric(p)
  q1 <- min(18L, p)
  v_type[seq_len(q1)] <- seq(1, 0.35, length.out = q1)
  v_type <- standardize_direction(v_type)

  v_nuis <- numeric(p)
  if (p >= 10L) {
    lo <- min(19L, p); hi <- min(30L, p)
    if (lo <= hi) v_nuis[lo:hi] <- seq(1, 0.4, length.out = hi - lo + 1L)
  } else {
    v_nuis[seq_len(p)] <- rev(seq(0.4, 1, length.out = p))
  }
  # Orthogonalize nuisance direction against type direction.
  v_nuis <- v_nuis - sum(v_nuis * v_type) * v_type
  v_nuis <- standardize_direction(v_nuis)

  type_score <- seq(-(K_type - 1) / 2, (K_type - 1) / 2, length.out = K_type)
  type_center <- lapply(seq_len(K_type), function(k) type_separation * type_score[k] * v_type)

  mode_center <- matrix(0, M, p)
  for (m in seq_len(M)) {
    t <- mode_to_type[m]
    within <- ((m - 1L) %% modes_per_type) + 1L
    nuisance_code <- if (modes_per_type == 2L) c(-1, 1)[within] else seq(-1, 1, length.out = modes_per_type)[within]
    mode_center[m, ] <- type_center[[t]] + nuisance_separation * nuisance_code * v_nuis
  }

  # Moderate within-mode spread; later coordinates are somewhat noisier.
  variances <- rep(0.35, p)
  if (p > 30L) variances[31:p] <- 0.55
  Sigma <- ar_cov(p, rho = 0.20, variances = variances)
  unit_jitter <- matrix(rnorm(R_unit * p, sd = unit_jitter_sd), R_unit, p)

  make_X <- function(C, offset) {
    X <- matrix(0, sum(C), p)
    unit_lab <- rep(seq_len(R_unit), C)
    start <- 0L
    for (u in seq_len(R_unit)) {
      set.seed(seed + offset + 1009L * u)
      idx <- (start + 1L):(start + C[u])
      mu <- mode_center[unit_to_mode[u], ] + unit_jitter[u, ]
      X[idx, ] <- draw_mvn(C[u], mu, Sigma)
      start <- start + C[u]
    }
    list(X = X, unit = unit_lab)
  }

  tr <- make_X(sz$train, 100000L)
  rf <- make_X(sz$ref,   200000L)
  te <- make_X(sz$test,  300000L)

  true_tr <- unit_to_type[tr$unit]
  true_te <- unit_to_type[te$unit]

  # Random effects belong to predictive TYPE, not to administrative unit or
  # geometric mode. This is the central regrouping feature of Case 10.
  set.seed(seed + 400001L)
  a <- rnorm(K_type, sd = sqrt(var_a))
  B <- if (var_b > 0) matrix(rnorm(K_type * p, sd = sqrt(var_b)), K_type, p) else matrix(0, K_type, p)
  e_tr <- rnorm(nrow(tr$X), sd = sqrt(var_e))
  e_te <- rnorm(nrow(te$X), sd = sqrt(var_e))
  re_tr <- a[true_tr] + rowSums(tr$X * B[true_tr, , drop = FALSE])
  re_te <- a[true_te] + rowSums(te$X * B[true_te, , drop = FALSE])
  y_tr <- 1 + drop(tr$X %*% beta) + re_tr + e_tr
  y_te <- 1 + drop(te$X %*% beta) + re_te + e_te

  list(
    scenario = "case10_nested_predictive_types",
    X_train = tr$X, X_ref = rf$X, X_test = te$X,
    y_train = y_tr, y_test = y_te, beta = beta, model_type = model_type,
    var_a_metric = var_a, var_b_metric = var_b, var_e = var_e,
    G_true = if (var_b > 0) diag(c(var_a, rep(var_b, p))) else NULL,
    oracle_train = true_tr, oracle_test = true_te, oracle_K = K_type,
    true_partition_test = true_te,
    observed_train = tr$unit, observed_test = te$unit,
    unit_train = tr$unit, unit_ref = rf$unit, unit_test = te$unit,
    mode_train = unit_to_mode[tr$unit], mode_test = unit_to_mode[te$unit],
    unit_to_type = unit_to_type, unit_to_mode = unit_to_mode,
    K_type = K_type, K_mode = M, R_unit = R_unit,
    a_true = a, B_true = B
  )
}

run_case10_comparison <- function(nloop = 50,
                                  p = 50,
                                  K_type = 6,
                                  modes_per_type = 2,
                                  units_per_mode = 5,
                                  n_test_per_unit = 30,
                                  Var.a = 2.25,
                                  Var.b = 0,
                                  Var.e = 9,
                                  tau = 5 * (p + 1L),
                                  K_grid = 2:16,
                                  initial_Cn = 4L,
                                  sa_max_iter = 80,
                                  em_tol = 1e-5,
                                  em_max_iter = 500,
                                  seed = 12345L,
                                  methods = c("ORACLE", "OBS", "LM", "KM", "GMM", "CPF", "BLM", "CIRG-U"),
                                  keep_extra = FALSE,
                                  verbose = FALSE) {
  raw_all <- list(); extra_all <- if (keep_extra) list() else NULL
  for (r in seq_len(nloop)) {
    s <- as.integer(seed + 1000003L * r)
    dat <- generate_case10(p, K_type, modes_per_type, units_per_mode,
                           n_test_per_unit, 3L, Var.a, Var.b, Var.e, seed = s)
    ans <- evaluate_generated_dataset(dat, methods, tau, K_grid,
                                      initial_Cn, sa_max_iter, em_tol, em_max_iter,
                                      s + 77L, verbose = verbose)
    d <- ans$metrics; d$replication <- r
    raw_all[[r]] <- d
    if (keep_extra) extra_all[[r]] <- ans$extra
    if (verbose) cat("case10 replication ", r, " done\n", sep = "")
  }
  raw <- do.call(rbind, raw_all)
  list(raw = raw, summary = summarize_structured_comparison(raw), extra = extra_all)
}

# ------------------------------- CASE 11 -------------------------------------

# Continuous heterogeneity. Unit i has z_i in [-1,1]. There is no true discrete
# partition. X and random effects vary smoothly with z. OBS is a coarse noisy
# 4-bin administrative categorization. CPF/BLM receive repeated-measure unit IDs
# because their principles estimate/discretize unit heterogeneity. CIRG sees only
# X and chooses its finite predictive partition through IMSPE.
generate_case11 <- function(p = 50,
                            R_unit = 100,
                            n_test_per_unit = 20,
                            train_multiplier = 3L,
                            obs_bins = 4L,
                            proxy_sd = 0.35,
                            var_a_target = 2.25,
                            var_b_target = 0,
                            var_e = 9,
                            within_var = 0.30,
                            seed = 12345L) {
  set.seed(seed)
  beta <- rep(1, p)
  model_type <- if (var_b_target > 0) "RS" else "RI"
  sz <- balanced_unit_sizes(R_unit, n_test_per_unit, train_multiplier)

  # Stratified random z values cover the support more evenly than iid Uniform,
  # while retaining replicate-to-replicate randomness.
  z <- ((seq_len(R_unit) - 0.5) / R_unit) * 2 - 1
  z <- pmin(1, pmax(-1, z + rnorm(R_unit, sd = 0.04)))
  z <- sample(z, R_unit, replace = FALSE)

  mu_fun <- function(zv) {
    out <- numeric(p)
    blocks <- split(seq_len(min(p, 30L)), ceiling(seq_len(min(p, 30L)) / 10))
    if (length(blocks) >= 1L) out[blocks[[1]]] <- seq(1.3, 0.5, length.out = length(blocks[[1]])) * zv
    if (length(blocks) >= 2L) out[blocks[[2]]] <- seq(1.0, 0.4, length.out = length(blocks[[2]])) * (zv^2 - 1/3)
    if (length(blocks) >= 3L) out[blocks[[3]]] <- seq(0.9, 0.35, length.out = length(blocks[[3]])) * sin(pi * zv)
    out
  }
  centers <- t(vapply(z, mu_fun, numeric(p)))

  variances <- rep(within_var, p)
  if (p > 30L) variances[31:p] <- within_var * 1.5
  Sigma <- ar_cov(p, rho = 0.15, variances = variances)

  make_X <- function(C, offset) {
    X <- matrix(0, sum(C), p); unit_lab <- rep(seq_len(R_unit), C)
    start <- 0L
    for (u in seq_len(R_unit)) {
      set.seed(seed + offset + 1013L * u)
      idx <- (start + 1L):(start + C[u])
      X[idx, ] <- draw_mvn(C[u], centers[u, ], Sigma)
      start <- start + C[u]
    }
    list(X = X, unit = unit_lab)
  }

  tr <- make_X(sz$train, 100000L)
  rf <- make_X(sz$ref,   200000L)
  te <- make_X(sz$test,  300000L)

  # Smooth continuous intercept heterogeneity, centered and scaled to a stable
  # empirical variance across replications.
  set.seed(seed + 400001L)
  a_raw <- 1.25 * z + 0.75 * sin(pi * z) + rnorm(R_unit, sd = 0.20)
  a <- a_raw - mean(a_raw)
  if (sd(a) > 1e-10) a <- a / sd(a) * sqrt(var_a_target)

  B <- matrix(0, R_unit, p)
  if (var_b_target > 0) {
    # Smooth slope heterogeneity is spread over all p coordinates so the
    # diagonal-common-variance RS fitter remains a reasonable approximation.
    load <- seq(1, 0.25, length.out = p)
    for (j in seq_len(p)) {
      base <- load[j] * (0.65 * z + 0.35 * sin(pi * z + 2 * pi * j / p))
      B[, j] <- base + rnorm(R_unit, sd = 0.20)
      B[, j] <- B[, j] - mean(B[, j])
    }
    current <- mean(apply(B, 2, var))
    if (is.finite(current) && current > 1e-12) B <- B * sqrt(var_b_target / current)
  }

  unit_tr <- tr$unit; unit_te <- te$unit
  e_tr <- rnorm(nrow(tr$X), sd = sqrt(var_e))
  e_te <- rnorm(nrow(te$X), sd = sqrt(var_e))
  re_tr <- a[unit_tr] + rowSums(tr$X * B[unit_tr, , drop = FALSE])
  re_te <- a[unit_te] + rowSums(te$X * B[unit_te, , drop = FALSE])
  y_tr <- 1 + drop(tr$X %*% beta) + re_tr + e_tr
  y_te <- 1 + drop(te$X %*% beta) + re_te + e_te

  # Coarse noisy administrative grouping based on an imperfect proxy for z.
  set.seed(seed + 500001L)
  proxy <- z + rnorm(R_unit, sd = proxy_sd)
  cuts <- unique(as.numeric(stats::quantile(proxy, probs = seq(0, 1, length.out = obs_bins + 1L), type = 8)))
  if (length(cuts) < 3L) cuts <- seq(min(proxy) - 1e-8, max(proxy) + 1e-8, length.out = obs_bins + 1L)
  obs_unit <- cut(proxy, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  obs_unit <- as.integer(factor(obs_unit))
  obs_tr <- obs_unit[unit_tr]; obs_te <- obs_unit[unit_te]

  var_a_emp <- stats::var(a)
  var_b_emp <- if (var_b_target > 0) mean(apply(B, 2, var)) else 0
  oracle_pred <- 1 + drop(te$X %*% beta) + re_te

  list(
    scenario = "case11_continuous_heterogeneity",
    X_train = tr$X, X_ref = rf$X, X_test = te$X,
    y_train = y_tr, y_test = y_te, beta = beta, model_type = model_type,
    var_a_metric = var_a_emp, var_b_metric = var_b_emp, var_e = var_e,
    G_true = NULL,
    oracle_pred = oracle_pred, oracle_K = NA_integer_,
    observed_train = obs_tr, observed_test = obs_te,
    unit_train = unit_tr, unit_ref = rf$unit, unit_test = unit_te,
    z_unit = z, z_train = z[unit_tr], z_test = z[unit_te],
    a_true = a, B_true = B, observed_unit_group = obs_unit,
    R_unit = R_unit
  )
}

run_case11_comparison <- function(nloop = 50,
                                  p = 50,
                                  R_unit = 100,
                                  n_test_per_unit = 20,
                                  obs_bins = 4L,
                                  proxy_sd = 0.35,
                                  Var.a = 2.25,
                                  Var.b = 0,
                                  Var.e = 9,
                                  tau = 5 * (p + 1L),
                                  K_grid = 2:15,
                                  initial_Cn = 4L,
                                  sa_max_iter = 80,
                                  em_tol = 1e-5,
                                  em_max_iter = 500,
                                  seed = 12345L,
                                  methods = c("ORACLE", "UNIT", "OBS", "LM", "KM", "GMM", "CPF", "BLM", "CIRG-U"),
                                  keep_extra = FALSE,
                                  verbose = FALSE) {
  raw_all <- list(); extra_all <- if (keep_extra) list() else NULL
  for (r in seq_len(nloop)) {
    s <- as.integer(seed + 1000003L * r)
    dat <- generate_case11(p, R_unit, n_test_per_unit, 3L, obs_bins, proxy_sd,
                           Var.a, Var.b, Var.e, seed = s)
    ans <- evaluate_generated_dataset(dat, methods, tau, K_grid,
                                      initial_Cn, sa_max_iter, em_tol, em_max_iter,
                                      s + 91L, verbose = verbose)
    d <- ans$metrics; d$replication <- r
    raw_all[[r]] <- d
    if (keep_extra) extra_all[[r]] <- ans$extra
    if (verbose) cat("case11 replication ", r, " done\n", sep = "")
  }
  raw <- do.call(rbind, raw_all)
  list(raw = raw, summary = summarize_structured_comparison(raw), extra = extra_all)
}

# No experiment runs automatically when sourced.
