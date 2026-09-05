# ROG_method_comparison.R
# -----------------------------------------------------------------------------
# Baseline-comparison extension for the corrected CIRG simulation framework.
# Source ROG_aligned_fixed.R first (done automatically below).
#
# Methods:
#   ORACLE : correctly specified LMM using the true data-generating groups.
#   OBS    : LMM using deliberately misspecified observed group labels.
#   LM     : pooled linear model, ignoring grouping.
#   KM     : X-only K-means + LMM; K chosen by a BIC-like clustering criterion.
#   GMM    : Gaussian-mixture clustering of X + LMM (mclust; optional package).
#   CPF    : Ma-Huang-style concave pairwise fusion adaptation for repeated
#            grouped observations. Pairwise MCP fusion is applied to observed-
#            group intercepts via ADMM, then the fused labels are used in the
#            same RI/RS LMM fitting/prediction pipeline as all other baselines.
#   BLM    : Bonhomme-Lamadon-Manresa-style two-step discretization. Informative
#            group-level moments are estimated in step 1, K-means discretizes
#            those moments in step 2, then the same RI/RS LMM is fitted.
#   CIRG   : proposed RASC + SGA, model-matched IMSPE search.
#
# IMPORTANT:
# CPF and BLM here are LMM-compatible adaptations, not literal replications of
# the original papers' model-specific software. The comparison isolates their
# grouping principles while holding the final LMM estimator and MSPE evaluation
# fixed across methods.
# -----------------------------------------------------------------------------

if (!exists("run_cirg_simulation", mode = "function")) {
  if (!file.exists("ROG_aligned_fixed.R")) stop("ROG_aligned_fixed.R not found.")
  source("ROG_aligned_fixed.R")
}

# ------------------------------ label helpers --------------------------------

compact_label_pair <- function(train_labels, test_labels) {
  train_labels <- as.integer(train_labels)
  test_labels <- as.integer(test_labels)
  lev <- sort(unique(train_labels))
  tr <- match(train_labels, lev)
  te <- match(test_labels, lev)
  if (anyNA(te)) stop("Test labels contain levels absent from training labels.")
  list(train = tr, test = te, K = length(lev), levels = lev)
}

labels_to_sizes <- function(labels) {
  as.integer(tabulate(as.integer(labels), nbins = max(labels)))
}

sort_by_labels <- function(X, y, labels) {
  o <- order(labels)
  ls <- as.integer(labels[o])
  list(X = X[o, , drop = FALSE], y = y[o], labels = ls,
       sizes = labels_to_sizes(ls), order = o)
}

adjusted_rand_index <- function(a, b) {
  a <- as.integer(factor(a))
  b <- as.integer(factor(b))
  tab <- table(a, b)
  choose2 <- function(x) x * (x - 1) / 2
  sum_nij <- sum(choose2(tab))
  sum_ai <- sum(choose2(rowSums(tab)))
  sum_bj <- sum(choose2(colSums(tab)))
  n <- length(a)
  den_pairs <- choose2(n)
  if (den_pairs <= 0) return(NA_real_)
  expected <- sum_ai * sum_bj / den_pairs
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (abs(denom) < 1e-14) return(1)
  (sum_nij - expected) / denom
}

# -------------------------- observed group errors -----------------------------

permute_selected_labels <- function(labels, rho, seed) {
  labels <- as.integer(labels)
  if (rho <= 0) return(labels)
  set.seed(seed)
  n <- length(labels)
  m <- max(2L, min(n, as.integer(round(rho * n))))
  idx <- sample.int(n, m)
  if (length(unique(labels[idx])) <= 1L) {
    # fall back to reassignment when the selected block contains one label only
    K <- max(labels)
    for (i in idx) {
      cand <- setdiff(seq_len(K), labels[i])
      labels[i] <- sample(cand, 1L)
    }
  } else {
    labels[idx] <- sample(labels[idx], length(idx), replace = FALSE)
  }
  labels
}

make_observed_groups <- function(true_train, true_test, X_train, X_test,
                                 type = c("correct", "contam", "merge", "split"),
                                 rho = 0.25, merge_factor = 2L, seed = 1L) {
  type <- match.arg(type)
  true_train <- as.integer(true_train)
  true_test <- as.integer(true_test)
  R <- max(true_train)

  if (type == "correct") {
    tr <- true_train; te <- true_test
  } else if (type == "contam") {
    # Observation-level label contamination. Within each split the selected
    # labels are permuted, preserving the marginal label counts as far as
    # possible while breaking the true covariance grouping.
    tr <- permute_selected_labels(true_train, rho, seed + 11L)
    te <- permute_selected_labels(true_test,  rho, seed + 29L)
  } else if (type == "merge") {
    merge_factor <- max(2L, as.integer(merge_factor))
    map <- ceiling(seq_len(R) / merge_factor)
    tr <- map[true_train]
    te <- map[true_test]
  } else {
    # Over-grouping: each true group is split into two observed labels using a
    # training-defined threshold in X_1; the same threshold is applied at test.
    tr <- integer(length(true_train))
    te <- integer(length(true_test))
    for (g in seq_len(R)) {
      itr <- which(true_train == g)
      ite <- which(true_test == g)
      cut <- stats::median(X_train[itr, 1])
      tr[itr] <- 2L * (g - 1L) + 1L + as.integer(X_train[itr, 1] > cut)
      te[ite] <- 2L * (g - 1L) + 1L + as.integer(X_test[ite, 1] > cut)
    }
  }

  z <- compact_label_pair(tr, te)
  list(train = z$train, test = z$test, K = z$K,
       type = type, rho = rho, merge_factor = merge_factor)
}

# -------------------------- common LMM prediction -----------------------------

fit_lmm_from_labels <- function(X_train, y_train, labels_train,
                                model_type = c("RI", "RS"),
                                em_tol = 1e-5, em_max_iter = 500) {
  model_type <- match.arg(model_type)
  s <- sort_by_labels(X_train, y_train, labels_train)
  fit <- if (model_type == "RI") {
    fit_ri_lmm_cpp(s$X, s$y, s$sizes, em_tol, em_max_iter)
  } else {
    fit_rs_lmm_cpp(s$X, s$y, s$sizes, em_tol, em_max_iter)
  }
  fit
}

predict_lmm_from_labels <- function(fit, X_test, labels_test,
                                    model_type = c("RI", "RS")) {
  model_type <- match.arg(model_type)
  if (model_type == "RI") {
    predict_ri_known_groups(fit, X_test, labels_test)
  } else {
    predict_rs_known_groups(fit, X_test, labels_test)
  }
}

fit_predict_labeled <- function(X_train, y_train, X_test,
                                labels_train, labels_test,
                                model_type, em_tol, em_max_iter) {
  z <- compact_label_pair(labels_train, labels_test)
  fit <- fit_lmm_from_labels(X_train, y_train, z$train,
                             model_type, em_tol, em_max_iter)
  pred <- predict_lmm_from_labels(fit, X_test, z$test, model_type)
  list(fit = fit, pred = pred, labels_train = z$train,
       labels_test = z$test, K = z$K)
}

# -------------------------------- pooled LM -----------------------------------

fit_pooled_lm <- function(X, y) {
  F <- cbind(1, X)
  beta <- tryCatch(qr.solve(F, y), error = function(e) drop(MASS::ginv(crossprod(F)) %*% crossprod(F, y)))
  resid <- drop(y - F %*% beta)
  ve <- mean(resid^2)
  list(beta = beta, Var.a = 0, Var.b = NA_real_, Var.e = ve,
       G = NULL, u_hat = numeric(0), converged = TRUE, iterations = 1L)
}

predict_pooled_lm <- function(fit, Xnew) drop(cbind(1, Xnew) %*% fit$beta)

# ------------------------------- KM baseline ---------------------------------

fit_kmeans_bic <- function(X, K_grid, seed = 1L, tau = 2L,
                           batch_size = 1024, max_iters = 50) {
  n <- nrow(X); p <- ncol(X)
  tau <- max(2L, as.integer(tau))
  K_grid <- sort(unique(as.integer(K_grid)))
  K_grid <- K_grid[K_grid >= 2L & K_grid <= floor(n / tau)]
  if (!length(K_grid)) stop("No admissible K in K_grid.")
  best <- NULL
  rows <- list()
  for (K in K_grid) {
    cl <- try(cluster_x(seed + 1009L * K, X, rep(0, n), K,
                        tau = tau, batch_size = batch_size,
                        max_iters = max_iters), silent = TRUE)
    if (inherits(cl, "try-error")) next
    lab <- cl$labels
    wss <- 0
    for (k in seq_len(cl$K)) {
      idx <- which(lab == k)
      cen <- colMeans(X[idx, , drop = FALSE])
      wss <- wss + sum((X[idx, , drop = FALSE] - rep(cen, each = length(idx)))^2)
    }
    crit <- n * log(max(wss / n, 1e-12)) + cl$K * p * log(n)
    rows[[length(rows) + 1L]] <- data.frame(K = cl$K, criterion = crit)
    if (is.null(best) || crit < best$criterion) best <- list(cluster = cl, criterion = crit)
  }
  if (is.null(best)) stop("All K-means candidates failed.")
  best$trace <- do.call(rbind, rows)
  best
}

# ------------------------------- GMM baseline --------------------------------

fit_gmm_bic <- function(X, K_grid, seed = 1L,
                        modelNames = c("EII", "VII", "EEI", "VEI")) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("GMM baseline requires package 'mclust'. Install with install.packages('mclust').")
  }
  set.seed(seed)
  fit <- mclust::Mclust(X, G = K_grid, modelNames = modelNames, verbose = FALSE)
  if (is.null(fit$classification)) stop("mclust failed to return a classification.")
  list(model = fit, labels = as.integer(fit$classification), K = fit$G)
}

predict_gmm_labels <- function(gmm_fit, Xnew) {
  as.integer(predict(gmm_fit$model, newdata = Xnew)$classification)
}

# ---------------------- CPF-style MCP pairwise fusion -------------------------

pair_index <- function(G) {
  if (G < 2L) return(matrix(integer(0), ncol = 2L))
  t(utils::combn(seq_len(G), 2L))
}

make_difference_matrix <- function(G) {
  pr <- pair_index(G)
  D <- matrix(0, nrow(pr), G)
  if (nrow(pr)) {
    D[cbind(seq_len(nrow(pr)), pr[, 1])] <- 1
    D[cbind(seq_len(nrow(pr)), pr[, 2])] <- -1
  }
  list(D = D, pairs = pr)
}

soft_threshold <- function(x, a) sign(x) * pmax(abs(x) - a, 0)

prox_mcp <- function(delta, lambda, theta = 1, gamma = 3) {
  if (gamma <= 1 / theta) stop("MCP requires gamma > 1/theta.")
  out <- delta
  ind <- abs(delta) <= gamma * lambda
  out[ind] <- soft_threshold(delta[ind], lambda / theta) / (1 - 1 / (gamma * theta))
  out
}

union_find_components <- function(G, pairs, fused) {
  parent <- seq_len(G)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  unite <- function(a, b) {
    ra <- find_root(a); rb <- find_root(b)
    if (ra != rb) parent[rb] <<- ra
  }
  if (length(fused)) {
    for (e in which(fused)) unite(pairs[e, 1], pairs[e, 2])
  }
  roots <- vapply(seq_len(G), find_root, integer(1))
  as.integer(factor(roots))
}

cpf_admm_once <- function(X, y, unit_labels, lambda,
                          theta = 1, gamma = 3,
                          tol = 1e-3, max_iter = 1000,
                          ridge = 1e-8) {
  unit_labels <- as.integer(factor(unit_labels))
  G <- max(unit_labels); p <- ncol(X); N <- nrow(X)
  Z <- matrix(0, N, G)
  Z[cbind(seq_len(N), unit_labels)] <- 1
  H <- cbind(X, Z)  # no global intercept: group intercepts play that role
  di <- make_difference_matrix(G)
  D <- di$D; pairs <- di$pairs; E <- nrow(D)

  # Initial OLS/ridge estimate.
  A0 <- crossprod(H) + diag(ridge, ncol(H))
  coef <- tryCatch(solve(A0, crossprod(H, y)), error = function(e) drop(MASS::ginv(A0) %*% crossprod(H, y)))
  beta <- coef[seq_len(p)]
  mu <- coef[p + seq_len(G)]
  eta <- if (E) drop(D %*% mu) else numeric(0)
  v <- numeric(E)

  penalty_block <- matrix(0, p + G, p + G)
  if (E) penalty_block[p + seq_len(G), p + seq_len(G)] <- theta * crossprod(D)
  A <- crossprod(H) + penalty_block + diag(ridge, p + G)
  Hty <- crossprod(H, y)

  converged <- FALSE
  used <- 0L
  for (it in seq_len(max_iter)) {
    rhs <- Hty
    if (E) rhs[p + seq_len(G)] <- rhs[p + seq_len(G)] + theta * crossprod(D, eta - v / theta)
    coef_new <- tryCatch(solve(A, rhs), error = function(e) drop(MASS::ginv(A) %*% rhs))
    beta_new <- coef_new[seq_len(p)]
    mu_new <- coef_new[p + seq_len(G)]

    if (E) {
      dmu <- drop(D %*% mu_new)
      delta <- dmu + v / theta
      eta_new <- prox_mcp(delta, lambda, theta, gamma)
      r <- dmu - eta_new
      s <- theta * drop(crossprod(D, eta_new - eta))
      v <- v + theta * r
      primal <- sqrt(mean(r^2))
      dual <- sqrt(mean(s^2))
      eta <- eta_new
    } else {
      primal <- dual <- 0
    }
    beta <- beta_new; mu <- mu_new
    used <- it
    if (max(primal, dual) < tol) { converged <- TRUE; break }
  }

  fitted <- drop(X %*% beta + mu[unit_labels])
  rss <- sum((y - fitted)^2)
  fused <- if (E) abs(eta) < max(1e-7, tol * 0.1) else logical(0)
  group_map <- union_find_components(G, pairs, fused)
  list(beta = beta, mu = mu, eta = eta, group_map = group_map,
       K = max(group_map), rss = rss, converged = converged,
       iterations = used)
}

fit_cpf_bic <- function(X, y, observed_labels,
                        lambda_grid = NULL,
                        theta = 1, gamma = 3, bic_c = 10,
                        tol = 1e-3, max_iter = 1000) {
  observed_labels <- as.integer(factor(observed_labels))
  G <- max(observed_labels); p <- ncol(X); N <- nrow(X)
  if (G == 1L) return(list(group_map = 1L, K = 1L, trace = NULL))

  # Scale lambda to the initial spread of group residual means.
  F <- cbind(1, X)
  b0 <- tryCatch(qr.solve(F, y), error = function(e) drop(MASS::ginv(crossprod(F)) %*% crossprod(F, y)))
  r0 <- drop(y - F %*% b0)
  mu0 <- tapply(r0, observed_labels, mean)
  maxdiff <- max(abs(outer(mu0, mu0, "-")))
  if (!is.finite(maxdiff) || maxdiff < 1e-8) maxdiff <- stats::sd(y)
  if (!is.finite(maxdiff) || maxdiff < 1e-8) maxdiff <- 1
  if (is.null(lambda_grid)) lambda_grid <- seq(0.03, 1.20, length.out = 16) * maxdiff

  Cn <- bic_c * log(log(max(G + p, 3)))
  best <- NULL; tr <- list()
  for (lam in lambda_grid) {
    fit <- cpf_admm_once(X, y, observed_labels, lam, theta, gamma, tol, max_iter)
    bic <- log(max(fit$rss / N, 1e-12)) + Cn * log(N) / N * (fit$K + p)
    tr[[length(tr) + 1L]] <- data.frame(lambda = lam, K = fit$K, BIC = bic,
                                         converged = fit$converged,
                                         iterations = fit$iterations)
    if (is.null(best) || bic < best$BIC) best <- c(fit, list(lambda = lam, BIC = bic))
  }
  best$trace <- do.call(rbind, tr)
  best
}

# ---------------- BLM-style informative-moment discretization ----------------

ridge_unit_effects <- function(X, y, unit_labels, pooled_beta, ridge = 1) {
  unit_labels <- as.integer(factor(unit_labels))
  G <- max(unit_labels); p <- ncol(X)
  F <- cbind(1, X)
  resid <- drop(y - F %*% pooled_beta)
  U <- matrix(0, G, p + 1L)
  for (g in seq_len(G)) {
    idx <- which(unit_labels == g)
    Fg <- F[idx, , drop = FALSE]
    rg <- resid[idx]
    A <- crossprod(Fg) + diag(ridge, p + 1L)
    U[g, ] <- tryCatch(solve(A, crossprod(Fg, rg)),
                        error = function(e) drop(MASS::ginv(A) %*% crossprod(Fg, rg)))
  }
  U
}

blm_informative_moments <- function(X, y, unit_labels,
                                    model_type = c("RI", "RS"), ridge = 1) {
  model_type <- match.arg(model_type)
  unit_labels <- as.integer(factor(unit_labels))
  F <- cbind(1, X)
  pooled_beta <- tryCatch(qr.solve(F, y), error = function(e) drop(MASS::ginv(crossprod(F)) %*% crossprod(F, y)))
  if (model_type == "RI") {
    resid <- drop(y - F %*% pooled_beta)
    M <- matrix(tapply(resid, unit_labels, mean), ncol = 1)
    colnames(M) <- "mean_residual"
  } else {
    M <- ridge_unit_effects(X, y, unit_labels, pooled_beta, ridge)
  }
  # Standardize informative moments before k-means; zero-variance coordinates
  # are retained as zeros rather than producing NaNs.
  mu <- colMeans(M)
  s <- apply(M, 2, sd)
  s[!is.finite(s) | s < 1e-8] <- 1
  Mz <- sweep(sweep(M, 2, mu, "-"), 2, s, "/")
  list(M = Mz, raw = M, pooled_beta = pooled_beta)
}

fit_blm_discretization <- function(X, y, observed_labels,
                                   model_type = c("RI", "RS"),
                                   K_grid = NULL, seed = 1L, ridge = 1) {
  model_type <- match.arg(model_type)
  observed_labels <- as.integer(factor(observed_labels))
  G <- max(observed_labels)
  mom <- blm_informative_moments(X, y, observed_labels, model_type, ridge)
  M <- mom$M; d <- ncol(M)
  if (G <= 2L) return(list(group_map = seq_len(G), K = G, moments = mom, trace = NULL))
  if (is.null(K_grid)) K_grid <- seq_len(min(G - 1L, max(2L, floor(sqrt(G)) + 3L)))
  K_grid <- sort(unique(K_grid[K_grid >= 1L & K_grid <= G]))
  best <- NULL; tr <- list()
  for (K in K_grid) {
    if (K == 1L) {
      lab <- rep(1L, G)
      wss <- sum(M^2)
    } else {
      set.seed(seed + 7919L * K)
      km <- stats::kmeans(M, centers = K, nstart = 25, iter.max = 100)
      lab <- km$cluster
      wss <- km$tot.withinss
    }
    # A transparent data-driven approximation criterion for the first-step
    # discretization. This is BLM-style rather than their model-specific rule.
    crit <- G * log(max(wss / max(G, 1), 1e-12)) + K * d * log(max(G, 2))
    tr[[length(tr) + 1L]] <- data.frame(K = K, criterion = crit)
    if (is.null(best) || crit < best$criterion) best <- list(group_map = as.integer(lab), K = K, criterion = crit)
  }
  best$moments <- mom
  best$trace <- do.call(rbind, tr)
  best
}

# ---------------------------- result utilities --------------------------------

metric_row <- function(method, fit, pred, y_test, beta, var_a, var_b, var_e,
                       G_true, selected_K, runtime, ari = NA_real_,
                       objective = NA_real_) {
  d <- method_metrics(method, fit, pred, y_test, beta, var_a, var_b, var_e,
                      G_true, selected_K, objective, runtime)
  d$group_ARI <- ari
  d
}

summarize_comparison <- function(raw) {
  num <- c("MSPE", "beta_MSE", "beta_SSE", "intercept_MSE",
           "Var_a_hat", "Var_b_hat", "Var_e_hat", "Var_a_MSE",
           "Var_b_MSE", "Var_e_MSE", "G_MSE", "selected_K", "runtime",
           "iterations", "group_ARI")
  scenario <- ifelse(raw$mis_type == "contam", paste0("contam_", sprintf("%g", raw$rho)), raw$mis_type)
  parts <- lapply(split(raw, interaction(raw$method, scenario, drop = TRUE)), function(d) {
    out <- data.frame(method = d$method[1], mis_type = d$mis_type[1], rho = d$rho[1],
                      n = nrow(d), convergence_rate = mean(d$converged))
    for (nm in num) {
      x <- d[[nm]]
      out[[paste0(nm, "_mean")]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      out[[paste0(nm, "_MCSE")]] <- if (sum(is.finite(x)) <= 1L) NA_real_ else sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
    }
    out
  })
  out <- do.call(rbind, parts); rownames(out) <- NULL; out
}

# --------------------------- one comparison run ------------------------------

run_comparison_replication <- function(N = 2500, p = 50, R = 20,
                                       dist_x = "case1", groupsize = "large",
                                       var_a = 2.25, var_b = 0, var_e = 9,
                                       mis_type = "contam", rho = 0.25,
                                       merge_factor = 2L,
                                       tau = 5 * (p + 1L),
                                       initial_Cn = 2L,
                                       sa_max_iter = 80,
                                       em_tol = 1e-5, em_max_iter = 500,
                                       seed = 12345L,
                                       methods = c("ORACLE", "OBS", "LM", "KM", "GMM", "CPF", "BLM", "CIRG"),
                                       gmm_models = c("EII", "VII", "EEI", "VEI"),
                                       verbose = FALSE) {
  set.seed(seed)
  methods <- unique(toupper(methods))
  beta <- rep(1, p)
  model_type <- if (var_b > 0) "RS" else "RI"
  m <- N / (10 * R)
  C_test <- generate_groups(R, m, N, groupsize)
  C_train <- as.integer(3L * C_test)
  C_ref <- C_test

  X_train <- generate_covariates(C_train, p, dist_x, "train", seed + 11L)
  X_ref <- generate_covariates(C_ref, p, dist_x, "reference", seed + 23L)
  X_test <- generate_covariates(C_test, p, dist_x, "test", seed + 37L)
  resp <- generate_responses(X_train, X_test, C_train, C_test, beta,
                             var_a, var_b, var_e, seed + 101L, "normal")
  true_tr <- resp$true_group_train; true_te <- resp$true_group_test
  obs <- make_observed_groups(true_tr, true_te, X_train, X_test,
                              mis_type, rho, merge_factor, seed + 313L)
  G_true <- if (var_b > 0) resp$G_true else NULL
  rows <- list(); extra <- list()

  add <- function(d) rows[[length(rows) + 1L]] <<- d

  if ("ORACLE" %in% methods) {
    t0 <- proc.time()[3]
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, true_tr, true_te,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    add(metric_row("ORACLE", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, rt, adjusted_rand_index(fp$labels_test, true_te)))
  }

  if ("OBS" %in% methods) {
    t0 <- proc.time()[3]
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, obs$train, obs$test,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    add(metric_row("OBS", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, rt, adjusted_rand_index(fp$labels_test, true_te)))
  }

  if ("LM" %in% methods) {
    t0 <- proc.time()[3]
    fit <- fit_pooled_lm(X_train, resp$y_train)
    pred <- predict_pooled_lm(fit, X_test)
    rt <- proc.time()[3] - t0
    add(metric_row("LM", fit, pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, 1L, rt, NA_real_))
  }

  maxK <- max(2L, floor(nrow(X_train) / max(2L, tau)))
  K_grid <- seq.int(2L, maxK)

  if ("KM" %in% methods) {
    t0 <- proc.time()[3]
    km <- fit_kmeans_bic(X_train, K_grid, seed + 401L, tau = tau)
    trlab <- km$cluster$labels
    telab <- assign_kmeans(X_test, km$cluster$centroids)
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, trlab, telab,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    add(metric_row("KM", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, rt, adjusted_rand_index(telab, true_te)))
    extra$KM <- km
  }

  if ("GMM" %in% methods) {
    if (requireNamespace("mclust", quietly = TRUE)) {
      t0 <- proc.time()[3]
      gm <- fit_gmm_bic(X_train, K_grid, seed + 503L, gmm_models)
      trlab <- gm$labels
      telab <- predict_gmm_labels(gm, X_test)
      # RS requires enough observations per group. If a GMM component is too
      # small, skip this replicate rather than silently alter the method.
      if (model_type == "RS" && min(tabulate(trlab, nbins = max(trlab))) < p + 1L) {
        warning("GMM produced a group smaller than p+1; GMM row omitted for this replication.")
      } else {
        fp <- fit_predict_labeled(X_train, resp$y_train, X_test, trlab, telab,
                                  model_type, em_tol, em_max_iter)
        rt <- proc.time()[3] - t0
        add(metric_row("GMM", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                       G_true, fp$K, rt, adjusted_rand_index(telab, true_te)))
      }
      extra$GMM <- gm
    } else {
      warning("Package mclust is unavailable; GMM baseline skipped.")
    }
  }

  if ("CPF" %in% methods) {
    t0 <- proc.time()[3]
    cpf <- fit_cpf_bic(X_train, resp$y_train, obs$train,
                       theta = 1, gamma = 3, bic_c = 10,
                       tol = 1e-3, max_iter = 1000)
    cpf_tr <- cpf$group_map[obs$train]
    cpf_te <- cpf$group_map[obs$test]
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, cpf_tr, cpf_te,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    add(metric_row("CPF", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, rt, adjusted_rand_index(cpf_te, true_te)))
    extra$CPF <- cpf
  }

  if ("BLM" %in% methods) {
    t0 <- proc.time()[3]
    blm <- fit_blm_discretization(X_train, resp$y_train, obs$train,
                                  model_type, seed = seed + 607L, ridge = 1)
    blm_tr <- blm$group_map[obs$train]
    blm_te <- blm$group_map[obs$test]
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, blm_tr, blm_te,
                              model_type, em_tol, em_max_iter)
    rt <- proc.time()[3] - t0
    add(metric_row("BLM", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, rt, adjusted_rand_index(blm_te, true_te)))
    extra$BLM <- blm
  }

  if ("CIRG" %in% methods) {
    t0 <- proc.time()[3]
    search <- cirg_search(X_train, resp$y_train, X_ref,
                          initial_Cn = initial_Cn, tau = tau,
                          model_type = model_type,
                          max_iter = sa_max_iter, seed = seed + 701L,
                          em_tol = em_tol, em_max_iter = em_max_iter,
                          verbose = verbose)
    best <- search$best
    fit <- best$imspe_fit
    if (model_type == "RI") {
      pred <- predict_ri_soft(fit, X_test, best$cluster$params)
    } else {
      pred <- predict_rs_soft(fit, X_test, best$cluster$params)
    }
    telab <- predict_soft_labels(X_test, best$cluster$params)
    rt <- proc.time()[3] - t0
    add(metric_row("CIRG", fit, pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, best$K, rt, adjusted_rand_index(telab, true_te), best$objective))
    extra$CIRG <- search
  }

  raw <- do.call(rbind, rows)
  raw$mis_type <- mis_type
  raw$rho <- if (mis_type == "contam") rho else NA_real_
  raw$merge_factor <- merge_factor
  raw$dist_x <- dist_x
  raw$var_b <- var_b
  raw$R_true <- R
  raw$N_test <- N
  list(metrics = raw, observed = obs, extra = extra)
}

run_method_comparison <- function(N_all = 2500, p = 50, R = 20,
                                  Var.e = 9, Var.a = 2.25, Var.b = 0,
                                  nloop = 50, dist_x = "case1",
                                  groupsize = "large",
                                  mis_type = "contam", rho = 0.25,
                                  merge_factor = 2L,
                                  tau = 5 * (p + 1L),
                                  initial_Cn = 2L,
                                  sa_max_iter = 80,
                                  em_tol = 1e-5, em_max_iter = 500,
                                  seed = 12345L,
                                  methods = c("ORACLE", "OBS", "LM", "KM", "GMM", "CPF", "BLM", "CIRG"),
                                  keep_extra = FALSE,
                                  verbose = FALSE) {
  raw_all <- list(); extra_all <- if (keep_extra) list() else NULL; pos <- 0L
  for (N in N_all) {
    for (r in seq_len(nloop)) {
      rep_seed <- as.integer(seed + 1000003L * r + 7919L * N)
      ans <- run_comparison_replication(N, p, R, dist_x, groupsize,
                                        Var.a, Var.b, Var.e,
                                        mis_type, rho, merge_factor,
                                        tau, initial_Cn, sa_max_iter,
                                        em_tol, em_max_iter, rep_seed,
                                        methods, verbose = verbose)
      pos <- pos + 1L
      d <- ans$metrics
      d$replication <- r
      raw_all[[pos]] <- d
      if (keep_extra) extra_all[[pos]] <- ans$extra
      if (verbose) cat("N=", N, " replication=", r, " done\n", sep = "")
    }
  }
  raw <- do.call(rbind, raw_all)
  summary <- summarize_comparison(raw)
  list(raw = raw, summary = summary, extra = extra_all)
}

# No automatic experiment is executed when this file is sourced.
