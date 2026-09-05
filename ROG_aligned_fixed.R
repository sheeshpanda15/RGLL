# ROG_aligned_fixed.R
# -----------------------------------------------------------------------------
# CIRG / RASC simulation code aligned with the paper description.
#
# Statistical alignment and fixes:
#   * CIRG uses Residual-Augmented Supervised Clustering (RASC): cluster on
#     [X, lambda * standardized residual].
#   * CIRG uses soft Gaussian posterior assignment (SGA) out of sample.
#   * Minimum group size tau defaults to 5*(p+1), matching the simulation setup.
#   * RI experiments optimize an RI-specific IMSPE; RS experiments optimize
#     the RS IMSPE.  The objective is now matched to the fitted model.
#   * An independent reference X sample estimates target moments; test X is
#     never used to tune K.
#   * Random-slope covariance is fitted as
#       G = diag(var_a, var_b, ..., var_b),
#     matching the simulation DGP and avoiding a poorly identified full
#     (p+1)x(p+1) covariance matrix when R is small relative to p.
#   * Train/test observations from the same true group share the same random
#     effects.
#   * Random-slope prediction includes both a_k and x' b_k.
#   * EM fitting iterates to convergence.
#
# Important statistical note:
# The C++ estimators are convergence-based ML-EM implementations, not exact
# REML optimizers.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Rcpp)
  library(ClusterR)
  library(MASS)
})

if (!exists("fit_rs_lmm_cpp", mode = "function") ||
    !exists("imspe_ri_cpp", mode = "function")) {
  if (!file.exists("lmm_fast_aligned.cpp")) {
    stop("Put lmm_fast_aligned.cpp in the working directory before sourcing ROG_aligned.R")
  }
  Rcpp::sourceCpp("lmm_fast_aligned.cpp")
}

# ------------------------------- utilities -----------------------------------

make_equicor_sigma <- function(p, rho = 0.5) {
  S <- matrix(rho, p, p)
  diag(S) <- 1
  S
}

make_ar_sigma <- function(p, rho = 0.8) {
  outer(seq_len(p), seq_len(p), function(i, j) rho^abs(i - j))
}

make_case2_sigma <- function(p) {
  make_ar_sigma(p, 0.5)
}

draw_mvn <- function(n, mu, Sigma) {
  z <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  matrix(z, nrow = n, ncol = length(mu))
}

generate_groups <- function(R, m, N, imbalance = c("large", "small")) {
  imbalance <- match.arg(imbalance)
  if (N <= R * m) stop("N must be greater than R*m.")
  Vu <- if (imbalance == "large") 300 * N / (5 * R) else N / (5 * R)
  adjusted_sum <- N - R * m
  raw <- abs(rnorm(R, mean = N / R, sd = Vu))
  raw <- raw / sum(raw) * adjusted_sum
  out <- ceiling(raw + m)
  difference <- N - sum(out)
  while (difference != 0) {
    if (difference > 0) {
      j <- sample.int(R, 1)
      out[j] <- out[j] + 1
      difference <- difference - 1
    } else {
      j <- sample.int(R, 1)
      if (out[j] - 1 > m) {
        out[j] <- out[j] - 1
        difference <- difference + 1
      }
    }
  }
  as.integer(out)
}

# Generate X for one of the nine paper scenarios.  "reference" is an
# independent draw from the target/test covariate distribution and is used only
# to estimate W and W_k for the optimization criterion.
generate_covariates <- function(C, p, dist_x, split = c("train", "reference", "test"),
                                seed = 1L,
                                case5_df = 10,
                                case6_logvar = 0.2) {
  split <- match.arg(split)
  R <- length(C)
  SC <- c(0L, cumsum(C))
  X <- matrix(0, nrow = sum(C), ncol = p)
  Sigma_case2 <- make_case2_sigma(p)
  Sigma_ar <- make_ar_sigma(p, 0.8)
  split_offset <- switch(split, train = 0L, reference = 10000000L, test = 20000000L)

  for (i in seq_len(R)) {
    set.seed(as.integer(seed + split_offset + i * 1009L))
    idx <- (SC[i] + 1L):SC[i + 1L]
    ni <- C[i]

    if (dist_x == "case1") {
      Xi <- matrix(runif(ni * p, -1, 1), ni, p)
    } else if (dist_x == "case2") {
      Xi <- draw_mvn(ni, rep(0, p), Sigma_case2)
    } else if (dist_x == "case3") {
      Xi <- matrix(runif(ni * p, -1.55 + i / 20, 0.45 + i / 20), ni, p)
    } else if (dist_x == "case4") {
      Xi <- draw_mvn(ni, rep(-2 + (i - 1) / 5, p), Sigma_case2)
    } else if (dist_x == "case5") {
      # Paper default: t_10.  Use case5_df=3 for the heavy-tail sensitivity run.
      Xi <- matrix(rt(ni * p, df = case5_df), ni, p)
    } else if (dist_x == "case6") {
      # If log(X) ~ N(0, 0.2), the sdlog is sqrt(0.2), not 1.
      Xi <- matrix(rlnorm(ni * p, meanlog = 0, sdlog = sqrt(case6_logvar)), ni, p)
    } else if (dist_x == "case7") {
      mix <- rbinom(ni, 1, 0.5)
      X1 <- matrix(rnorm(ni * p, -2, 1), ni, p)
      X2 <- matrix(rnorm(ni * p,  2, 1), ni, p)
      Xi <- X1 * (1 - mix) + X2 * mix
    } else if (dist_x == "case8") {
      mu <- if (split == "train") rep(0, p) else rep(1.5, p)
      Xi <- draw_mvn(ni, mu, Sigma_case2)
    } else if (dist_x == "case9") {
      Sigma <- if (split == "train") diag(p) else 2 * Sigma_ar
      Xi <- draw_mvn(ni, rep(0, p), Sigma)
    } else {
      stop("Unknown dist_x: ", dist_x)
    }
    X[idx, ] <- Xi
  }
  X
}

# ------------------------- response generation -------------------------------

draw_group_intercepts <- function(R, var_a, distribution = c("normal", "t3")) {
  distribution <- match.arg(distribution)
  if (var_a <= 0) return(rep(0, R))
  if (distribution == "normal") return(rnorm(R, sd = sqrt(var_a)))
  # t_3 has variance 3, so rescale it to the requested variance.
  rt(R, df = 3) * sqrt(var_a / 3)
}

generate_responses <- function(X_train, X_test, C_train, C_test, beta,
                               var_a = 2.25, var_b = 0, var_e = 9,
                               seed = 1L,
                               intercept_distribution = c("normal", "t3")) {
  intercept_distribution <- match.arg(intercept_distribution)
  R <- length(C_train)
  p <- ncol(X_train)
  stopifnot(length(C_test) == R)
  set.seed(seed)

  a <- draw_group_intercepts(R, var_a, intercept_distribution)
  B <- if (var_b > 0) matrix(rnorm(R * p, sd = sqrt(var_b)), R, p) else matrix(0, R, p)

  g_train <- rep(seq_len(R), C_train)
  g_test  <- rep(seq_len(R), C_test)
  e_train <- rnorm(nrow(X_train), sd = sqrt(var_e))
  e_test  <- rnorm(nrow(X_test),  sd = sqrt(var_e))

  random_train <- a[g_train] + rowSums(X_train * B[g_train, , drop = FALSE])
  random_test  <- a[g_test]  + rowSums(X_test  * B[g_test,  , drop = FALSE])

  y_train <- 1 + drop(X_train %*% beta) + random_train + e_train
  y_test  <- 1 + drop(X_test  %*% beta) + random_test  + e_test

  list(
    y_train = y_train,
    y_test = y_test,
    true_group_train = g_train,
    true_group_test = g_test,
    a = a,
    B = B,
    G_true = diag(c(var_a, rep(var_b, p)))
  )
}

# -------------------------- X-only clustering -------------------------------

normalize_cluster_labels <- function(labels, K) {
  labels <- as.integer(labels)
  if (length(labels) == 0L) stop("Empty clustering result.")
  if (min(labels) == 0L && max(labels) <= K - 1L) labels <- labels + 1L
  if (any(labels < 1L | labels > K)) {
    u <- sort(unique(labels))
    if (length(u) != K) stop("Unexpected cluster labels returned by ClusterR.")
    labels <- match(labels, u)
  }
  labels
}

assign_kmeans <- function(Xnew, centroids) {
  K <- nrow(centroids)
  lab <- ClusterR::predict_KMeans(Xnew, centroids)
  normalize_cluster_labels(lab, K)
}

regularize_covariance <- function(S, ridge = 1e-6) {
  S <- (S + t(S)) / 2
  rr <- ridge
  for (j in 1:8) {
    ans <- try(chol(S), silent = TRUE)
    if (!inherits(ans, "try-error")) return(list(S = S, chol = ans))
    diag(S) <- diag(S) + rr
    rr <- rr * 10
  }
  ee <- eigen(S, symmetric = TRUE)
  ee$values[ee$values < ridge] <- ridge
  S <- ee$vectors %*% diag(ee$values, nrow = length(ee$values)) %*% t(ee$vectors)
  list(S = S, chol = chol((S + t(S)) / 2))
}

fit_group_params <- function(X, labels, ridge = 1e-6) {
  K <- max(labels)
  p <- ncol(X)
  n <- nrow(X)
  means <- vector("list", K)
  covs <- vector("list", K)
  chols <- vector("list", K)
  priors <- numeric(K)

  for (k in seq_len(K)) {
    idx <- which(labels == k)
    Xk <- X[idx, , drop = FALSE]
    priors[k] <- length(idx) / n
    means[[k]] <- colMeans(Xk)
    S <- if (nrow(Xk) <= p + 1L) diag(p) else stats::cov(Xk)
    reg <- regularize_covariance(S + diag(ridge, p), ridge)
    covs[[k]] <- reg$S
    chols[[k]] <- reg$chol
  }
  list(mu = means, Sigma = covs, chol = chols, pi = priors, K = K)
}

soft_group_weights <- function(Xnew, params) {
  n <- nrow(Xnew)
  p <- ncol(Xnew)
  K <- params$K
  logw <- matrix(NA_real_, n, K)

  for (k in seq_len(K)) {
    Rchol <- params$chol[[k]]
    centered <- sweep(Xnew, 2, params$mu[[k]], "-")
    z <- forwardsolve(t(Rchol), t(centered))
    md2 <- colSums(z^2)
    logdet <- 2 * sum(log(diag(Rchol)))
    logw[, k] <- log(max(params$pi[k], .Machine$double.xmin)) -
      0.5 * (p * log(2 * pi) + logdet + md2)
  }

  mx <- apply(logw, 1, max)
  w <- exp(logw - mx)
  w / rowSums(w)
}

cluster_x <- function(setseed, X, y, Cn,
                      lambda = 1,  # retained only for backward compatibility
                      tau = ncol(X) + 1L,
                      batch_size = 1024,
                      num_init = 3,
                      max_iters = 50) {
  n <- nrow(X)
  tau <- max(2L, as.integer(tau))
  maxK <- floor(n / tau)
  if (maxK < 2L) {
    stop("Training sample is too small for two groups under tau=", tau,
         ". Increase N or reduce tau.")
  }
  current_Cn <- min(max(2L, as.integer(Cn)), maxK)

  repeat {
    set.seed(as.integer(setseed + current_Cn * 17L))
    km <- ClusterR::MiniBatchKmeans(
      X, clusters = current_Cn,
      batch_size = min(batch_size, n),
      num_init = num_init,
      max_iters = max_iters,
      initializer = "kmeans++"
    )
    labels <- assign_kmeans(X, km$centroids)
    sizes <- tabulate(labels, nbins = current_Cn)
    if (all(sizes >= tau)) break
    current_Cn <- current_Cn - 1L
    if (current_Cn < 2L) {
      stop("K-means could not produce >=2 clusters satisfying tau=", tau)
    }
  }

  ord <- order(labels)
  list(
    K = current_Cn,
    labels = labels,
    sizes = as.integer(tabulate(labels, nbins = current_Cn)),
    X_sorted = X[ord, , drop = FALSE],
    y_sorted = y[ord],
    order = ord,
    centroids = km$centroids,
    cluster_mode = "x_only"
  )
}

rasc_cluster <- function(setseed, X, y, Cn,
                         lambda = 1,
                         tau = 5 * (ncol(X) + 1L),
                         batch_size = 1024,
                         num_init = 3,
                         max_iters = 50) {
  n <- nrow(X)
  p <- ncol(X)
  tau <- max(2L, as.integer(tau))
  if (tau < p + 1L) warning("tau < p+1: random-slope group estimates may be unstable.")
  maxK <- floor(n / tau)
  if (maxK < 2L) {
    stop("Training sample is too small for two groups under tau=", tau,
         ". Increase N or reduce tau.")
  }
  current_Cn <- min(max(2L, as.integer(Cn)), maxK)

  Xfix <- cbind(1, X)
  beta_ols <- tryCatch(
    qr.solve(Xfix, y),
    error = function(e) drop(MASS::ginv(crossprod(Xfix)) %*% crossprod(Xfix, y))
  )
  resid <- drop(y - Xfix %*% beta_ols)
  resid_sd <- sd(resid)
  rstd <- if (!is.finite(resid_sd) || resid_sd < 1e-12) {
    rep(0, n)
  } else {
    (resid - mean(resid)) / resid_sd
  }
  X_aug <- cbind(X, lambda * rstd)

  repeat {
    set.seed(as.integer(setseed + current_Cn * 17L))
    km <- ClusterR::MiniBatchKmeans(
      X_aug, clusters = current_Cn,
      batch_size = min(batch_size, n),
      num_init = num_init,
      max_iters = max_iters,
      initializer = "kmeans++"
    )
    labels <- ClusterR::predict_KMeans(X_aug, km$centroids)
    labels <- normalize_cluster_labels(labels, current_Cn)
    sizes <- tabulate(labels, nbins = current_Cn)
    if (all(sizes >= tau)) break
    current_Cn <- current_Cn - 1L
    if (current_Cn < 2L) {
      stop("RASC could not produce >=2 clusters satisfying tau=", tau)
    }
  }

  ord <- order(labels)
  list(
    K = current_Cn,
    labels = labels,
    sizes = as.integer(tabulate(labels, nbins = current_Cn)),
    X_sorted = X[ord, , drop = FALSE],
    y_sorted = y[ord],
    order = ord,
    augmented_centroids = km$centroids,
    centroids_x = km$centroids[, seq_len(p), drop = FALSE],
    params = fit_group_params(X, labels),
    beta_ols = beta_ols,
    residual = resid,
    cluster_mode = "rasc"
  )
}

# -------------------------- target moments -----------------------------------

empirical_target_moments <- function(X_reference, assignment) {
  if (is.list(assignment) && !is.null(assignment$K) && !is.null(assignment$mu)) {
    prob <- soft_group_weights(X_reference, assignment)
    hard <- max.col(prob, ties.method = "first")
    K <- assignment$K
  } else {
    prob <- NULL
    hard <- assign_kmeans(X_reference, assignment)
    K <- nrow(assignment)
  }
  Fref <- cbind(1, X_reference)
  nref <- nrow(Fref)

  W <- crossprod(Fref) / nref
  Wk <- vector("list", K)
  mk <- vector("list", K)
  pk <- numeric(K)

  for (k in seq_len(K)) {
    idx <- which(hard == k)
    pk[k] <- length(idx) / nref
    if (length(idx) == 0L) {
      Wk[[k]] <- matrix(0, ncol(Fref), ncol(Fref))
      mk[[k]] <- rep(0, ncol(Fref))
    } else {
      Fk <- Fref[idx, , drop = FALSE]
      Wk[[k]] <- crossprod(Fk) / nref
      mk[[k]] <- colSums(Fk) / nref
    }
  }
  list(W = W, Wk = Wk, mk = mk, pk = pk, hard = hard, prob = prob)
}

# --------------------------- CIRG objective/search ----------------------------

evaluate_cirg_candidate <- function(setseed, X_train, y_train, X_reference,
                                    Cn, lambda, tau,
                                    model_type = c("RI", "RS"),
                                    em_tol = 1e-6, em_max_iter = 100) {
  model_type <- match.arg(model_type)
  cl <- rasc_cluster(setseed, X_train, y_train, Cn, lambda = lambda, tau = tau)
  moments <- empirical_target_moments(X_reference, cl$params)

  fit_obj <- if (model_type == "RI") {
    imspe_ri_cpp(
      cl$X_sorted, cl$y_sorted, as.integer(cl$sizes),
      moments$W, moments$mk, moments$pk,
      tol = em_tol, max_iter = em_max_iter
    )
  } else {
    imspe_rs_cpp(
      cl$X_sorted, cl$y_sorted, as.integer(cl$sizes),
      moments$W, moments$Wk,
      tol = em_tol, max_iter = em_max_iter
    )
  }

  list(
    requested_K = Cn,
    K = cl$K,
    cluster = cl,
    moments = moments,
    objective = fit_obj$objective,
    IMSPE = fit_obj$IMSPE,
    imspe_fit = fit_obj,
    model_type = model_type
  )
}

cirg_search <- function(X_train, y_train, X_reference,
                        initial_Cn = 2,
                        lambda = 1,
                        tau = 5 * (ncol(X_train) + 1L),
                        model_type = c("RI", "RS"),
                        T0 = 50,
                        alpha = 0.95,
                        max_iter = 80,
                        seed = 1L,
                        em_tol = 1e-6,
                        em_max_iter = 100,
                        verbose = FALSE) {
  model_type <- match.arg(model_type)
  tau <- max(2L, as.integer(tau))
  maxK <- floor(nrow(X_train) / tau)
  if (maxK < 2L) stop("No admissible K>=2 under the chosen tau.")
  set.seed(seed)

  current <- evaluate_cirg_candidate(
    seed, X_train, y_train, X_reference,
    min(max(2L, initial_Cn), maxK), lambda, tau, model_type,
    em_tol, em_max_iter
  )
  best <- current
  Tcur <- T0
  trace <- vector("list", max_iter + 1L)
  trace[[1L]] <- data.frame(
    iter = 0L, requested_K = current$requested_K, K = current$K,
    objective = current$objective, IMSPE = current$IMSPE,
    accepted = TRUE, temperature = Tcur
  )

  steps <- c(-3L, -2L, -1L, 1L, 2L, 3L)
  used <- 0L
  for (iter in seq_len(max_iter)) {
    if (Tcur < 1e-4) break
    step <- sample(steps, 1L)
    Ccand <- min(maxK, max(2L, current$K + step))

    cand <- evaluate_cirg_candidate(
      seed + iter * 10007L, X_train, y_train, X_reference,
      Ccand, lambda, tau, model_type, em_tol, em_max_iter
    )

    delta <- cand$objective - current$objective
    accept_prob <- if (delta >= 0) 1 else exp(delta / Tcur)
    accepted <- runif(1) < accept_prob
    if (accepted) current <- cand
    if (cand$objective > best$objective) best <- cand

    Tcur <- Tcur * alpha
    used <- iter
    trace[[iter + 1L]] <- data.frame(
      iter = iter, requested_K = Ccand, K = cand$K,
      objective = cand$objective, IMSPE = cand$IMSPE,
      accepted = accepted, temperature = Tcur
    )
    if (verbose) {
      cat(sprintf("iter=%d model=%s K=%d candK=%d IMSPE=%.6f best=%.6f accepted=%s\n",
                  iter, model_type, current$K, cand$K,
                  cand$IMSPE, best$IMSPE, accepted))
    }
  }

  trace <- do.call(rbind, trace[seq_len(used + 1L)])
  list(best = best, current = current, trace = trace,
       lambda = lambda, tau = tau, maxK = maxK,
       model_type = model_type, cluster_mode = "rasc")
}

# ------------------------------- prediction ----------------------------------

predict_ri_soft <- function(fit, Xnew, group_params) {
  F <- cbind(1, Xnew)
  w <- soft_group_weights(Xnew, group_params)
  fixed <- drop(F %*% fit$beta)
  random <- drop(w %*% as.numeric(fit$u_hat))
  fixed + random
}

predict_rs_soft <- function(fit, Xnew, group_params) {
  F <- cbind(1, Xnew)
  w <- soft_group_weights(Xnew, group_params)
  U <- as.matrix(fit$u_hat)
  usoft <- w %*% U
  drop(F %*% fit$beta) + rowSums(F * usoft)
}

predict_soft_labels <- function(Xnew, group_params) {
  prob <- soft_group_weights(Xnew, group_params)
  max.col(prob, ties.method = "first")
}

predict_ri_kmeans <- function(fit, Xnew, centroids) {
  F <- cbind(1, Xnew)
  g <- assign_kmeans(Xnew, centroids)
  drop(F %*% fit$beta) + as.numeric(fit$u_hat)[g]
}

predict_rs_kmeans <- function(fit, Xnew, centroids) {
  F <- cbind(1, Xnew)
  g <- assign_kmeans(Xnew, centroids)
  U <- as.matrix(fit$u_hat)
  drop(F %*% fit$beta) + rowSums(F * U[g, , drop = FALSE])
}

predict_ri_known_groups <- function(fit, Xnew, group_labels) {
  F <- cbind(1, Xnew)
  drop(F %*% fit$beta) + as.numeric(fit$u_hat)[group_labels]
}

predict_rs_known_groups <- function(fit, Xnew, group_labels) {
  F <- cbind(1, Xnew)
  U <- as.matrix(fit$u_hat)
  drop(F %*% fit$beta) + rowSums(F * U[group_labels, , drop = FALSE])
}

# ----------------------------- metrics/reporting -----------------------------

method_metrics <- function(method, fit, pred, y_test, beta_true,
                           var_a_true, var_b_true, var_e_true, G_true = NULL,
                           selected_K = NA_integer_, objective = NA_real_,
                           runtime = NA_real_) {
  b <- as.numeric(fit$beta)
  slope_err <- b[-1] - beta_true
  G_err <- NA_real_
  if (!is.null(G_true) && !is.null(fit$G)) {
    Ghat <- as.matrix(fit$G)
    G_err <- mean((Ghat - G_true)^2)
  }
  vb_hat <- if (!is.null(fit$Var.b)) as.numeric(fit$Var.b) else NA_real_

  data.frame(
    method = method,
    MSPE = mean((y_test - pred)^2),
    beta_MSE = mean(slope_err^2),
    beta_SSE = sum(slope_err^2),
    intercept_MSE = (b[1] - 1)^2,
    Var_a_hat = as.numeric(fit$Var.a),
    Var_b_hat = vb_hat,
    Var_e_hat = as.numeric(fit$Var.e),
    Var_a_MSE = (as.numeric(fit$Var.a) - var_a_true)^2,
    Var_b_MSE = if (is.finite(vb_hat)) (vb_hat - var_b_true)^2 else NA_real_,
    Var_e_MSE = (as.numeric(fit$Var.e) - var_e_true)^2,
    G_MSE = G_err,
    selected_K = selected_K,
    objective = objective,
    runtime = runtime,
    converged = isTRUE(fit$converged),
    iterations = as.integer(fit$iterations)
  )
}

summarize_results <- function(raw) {
  numeric_cols <- c("MSPE", "beta_MSE", "beta_SSE", "intercept_MSE",
                    "Var_a_hat", "Var_b_hat", "Var_e_hat",
                    "Var_a_MSE", "Var_b_MSE", "Var_e_MSE",
                    "G_MSE", "selected_K", "objective", "runtime", "iterations")
  pieces <- lapply(split(raw, raw$method), function(d) {
    out <- data.frame(method = d$method[1], n = nrow(d),
                      convergence_rate = mean(d$converged))
    for (nm in numeric_cols) {
      x <- d[[nm]]
      out[[paste0(nm, "_mean")]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      out[[paste0(nm, "_MCSE")]] <- if (sum(is.finite(x)) <= 1L) NA_real_ else
        sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
    }
    out
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

# ---------------------------- one replication --------------------------------

run_one_replication <- function(N, p, R, dist_x, groupsize,
                                var_a, var_b, var_e,
                                lambda, tau,
                                seed,
                                initial_Cn = 2,
                                T0 = 50, alpha = 0.95, sa_max_iter = 80,
                                em_tol = 1e-6, em_max_iter = 100,
                                case5_df = 10, case6_logvar = 0.2,
                                intercept_distribution = "normal",
                                verbose = FALSE) {
  set.seed(seed)
  beta <- rep(1, p)
  m <- N / (10 * R)
  C_test <- generate_groups(R, m, N, groupsize)
  C_train <- as.integer(3L * C_test)
  C_ref <- C_test

  X_train <- generate_covariates(C_train, p, dist_x, "train", seed + 11L,
                                 case5_df, case6_logvar)
  X_ref <- generate_covariates(C_ref, p, dist_x, "reference", seed + 23L,
                               case5_df, case6_logvar)
  X_test <- generate_covariates(C_test, p, dist_x, "test", seed + 37L,
                                case5_df, case6_logvar)

  resp <- generate_responses(
    X_train, X_test, C_train, C_test, beta,
    var_a = var_a, var_b = var_b, var_e = var_e,
    seed = seed + 101L,
    intercept_distribution = intercept_distribution
  )

  model_type <- if (var_b > 0) "RS" else "RI"

  t0 <- proc.time()[3]
  search <- cirg_search(
    X_train, resp$y_train, X_ref,
    initial_Cn = initial_Cn, lambda = lambda, tau = tau,
    model_type = model_type,
    T0 = T0, alpha = alpha, max_iter = sa_max_iter,
    seed = seed + 1001L,
    em_tol = em_tol, em_max_iter = em_max_iter,
    verbose = verbose
  )
  regroup_runtime <- proc.time()[3] - t0
  best <- search$best

  # Same selected RASC grouping for both regrouped models.
  # The model matching the DGP reuses the fit already computed by its IMSPE.
  if (model_type == "RI") {
    fit_gall <- best$imspe_fit
    fit_gallrs <- fit_rs_lmm_cpp(best$cluster$X_sorted, best$cluster$y_sorted,
                                 as.integer(best$cluster$sizes), em_tol, em_max_iter)
  } else {
    fit_gall <- fit_ri_lmm_cpp(best$cluster$X_sorted, best$cluster$y_sorted,
                               as.integer(best$cluster$sizes), em_tol, em_max_iter)
    fit_gallrs <- best$imspe_fit
  }

  pred_gall <- predict_ri_soft(fit_gall, X_test, best$cluster$params)
  pred_gallrs <- predict_rs_soft(fit_gallrs, X_test, best$cluster$params)

  # Oracle benchmark: true train/test group identities share the same random
  # effects, and ALL uses the correctly specified RI/RS model.
  if (var_b > 0) {
    fit_all <- fit_rs_lmm_cpp(X_train, resp$y_train, as.integer(C_train),
                              em_tol, em_max_iter)
    pred_all <- predict_rs_known_groups(fit_all, X_test, resp$true_group_test)
  } else {
    fit_all <- fit_ri_lmm_cpp(X_train, resp$y_train, as.integer(C_train),
                              em_tol, em_max_iter)
    pred_all <- predict_ri_known_groups(fit_all, X_test, resp$true_group_test)
  }

  out <- rbind(
    method_metrics("ALL", fit_all, pred_all, resp$y_test, beta,
                   var_a, var_b, var_e,
                   if (var_b > 0) resp$G_true else NULL,
                   selected_K = R, runtime = NA_real_),
    method_metrics("GALL", fit_gall, pred_gall, resp$y_test, beta,
                   var_a, var_b, var_e,
                   selected_K = best$K, objective = best$objective,
                   runtime = regroup_runtime),
    method_metrics("GALLRS", fit_gallrs, pred_gallrs, resp$y_test, beta,
                   var_a, var_b, var_e, resp$G_true,
                   selected_K = best$K, objective = best$objective,
                   runtime = regroup_runtime)
  )
  out$IMSPE_selected <- c(NA_real_, best$IMSPE, best$IMSPE)
  out$lambda <- lambda
  out$tau <- tau
  out$dist_x <- dist_x
  out$var_a <- var_a
  out$var_b <- var_b
  out$var_e <- var_e
  out$R_true <- R
  out$N_test <- N
  out$objective_model <- model_type
  out$cluster_mode <- "rasc"
  list(metrics = out, search = search)
}

# ---------------------------- full experiment --------------------------------

run_cirg_simulation <- function(N_all = 2500,
                                p = 50,
                                R = 20,
                                Var.e = 9,
                                Var.a = 2.25,
                                Var.b = 0,
                                nloop = 50,
                                dist_x = "case1",
                                groupsize = "large",
                                lambda = 1,
                                tau = 5 * (p + 1L),
                                initial_Cn = 2,
                                T0 = 50,
                                alpha = 0.95,
                                sa_max_iter = 80,
                                em_tol = 1e-6,
                                em_max_iter = 100,
                                case5_df = 10,
                                case6_logvar = 0.2,
                                intercept_distribution = "normal",
                                seed = 12345,
                                keep_search = FALSE,
                                verbose = FALSE) {
  raw_all <- list()
  searches <- if (keep_search) list() else NULL
  pos <- 0L

  for (N in N_all) {
    for (k in seq_len(nloop)) {
      rep_seed <- as.integer(seed + 1000003L * k + 7919L * N)
      ans <- run_one_replication(
        N, p, R, dist_x, groupsize,
        Var.a, Var.b, Var.e,
        lambda, tau, rep_seed,
        initial_Cn, T0, alpha, sa_max_iter,
        em_tol, em_max_iter,
        case5_df, case6_logvar,
        intercept_distribution, verbose
      )
      pos <- pos + 1L
      d <- ans$metrics
      d$replication <- k
      raw_all[[pos]] <- d
      if (keep_search) searches[[pos]] <- ans$search
      if (verbose) cat("N=", N, " replication=", k, " done\n", sep = "")
    }
  }

  raw <- do.call(rbind, raw_all)
  summary <- do.call(rbind, lapply(split(raw, interaction(raw$N_test, raw$method, drop = TRUE)), function(d) {
    ss <- summarize_results(d)
    ss$N_test <- d$N_test[1]
    ss
  }))
  rownames(summary) <- NULL
  list(raw = raw, summary = summary, searches = searches)
}

# ---------------------- compatibility with legacy names ----------------------

legacy_recaps <- function(ans, N_all) {
  raw <- ans$raw
  byN <- split(raw, raw$N_test)
  getmat <- function(field) {
    do.call(rbind, lapply(as.character(N_all), function(nn) {
      d <- byN[[nn]]
      sapply(c("ALL", "GALL", "GALLRS"), function(m) mean(d[d$method == m, field], na.rm = TRUE))
    }))
  }
  rec1 <- getmat("beta_MSE")
  rec2 <- getmat("MSPE")
  rec3 <- getmat("intercept_MSE")
  rec4 <- getmat("Var_a_hat")
  rec5 <- getmat("Var_e_hat")
  colnames(rec1) <- colnames(rec2) <- colnames(rec3) <- colnames(rec4) <- colnames(rec5) <- c("ALL", "GALL", "GALLRS")
  list(rec1, rec2, rec3, rec4, rec5, raw = raw, summary = ans$summary)
}

Comp <- function(N_all, p, R, Var.e, nloop, n = NULL,
                 dist_x = "case1", dist_a = "N.ori", groupsize = "large",
                 obj.c = NULL, lambda = 1, tau = 5 * (p + 1L),
                 Var.a = NULL, ...) {
  var_a <- if (!is.null(Var.a)) Var.a else if (dist_a == "N.ML") 0 else if (dist_a == "T") 3 else 2.25
  int_dist <- if (dist_a == "T") "t3" else "normal"
  ans <- run_cirg_simulation(
    N_all = N_all, p = p, R = R, Var.e = Var.e,
    Var.a = var_a, Var.b = 0, nloop = nloop,
    dist_x = dist_x, groupsize = groupsize,
    lambda = lambda, tau = tau,
    intercept_distribution = int_dist, ...
  )
  legacy_recaps(ans, N_all)
}

Comp_RS <- function(N_all, p, R, Var.e, nloop, n = NULL,
                    dist_x = "case1", dist_a = "N.ori", groupsize = "large",
                    obj.c = NULL, lambda = 1, tau = 5 * (p + 1L),
                    Var.a = NULL, Var.b = 0.1, ...) {
  var_a <- if (!is.null(Var.a)) Var.a else if (dist_a == "N.ML") 0 else if (dist_a == "T") 3 else 2.25
  int_dist <- if (dist_a == "T") "t3" else "normal"
  ans <- run_cirg_simulation(
    N_all = N_all, p = p, R = R, Var.e = Var.e,
    Var.a = var_a, Var.b = Var.b, nloop = nloop,
    dist_x = dist_x, groupsize = groupsize,
    lambda = lambda, tau = tau,
    intercept_distribution = int_dist, ...
  )
  legacy_recaps(ans, N_all)
}


# -----------------------------------------------------------------------------
# HPC note: definitions only. No experiment is run automatically when sourced.
# Use run_RI_case1.R ... run_RS_case9.R to launch individual experiments.
# -----------------------------------------------------------------------------
