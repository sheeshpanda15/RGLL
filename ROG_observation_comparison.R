# ROG_observation_comparison.R
# -----------------------------------------------------------------------------
# Observation-level method-comparison framework.
# Every clustering/regrouping method produces labels for INDIVIDUAL OBSERVATIONS.
# No method is allowed to use an original repeated-measure/unit ID as its
# indivisible clustering object.
#
# Main methods:
#   ORACLE : true observation labels (simulation upper benchmark only)
#   OBS    : misspecified observed observation labels
#   LM     : no grouping
#   KM     : observation-level K-means on X
#   GMM    : observation-level Gaussian mixture on X
#   CPF    : observation-level Ma-Huang-style concave pairwise fusion.
#            Each training observation has its own latent intercept. For large N
#            a sparse residual-neighbour graph approximates the all-pairs fusion
#            graph; for small N the full pairwise graph is used. Test labels are
#            transported using an X-only classifier fitted to the training labels.
#   MIXREG : observation-level latent-intercept finite mixture regression.
#            Train labels use (X,y); test labels are transported by an X-only
#            classifier. This replaces BLM, whose original object is a panel
#            individual rather than an observation.
#   CIRG   : proposed observation-level RASC + SGA, model-matched IMSPE regrouping.
#
# BLM is deliberately NOT a main baseline here: its defining first-step moments
# are individual/panel-level informative moments. Turning a single row into that
# object would no longer be the Bonhomme-Lamadon-Manresa method.
# -----------------------------------------------------------------------------

if (!exists("run_comparison_replication", mode = "function")) {
  if (!file.exists("ROG_method_comparison.R")) stop("ROG_method_comparison.R not found.")
  source("ROG_method_comparison.R")
}

# ----------------------- generic X -> label transport -------------------------

fit_x_label_classifier <- function(X, labels, seed = 1L, pca_dim = 10L,
                                   multinom_max_k = 20L, maxit = 250L) {
  labels <- as.integer(factor(labels))
  K <- max(labels)
  if (K <= 1L) return(list(type = "constant", K = 1L, label = 1L))

  X <- as.matrix(X)
  cen <- colMeans(X)
  sc <- apply(X, 2, stats::sd)
  sc[!is.finite(sc) | sc < 1e-10] <- 1
  Z0 <- sweep(sweep(X, 2, cen, "-"), 2, sc, "/")
  q <- min(as.integer(pca_dim), ncol(Z0), max(1L, nrow(Z0) - 1L))
  pc <- try(stats::prcomp(Z0, center = FALSE, scale. = FALSE, rank. = q), silent = TRUE)
  if (inherits(pc, "try-error")) {
    rot <- diag(ncol(Z0))[, seq_len(min(q, ncol(Z0))), drop = FALSE]
    Z <- Z0 %*% rot
  } else {
    rot <- pc$rotation[, seq_len(min(q, ncol(pc$rotation))), drop = FALSE]
    Z <- Z0 %*% rot
  }

  # Multinomial transport is preferred when it is numerically feasible.
  if (K <= multinom_max_k && requireNamespace("nnet", quietly = TRUE)) {
    set.seed(seed)
    dd <- as.data.frame(Z)
    dd$.g <- factor(labels, levels = seq_len(K))
    fm <- stats::as.formula(paste(".g ~", paste(names(dd)[names(dd) != ".g"], collapse = "+")))
    fit <- try(nnet::multinom(fm, data = dd, trace = FALSE, maxit = maxit,
                              MaxNWts = max(10000L, (ncol(Z) + 2L) * K * 10L)), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      return(list(type = "multinom", model = fit, center = cen, scale = sc,
                  rotation = rot, K = K))
    }
  }

  # Robust fallback: nearest centroid in a low-dimensional X representation.
  cents <- matrix(0, K, ncol(Z))
  for (k in seq_len(K)) cents[k, ] <- colMeans(Z[labels == k, , drop = FALSE])
  list(type = "centroid", centroids = cents, center = cen, scale = sc,
       rotation = rot, K = K)
}

predict_x_label_classifier <- function(object, Xnew) {
  if (object$type == "constant") return(rep(object$label, nrow(Xnew)))
  Xnew <- as.matrix(Xnew)
  Z0 <- sweep(sweep(Xnew, 2, object$center, "-"), 2, object$scale, "/")
  Z <- Z0 %*% object$rotation
  if (object$type == "multinom") {
    pr <- predict(object$model, newdata = as.data.frame(Z), type = "class")
    return(as.integer(pr))
  }
  C <- object$centroids
  d2 <- sapply(seq_len(nrow(C)), function(k) rowSums((Z - rep(C[k, ], each = nrow(Z)))^2))
  if (is.vector(d2)) return(rep(1L, nrow(Z)))
  max.col(-d2, ties.method = "first")
}

# ---------------- observation-level CPF: sparse/full pairwise fusion ----------

cpf_obs_edges <- function(residual, full_pair_max_n = 500L, graph_k = 8L) {
  n <- length(residual)
  if (n <= full_pair_max_n) return(t(utils::combn(seq_len(n), 2L)))
  graph_k <- max(1L, min(as.integer(graph_k), n - 1L))
  ord <- order(residual)
  E <- vector("list", graph_k)
  for (d in seq_len(graph_k)) {
    E[[d]] <- cbind(ord[seq_len(n - d)], ord[(1L + d):n])
  }
  ed <- do.call(rbind, E)
  ed <- t(apply(ed, 1, sort))
  unique(ed)
}

cpf_obs_prepare <- function(X, y, theta = 1, full_pair_max_n = 500L,
                            graph_k = 8L, ridge_beta = 1e-6) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("CPF requires package 'Matrix'.")
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  b0 <- tryCatch(qr.solve(cbind(1, X), y),
                 error = function(e) drop(MASS::ginv(crossprod(cbind(1, X))) %*% crossprod(cbind(1, X), y)))
  residual <- drop(y - cbind(1, X) %*% b0)
  edges <- cpf_obs_edges(residual, full_pair_max_n, graph_k)
  E <- nrow(edges)
  D <- Matrix::sparseMatrix(i = c(seq_len(E), seq_len(E)),
                            j = c(edges[, 1], edges[, 2]),
                            x = c(rep(1, E), rep(-1, E)),
                            dims = c(E, n))
  A <- Matrix::Diagonal(n) + theta * Matrix::crossprod(D)
  fac <- Matrix::Cholesky(A, LDL = FALSE, perm = TRUE)
  solveA <- function(B) as.matrix(Matrix::solve(fac, B))
  AX <- solveA(X)
  S <- crossprod(X) - crossprod(X, AX) + diag(ridge_beta, p)
  Sinv_rhs <- function(b) tryCatch(solve(S, b), error = function(e) drop(MASS::ginv(S) %*% b))
  list(X = X, y = y, n = n, p = p, edges = edges, D = D, theta = theta,
       solveA = solveA, AX = AX, Sinv_rhs = Sinv_rhs,
       Xty = crossprod(X, y), residual0 = residual)
}

cpf_obs_admm_once <- function(prep, lambda, gamma = 3,
                              tol = 1e-3, max_iter = 400L) {
  X <- prep$X; y <- prep$y; D <- prep$D; theta <- prep$theta
  n <- prep$n; edges <- prep$edges
  beta <- tryCatch(qr.solve(X, y), error = function(e) rep(0, ncol(X)))
  alpha <- drop(y - X %*% beta)
  eta <- as.numeric(D %*% alpha)
  v <- numeric(length(eta))
  converged <- FALSE; used <- 0L

  for (it in seq_len(max_iter)) {
    cvec <- as.numeric(Matrix::crossprod(D, theta * eta - v))
    qv <- as.numeric(prep$solveA(y + cvec))
    beta <- drop(prep$Sinv_rhs(prep$Xty - crossprod(X, qv)))
    alpha <- as.numeric(qv - drop(prep$AX %*% beta))
    d <- as.numeric(D %*% alpha)
    eta_old <- eta
    eta <- prox_mcp(d + v / theta, lambda, theta, gamma)
    r <- d - eta
    v <- v + theta * r
    primal <- sqrt(mean(r^2))
    dual <- sqrt(mean((eta - eta_old)^2)) * theta
    used <- it
    if (max(primal, dual) < tol) { converged <- TRUE; break }
  }

  fused <- abs(eta) < max(1e-7, tol * 0.25)
  labels <- union_find_components(n, edges, fused)
  list(beta = beta, alpha = alpha, eta = eta, labels = labels,
       K = max(labels), converged = converged, iterations = used)
}

merge_small_groups_by_alpha <- function(labels, alpha, min_size = 2L) {
  labels <- as.integer(factor(labels)); min_size <- max(1L, as.integer(min_size))
  repeat {
    K <- max(labels); sz <- tabulate(labels, K)
    bad <- which(sz < min_size)
    if (!length(bad) || K <= 1L) break
    g <- bad[which.min(sz[bad])]
    mus <- vapply(seq_len(K), function(k) mean(alpha[labels == k]), numeric(1))
    cand <- setdiff(seq_len(K), g)
    target <- cand[which.min(abs(mus[cand] - mus[g]))]
    labels[labels == g] <- target
    labels <- as.integer(factor(labels))
  }
  labels
}

rss_group_intercepts <- function(X, y, labels) {
  labels <- as.integer(factor(labels)); K <- max(labels)
  Xbar <- matrix(0, K, ncol(X)); ybar <- numeric(K)
  for (k in seq_len(K)) {
    ii <- labels == k
    Xbar[k, ] <- colMeans(X[ii, , drop = FALSE]); ybar[k] <- mean(y[ii])
  }
  Xc <- X - Xbar[labels, , drop = FALSE]
  yc <- y - ybar[labels]
  A <- crossprod(Xc) + diag(1e-8, ncol(X))
  beta <- tryCatch(solve(A, crossprod(Xc, yc)), error = function(e) drop(MASS::ginv(A) %*% crossprod(Xc, yc)))
  r <- drop(y - X %*% beta)
  mu <- vapply(seq_len(K), function(k) mean(r[labels == k]), numeric(1))
  fitted <- drop(X %*% beta + mu[labels])
  list(rss = sum((y - fitted)^2), beta = beta, mu = mu)
}

fit_cpf_observation <- function(X, y, min_group_size = 2L,
                                lambda_grid = NULL, theta = 1, gamma = 3,
                                bic_c = 10, tol = 1e-3, max_iter = 400L,
                                full_pair_max_n = 500L, graph_k = 8L) {
  prep <- cpf_obs_prepare(X, y, theta, full_pair_max_n, graph_k)
  scale0 <- stats::sd(prep$residual0)
  if (!is.finite(scale0) || scale0 < 1e-8) scale0 <- 1
  if (is.null(lambda_grid)) lambda_grid <- seq(0.05, 1.10, length.out = 10L) * scale0
  n <- nrow(X); p <- ncol(X)
  Cn <- bic_c * log(log(max(n + p, 3)))
  best <- NULL; trace <- list()
  for (lam in lambda_grid) {
    z <- cpf_obs_admm_once(prep, lam, gamma, tol, max_iter)
    lab <- merge_small_groups_by_alpha(z$labels, z$alpha, min_group_size)
    ss <- rss_group_intercepts(X, y, lab)
    K <- max(lab)
    bic <- log(max(ss$rss / n, 1e-12)) + Cn * log(n) / n * (K + p)
    trace[[length(trace) + 1L]] <- data.frame(lambda = lam, K = K, BIC = bic,
                                               converged = z$converged,
                                               iterations = z$iterations)
    if (is.null(best) || bic < best$BIC) {
      best <- list(labels = lab, K = K, BIC = bic, lambda = lam,
                   alpha = z$alpha, beta = ss$beta, mu = ss$mu,
                   converged = z$converged, iterations = z$iterations,
                   graph_edges = nrow(prep$edges),
                   graph_type = if (n <= full_pair_max_n) "full" else "sparse-residual-neighbour")
    }
  }
  best$trace <- do.call(rbind, trace)
  best
}

# ---------------- observation-level latent-intercept mixture regression -------

logsumexp_rows <- function(M) {
  mx <- apply(M, 1, max)
  mx + log(rowSums(exp(M - mx)))
}

fit_mixreg_one <- function(X, y, K, seed = 1L, max_iter = 250L,
                           tol = 1e-6, ridge = 1e-6) {
  X <- as.matrix(X); y <- as.numeric(y); n <- nrow(X); p <- ncol(X)
  set.seed(seed)
  b0 <- tryCatch(qr.solve(cbind(1, X), y), error = function(e) rep(0, p + 1L))
  beta <- b0[-1L]
  r0 <- drop(y - X %*% beta)
  init <- stats::kmeans(r0, centers = K, nstart = 5, iter.max = 50)$cluster
  mu <- vapply(seq_len(K), function(k) mean(r0[init == k]), numeric(1))
  mixpi <- pmax(tabulate(init, K) / n, 1e-6); mixpi <- mixpi / sum(mixpi)
  sig2 <- max(stats::var(r0), 1e-4)
  ll_old <- -Inf; converged <- FALSE; used <- 0L

  for (it in seq_len(max_iter)) {
    r <- drop(y - X %*% beta)
    logw <- sapply(seq_len(K), function(k) log(mixpi[k]) - 0.5 * (log(2 * base::pi * sig2) + (r - mu[k])^2 / sig2))
    den <- logsumexp_rows(logw)
    tau <- exp(logw - den)
    Nk <- colSums(tau)
    mixpi <- pmax(Nk / n, 1e-8); mixpi <- mixpi / sum(mixpi)
    for (mm in 1:2) {
      mu_exp <- drop(tau %*% mu)
      A <- crossprod(X) + diag(ridge, p)
      beta <- tryCatch(solve(A, crossprod(X, y - mu_exp)), error = function(e) drop(MASS::ginv(A) %*% crossprod(X, y - mu_exp)))
      r <- drop(y - X %*% beta)
      mu <- colSums(tau * r) / pmax(Nk, 1e-8)
    }
    err <- sapply(seq_len(K), function(k) (r - mu[k])^2)
    sig2 <- max(sum(tau * err) / n, 1e-6)
    logw <- sapply(seq_len(K), function(k) log(mixpi[k]) - 0.5 * (log(2 * base::pi * sig2) + (r - mu[k])^2 / sig2))
    ll <- sum(logsumexp_rows(logw))
    used <- it
    if (is.finite(ll_old) && abs(ll - ll_old) / (1 + abs(ll_old)) < tol) { ll_old <- ll; converged <- TRUE; break }
    ll_old <- ll
  }
  r <- drop(y - X %*% beta)
  logw <- sapply(seq_len(K), function(k) log(mixpi[k]) - 0.5 * (log(2 * base::pi * sig2) + (r - mu[k])^2 / sig2))
  den <- logsumexp_rows(logw)
  tau <- exp(logw - den)
  labels <- max.col(tau, ties.method = "first")
  list(beta = beta, mu = mu, pi = mixpi, sigma2 = sig2, labels = labels,
       logLik = ll_old, converged = converged, iterations = used)
}

fit_mixreg_bic <- function(X, y, K_grid = 2:10, min_group_size = 2L,
                           seed = 1L, max_iter = 250L) {
  n <- nrow(X); p <- ncol(X)
  K_grid <- sort(unique(as.integer(K_grid)))
  best <- NULL; trace <- list()
  for (K in K_grid) {
    z <- try(fit_mixreg_one(X, y, K, seed + 1009L * K, max_iter = max_iter), silent = TRUE)
    if (inherits(z, "try-error")) next
    sz <- tabulate(z$labels, K)
    if (any(sz < min_group_size)) next
    # shared slope beta, K latent intercepts, K-1 mixing weights, one variance
    df <- p + K + (K - 1L) + 1L
    bic <- -2 * z$logLik + df * log(n)
    trace[[length(trace) + 1L]] <- data.frame(K = K, BIC = bic,
                                               converged = z$converged,
                                               iterations = z$iterations)
    if (is.null(best) || bic < best$BIC) best <- c(z, list(K = K, BIC = bic))
  }
  if (is.null(best)) stop("All MIXREG candidates failed/minimum group size violated.")
  best$trace <- do.call(rbind, trace)
  best
}

# ----------------------- common observation-level evaluator ------------------

run_observation_replication <- function(N = 2500, p = 50, R = 20,
                                        dist_x = "case1", groupsize = "large",
                                        var_a = 2.25, var_b = 0, var_e = 9,
                                        mis_type = "contam", rho = 0.25,
                                        merge_factor = 2L,
                                        tau = 5 * (p + 1L),
                                        initial_Cn = 2L,
                                        sa_max_iter = 80,
                                        em_tol = 1e-5, em_max_iter = 500,
                                        seed = 12345L,
                                        methods = c("ORACLE", "OBS", "LM", "KM", "GMM", "CPF", "MIXREG", "CIRG"),
                                        gmm_models = c("EII", "VII", "EEI", "VEI"),
                                        cpf_graph_k = 8L,
                                        verbose = FALSE) {
  set.seed(seed)
  methods <- unique(toupper(methods))
  beta <- rep(1, p)
  model_type <- if (var_b > 0) "RS" else "RI"
  m <- N / (10 * R)
  C_test <- generate_groups(R, m, N, groupsize)
  C_train <- as.integer(3L * C_test); C_ref <- C_test
  X_train <- generate_covariates(C_train, p, dist_x, "train", seed + 11L)
  X_ref <- generate_covariates(C_ref, p, dist_x, "reference", seed + 23L)
  X_test <- generate_covariates(C_test, p, dist_x, "test", seed + 37L)
  resp <- generate_responses(X_train, X_test, C_train, C_test, beta,
                             var_a, var_b, var_e, seed + 101L, "normal")
  true_tr <- resp$true_group_train; true_te <- resp$true_group_test
  obs <- make_observed_groups(true_tr, true_te, X_train, X_test,
                              mis_type, rho, merge_factor, seed + 313L)
  G_true <- if (var_b > 0) resp$G_true else NULL
  rows <- list(); extra <- list(); add <- function(d) rows[[length(rows) + 1L]] <<- d
  min_group <- if (model_type == "RS") p + 1L else 2L

  if ("ORACLE" %in% methods) {
    t0 <- proc.time()[3]; fp <- fit_predict_labeled(X_train, resp$y_train, X_test, true_tr, true_te, model_type, em_tol, em_max_iter)
    add(metric_row("ORACLE", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, proc.time()[3] - t0, adjusted_rand_index(fp$labels_test, true_te)))
  }
  if ("OBS" %in% methods) {
    t0 <- proc.time()[3]; fp <- fit_predict_labeled(X_train, resp$y_train, X_test, obs$train, obs$test, model_type, em_tol, em_max_iter)
    add(metric_row("OBS", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, proc.time()[3] - t0, adjusted_rand_index(fp$labels_test, true_te)))
  }
  if ("LM" %in% methods) {
    t0 <- proc.time()[3]; fit <- fit_pooled_lm(X_train, resp$y_train); pred <- predict_pooled_lm(fit, X_test)
    add(metric_row("LM", fit, pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, 1L, proc.time()[3] - t0, NA_real_))
  }

  maxK <- max(2L, floor(nrow(X_train) / max(2L, tau)))
  K_grid <- seq.int(2L, maxK)

  if ("KM" %in% methods) {
    t0 <- proc.time()[3]; km <- fit_kmeans_bic(X_train, K_grid, seed + 401L, tau = tau)
    trlab <- km$cluster$labels; telab <- assign_kmeans(X_test, km$cluster$centroids)
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, trlab, telab, model_type, em_tol, em_max_iter)
    add(metric_row("KM", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, proc.time()[3] - t0, adjusted_rand_index(telab, true_te))); extra$KM <- km
  }
  if ("GMM" %in% methods && requireNamespace("mclust", quietly = TRUE)) {
    t0 <- proc.time()[3]; gm <- fit_gmm_bic(X_train, K_grid, seed + 503L, gmm_models)
    trlab <- gm$labels; telab <- predict_gmm_labels(gm, X_test)
    if (model_type == "RI" || min(tabulate(trlab, max(trlab))) >= p + 1L) {
      fp <- fit_predict_labeled(X_train, resp$y_train, X_test, trlab, telab, model_type, em_tol, em_max_iter)
      add(metric_row("GMM", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                     G_true, fp$K, proc.time()[3] - t0, adjusted_rand_index(telab, true_te)))
    }
    extra$GMM <- gm
  }
  if ("CPF" %in% methods) {
    t0 <- proc.time()[3]
    cpf <- fit_cpf_observation(X_train, resp$y_train, min_group_size = min_group,
                               graph_k = cpf_graph_k)
    clf <- fit_x_label_classifier(X_train, cpf$labels, seed + 613L)
    telab <- predict_x_label_classifier(clf, X_test)
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, cpf$labels, telab, model_type, em_tol, em_max_iter)
    add(metric_row("CPF", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, proc.time()[3] - t0, adjusted_rand_index(telab, true_te)))
    extra$CPF <- list(fit = cpf, classifier = clf)
  }
  if ("MIXREG" %in% methods) {
    t0 <- proc.time()[3]
    mixK <- K_grid[K_grid <= min(10L, max(K_grid))]
    mr <- fit_mixreg_bic(X_train, resp$y_train, mixK, min_group_size = min_group, seed = seed + 641L)
    clf <- fit_x_label_classifier(X_train, mr$labels, seed + 647L)
    telab <- predict_x_label_classifier(clf, X_test)
    fp <- fit_predict_labeled(X_train, resp$y_train, X_test, mr$labels, telab, model_type, em_tol, em_max_iter)
    add(metric_row("MIXREG", fp$fit, fp$pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, fp$K, proc.time()[3] - t0, adjusted_rand_index(telab, true_te)))
    extra$MIXREG <- list(fit = mr, classifier = clf)
  }
  if ("CIRG" %in% methods) {
    t0 <- proc.time()[3]
    search <- cirg_search(X_train, resp$y_train, X_ref,
                          initial_Cn = initial_Cn, tau = tau,
                          model_type = model_type, max_iter = sa_max_iter,
                          seed = seed + 701L, em_tol = em_tol,
                          em_max_iter = em_max_iter, verbose = verbose)
    best <- search$best; fit <- best$imspe_fit
    pred <- if (model_type == "RI") {
      predict_ri_soft(fit, X_test, best$cluster$params)
    } else {
      predict_rs_soft(fit, X_test, best$cluster$params)
    }
    telab <- predict_soft_labels(X_test, best$cluster$params)
    add(metric_row("CIRG", fit, pred, resp$y_test, beta, var_a, var_b, var_e,
                   G_true, best$K, proc.time()[3] - t0, adjusted_rand_index(telab, true_te), best$objective))
    extra$CIRG <- search
  }

  raw <- do.call(rbind, rows)
  raw$mis_type <- mis_type; raw$rho <- if (mis_type == "contam") rho else NA_real_
  raw$merge_factor <- merge_factor; raw$dist_x <- dist_x; raw$var_b <- var_b
  raw$R_true <- R; raw$N_test <- N
  list(metrics = raw, observed = obs, extra = extra)
}

run_observation_comparison <- function(N_all = 2500, p = 50, R = 20,
                                       Var.e = 9, Var.a = 2.25, Var.b = 0,
                                       nloop = 50, dist_x = "case1",
                                       groupsize = "large", mis_type = "contam",
                                       rho = 0.25, merge_factor = 2L,
                                       tau = 5 * (p + 1L), initial_Cn = 2L,
                                       sa_max_iter = 80, em_tol = 1e-5,
                                       em_max_iter = 500, seed = 12345L,
                                       methods = c("ORACLE", "OBS", "LM", "KM", "GMM", "CPF", "MIXREG", "CIRG"),
                                       keep_extra = FALSE, verbose = FALSE) {
  raw_all <- list(); extra_all <- if (keep_extra) list() else NULL; pos <- 0L
  for (N in N_all) for (r in seq_len(nloop)) {
    s <- as.integer(seed + 1000003L * r + 7919L * N)
    ans <- run_observation_replication(N, p, R, dist_x, groupsize,
                                       Var.a, Var.b, Var.e, mis_type, rho,
                                       merge_factor, tau, initial_Cn,
                                       sa_max_iter, em_tol, em_max_iter,
                                       s, methods, verbose = verbose)
    pos <- pos + 1L; d <- ans$metrics; d$replication <- r; raw_all[[pos]] <- d
    if (keep_extra) extra_all[[pos]] <- ans$extra
  }
  raw <- do.call(rbind, raw_all)
  list(raw = raw, summary = summarize_comparison(raw), extra = extra_all)
}

# No experiment runs automatically.
