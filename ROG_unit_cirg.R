# ROG_unit_cirg.R
# -----------------------------------------------------------------------------
# Unit-level CIRG (CIRG-U)
#
# The basic object being regrouped is an ORIGINAL REPEATED-MEASURE UNIT, not an
# individual observation.  All observations from one unit are therefore forced
# to remain together throughout the search.  This aligns the regrouping layer
# with the level at which LMM random effects are defined.
#
# Candidate-generation layer:
#   RI: unit summaries built from X means + pooled residual mean.
#   RS: unit summaries built from X means + ridge-stabilized unit effect proxies.
#
# Selection layer:
#   candidate unit partitions are evaluated by the same model-matched empirical
#   IMSPE used by CIRG.  K-means is only an initializer; simulated annealing can
#   move/swap/merge/split whole units and directly improve the partition.
#
# No test y values are used.  Reference observations enter only through target
# moments and are mapped to candidate super-groups by their known unit IDs.
# -----------------------------------------------------------------------------

# ----------------------------- basic helpers ---------------------------------

uc_compact_units <- function(unit_train, unit_ref = NULL, unit_test = NULL) {
  tr_raw <- as.character(unit_train)
  lev <- unique(tr_raw)
  tr <- match(tr_raw, lev)
  out <- list(train = as.integer(tr), levels = lev, R = length(lev))
  map_one <- function(x, nm) {
    if (is.null(x)) return(NULL)
    z <- match(as.character(x), lev)
    if (anyNA(z)) stop(nm, " contains unit IDs absent from training data.")
    as.integer(z)
  }
  out$ref <- map_one(unit_ref, "unit_ref")
  out$test <- map_one(unit_test, "unit_test")
  out
}

uc_normalize_partition <- function(group_map) {
  group_map <- as.integer(group_map)
  if (!length(group_map) || anyNA(group_map)) stop("Invalid unit partition.")
  as.integer(match(group_map, sort(unique(group_map))))
}

uc_standardize <- function(M) {
  M <- as.matrix(M)
  if (!ncol(M)) return(M)
  s <- apply(M, 2, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  if (!any(keep)) return(matrix(0, nrow(M), 1L))
  M <- M[, keep, drop = FALSE]
  mu <- colMeans(M)
  s <- apply(M, 2, stats::sd)
  sweep(sweep(M, 2, mu, "-"), 2, s, "/")
}

uc_pc_scores <- function(M, ncomp = 5L) {
  M <- as.matrix(M)
  if (nrow(M) <= 1L || ncol(M) == 0L || ncomp <= 0L) return(matrix(numeric(0), nrow(M), 0L))
  s <- apply(M, 2, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  if (!any(keep)) return(matrix(numeric(0), nrow(M), 0L))
  Z <- scale(M[, keep, drop = FALSE], center = TRUE, scale = TRUE)
  npc <- min(as.integer(ncomp), ncol(Z), nrow(Z) - 1L)
  if (npc <= 0L) return(matrix(numeric(0), nrow(M), 0L))
  pc <- try(stats::prcomp(Z, center = FALSE, scale. = FALSE, rank. = npc), silent = TRUE)
  if (inherits(pc, "try-error")) return(matrix(numeric(0), nrow(M), 0L))
  as.matrix(pc$x[, seq_len(npc), drop = FALSE])
}

uc_pooled_fit <- function(X, y) {
  F <- cbind(1, X)
  b <- try(qr.solve(F, y), silent = TRUE)
  if (inherits(b, "try-error")) {
    XtX <- crossprod(F) + diag(1e-8, ncol(F))
    b <- solve(XtX, crossprod(F, y))
  }
  b <- as.numeric(b)
  list(beta = b, residuals = as.numeric(y - F %*% b))
}

uc_ridge_unit_effects <- function(X, y, units, beta, ridge = 1) {
  units <- as.integer(units)
  R <- max(units); p <- ncol(X)
  F <- cbind(1, X)
  residual <- as.numeric(y - F %*% beta)
  U <- matrix(0, R, p + 1L)
  for (u in seq_len(R)) {
    idx <- which(units == u)
    Fu <- F[idx, , drop = FALSE]
    ru <- residual[idx]
    pen <- diag(ridge, p + 1L)
    pen[1, 1] <- ridge * 0.25
    U[u, ] <- as.numeric(solve(crossprod(Fu) + pen, crossprod(Fu, ru)))
  }
  U
}

# -------------------------- unit-level embedding ------------------------------

build_unit_embedding <- function(X, y, unit_labels,
                                 model_type = c("RI", "RS"),
                                 x_pc = 5L,
                                 effect_pc = 5L,
                                 ridge = 1,
                                 effect_weight = 1,
                                 include_size = FALSE) {
  model_type <- match.arg(model_type)
  uc <- uc_compact_units(unit_labels)
  u <- uc$train; R <- uc$R; p <- ncol(X)
  pooled <- uc_pooled_fit(X, y)
  r <- pooled$residuals

  xmean <- matrix(0, R, p)
  rmean <- numeric(R)
  nunit <- integer(R)
  for (i in seq_len(R)) {
    idx <- which(u == i)
    nunit[i] <- length(idx)
    xmean[i, ] <- colMeans(X[idx, , drop = FALSE])
    rmean[i] <- mean(r[idx])
  }

  Xscore <- uc_pc_scores(xmean, x_pc)
  pieces <- list()

  if (model_type == "RI") {
    # rbar_i is a noisy proxy for a_i under a random-intercept LMM.
    pieces[[length(pieces) + 1L]] <- matrix(effect_weight * rmean, ncol = 1L,
                                             dimnames = list(NULL, "rmean"))
  } else {
    # For RS, use a ridge-stabilized unit-specific intercept/slope proxy.  Keep
    # the intercept explicitly and compress the p slope coordinates by PCA.
    Uproxy <- uc_ridge_unit_effects(X, y, u, pooled$beta, ridge)
    pieces[[length(pieces) + 1L]] <- matrix(effect_weight * Uproxy[, 1], ncol = 1L,
                                             dimnames = list(NULL, "u0_proxy"))
    slope_pc <- uc_pc_scores(Uproxy[, -1, drop = FALSE], effect_pc)
    if (ncol(slope_pc)) pieces[[length(pieces) + 1L]] <- effect_weight * slope_pc
  }

  if (ncol(Xscore)) pieces[[length(pieces) + 1L]] <- Xscore
  if (include_size) pieces[[length(pieces) + 1L]] <- matrix(log1p(nunit), ncol = 1L)

  H <- do.call(cbind, pieces)
  H <- uc_standardize(H)
  rownames(H) <- uc$levels

  list(
    H = H,
    unit_levels = uc$levels,
    unit_labels = u,
    pooled_beta = pooled$beta,
    residual_mean = rmean,
    X_mean = xmean,
    n_unit = nunit,
    model_type = model_type
  )
}

# ----------------------- candidate partition helpers --------------------------

uc_partition_labels <- function(unit_labels, group_map) {
  group_map <- uc_normalize_partition(group_map)
  group_map[as.integer(unit_labels)]
}

uc_partition_valid <- function(group_map, unit_train, tau = 2L, min_units = 1L) {
  group_map <- uc_normalize_partition(group_map)
  K <- max(group_map)
  unit_counts <- tabulate(group_map, nbins = K)
  if (any(unit_counts < as.integer(min_units))) return(FALSE)
  row_groups <- group_map[as.integer(unit_train)]
  obs_counts <- tabulate(row_groups, nbins = K)
  if (any(obs_counts < as.integer(tau))) return(FALSE)
  TRUE
}

uc_kmeans_partition <- function(H, K, seed = 1L) {
  R <- nrow(H)
  K <- max(2L, min(as.integer(K), R))
  set.seed(seed)
  if (K == R) return(seq_len(R))
  km <- stats::kmeans(H, centers = K, nstart = min(30L, max(5L, R)), iter.max = 100)
  uc_normalize_partition(km$cluster)
}

uc_grouped_training <- function(X, y, unit_train, group_map) {
  labels <- uc_partition_labels(unit_train, group_map)
  ord <- order(labels)
  K <- max(labels)
  list(
    labels = labels,
    order = ord,
    X_sorted = X[ord, , drop = FALSE],
    y_sorted = y[ord],
    sizes = as.integer(tabulate(labels, nbins = K)),
    K = K
  )
}

unit_target_moments <- function(X_reference, unit_ref, group_map) {
  glab <- uc_partition_labels(unit_ref, group_map)
  K <- max(glab)
  F <- cbind(1, X_reference)
  nref <- nrow(F)
  W <- crossprod(F) / nref
  Wk <- vector("list", K)
  mk <- vector("list", K)
  pk <- numeric(K)
  for (k in seq_len(K)) {
    idx <- which(glab == k)
    pk[k] <- length(idx) / nref
    if (!length(idx)) {
      Wk[[k]] <- matrix(0, ncol(F), ncol(F))
      mk[[k]] <- rep(0, ncol(F))
    } else {
      Fk <- F[idx, , drop = FALSE]
      Wk[[k]] <- crossprod(Fk) / nref
      mk[[k]] <- colSums(Fk) / nref
    }
  }
  list(W = W, Wk = Wk, mk = mk, pk = pk, labels = glab)
}

# ---------------------------- IMSPE evaluator --------------------------------

evaluate_unit_cirg_partition <- function(group_map,
                                         X_train, y_train, unit_train,
                                         X_reference, unit_ref,
                                         model_type = c("RI", "RS"),
                                         tau = ncol(X_train) + 1L,
                                         min_units = 1L,
                                         em_tol = 1e-5,
                                         em_max_iter = 500) {
  model_type <- match.arg(model_type)
  group_map <- uc_normalize_partition(group_map)
  if (!uc_partition_valid(group_map, unit_train, tau, min_units)) return(NULL)

  gd <- uc_grouped_training(X_train, y_train, unit_train, group_map)
  tm <- unit_target_moments(X_reference, unit_ref, group_map)

  fit <- if (model_type == "RI") {
    imspe_ri_cpp(gd$X_sorted, gd$y_sorted, gd$sizes,
                 tm$W, tm$mk, tm$pk,
                 tol = em_tol, max_iter = em_max_iter)
  } else {
    imspe_rs_cpp(gd$X_sorted, gd$y_sorted, gd$sizes,
                 tm$W, tm$Wk,
                 tol = em_tol, max_iter = em_max_iter)
  }

  list(
    group_map = group_map,
    K = max(group_map),
    grouped = gd,
    moments = tm,
    objective = as.numeric(fit$objective),
    IMSPE = as.numeric(fit$IMSPE),
    imspe_fit = fit,
    model_type = model_type
  )
}

# -------------------------- whole-unit neighbors ------------------------------

uc_neighbor <- function(group_map, H, min_units = 1L, maxK = nrow(H), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g <- uc_normalize_partition(group_map)
  K <- max(g); R <- length(g)
  counts <- tabulate(g, nbins = K)

  ops <- character(0)
  if (K >= 2L && any(counts > min_units)) ops <- c(ops, "move")
  if (K >= 2L) ops <- c(ops, "swap")
  if (K > 2L) ops <- c(ops, "merge")
  if (K < maxK && any(counts >= 2L * min_units)) ops <- c(ops, "split")
  if (!length(ops)) return(list(map = g, op = "none"))

  # Favor local moves while retaining explicit K-changing proposals.
  w <- c(move = 0.45, swap = 0.20, merge = 0.15, split = 0.20)
  p <- w[ops] / sum(w[ops])
  op <- sample(ops, 1L, prob = p)
  out <- g

  if (op == "move") {
    source_groups <- which(counts > min_units)
    gs <- sample(source_groups, 1L)
    u <- sample(which(g == gs), 1L)
    gt <- sample(setdiff(seq_len(K), gs), 1L)
    out[u] <- gt
  } else if (op == "swap") {
    pairg <- sample(seq_len(K), 2L)
    u1 <- sample(which(g == pairg[1]), 1L)
    u2 <- sample(which(g == pairg[2]), 1L)
    out[c(u1, u2)] <- c(pairg[2], pairg[1])
  } else if (op == "merge") {
    pairg <- sample(seq_len(K), 2L)
    out[out == pairg[2]] <- pairg[1]
  } else if (op == "split") {
    splittable <- which(counts >= 2L * min_units)
    gs <- sample(splittable, 1L)
    idx <- which(g == gs)
    if (length(idx) < 2L) return(list(map = g, op = "none"))
    set.seed(sample.int(.Machine$integer.max, 1L))
    km <- try(stats::kmeans(H[idx, , drop = FALSE], centers = 2L,
                            nstart = min(20L, length(idx)), iter.max = 100), silent = TRUE)
    if (inherits(km, "try-error")) return(list(map = g, op = "none"))
    if (any(tabulate(km$cluster, nbins = 2L) < min_units)) return(list(map = g, op = "none"))
    newg <- K + 1L
    out[idx[km$cluster == 2L]] <- newg
  }

  list(map = uc_normalize_partition(out), op = op)
}

# ---------------------------- CIRG-U search ----------------------------------

unit_cirg_search <- function(X_train, y_train, unit_train,
                             X_reference, unit_ref,
                             model_type = c("RI", "RS"),
                             K_grid = NULL,
                             initial_Cn = 2L,
                             tau = ncol(X_train) + 1L,
                             min_units = 2L,
                             x_pc = 5L,
                             effect_pc = 5L,
                             ridge = 1,
                             effect_weight = 1,
                             max_iter = 80L,
                             alpha = 0.95,
                             T0 = NULL,
                             seed = 1L,
                             em_tol = 1e-5,
                             em_max_iter = 500,
                             verbose = FALSE) {
  model_type <- match.arg(model_type)
  units <- uc_compact_units(unit_train, unit_ref)
  u_train <- units$train; u_ref <- units$ref; R <- units$R
  min_units <- max(1L, as.integer(min_units))
  tau <- max(2L, as.integer(tau))

  emb <- build_unit_embedding(X_train, y_train, u_train, model_type,
                              x_pc, effect_pc, ridge, effect_weight)
  H <- emb$H

  max_by_units <- floor(R / min_units)
  if (max_by_units < 2L) stop("Too few units for K>=2 under min_units.")
  # Observation-count constraints can make this smaller, but are checked exactly
  # in the evaluator because unit sizes need not be balanced.
  maxK <- min(R, max_by_units)
  if (is.null(K_grid)) K_grid <- seq.int(2L, maxK)
  K_grid <- sort(unique(as.integer(K_grid)))
  K_grid <- K_grid[K_grid >= 2L & K_grid <= maxK]
  if (!length(K_grid)) stop("No admissible K for unit-level CIRG.")

  cache <- new.env(parent = emptyenv())
  eval_cached <- function(map) {
    map <- uc_normalize_partition(map)
    key <- paste(map, collapse = ",")
    if (exists(key, envir = cache, inherits = FALSE)) return(get(key, envir = cache))
    ans <- evaluate_unit_cirg_partition(map, X_train, y_train, u_train,
                                        X_reference, u_ref, model_type,
                                        tau, min_units, em_tol, em_max_iter)
    assign(key, ans, envir = cache)
    ans
  }

  init_rows <- list(); init_candidates <- list()
  for (K in K_grid) {
    mp <- try(uc_kmeans_partition(H, K, seed + 1009L * K), silent = TRUE)
    if (inherits(mp, "try-error")) next
    ev <- eval_cached(mp)
    if (is.null(ev)) next
    init_candidates[[length(init_candidates) + 1L]] <- ev
    init_rows[[length(init_rows) + 1L]] <- data.frame(
      stage = "init", iter = 0L, op = "kmeans_init", K = ev$K,
      objective = ev$objective, IMSPE = ev$IMSPE,
      accepted = TRUE, temperature = NA_real_
    )
  }
  if (!length(init_candidates)) stop("All unit-level CIRG initial partitions failed constraints.")
  obj <- vapply(init_candidates, function(z) z$objective, numeric(1))
  current <- init_candidates[[which.max(obj)]]
  best <- current

  if (is.null(T0)) {
    s <- stats::sd(obj)
    if (!is.finite(s) || s < 1e-4) s <- 0.10
    Tcur <- max(0.05, s)
  } else {
    Tcur <- max(1e-6, as.numeric(T0))
  }

  trace <- init_rows
  set.seed(seed + 700001L)
  for (iter in seq_len(as.integer(max_iter))) {
    nb <- uc_neighbor(current$group_map, H, min_units, maxK,
                      seed = seed + 5003L * iter)
    if (identical(nb$op, "none")) next
    cand <- eval_cached(nb$map)
    if (is.null(cand)) {
      trace[[length(trace) + 1L]] <- data.frame(
        stage = "sa", iter = iter, op = nb$op, K = max(nb$map),
        objective = NA_real_, IMSPE = NA_real_, accepted = FALSE,
        temperature = Tcur
      )
      Tcur <- Tcur * alpha
      next
    }

    delta <- cand$objective - current$objective
    accp <- if (delta >= 0) 1 else exp(delta / Tcur)
    accepted <- stats::runif(1) < accp
    if (accepted) current <- cand
    if (cand$objective > best$objective) best <- cand

    trace[[length(trace) + 1L]] <- data.frame(
      stage = "sa", iter = iter, op = nb$op, K = cand$K,
      objective = cand$objective, IMSPE = cand$IMSPE,
      accepted = accepted, temperature = Tcur
    )
    if (verbose) {
      cat(sprintf("CIRG-U iter=%d op=%s candK=%d IMSPE=%.6f best=%.6f accepted=%s\n",
                  iter, nb$op, cand$K, cand$IMSPE, best$IMSPE, accepted))
    }
    Tcur <- Tcur * alpha
  }

  list(
    best = best,
    current = current,
    trace = do.call(rbind, trace),
    embedding = emb,
    K_grid = K_grid,
    min_units = min_units,
    tau = tau,
    model_type = model_type,
    cluster_mode = "unit_level_hybrid",
    n_evaluated = length(ls(cache))
  )
}

predict_unit_cirg <- function(search, X_test, unit_test,
                              model_type = c("RI", "RS")) {
  model_type <- match.arg(model_type)
  unit_test <- as.integer(unit_test)
  mp <- search$best$group_map
  if (any(unit_test < 1L | unit_test > length(mp))) stop("unit_test outside fitted unit map.")
  g <- mp[unit_test]
  fit <- search$best$imspe_fit
  pred <- if (model_type == "RI") {
    predict_ri_known_groups(fit, X_test, g)
  } else {
    predict_rs_known_groups(fit, X_test, g)
  }
  list(pred = pred, labels = g, K = max(g), fit = fit)
}
