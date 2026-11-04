# Source after ROG_aligned_fixed.R. BLM is the Gaussian, time-invariant specialization.
atomic_save_rds <- function(object, path, compress = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".result-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp, compress = compress)
  if (!file.rename(tmp, path)) stop("Cannot replace result: ", path)
}

validate_literature_xy <- function(X, y) {
  if (!is.matrix(X) || !is.numeric(X) || nrow(X) != length(y) ||
      any(!is.finite(X)) || any(!is.finite(y)) || ncol(X) < 1L ||
      nrow(X) <= ncol(X) + 1L) stop("Supply finite X and y with n > p + 1.")
  if (qr(cbind(1, X))$rank < ncol(X) + 1L)
    stop("The common-intercept design must have full column rank.")
}

# BLM eq. (4) and eq. (6): first K with Q(K) <= gamma * V_h.
# V_h must describe measurement error in these exact moment coordinates.
fit_blm_external_moments <- function(moments, moment_variance = NULL,
                                     K = NULL, K_max = nrow(moments),
                                     gamma = 1, seed = 1L, nstart = 20L,
                                     iter.max = 200L) {
  M <- as.matrix(moments)
  if (!is.numeric(M) || !nrow(M) || !ncol(M) || any(!is.finite(M)))
    stop("BLM requires a finite individual-moment matrix.")
  if (!is.finite(gamma) || gamma <= 0 || gamma > 1)
    stop("BLM gamma must lie in (0,1].")
  if (is.null(K) && (is.null(moment_variance) || length(moment_variance) != 1L ||
                    !is.finite(moment_variance) || moment_variance < 0))
    stop("BLM automatic K requires the sampling-error variance V_h of informative moments.")
  n <- nrow(M)
  distinct <- nrow(unique(M))
  K_max <- min(as.integer(K_max), distinct)
  if (length(K_max) != 1L || !is.finite(K_max) || K_max < 1L)
    stop("K_max must be positive.")
  if (!is.null(K) && (length(K) != 1L || !is.finite(K) || K != as.integer(K) ||
                      K < 1L || K > K_max)) stop("Invalid fixed BLM K.")
  candidates <- if (is.null(K)) seq_len(K_max) else as.integer(K)
  trace <- list(); chosen <- NULL
  for (k in candidates) {
    if (k == 1L) {
      centers <- matrix(colMeans(M), nrow = 1L)
      labels <- rep(1L, n)
      wss <- sum(sweep(M, 2, centers[1, ], "-")^2)
      ok <- TRUE
    } else if (k == distinct) {
      centers <- unique(M)
      labels <- nearest_moment_centroid(M, centers)
      wss <- 0
      ok <- TRUE
    } else {
      set.seed(seed + 7919L * k)
      km <- stats::kmeans(M, k, nstart = nstart, iter.max = iter.max)
      centers <- km$centers; labels <- km$cluster; wss <- km$tot.withinss
      ok <- is.null(km$ifault) || km$ifault == 0L
    }
    criterion <- wss / n
    threshold <- if (is.null(moment_variance)) NA_real_ else gamma * moment_variance
    accepted <- ok && (!is.null(K) || criterion <= threshold)
    trace[[length(trace) + 1L]] <- data.frame(
      K = k, Q = criterion, threshold = threshold, converged = ok, accepted = accepted)
    # A failed smaller K prevents claiming the smallest admissible K.
    if (!ok) return(list(converged = FALSE, status = "kmeans_failed",
                         trace = do.call(rbind, trace)))
    if (accepted) {
      chosen <- list(group_map = as.integer(labels), centers = centers,
                     K = k, criterion = criterion, converged = TRUE,
                     status = "ok", moment_variance = moment_variance,
                     selection = if (is.null(K)) "BLM_eq6" else "fixed_K")
      break
    }
  }
  if (is.null(chosen)) chosen <- list(converged = FALSE, status = "threshold_not_reached")
  chosen$trace <- do.call(rbind, trace)
  chosen
}

nearest_moment_centroid <- function(M, centers) {
  M <- as.matrix(M)
  if (ncol(M) != ncol(centers) || any(!is.finite(M))) stop("Invalid new moments.")
  d <- outer(rowSums(M^2), rowSums(centers^2), "+") - 2 * tcrossprod(M, centers)
  max.col(-d, ties.method = "first")
}

assign_blm_external_moments <- function(moments_new, fit) {
  if (!isTRUE(fit$converged)) stop("Cannot assign with a failed BLM fit.")
  nearest_moment_centroid(moments_new, fit$centers)
}

blm_panel_moments <- function(moment_observations, unit_id,
                              independent_repeats = FALSE) {
  H <- as.matrix(moment_observations)
  if (is.null(unit_id) || length(unit_id) != nrow(H) || anyNA(unit_id) ||
      !is.numeric(H) || any(!is.finite(H)))
    stop("BLM panel moments require observed unit identities and finite measurements.")
  if (!isTRUE(independent_repeats))
    stop("Supply external V_h for dependent repeats, or explicitly assume independent repeats.")
  units <- unique(unit_id)
  row_unit <- match(unit_id, units)
  counts <- tabulate(row_unit, length(units))
  if (any(counts < 2L)) stop("Repeated observations per unit are needed to estimate V_h.")
  if (length(unique(counts)) != 1L)
    stop("This paper-aligned panel interface requires a balanced panel.")
  moments <- rowsum(H, row_unit, reorder = TRUE) / counts
  residuals <- H - moments[row_unit, , drop = FALSE]
  # BLM footnote 7, the stated estimator (no finite-T correction).
  Vh <- sum(residuals^2) / (length(units) * counts[1]^2)
  list(moments = moments, moment_variance = Vh, unit_id = units,
       row_unit = row_unit, T = counts[1])
}

fit_grouped_gaussian <- function(X, y, labels) {
  validate_literature_xy(X, y)
  labels <- match(labels, sort(unique(labels)))
  K <- max(labels); p <- ncol(X)
  counts <- tabulate(labels, K)
  xbar <- rowsum(X, labels, reorder = TRUE) / counts
  ybar <- drop(rowsum(matrix(y), labels, reorder = TRUE)) / counts
  xc <- X - xbar[labels, , drop = FALSE]
  yc <- y - ybar[labels]
  fit <- lm.fit(xc, yc)
  if (fit$rank != p || nrow(X) <= K + p)
    stop("BLM grouped regression is unidentified or has no residual degrees of freedom.")
  beta <- fit$coefficients
  alpha <- ybar - drop(xbar %*% beta)
  residual <- y - alpha[labels] - drop(X %*% beta)
  list(beta = beta, alpha = alpha, rss = sum(residual^2),
       residual_variance = mean(residual^2), K = K, converged = TRUE)
}

# Gaussian time-invariant version of BLM eq. (5); caller supplies individual
# moments and a row-to-unit mapping. No proxy identities are manufactured.
fit_blm_original <- function(X, y, moments, row_unit = seq_len(nrow(X)),
                             moment_variance = NULL, ...) {
  if (length(row_unit) != nrow(X) || anyNA(row_unit) ||
      any(row_unit != as.integer(row_unit)) || any(row_unit < 1L) ||
      any(row_unit > nrow(moments)) || length(unique(row_unit)) != nrow(moments))
    stop("row_unit must map the training rows onto the supplied individual moments.")
  classification <- fit_blm_external_moments(moments, moment_variance, ...)
  if (!classification$converged) return(classification)
  labels <- classification$group_map[row_unit]
  estimate <- fit_grouped_gaussian(X, y, labels)
  estimate$classification <- classification
  estimate$labels <- labels
  estimate$status <- "ok"
  estimate
}

predict_blm_original <- function(fit, Xnew, moments_new = NULL,
                                 known_unit = NULL) {
  if (!isTRUE(fit$converged)) stop("Cannot predict with a failed BLM fit.")
  if (is.null(known_unit)) {
    if (is.null(moments_new) || nrow(moments_new) != nrow(Xnew))
      stop("New units require independent informative moments, not their held-out y.")
    # Nearest-centroid extension for new units, explicitly outside eq. (5).
    labels <- assign_blm_external_moments(moments_new, fit$classification)
  } else {
    if (length(known_unit) != nrow(Xnew) || anyNA(known_unit) ||
        any(known_unit != as.integer(known_unit)) || any(known_unit < 1L) ||
        any(known_unit > length(fit$classification$group_map)))
      stop("Invalid known-unit indices.")
    labels <- fit$classification$group_map[known_unit]
  }
  drop(Xnew %*% fit$beta) + fit$alpha[labels]
}
