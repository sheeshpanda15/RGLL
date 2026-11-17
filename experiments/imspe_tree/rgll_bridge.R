args <- commandArgs(trailingOnly = TRUE)
setwd(args[[1L]])
source("ROG_method_comparison.R")

generate_rgll_data <- function(old_case, R, var_a, var_b, seed) {
  # Mirror the data-only prefix of RGLL's run_comparison_replication.
  # Call its original generators so that its R RNG stream is preserved.
  N <- 2500L
  p <- 50L
  beta <- rep(1, p)
  set.seed(seed)
  C_test <- generate_groups(R, N / (10 * R), N, "large")
  C_train <- 3L * C_test
  if (old_case >= 10L) {
    generator <- get(sprintf("generate_case%d_comparison_data", old_case))
    return(generator(C_train, C_test, p, beta, var_a, var_b, 9, seed + 101L))
  }
  dist_x <- paste0("case", old_case)
  X_train <- generate_covariates(C_train, p, dist_x, "train", seed + 11L)
  X_ref <- generate_covariates(C_test, p, dist_x, "reference", seed + 23L)
  X_test <- generate_covariates(C_test, p, dist_x, "test", seed + 37L)
  resp <- generate_responses(X_train, X_test, C_train, C_test, beta,
                             var_a, var_b, 9, seed + 101L, "normal")
  list(X_train = X_train, X_ref = X_ref, X_test = X_test, resp = resp)
}

request <- function(parts) {
  mode <- parts[[1L]]
  if (mode == "quit") return(FALSE)
  old_case <- as.integer(parts[[2L]])
  R <- as.integer(parts[[3L]])
  var_a <- as.numeric(parts[[4L]])
  var_b <- as.numeric(parts[[5L]])
  rep_id <- as.integer(parts[[6L]])
  output <- parts[[7L]]
  seed <- as.integer(12345L + 1000003L * rep_id + 7919L * 2500L)
  if (mode == "data") {
    d <- generate_rgll_data(old_case, R, var_a, var_b, seed)
    values <- c(as.vector(t(d$X_train)), as.vector(t(d$X_ref)), as.vector(t(d$X_test)),
                d$resp$y_train, d$resp$y_test,
                d$resp$true_group_train, d$resp$true_group_test)
    writeBin(as.double(values), output, size = 8L, endian = "little")
  } else if (mode == "reference") {
    # Test fixture for one fixed two-leaf tree, matching Python's threshold.
    d <- generate_rgll_data(old_case, R, var_a, var_b, seed)
    feature <- as.integer(parts[[8L]]) + 1L
    threshold <- as.numeric(parts[[9L]])
    train_labels <- 1L + (d$X_train[, feature] > threshold)
    ref_labels <- 1L + (d$X_ref[, feature] > threshold)
    o <- order(train_labels)
    f <- cbind(1, d$X_ref)
    groups <- lapply(1:2, function(k) f[ref_labels == k, , drop = FALSE])
    W <- crossprod(f) / nrow(f)
    if (var_b == 0) {
      score <- imspe_ri_cpp(d$X_train[o, ], d$resp$y_train[o], tabulate(train_labels), W,
                            lapply(groups, function(g) colSums(g) / nrow(f)),
                            tabulate(ref_labels, 2L) / nrow(f), tol = 1e-5, max_iter = 300L)
    } else {
      score <- imspe_rs_cpp(d$X_train[o, ], d$resp$y_train[o], tabulate(train_labels), W,
                            lapply(groups, function(g) crossprod(g) / nrow(f)),
                            tol = 1e-5, max_iter = 300L)
    }
    values <- c(score$IMSPE, score$Var.a, if (var_b > 0) score$Var.b else 0,
                score$Var.e, score$iterations, score$converged, score$beta,
                if (var_b > 0) as.vector(t(score$u_hat)) else score$u_hat)
    writeBin(as.double(values), output, size = 8L, endian = "little")
  } else if (mode == "live") {
    methods <- strsplit(parts[[8L]], ",", fixed = TRUE)[[1L]]
    rows <- lapply(methods, function(method) {
      tryCatch({
        d <- run_comparison_replication(N = 2500L, p = 50L, R = R,
                 dist_x = paste0("case", old_case), var_a = var_a, var_b = var_b,
                 tau = 255L, label_policy = "unknown", seed = seed,
                 methods = method)$metrics
        d$error <- ""
        d[, c("method", "MSPE", "selected_K", "runtime", "group_ARI", "IMSPE_selected", "iterations", "converged", "error")]
      }, error = function(e) {
        data.frame(method = method, MSPE = NA_real_, selected_K = NA_integer_, runtime = NA_real_,
                   group_ARI = NA_real_, IMSPE_selected = NA_real_, iterations = NA_integer_,
                   converged = NA, error = conditionMessage(e))
      })
    })
    write.csv(do.call(rbind, rows), output, row.names = FALSE)
  } else {
    stop("Unknown request mode: ", mode)
  }
  TRUE
}

cat("READY\n")
flush.console()
input <- file("stdin", open = "r")
repeat {
  line <- readLines(input, n = 1L)
  if (!length(line)) break
  answer <- tryCatch(request(strsplit(line, "\t", fixed = TRUE)[[1L]]),
                     error = function(e) gsub("[\r\n]", " ", conditionMessage(e)))
  if (identical(answer, FALSE)) break
  cat(if (isTRUE(answer)) "OK" else paste0("ERROR ", answer), "\n", sep = "")
  flush.console()
}
