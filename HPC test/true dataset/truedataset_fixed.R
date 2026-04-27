#######################
rm(list=ls())
library(devtools)
library(ClusterR)
library(MASS)
Rcpp::sourceCpp("lmm_fast.cpp")

calculate_K_dynamic <- function(FXX_test) {
  X_test_with_intercept <- cbind(1, FXX_test)
  K <- (t(X_test_with_intercept) %*% X_test_with_intercept) / nrow(FXX_test)
  return(K)
}
MSPE_fn = function(fy_test, fx_test, sx_train, sy_train, beta_hat, Var.a, Var.e, nc_train, centroids){
  R <- length(nc_train)
  mv_hat <- numeric(R)
  index <- 1
  fixed_pred_train <- cbind(1, sx_train) %*% beta_hat
  for (i in 1:R) {
    current_indices <- index:(index + nc_train[i] - 1)
    residuals_i <- sy_train[current_indices] - fixed_pred_train[current_indices]
    mv_hat[i] <- (Var.a / (Var.e + nc_train[i] * Var.a)) * sum(residuals_i)
    index <- index + nc_train[i]
  }
  pred_labels <- ClusterR::predict_KMeans(fx_test, centroids)
  u_assigned <- mv_hat[pred_labels]
  y_hat <- cbind(1, fx_test) %*% beta_hat + u_assigned
  return(mean((fy_test - y_hat)^2))
}
MSPE_fn_RS = function(fy_test, fx_test, sx_train, sy_train, beta_hat, G_hat, Var.e, nc_train, centroids){
  R <- length(nc_train)
  q <- ncol(fx_test) + 1
  u_vecs <- matrix(0, nrow = R, ncol = q)
  index <- 1
  X_train_full <- cbind(1, sx_train)
  G_inv <- solve(G_hat + diag(1e-7, q))
  for (i in 1:R) {
    curr_idx <- index:(index + nc_train[i] - 1)
    Xg <- X_train_full[curr_idx, , drop = FALSE]
    yg <- sy_train[curr_idx]
    res_g <- yg - Xg %*% beta_hat
    M <- solve(Var.e * G_inv + t(Xg) %*% Xg)
    u_vecs[i, ] <- M %*% (t(Xg) %*% res_g)
    index <- index + nc_train[i]
  }
  pred_labels <- ClusterR::predict_KMeans(fx_test, centroids)
  u_assigned <- u_vecs[pred_labels, ]
  X_test_full <- cbind(1, fx_test)
  y_hat <- X_test_full %*% beta_hat + rowSums(X_test_full * u_assigned)
  return(mean((fy_test - y_hat)^2))
}
MSPE_tru_RS = function(fy_test, fx_test, sx_train, sy_train, beta_hat, G_hat, Var.e, nc_train, nc_test, R){
  q <- ncol(fx_test) + 1
  index_train <- 1
  index_test <- 1
  X_train_full <- cbind(1, sx_train)
  X_test_full <- cbind(1, fx_test)
  G_inv <- solve(G_hat + diag(1e-7, q))
  y_hat <- numeric(length(fy_test))
  for (i in 1:R) {
    curr_idx_train <- index_train:(index_train + nc_train[i] - 1)
    curr_idx_test  <- index_test:(index_test + nc_test[i] - 1)
    Xg_train <- X_train_full[curr_idx_train, , drop = FALSE]
    yg_train <- sy_train[curr_idx_train]
    res_g <- yg_train - Xg_train %*% beta_hat
    M <- solve(Var.e * G_inv + t(Xg_train) %*% Xg_train)
    u_g <- M %*% (t(Xg_train) %*% res_g)
    Xg_test <- X_test_full[curr_idx_test, , drop = FALSE]
    y_hat[curr_idx_test] <- Xg_test %*% beta_hat + Xg_test %*% u_g
    index_train <- index_train + nc_train[i]
    index_test  <- index_test  + nc_test[i]
  }
  return(mean((fy_test - y_hat)^2))
}
mbky <- function(setseed, FXX, y, Cn, threshold = 2) {
  current_Cn <- Cn
  repeat {
    if (current_Cn < 1) current_Cn <- 1
    set.seed(setseed)
    mini_batch_kmeans <- ClusterR::MiniBatchKmeans(FXX, clusters = current_Cn, batch_size = 1024,
                                                   num_init = 3, max_iters = 5, initializer = 'kmeans++')
    batchs <- ClusterR::predict_KMeans(FXX, mini_batch_kmeans$centroids)
    cluster_sizes <- table(batchs)
    if (min(cluster_sizes) >= threshold || current_Cn == 1) break
    current_Cn <- current_Cn - 1
    setseed <- setseed + 1
  }
  sort_idx <- order(batchs)
  return(list(R_CGOSS = length(cluster_sizes),
              data_matrix_sorted = FXX[sort_idx, , drop = FALSE],
              sorted_y = y[sort_idx],
              cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
              sorted_indices = (1:nrow(FXX))[sort_idx],
              centroids = mini_batch_kmeans$centroids,
              final_Cn = current_Cn))
}

Comp = function(dataset, obj.c = 0.5) {

  FXX.train <- dataset$FXX.train
  FXX.test  <- dataset$FXX.test
  FY.train  <- dataset$FY.train
  FY.test   <- dataset$FY.test
  C.train   <- dataset$nc.train
  C.test    <- dataset$nc.test
  R         <- dataset$R
  p         <- dataset$p
  setseed   <- 42

  # ── SA 初始化 ──────────────────────────────────────────────
  T.initial <- 50
  K_mat <- calculate_K_dynamic(FXX.test)
  Cn    <- 2

  time2.start  <- Sys.time()
  cluster.curr <- mbky(setseed, FXX.train, FY.train, Cn)
  R_CGOSS.curr <- cluster.curr$R_CGOSS
  FXXXX.curr   <- cluster.curr$data_matrix_sorted
  FYYY.curr    <- cluster.curr$sorted_y
  C.curr       <- cluster.curr$cluster_sizes_vector
  centroids.curr <- cluster.curr$centroids

  info_res_curr <- count_info_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)
  D.curr   <- info_res_curr$D
  A.curr   <- info_res_curr$A
  obj.curr <- 0.7 * log(D.curr) / p + 0.3 * log(A.curr / p)

  obj.best       <- obj.curr
  FXX.best       <- FXXXX.curr
  FY.best        <- FYYY.curr
  C.best         <- C.curr
  R.best         <- R_CGOSS.curr
  Cn.best        <- Cn
  centroids.best <- centroids.curr

  T.curr   <- T.initial
  alpha    <- 0.95
  iter     <- 0
  max_iter <- 50

  # ── SA 主循环 ──────────────────────────────────────────────
  repeat {
    iter <- iter + 1
    if (T.curr < 1e-4 || iter > max_iter) break

    step     <- sample(c(-1, 1, 2), 1)
    Cn.candi <- max(1, Cn + step)
    if (Cn.candi == Cn) Cn.candi <- Cn + 1

    cluster.candi  <- mbky(setseed + iter, FXX.train, FY.train, Cn.candi)
    R.candi  <- cluster.candi$R_CGOSS
    F.candi  <- cluster.candi$data_matrix_sorted
    Y.candi  <- cluster.candi$sorted_y
    C.candi  <- cluster.candi$cluster_sizes_vector

    info_res_candi <- count_info_cpp(F.candi, Y.candi, C.candi, R.candi, p)
    D.candi   <- info_res_candi$D
    A.candi   <- info_res_candi$A
    obj.candi <- 0.7 * log(D.candi) / p + 0.3 * log(A.candi / p)

    delta  <- obj.candi - obj.curr
    accept <- delta > 0 || runif(1) < exp(delta / T.curr)

    if (accept) {
      Cn           <- Cn.candi
      obj.curr     <- obj.candi
      FXXXX.curr   <- F.candi
      FYYY.curr    <- Y.candi
      C.curr       <- C.candi
      R_CGOSS.curr <- R.candi
      centroids.curr <- cluster.candi$centroids
      if (obj.candi > obj.best) {
        obj.best       <- obj.candi
        FXX.best       <- F.candi
        FY.best        <- Y.candi
        C.best         <- C.candi
        R.best         <- R.candi
        Cn.best        <- Cn.candi
        centroids.best <- centroids.curr
      }
    }
    T.curr <- T.curr * alpha
    cat(sprintf("Iter: %d | T: %.4f | Cn: %d | Obj: %.4f | Best: %.4f\n",
                iter, T.curr, Cn, obj.curr, obj.best))
  }

  time.CGOSS <- as.numeric(difftime(Sys.time(), time2.start, units = "secs"))
  cat("SA 用时:", time.CGOSS, "秒\n")

  # ── 估计与预测 ─────────────────────────────────────────────
  GALL.Est   <- Est_hat_cpp(xx = FXX.best, yy = FY.best, C.best, R.best, p)
  GALLRS.Est <- Est_hat_RS_cpp(xx = FXX.best, yy = FY.best, C.best, R.best, p)
  ALL.Est    <- Est_hat_RS_cpp(xx = FXX.train, yy = FY.train, C.train, R, p)

  rec1 <- cbind(
    ALL   = ALL.Est[[1]],
    GALL  = GALL.Est[[1]],
    GALLRS = GALLRS.Est[[1]]
  )
  rec2 <- cbind(
    ALL    = MSPE_tru_RS(FY.test, FXX.test, FXX.train, FY.train,
                         ALL.Est$beta2, ALL.Est$G.hat, ALL.Est$Var.e, C.train, C.test, R),
    GALL   = MSPE_fn(FY.test, FXX.test, FXX.best, FY.best,
                     GALL.Est[[5]], GALL.Est[[6]], GALL.Est[[7]], C.best, centroids.best),
    GALLRS = MSPE_fn_RS(FY.test, FXX.test, FXX.best, FY.best,
                        GALLRS.Est$beta2, GALLRS.Est$G.hat, GALLRS.Est$Var.e, C.best, centroids.best)
  )
  rec3 <- cbind(ALL = ALL.Est[[4]],   GALL = GALL.Est[[4]],   GALLRS = GALLRS.Est[[4]])
  rec4 <- cbind(ALL = ALL.Est[[6]],   GALL = GALL.Est[[6]],   GALLRS = GALLRS.Est[[6]])
  rec5 <- cbind(ALL = ALL.Est[[7]],   GALL = GALL.Est[[7]],   GALLRS = GALLRS.Est[[7]])

  return(list(MSE_beta = rec1, MSPE = rec2, MSE_beta0 = rec3, Var_a = rec4, Var_e = rec5))
}

# ── 读取数据 ───────────────────────────────────────────────
dat <- read.csv("default_clean.csv", header = TRUE)
y   <- as.numeric(dat$S2019400)

exclude <- c("S2019400", "R1208500", "R0000100", "E8033100",
             "R0536401", "R0536402", "S3588900",
             "S2261000", "S2005400", "S2011401")
X <- as.matrix(dat[, setdiff(colnames(dat), exclude)])
for (j in 1:ncol(X)) {
  lo <- min(X[, j]); hi <- max(X[, j])
  X[, j] <- if (hi > lo) 2 * (X[, j] - lo) / (hi - lo) - 1 else 0
}

# ── 分组 + 70/30 分割 ──────────────────────────────────────
make_group_data <- function(group_col, train_ratio = 0.7, seed = 42) {
  set.seed(seed)
  g        <- as.integer(as.factor(dat[[group_col]]))
  # 剔除样本数 < 2 的组
  keep     <- g %in% as.integer(names(table(g)[table(g) >= 2]))
  g        <- as.integer(as.factor(g[keep]))   # 重新编号
  X_k      <- X[keep, ]
  y_k      <- y[keep]
  ord      <- order(g)
  X_s      <- X_k[ord, ]; y_s <- y_k[ord]; g_s <- g[ord]
  groups   <- unique(g_s)
  train_idx <- c(); test_idx <- c()
  for (grp in groups) {
    idx     <- which(g_s == grp)
    n_tr    <- min(max(1, round(length(idx) * train_ratio)), length(idx) - 1)
    sel     <- sample(idx, n_tr)
    train_idx <- c(train_idx, sort(sel))
    test_idx  <- c(test_idx,  sort(setdiff(idx, sel)))
  }
  list(
    FXX.train = X_s[train_idx, ],
    FY.train  = as.matrix(y_s[train_idx]),
    nc.train  = as.integer(table(factor(g_s[train_idx], levels = groups))),
    FXX.test  = X_s[test_idx, ],
    FY.test   = as.matrix(y_s[test_idx]),
    nc.test   = as.integer(table(factor(g_s[test_idx],  levels = groups))),
    R = length(groups),
    p = ncol(X)
  )
}

d_S2261000 <- make_group_data("S2261000")
d_S2005400 <- make_group_data("S2005400")
d_S2011401 <- make_group_data("S2011401")

for (name in c("d_S2261000", "d_S2005400", "d_S2011401")) {
  d <- get(name)
  cat(name, "=> R:", d$R, "| p:", d$p,
      "| Train N:", sum(d$nc.train), "| Test N:", sum(d$nc.test), "\n")
}

# ── 运行 ───────────────────────────────────────────────────
result <- Comp(d_S2005400, obj.c = 0.1)
result
