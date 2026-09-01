
#######################
rm(list=ls())
library(devtools)
library(ClusterR)
library(MASS)
library(mvtnorm)
Rcpp::sourceCpp("lmm_fast.cpp")

filename<-CASE<-"case1"
calculate_K <- function(p) {
  
  k_diag <- c(1, rep(1/3, p - 1))
  K <- diag(k_diag)
  return(K)
}

iboss=function(x,k){
  ind=NULL
  m = ncol(x)
  r = rep(floor(k/2/m),m)
  if(sum(r)<k/2) r[1:((k-2*sum(r))/2)] = r[1:((k-2*sum(r))/2)]+1
  candi=1:nrow(x)
  for (i in 1:m)
  {
    xi = x[candi,i]
    j1 = top_k(xi,r[i])+1
    j2 = bottom_k(xi,r[i])+1
    j = unique(c(j1,j2))
    if(length(j)<2*r[i]) {jj=(1:length(candi))[-j];j=c(j,jj[1:(2*r[i]-length(j))])}
    ind = c(ind,candi[j])
    candi=setdiff(candi,ind)
  }
  return(ind)
}
scalex=function(a){
  2*(a-min(a))/(max(a)-min(a))-1
}
calculate_K_dynamic <- function(FXX_test) {
  X_test_with_intercept <- cbind(1, FXX_test)
  K <- (t(X_test_with_intercept) %*% X_test_with_intercept) / nrow(FXX_test)
  return(K)
}
MSPE_tru=function(fy,fx, sx, sy, beta, Var.a, Var.e, nc,C, R){
  index <- 1
  mv_hat <- c()
  for (i in 1:R) {
    mv_hat[i] <- (Var.a/(Var.e+nc[i]*Var.a)) * sum((sy - cbind(1, sx)%*%beta)[index:(index + nc[i] - 1)]) 
    index <- index + nc[i]
  }
  
  y_hat <- cbind(1, fx)%*%beta + rep(mv_hat, C)
  mspe <- mean((fy - y_hat)^2)
  
  return(mspe)
}
MSPE_fn = function(fy_test, fx_test, sx_train, sy_train, beta_hat, Var.a, Var.e, nc_train, centroids){
  R <- length(nc_train)
  mv_hat <- numeric(R)
  index <- 1
  fixed_pred_train <- cbind(1, sx_train) %*% beta_hat
  for (i in 1:R) {
    current_indices <- index:(index + nc_train[i] - 1)
    residuals_i <- sy_train[current_indices] - fixed_pred_train[current_indices]
    term1 <- Var.a / (Var.e + nc_train[i] * Var.a)
    term2 <- sum(residuals_i)
    mv_hat[i] <- term1 * term2
    index <- index + nc_train[i]
  }
  pred_labels <- ClusterR::predict_KMeans(fx_test, centroids)
  u_assigned <- mv_hat[pred_labels]
  y_hat <- cbind(1, fx_test) %*% beta_hat + u_assigned
  mspe <- mean((fy_test - y_hat)^2)
  return(mspe)
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
  
  mspe <- mean((fy_test - y_hat)^2)
  return(mspe)
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
    curr_idx_test <- index_test:(index_test + nc_test[i] - 1)
    
    Xg_train <- X_train_full[curr_idx_train, , drop = FALSE]
    yg_train <- sy_train[curr_idx_train]
    
    res_g <- yg_train - Xg_train %*% beta_hat
    M <- solve(Var.e * G_inv + t(Xg_train) %*% Xg_train)
    u_g <- M %*% (t(Xg_train) %*% res_g)
    
    Xg_test <- X_test_full[curr_idx_test, , drop = FALSE]
    y_hat[curr_idx_test] <- Xg_test %*% beta_hat + Xg_test %*% u_g
    
    index_train <- index_train + nc_train[i]
    index_test <- index_test + nc_test[i]
  }
  
  mspe <- mean((fy_test - y_hat)^2)
  return(mspe)
}
generate_groups <- function(R, m, N,V) {
  if (N <= R * m) {
    stop("N must be greater than R * m to ensure all integers are greater than m")
  }
  if(V=="large"){Vu<-300*N/(5*R)}
  if(V=="small"){Vu<-N/(5*R)}
  adjusted_sum <- N - R * m
  random_numbers <- abs(rnorm(R, mean = N/R, sd = Vu))
  random_numbers <- random_numbers / sum(random_numbers) * adjusted_sum
  result <- ceiling(random_numbers + m)
  current_sum <- sum(result)
  difference <- N - current_sum
  
  while (difference != 0) {
    if (difference > 0) {
      idx <- sample(1:R, 1)
      result[idx] <- result[idx] + 1
      difference <- difference - 1
    } else {
      idx <- sample(1:R, 1)
      if (result[idx] - 1 > m) {
        result[idx] <- result[idx] - 1
        difference <- difference + 1
      }
    }
  }
  
  return(result)
}
mbky <- function(setseed, FXX, y, Cn, threshold = 2) {
  # 1. 初始化 current_Cn
  current_Cn <- Cn 
  
  repeat {
    if (current_Cn < 1) current_Cn <- 1
    
    # 2. 这里的 set.seed 会根据上面或下面更新的 setseed 重新生效
    set.seed(setseed) 
    mini_batch_kmeans <- ClusterR::MiniBatchKmeans(FXX, clusters = current_Cn, batch_size = 1024, 
                                                   num_init = 3, max_iters = 5, 
                                                   initializer = 'kmeans++')
    
    batchs <- ClusterR::predict_KMeans(FXX, mini_batch_kmeans$centroids)
    cluster_sizes <- table(batchs)
    
    # 检测是否所有簇的大小都满足阈值，或者已经降无可降
    if (min(cluster_sizes) >= threshold || current_Cn == 1) {
      break
    } else {
      # 如果不满足，聚类数减 1，改变种子以避免陷入相同的局部特征，重新聚类
      current_Cn <- current_Cn - 1
      setseed <- setseed + 1 
    }
  }
  
  R_CGOSS <- length(cluster_sizes)
  sort_idx <- order(batchs)
  
  data_matrix_sorted <- FXX[sort_idx, , drop = FALSE]
  sorted_y <- y[sort_idx]
  sorted_indices <- (1:nrow(FXX))[sort_idx]
  cluster_sizes_vector <- as.vector(table(batchs[sort_idx]))
  
  return(list(R_CGOSS = R_CGOSS, 
              data_matrix_sorted = data_matrix_sorted, 
              sorted_y = sorted_y, 
              cluster_sizes_vector = cluster_sizes_vector, 
              sorted_indices = sorted_indices,
              centroids = mini_batch_kmeans$centroids,
              final_Cn = current_Cn)) # 3. 把最新的聚类数传回给外面的主循环
}

# ================================================================
# supervised_kmy：残差增广聚类 + 保存组参数供软分配使用
# ================================================================
supervised_kmy <- function(setseed, FXX, FY, Cn,
                            lambda    = 1.0,
                            threshold = 5L * (ncol(FXX) + 1L)) {
  p          <- ncol(FXX)
  current_Cn <- Cn
  km_centers_X <- NULL
  batchs       <- NULL
  cluster_sizes <- NULL

  # 1. OLS 残差（只算一次，零额外开销）
  ols_resid    <- as.vector(residuals(lm(FY ~ FXX)))
  resid_scaled <- scale(ols_resid) * lambda
  FXX_aug      <- cbind(FXX, resid_scaled)   # N × (p+1)

  repeat {
    if (current_Cn < 1L) current_Cn <- 1L

    if (current_Cn == 1L) {
      batchs        <- rep(1L, nrow(FXX))
      cluster_sizes <- table(batchs)
      km_centers_X  <- matrix(colMeans(FXX), nrow = 1L)
      break
    }

    set.seed(setseed)
    km            <- stats::kmeans(FXX_aug, centers = current_Cn,
                                   nstart = 25, iter.max = 300)
    batchs        <- km$cluster
    cluster_sizes <- table(batchs)
    km_centers_X  <- km$centers[, seq_len(p), drop = FALSE]  # 只保留X列

    if (min(cluster_sizes) >= threshold) break
    current_Cn <- current_Cn - 1L
    setseed    <- setseed + 1L
  }

  R_now    <- max(batchs)
  sort_idx <- order(batchs)

  # 2. 每组 X 的均值、协方差和先验概率，供软分配使用
  group_params <- vector("list", R_now)
  for (g in seq_len(R_now)) {
    idx_g <- which(batchs == g)
    Xg    <- FXX[idx_g, , drop = FALSE]
    mu_g  <- colMeans(Xg)
    Sigma_g <- if (nrow(Xg) > p + 1L) {
      cov(Xg) + diag(1e-6, p)
    } else {
      diag(1.0, p)   # 组太小时退化为球形协方差
    }
    group_params[[g]] <- list(mu    = mu_g,
                              Sigma = Sigma_g,
                              prior = nrow(Xg) / nrow(FXX))
  }

  return(list(
    R_CGOSS              = R_now,
    data_matrix_sorted   = FXX[sort_idx, , drop = FALSE],
    sorted_y             = FY[sort_idx],
    cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
    sorted_indices       = (1:nrow(FXX))[sort_idx],
    centroids            = km_centers_X,  # X 质心（备用硬分配）
    group_params         = group_params,  # EM 软分配所需
    final_Cn             = current_Cn
  ))
}

# ================================================================
# soft_assign：给单个测试点计算属于各组的后验概率权重
# ================================================================
soft_assign <- function(x_new, group_params) {
  log_probs <- sapply(group_params, function(par) {
    log(par$prior + 1e-12) +
      mvtnorm::dmvnorm(x_new, mean = par$mu,
                       sigma = par$Sigma, log = TRUE)
  })
  log_probs <- log_probs - max(log_probs)   # 数值稳定
  w <- exp(log_probs)
  w / sum(w)
}

# ================================================================
# MSPE_fn_soft：随机截距版预测，使用软分配
# ================================================================
MSPE_fn_soft <- function(fy_test, fx_test,
                          sx_train, sy_train,
                          beta_hat, Var.a, Var.e,
                          nc_train, group_params) {
  R <- length(nc_train)

  # 训练集估计每组随机截距
  mv_hat           <- numeric(R)
  index            <- 1L
  fixed_pred_train <- cbind(1, sx_train) %*% beta_hat
  for (i in seq_len(R)) {
    idx       <- index:(index + nc_train[i] - 1L)
    res_i     <- sy_train[idx] - fixed_pred_train[idx]
    mv_hat[i] <- (Var.a / (Var.e + nc_train[i] * Var.a)) * sum(res_i)
    index     <- index + nc_train[i]
  }

  # 测试集软分配预测
  n_test          <- nrow(fx_test)
  fixed_pred_test <- cbind(1, fx_test) %*% beta_hat
  y_hat           <- numeric(n_test)
  for (j in seq_len(n_test)) {
    w         <- soft_assign(as.numeric(fx_test[j, ]), group_params)
    y_hat[j]  <- fixed_pred_test[j] + sum(w * mv_hat)
  }

  mean((fy_test - y_hat)^2)
}

# ================================================================
# MSPE_fn_RS_soft：随机斜率版预测，使用软分配
# ================================================================
MSPE_fn_RS_soft <- function(fy_test, fx_test,
                             sx_train, sy_train,
                             beta_hat, G_hat, Var.e,
                             nc_train, group_params) {
  R            <- length(nc_train)
  q            <- ncol(fx_test) + 1L
  u_vecs       <- matrix(0, nrow = R, ncol = q)
  X_train_full <- cbind(1, sx_train)
  G_inv        <- solve(G_hat + diag(1e-7, q))
  index        <- 1L

  # 训练集估计每组随机斜率向量
  for (i in seq_len(R)) {
    idx         <- index:(index + nc_train[i] - 1L)
    Xg          <- X_train_full[idx, , drop = FALSE]
    yg          <- sy_train[idx]
    res_g       <- yg - Xg %*% beta_hat
    M           <- solve(Var.e * G_inv + crossprod(Xg))
    u_vecs[i, ] <- M %*% (t(Xg) %*% res_g)
    index       <- index + nc_train[i]
  }

  # 测试集软分配预测
  n_test       <- nrow(fx_test)
  X_test_full  <- cbind(1, fx_test)
  y_hat        <- numeric(n_test)
  for (j in seq_len(n_test)) {
    w        <- soft_assign(as.numeric(fx_test[j, ]), group_params)
    u_soft   <- colSums(w * u_vecs)   # q 维加权随机斜率
    y_hat[j] <- sum(X_test_full[j, ] * (beta_hat + u_soft))
  }

  mean((fy_test - y_hat)^2)
}
kmy <- function(setseed, FXX, y, Cn, threshold = 5 * (ncol(FXX) + 1)) {
  current_Cn <- Cn
  km_centers <- NULL
  
  repeat {
    if (current_Cn < 1) current_Cn <- 1
    
    if (current_Cn == 1) {
      batchs       <- rep(1, nrow(FXX))
      cluster_sizes <- table(batchs)
      km_centers   <- matrix(colMeans(FXX), nrow = 1)
      break
    }
    
    set.seed(setseed)
    km            <- stats::kmeans(FXX, centers = current_Cn, nstart = 25, iter.max = 300)
    batchs        <- km$cluster
    cluster_sizes <- table(batchs)
    km_centers    <- km$centers
    
    if (min(cluster_sizes) >= threshold) break
    current_Cn <- current_Cn - 1
    setseed    <- setseed + 1
  }
  
  sort_idx <- order(batchs)
  
  return(list(
    R_CGOSS              = length(cluster_sizes),
    data_matrix_sorted   = FXX[sort_idx, , drop = FALSE],
    sorted_y             = y[sort_idx],
    cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
    sorted_indices       = (1:nrow(FXX))[sort_idx],
    centroids            = km_centers,
    final_Cn             = current_Cn
  ))
}
findsubforCGOSS<-function(n,R){
  if (n %% R != 0) {
    me=floor(n/R)
    loss=n-me*R
    
    mCGOSS=c(rep(me+1,loss),rep(me,R-loss))
  }
  else{
    mCGOSS=c(rep(n/R,R))
  }
  return(mCGOSS)
}
GOSS<-function(setseed,FXX,FY,n,Cn,p){
  cluster=mbky(setseed,FXX,FY,n,Cn)
  R_CGOSS=cluster$R_CGOSS
  FXXXX <- cluster$data_matrix_sorted
  FYYY<- as.matrix(cluster$sorted_y)
  SCC <- c(0, cumsum(cluster$cluster_sizes_vector))
  mcgoss<-findsubforCGOSS(n,R_CGOSS)
  index.CGOSS <- integer(0) 
  for (i in 1:(length(SCC) - 1)) {
    current_indices <- (SCC[i] + 1):SCC[i + 1]
    index.CGOSS <- c(index.CGOSS, OAJ2_cpp(apply(FXXXX[current_indices, ], 2, scalex), mcgoss[i], tPow=2) + SCC[i])
  }
  index_CGOSS_interation <- cluster$sorted_indices[index.CGOSS]
  ncCGOSS <- mcgoss
  D.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)$D
  A.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)$A
  return(list(index = index_CGOSS_interation,D = D.after,A = A.after,R = R_CGOSS,nc = ncCGOSS,C=cluster$cluster_sizes_vector,FX=FXXXX,FY=FYYY))
}
MSE_LM<-function(xx,yy,beta){
  p<-ncol(xx)
  beta0 <- as.matrix(lm(yy ~ xx)$coefficients)
  mse<-sum((beta0[-1]-beta)^2)
  bt0<-(beta0[1]-1)^2
  return(list(bt0,mse,beta0))
}
MSPE_LM<-function(xx,yy,beta){
  n<-nrow(xx)
  y.est<-cbind(1,xx)%*%beta
  mspe<- mean((yy-y.est)^2)
}



Comp=function(dataset,obj.c=0.5){
  big_column_vector<-c()
  FXX.train = dataset$FXX.train
  FXX.test = dataset$FXX.test
  FY.test <- dataset$FY.test
  FY.train <- dataset$FY.train
  R<-dataset$R
  p<-dataset$p
  C.train<-dataset$nc.train
  C.test<-dataset$nc.test
  
  
  
  
  
  
  
  setseed <- 42
  names=c("ALL.bt.mat","GALL.bt.mat", "GALLRS.bt.mat","ALL.var_a","GALL.var_a","GALLRS.var_a",
          "ALL.pred","GALL.pred", "GALLRS.pred","ALL.var_e","GALL.var_e","GALLRS.var_e",
          "ALL.bt0.dif","GALL.bt0.dif","GALLRS.bt0.dif"
  )
  mat_names=c( "ALL.bt", "GALL.bt", "GALLRS.bt")
  for(name in names) {
    assign(name, matrix(NA, 1, 1), envir = .GlobalEnv)
  }
  for(name in mat_names) {
    assign(name, matrix(NA, p+1, 1), envir = .GlobalEnv)
  }
  
  #######
    time.CGOSS=0
    meanR=0
      
      
      
      index.knowGOSS <-index.CGOSS<- index.GOSS<- index.GIBOSS<- index.knowGIBOSS <- c()
      cpu_time_index_goss<-0
      
       
      
      
      ##############
      
      ####################################################################################################
      
      K_mat <- calculate_K_dynamic(FXX.test)
      Cn=2
      time2.start<-Sys.time()
      
             # 更多迭代
      ################### 标准 SA 初始化 (必须在循环外) ###################
      # 1. 计算初始状态 (Current State)
      cluster.curr <- supervised_kmy(setseed, FXX.train, FY.train, Cn)
      R_CGOSS.curr <- cluster.curr$R_CGOSS
      FXXXX.curr   <- cluster.curr$data_matrix_sorted
      FYYY.curr    <- cluster.curr$sorted_y
      C.curr       <- cluster.curr$cluster_sizes_vector
      centroids.curr    <- cluster.curr$centroids
      group_params.curr <- cluster.curr$group_params
      # 计算初始目标函数值 (避免重复调用 C++ 函数)
      info_res_curr <- count_info_rs_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)
      I.curr <- -sum(diag( solve( info_res_curr$Information) %*% K_mat ))  
      #D.curr<- info_res_curr$D
      #A.curr<- info_res_curr$A
      #obj.curr <- 0.7*log(D.curr)/p+0.3*log(A.curr/p)
      obj.curr <- I.curr
      
      # 2. 初始化全局最优记录 (Global Best)
      obj.best          <- obj.curr
      FXX.best          <- FXXXX.curr
      FY.best           <- FYYY.curr
      C.best            <- C.curr
      R.best            <- R_CGOSS.curr
      Cn.best           <- cluster.curr$final_Cn
      centroids.best    <- centroids.curr
      group_params.best <- group_params.curr
      # 3. SA 参数设置
      T.curr <- 500        
      alpha     <- 0.97            
      iter   <- 0
      max_iter  <- 100            

      
      ################### SA 主循环 ###################
      repeat {
        iter <- iter + 1
        if (T.curr < 1e-4 || iter > max_iter) break
        
        set.seed(setseed + iter * 7)
        step     <- sample(c(-3, -2, -1, 1, 2, 3), 1)
        Cn.candi <- max(2, Cn + step)
        
        cluster.candi <- supervised_kmy(setseed + iter, FXX.train, FY.train, Cn.candi)
        
        R.candi <- cluster.candi$R_CGOSS
        F.candi <- cluster.candi$data_matrix_sorted
        Y.candi <- cluster.candi$sorted_y
        C.candi <- cluster.candi$cluster_sizes_vector
        # q = p+1 个随机效应，至少需要 min_factor * q 个观测才能识别
        min_factor    <- 5L          # 可调，保守用5，激进用3
        min_per_group <- min_factor * (p + 1L)
        
        if (min(C.candi) < min_per_group) {
          # 直接跳过，不计算 info，不参与 accept，但照常降温
          T.curr <- T.curr * alpha
          cat(sprintf("Iter: %d, T: %.4f, Cn: %d, REJECTED (min_size=%d < %d), Best: %.4f, Cn.best: %d\n",
                      iter, T.curr, Cn.candi, min(C.candi), min_per_group, obj.best, Cn.best))
          next
        }
        # ──────────────────────────────────────────────────────────────
        
        info_res_candi <- count_info_rs_cpp(F.candi, Y.candi, C.candi, R.candi, p)
        
        # 目标函数：I-optimality 或 DA-optimality，二选一取消注释
        obj.candi <- -sum(diag(solve(info_res_candi$Information) %*% K_mat))  # I-optimality
        # D.candi <- info_res_candi$D
        # A.candi <- info_res_candi$A
        # obj.candi <- 0.7 * log(D.candi) / p + 0.3 * log(A.candi / p)       # DA-optimality
        
        delta  <- obj.candi - obj.curr
        print(delta)
        accept <- FALSE
        
        if (delta > 0) {
          accept <- TRUE
        } else {
          prob <- exp(delta / T.curr)
          if (runif(1) < prob) accept <- TRUE
        }
        
        if (accept) {
          Cn                <- Cn.candi
          obj.curr          <- obj.candi
          FXXXX.curr        <- F.candi
          FYYY.curr         <- Y.candi
          C.curr            <- C.candi
          R_CGOSS.curr      <- R.candi
          centroids.curr    <- cluster.candi$centroids
          group_params.curr <- cluster.candi$group_params
          if (obj.candi > obj.best) {
            obj.best          <- obj.candi
            FXX.best          <- F.candi
            FY.best           <- Y.candi
            C.best            <- C.candi
            R.best            <- R.candi
            Cn.best           <- cluster.candi$final_Cn
            centroids.best    <- centroids.curr
            group_params.best <- group_params.curr
          }
        }
        
        T.curr <- T.curr * alpha
        cat(sprintf("Iter: %d, T: %.4f, Cn: %d, Obj: %.4f, Best: %.4f, Cn.best: %d\n",
                    iter, T.curr, Cn, obj.curr, obj.best, Cn.best))
      }
      
      meanR <- meanR + R.best
      time2.end <- Sys.time()
      time.CGOSS <- time.CGOSS + as.numeric(difftime(time2.end, time2.start, units = "secs"))
      
      ##############GALLL##############
      GALL.Est  <- Est_hat_cpp(xx=FXX.best, yy=FY.best, C.best, R.best, p)
      GALL.pred <- MSPE_fn_soft(FY.test, FXX.test, FXX.best, FY.best,
                                GALL.Est[[5]], GALL.Est[[6]], GALL.Est[[7]],
                                C.best, group_params.best)
      GALL.bt.mat  <- GALL.Est[[1]]
      GALL.bt0.dif <- GALL.Est[[4]]
      GALL.bt      <- GALL.Est[[5]]
      GALL.var_a   <- GALL.Est[[6]]
      GALL.var_e   <- GALL.Est[[7]]
      ##############GALLLRS##############
      GALLRS.Est  <- Est_hat_RS_cpp(xx=FXX.best, yy=FY.best, C.best, R.best, p)
      GALLRS.pred <- MSPE_fn_RS_soft(FY.test, FXX.test, FXX.best, FY.best,
                                     GALLRS.Est$beta2, GALLRS.Est$G.hat, GALLRS.Est$Var.e,
                                     C.best, group_params.best)
      GALLRS.bt.mat <- GALLRS.Est[[1]]
      GALLRS.bt0.dif <- GALLRS.Est[[4]]
      GALLRS.bt <- GALLRS.Est[[5]]
      GALLRS.var_a <- GALLRS.Est[[6]]
      GALLRS.var_e <- GALLRS.Est[[7]]
      ##############ALLL##############
      ALL.Est <- Est_hat_RS_cpp(xx=FXX.train, yy=FY.train, C.train, R, p)
      ALL.pred <- MSPE_tru_RS(FY.test, FXX.test, FXX.train, FY.train, 
                                    ALL.Est$beta2, ALL.Est$G.hat, ALL.Est$Var.e, 
                                    C.train, C.test, R)
      ALL.bt.mat <- ALL.Est[[1]]
      ALL.bt0.dif <- ALL.Est[[4]]
      ALL.bt <- ALL.Est[[5]]
      ALL.var_a <- ALL.Est[[6]]
      ALL.var_e <- ALL.Est[[7]]
      
      

  ##########################################################
  mse.ALL<-mse.GALL<-mse.GALLRS<-c()
  mse.ALL <- ALL.bt.mat
  mse.GALL <- GALL.bt.mat
  mse.GALLRS <- GALLRS.bt.mat
    
  rec1<-cbind(mse.ALL,mse.GALL,mse.GALLRS)
  
  
  
  ##################################################
  mspe.GALL<-  mspe.ALL<-  mspe.GALLRS <- c()

    mspe.ALL <- ALL.pred
    mspe.GALL <- GALL.pred
    mspe.GALLRS <- GALLRS.pred
  
  
  rec2 <- cbind(mspe.ALL,mspe.GALL,mspe.GALLRS)
  
  ################################################
  mse.bt0.ALL <- mse.bt0.GALL<- mse.bt0.GALLRS <- c()
  
    
    mse.bt0.ALL <- ALL.bt0.dif
    mse.bt0.GALL <- GALL.bt0.dif
    mse.bt0.GALLRS <- GALLRS.bt0.dif
    
  
  
  rec3 <- cbind(mse.bt0.ALL,mse.bt0.GALL,mse.bt0.GALLRS)
  
  
  ##################################################
  Vara.GALL<-  Vara.ALL <-Vara.GALLRS<- c()
 
    Vara.ALL <- ALL.var_a
    Vara.GALL <- GALL.var_a
    Vara.GALLRS <- GALLRS.var_a
  
  rec4 <- cbind(Vara.ALL,Vara.GALL,Vara.GALLRS)
  
  ##################################################
  Vare.GALL<-  Vare.ALL <-Vare.GALLRS<- c()
    Vare.ALL <- ALL.var_e
    Vare.GALL <- GALL.var_e
    Vare.GALLRS <- GALLRS.var_e
  
  
  rec5 <- cbind(Vare.ALL,Vare.GALL,Vare.GALLRS)
  
  
  
  ##################################################
  
  
  
  
  
  return(list(rec1,rec2,rec3,rec4,rec5))
}
##################################


# ── 读取数据 ───────────────────────────────────────────────
dat <- read.csv("2004.csv", header = TRUE)
dat[dat < 0] <- NA
dat$R8378703 <- NULL
dat <- dat[!is.na(dat$R8495200), ]
dat <- dat[complete.cases(dat), ]
cat("初始样本量:", nrow(dat), "\n")

# ── 提取 y 和 X ────────────────────────────────────────────
y <- as.numeric(dat$R8495200)
y <- (y - mean(y)) / sd(y)

exclude <- c("R8495200", "R0000100", "R0173600")
X_cols  <- setdiff(colnames(dat), exclude)
X <- as.matrix(dat[, X_cols])

for (j in 1:ncol(X)) {
  lo <- min(X[, j]); hi <- max(X[, j])
  X[, j] <- if (hi > lo) 2 * (X[, j] - lo) / (hi - lo) - 1 else 0
}

# ── 固定全局 70/30 划分 ────────────────────────────────────
set.seed(42)
n <- nrow(dat)
train_global <- sort(sample(1:n, round(n * 0.7)))
test_global  <- setdiff(1:n, train_global)

# ── 过滤：保证所有分组方案下 train 每组 >= 2，test 每组 >= 1
group_cols <- c("R0214700", "R0214800", "R8496500",
                "R7898500", "R8498600", "R8497000")

for (col in group_cols) {
  repeat {
    valid_train <- names(table(dat[[col]][train_global])[table(dat[[col]][train_global]) >= 2])
    valid_test  <- names(table(dat[[col]][test_global]) [table(dat[[col]][test_global])  >= 1])
    valid_both  <- intersect(valid_train, valid_test)
    
    keep_train <- dat[[col]][train_global] %in% valid_both
    keep_test  <- dat[[col]][test_global]  %in% valid_both
    
    # 如果已经满足条件就退出
    if (all(keep_train) && all(keep_test)) break
    
    train_global <- train_global[keep_train]
    test_global  <- test_global[keep_test]
  }
}

cat("过滤后 Train N:", length(train_global), 
    "| Test N:", length(test_global), "\n\n")

X_train <- X[train_global, ]
X_test  <- X[test_global,  ]
y_train <- y[train_global]
y_test  <- y[test_global]
# ================================================================
# 替换原有的 make_group_data，支持单列或多列组合分组
# 用法与原函数完全兼容：
#   单列：make_group_data_multi("R0214700")
#   多列：make_group_data_multi(c("R0214700", "R8496500"))
# ================================================================

make_group_data_multi <- function(group_cols,
                                  min_train = 2L,
                                  min_test  = 1L) {
  
  # ── 1. 构造组合标签（单列时退化为原始行为）─────────────────
  if (length(group_cols) == 1L) {
    raw_label_train <- as.character(dat[[group_cols]][train_global])
    raw_label_test  <- as.character(dat[[group_cols]][test_global])
  } else {
    # 多列：paste 拼接，例如 "1_2"、"1_3" ...
    raw_label_train <- do.call(paste, c(
      lapply(group_cols, function(col) as.character(dat[[col]][train_global])),
      sep = "_"
    ))
    raw_label_test  <- do.call(paste, c(
      lapply(group_cols, function(col) as.character(dat[[col]][test_global])),
      sep = "_"
    ))
  }
  
  # ── 2. 迭代过滤：保证 train 每组 >= min_train，test 每组 >= min_test
  #       且两侧的组标签完全匹配（train 有的 test 也有，反之亦然）
  keep_train <- rep(TRUE, length(train_global))
  keep_test  <- rep(TRUE, length(test_global))
  
  repeat {
    lbl_tr <- raw_label_train[keep_train]
    lbl_te <- raw_label_test[keep_test]
    
    valid_tr   <- names(which(table(lbl_tr) >= min_train))
    valid_te   <- names(which(table(lbl_te) >= min_test))
    valid_both <- intersect(valid_tr, valid_te)
    
    new_keep_train <- keep_train
    new_keep_test  <- keep_test
    new_keep_train[keep_train]  <- lbl_tr %in% valid_both
    new_keep_test[keep_test]    <- lbl_te %in% valid_both
    
    if (identical(new_keep_train, keep_train) &&
        identical(new_keep_test,  keep_test)) break
    
    keep_train <- new_keep_train
    keep_test  <- new_keep_test
  }
  
  # 过滤后的索引
  idx_train <- train_global[keep_train]
  idx_test  <- test_global[keep_test]
  lbl_tr    <- raw_label_train[keep_train]
  lbl_te    <- raw_label_test[keep_test]
  
  # ── 3. 动态排除分组列（从 X 中去掉参与分组的列）─────────────
  drop_idx <- which(colnames(X_train) %in% group_cols)
  if (length(drop_idx) > 0) {
    X_tr <- X[idx_train, -drop_idx, drop = FALSE]
    X_te <- X[idx_test,  -drop_idx, drop = FALSE]
  } else {
    X_tr <- X[idx_train, , drop = FALSE]
    X_te <- X[idx_test,  , drop = FALSE]
  }
  p_use <- ncol(X_tr)
  
  # ── 4. 按组标签排序 ───────────────────────────────────────
  g_train   <- as.integer(as.factor(lbl_tr))
  ord_train <- order(g_train)
  g_train_s <- g_train[ord_train]
  
  g_test    <- as.integer(as.factor(lbl_te))
  ord_test  <- order(g_test)
  g_test_s  <- g_test[ord_test]
  
  y_tr <- y[idx_train]
  y_te <- y[idx_test]
  
  # ── 5. 返回与原函数完全相同的结构 ────────────────────────
  list(
    FXX.train  = X_tr[ord_train, , drop = FALSE],
    FY.train   = as.matrix(y_tr[ord_train]),
    nc.train   = as.integer(table(g_train_s)),
    FXX.test   = X_te[ord_test,  , drop = FALSE],
    FY.test    = as.matrix(y_te[ord_test]),
    nc.test    = as.integer(table(g_test_s)),
    R          = length(unique(g_train_s)),
    p          = p_use,
    group_cols = group_cols,          # 记录用了哪些列（方便调试）
    group_labels = sort(unique(lbl_tr)) # 记录实际存在的组标签
  )
}


# ================================================================
# 打印分组信息的辅助函数（兼容新旧两种）
# ================================================================
print_group_info <- function(name, d) {
  cat(sprintf(
    "%-25s R=%-4d p=%-4d Train=%-6d Test=%-6d train组大小[%d-%d] test组大小[%d-%d]\n",
    name, d$R, d$p,
    sum(d$nc.train), sum(d$nc.test),
    min(d$nc.train),  max(d$nc.train),
    min(d$nc.test),   max(d$nc.test)
  ))
}


# ================================================================
# 使用示例（替换原来第 791-805 行）
# ================================================================

# ── 原有单列分组（用法不变）──────────────────────────────────
d_race   <- make_group_data_multi("R0214700")          # 种族,    3 组
d_region <- make_group_data_multi("R8496500")          # 地区,    4 组
d_worker <- make_group_data_multi("R7898500")          # 雇主,    5 组
d_urban  <- make_group_data_multi("R8498600")          # 城乡,    3 组
d_edu    <- make_group_data_multi("R8497000")          # 教育,  ~18 组

# ── 新增多列组合分组 ─────────────────────────────────────────
d_race_region  <- make_group_data_multi(c("R0214700", "R8496500"))          # 种族×地区,  最多 3×4=12 组
d_race_urban   <- make_group_data_multi(c("R0214700", "R8498600"))          # 种族×城乡,  最多 3×3=9  组
d_region_urban <- make_group_data_multi(c("R8496500", "R8498600"))          # 地区×城乡,  最多 4×3=12 组
d_edu_region   <- make_group_data_multi(c("R8497000", "R8496500"))          # 教育×地区,  最多 18×4  组
d_race_edu     <- make_group_data_multi(c("R0214700", "R8497000"))          # 种族×教育,  最多 3×18  组

# 三列组合（组数较多，需确认每组样本量足够）
d_race_region_urban <- make_group_data_multi(c("R0214700", "R8496500", "R8498600"))  # 最多 3×4×3=36 组


# ================================================================
# 打印所有分组信息
# ================================================================
cat("=== 分组信息 ===\n")
all_datasets <- list(
  # 单列
  d_race         = d_race,
  d_region       = d_region,
  d_worker       = d_worker,
  d_urban        = d_urban,
  d_edu          = d_edu,
  # 双列组合
  d_race_region  = d_race_region,
  d_race_urban   = d_race_urban,
  d_region_urban = d_region_urban,
  d_edu_region   = d_edu_region,
  d_race_edu     = d_race_edu,
  # 三列组合
  d_race_region_urban = d_race_region_urban
)

for (nm in names(all_datasets)) {
  print_group_info(nm, all_datasets[[nm]])
}


# ================================================================
# 批量跑 Comp 并收集结果
# ================================================================
results <- list()
for (nm in names(all_datasets)) {
  cat("\n>>> 正在运行:", nm, "\n")
  results[[nm]] <- tryCatch(
    Comp(all_datasets[[nm]], obj.c = 0.1),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
}

# ── 汇总 MSPE（rec2）到一张表 ────────────────────────────────
mspe_summary <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  if (is.null(r)) return(NULL)
  rec2 <- r[[2]]   # cbind(mspe.ALL, mspe.GALL, mspe.GALLRS)
  data.frame(
    dataset     = nm,
    R           = all_datasets[[nm]]$R,
    p           = all_datasets[[nm]]$p,
    mspe.ALL    = as.numeric(rec2[1]),
    mspe.GALL   = as.numeric(rec2[2]),
    mspe.GALLRS = as.numeric(rec2[3])
  )
}))

print(mspe_summary)









