
#######################
rm(list=ls())
library(devtools)
library(ClusterR)
library(MASS)
Rcpp::sourceCpp("lmm_fast.cpp")

filename<-CASE<-"case4"
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

Comp=function(N,p, R_all, Var.e, nloop, n, dist_x="case1", dist_a="N.ori",groupsize,obj.c=0.5){
  big_column_vector<-c()
  beta=rep(1, p)
  sigma=diag(0.5,p,p)+matrix(0.5,p,p)
  #sigma=diag(1,p,p)
  lrs=length(R_all)
  names=c("ALL.bt.mat","GALL.bt.mat", "GALLRS.bt.mat","ALL.var_a","GALL.var_a","GALLRS.var_a",
          "ALL.pred","GALL.pred", "GALLRS.pred","ALL.var_e","GALL.var_e","GALLRS.var_e",
          "ALL.bt0.dif","GALL.bt0.dif","GALLRS.bt0.dif"
  )
  mat_names=c( "ALL.bt", "GALL.bt", "GALLRS.bt")
  for(name in names) {
    assign(name, matrix(NA, 1, nloop*lrs), envir = .GlobalEnv)
  }
  for(name in mat_names) {
    assign(name, matrix(NA, p+1, nloop*lrs), envir = .GlobalEnv)
  }
  
  itr = 0
  
  
  #######
  for (j in 1:lrs) {
    time.CGOSS=0
    meanR=0
    for (k in 1:nloop) {
      R <- R_all[j]
      m <- N/(10*R)
      random_numbers <- generate_groups(R,m,N,groupsize)
      C.test <- round(random_numbers)
      C.train<- 3*C.test
      SC.test = c(0, cumsum(C.test))
      SC.train  = c(0, cumsum(C.train))
      
      
      if (k%/%100 == k/100) cat(k, "-")
      itr <- itr+1
      set.seed(k* 100000)
      if(dist_a == "N.ori") {Var.a = 2.25; Fa.test = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C.test)
      Var.a = 2.25; Fa.train = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C.train)
      }
      if(dist_a == "N.ML") { Var.a <- 0; Fa.train<-Fa.test <- 0 }
      if(dist_a=="T"){Var.a = 3;Fa.test = rep(rt(R,3), C.test)
      Fa.train = rep(rt(R,3), C.train)}
      
      Fe.train = rnorm(max(SC.train),mean = 0,sd = sqrt(Var.e))
      Fe.test = rnorm(max(SC.test),mean = 0,sd = sqrt(Var.e))
      FXX.train = matrix(0, nrow = max(SC.train), ncol = p)
      FXX.test = matrix(0, nrow = max(SC.test), ncol = p)
      
      
      index.knowGOSS <-index.CGOSS<- index.GOSS<- index.GIBOSS<- index.knowGIBOSS <- c()
      cpu_time_index_goss<-0
      
      
      
      sigma2 = matrix(0, p, p)
      for(a in 1:p){
        for(b in 1:p){
          sigma2[a, b] = 0.8^abs(a - b)
        }
      }
      
      
      ##############
      for (i in 1:R) {
        
        setseed =  k * 100000 + i * 100
        set.seed(setseed)
        if(dist_x=="case1") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=matrix(runif(C.train[i]*p, -1, 1),C.train[i],p)
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=matrix(runif(C.test[i]*p, -1, 1),C.test[i],p) }
        
        
        if(dist_x=="case2") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=mvrnorm(C.train[i], rep(0, p), sigma)
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=mvrnorm(C.test[i], rep(0, p), sigma)}
        
        
        if(dist_x=="case3") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=matrix(runif(C.train[i]*p, -1.55+i/20, 0.45+i/20),C.train[i],p)
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=matrix(runif(C.test[i]*p, -1.55+i/20, 0.45+i/20),C.test[i],p)}
        
        
        if(dist_x=="case4") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=mvrnorm(C.train[i], rep(-2+(i-1)/5, p), sigma) 
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=mvrnorm(C.test[i], rep(-2+(i-1)/5, p), sigma)}
        
        if(dist_x=="case5") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            matrix(rt(C.train[i] * p, df = 3), C.train[i], p)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            matrix(rt(C.test[i] * p, df = 3), C.test[i], p)
        }
        
        if(dist_x=="case6") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            matrix(rlnorm(C.train[i] * p, meanlog = 0, sdlog = 1), C.train[i], p)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            matrix(rlnorm(C.test[i] * p, meanlog = 0, sdlog = 1), C.test[i], p)
        }
        
        if(dist_x=="case7") {
          id.train = rbinom(C.train[i], 1, 0.5)
          id.test  = rbinom(C.test[i], 1, 0.5)
          
          X1.train = matrix(rnorm(C.train[i] * p, mean = -2, sd = 1), C.train[i], p)
          X2.train = matrix(rnorm(C.train[i] * p, mean =  2, sd = 1), C.train[i], p)
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            X1.train * (1 - id.train) + X2.train * id.train
          
          X1.test = matrix(rnorm(C.test[i] * p, mean = -2, sd = 1), C.test[i], p)
          X2.test = matrix(rnorm(C.test[i] * p, mean =  2, sd = 1), C.test[i], p)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            X1.test * (1 - id.test) + X2.test * id.test
        }
        
        if(dist_x=="case8") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            mvrnorm(C.train[i], rep(0, p), sigma)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            mvrnorm(C.test[i], rep(1.5, p), sigma)
        }
        
        if(dist_x=="case9") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            mvrnorm(C.train[i], rep(0, p), diag(p))
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            mvrnorm(C.test[i], rep(0, p), 2*sigma2)
        }
        
      }
      
      for (col_idx in 1:p) {
        # 联合训练集和测试集，寻找当前维度的全局最值
        global_min <- min(FXX.train[, col_idx], FXX.test[, col_idx])
        global_max <- max(FXX.train[, col_idx], FXX.test[, col_idx])
        
        # 避免分母为零的极端情况
        if (global_max > global_min) {
          # 使用极大极小线性缩放公式
          FXX.train[, col_idx] <- 2 * (FXX.train[, col_idx] - global_min) / (global_max - global_min) - 1
          FXX.test[, col_idx]  <- 2 * (FXX.test[, col_idx] - global_min) / (global_max - global_min) - 1
        } else {
          # 如果该列所有值都一样，统一置为 0
          FXX.train[, col_idx] <- 0
          FXX.test[, col_idx]  <- 0
        }
      }
      
      
      
      
      FY.test <- 1 + FXX.test%*%beta + Fa.test + Fe.test
      FY.train <- 1 + FXX.train%*%beta + Fa.train + Fe.train
      
      ####################################################################################################
      
      T.initial<-50
      K_mat <- calculate_K_dynamic(FXX.test)
      Cn=2
      time2.start<-Sys.time()
      ################### 标准 SA 初始化 (必须在循环外) ###################
      # 1. 计算初始状态 (Current State)
      cluster.curr <- mbky(setseed, FXX.train, FY.train, Cn)
      R_CGOSS.curr <- cluster.curr$R_CGOSS
      FXXXX.curr   <- cluster.curr$data_matrix_sorted
      FYYY.curr    <- cluster.curr$sorted_y
      C.curr       <- cluster.curr$cluster_sizes_vector
      centroids.curr <- cluster.curr$centroids
      # 计算初始目标函数值 (避免重复调用 C++ 函数)
      info_res_curr <- count_info_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)
      I.curr <- -sum(diag( solve( info_res_curr$Information) %*% K_mat ))  
      D.curr<- info_res_curr$D
      A.curr<- info_res_curr$A
      obj.curr <- 0.7*log(D.curr)/p+0.3*log(A.curr/p)
      
      # 2. 初始化全局最优记录 (Global Best)
      obj.best <- obj.curr
      FXX.best <- FXXXX.curr
      FY.best <- FYYY.curr  # 注意：这里你原本写的是 FY.bestM，但下面更新是 FY.best，建议统一命名为 FY.best
      C.best   <- C.curr
      R.best   <- R_CGOSS.curr
      Cn.best  <- Cn
      centroids.best <- centroids.curr
      # 3. SA 参数设置
      T.curr <- T.initial        
      alpha  <- 0.95            
      iter   <- 0
      max_iter <- 50           
      
      ################### SA 主循环 ###################
      repeat {
        iter <- iter + 1
        if (T.curr < 1e-4 || iter > max_iter) break 
        
        # 修正步长逻辑：允许加 1 或减 1 的随机游走
        step <- sample(c(-1,1, 2), 1) 
        Cn.candi <- Cn + step
        
        # 防止聚类数异常 (不能小于 1)
        if (Cn.candi < 1) Cn.candi <- 1
        # 如果因为小于 1 被拉回导致 Cn.candi 等于 Cn，强制让它向右走
        if (Cn.candi == Cn) { 
          Cn.candi <- Cn + 1 
        }
        
        cluster.candi <- mbky(setseed + iter, FXX.train, FY.train, Cn.candi)
        
        R.candi <- cluster.candi$R_CGOSS
        F.candi <- cluster.candi$data_matrix_sorted
        Y.candi <- cluster.candi$sorted_y
        C.candi <- cluster.candi$cluster_sizes_vector
        
        # 计算候选目标函数值 (同样提取出来避免重复调用)
        info_res_candi <- count_info_cpp(F.candi, Y.candi, C.candi, R.candi, p)
        I.candi <- -sum(diag( solve(info_res_candi$Information) %*% K_mat ))
        D.candi<- info_res_candi$D
        A.candi<- info_res_candi$A
        obj.candi <- 0.7*log(D.candi)/p+0.3*log(A.candi/p)
        
        
        
        
        
        
        delta <- obj.candi - obj.curr
        print(delta)
        accept <- FALSE
        
        if (delta > 0) {
          accept <- TRUE
        } else {
          prob <- exp(delta / T.curr)
          if (runif(1) < prob) {
            accept <- TRUE
          }
        }
        
        if (accept) {
          # 补全所有的状态同步，防止新旧状态错乱
          Cn           <- Cn.candi
          obj.curr     <- obj.candi
          FXXXX.curr   <- F.candi
          FYYY.curr    <- Y.candi      # 新增补全
          C.curr       <- C.candi      # 新增补全
          R_CGOSS.curr <- R.candi      # 新增补全
          centroids.curr <- cluster.candi$centroids
          if (obj.candi > obj.best) {
            obj.best <- obj.candi
            FXX.best <- F.candi
            FY.best  <- Y.candi      # 与上面的 FY.bestM 统一一下比较好
            C.best   <- C.candi
            R.best   <- R.candi
            Cn.best  <- Cn.candi
            centroids.best <- centroids.curr
          }
        }
        
        # -------------------------------------------------------
        # E. 降温 (Cooling)
        # -------------------------------------------------------
        T.curr <- T.curr * alpha
        
        # 可选：打印进度
        cat(sprintf("Iter: %d, T: %.4f, Cn: %d, Obj: %.4f, Best: %.4f\n", iter, T.curr, Cn, obj.curr, obj.best))
      }
      
      meanR <- meanR + R.best
      time2.end <- Sys.time()
      time.CGOSS <- time.CGOSS + as.numeric(difftime(time2.end, time2.start, units = "secs"))
      
      print(time.CGOSS)
      
      
      
      
      
      ##############GALLL##############
      
      GALL.Est <- Est_hat_cpp(xx=FXX.best, yy=FY.best, 
                              beta, Var.a, Var.e, C.best, R.best, p)
      GALL.pred[,itr] <- MSPE_fn(FY.test, FXX.test, FXX.best, FY.best, 
                                 GALL.Est[[5]], GALL.Est[[6]], GALL.Est[[7]], 
                                 C.best, centroids.best)
      GALL.bt.mat[,itr] <- GALL.Est[[1]]
      GALL.bt0.dif[,itr] <- GALL.Est[[4]]
      GALL.bt[,itr] <- GALL.Est[[5]]
      GALL.var_a[,itr] <- GALL.Est[[6]]
      GALL.var_e[,itr] <- GALL.Est[[7]]
      ##############GALLLRS##############
      GALLRS.Est <- Est_hat_RS_cpp(xx=FXX.best, yy=FY.best, 
                                   beta, Var.a, Var.e, C.best, R.best, p)
      GALLRS.pred[,itr] <- MSPE_fn_RS(FY.test, FXX.test, FXX.best, FY.best, 
                                      GALLRS.Est$beta2, GALLRS.Est$G.hat, GALLRS.Est$Var.e, 
                                      C.best, centroids.best)
      GALLRS.bt.mat[,itr] <- GALLRS.Est[[1]]
      GALLRS.bt0.dif[,itr] <- GALLRS.Est[[4]]
      GALLRS.bt[,itr] <- GALLRS.Est[[5]]
      GALLRS.var_a[,itr] <- GALLRS.Est[[6]]
      GALLRS.var_e[,itr] <- GALLRS.Est[[7]]
      ##############ALLL##############
      ALL.Est <- Est_hat_RS_cpp(xx=FXX.train, yy=FY.train, 
                                beta, Var.a, Var.e, C.train, R, p)
      ALL.pred[,itr] <- MSPE_tru_RS(FY.test, FXX.test, FXX.train, FY.train, 
                                    ALL.Est$beta2, ALL.Est$G.hat, ALL.Est$Var.e, 
                                    C.train, C.test, R)
      ALL.bt.mat[,itr] <- ALL.Est[[1]]
      ALL.bt0.dif[,itr] <- ALL.Est[[4]]
      ALL.bt[,itr] <- ALL.Est[[5]]
      ALL.var_a[,itr] <- ALL.Est[[6]]
      ALL.var_e[,itr] <- ALL.Est[[7]]
      
      
      cat(j,"-",k,"\n")
    }
    
    
    cat("mean time is",time.CGOSS/nloop,"\n")
    cat("mean CGOSS R is",meanR/nloop,"\n")
    
    
    cat("\n\n")
  }
  
  
  ##########################################################
  mse.ALL<-mse.GALL<-mse.GALLRS<-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.ALL <- c(mse.ALL, mean(ALL.bt.mat[,loc]))
    mse.GALL <- c(mse.GALL, mean(GALL.bt.mat[,loc]))
    mse.GALLRS <- c(mse.GALLRS, mean(GALLRS.bt.mat[,loc]))
    
  }
  
  rec1<-cbind(mse.ALL,mse.GALL,mse.GALLRS)
  
  
  
  ##################################################
  mspe.GALL<-  mspe.ALL<-  mspe.GALLRS <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mspe.ALL <- c(mspe.ALL, mean(ALL.pred[,loc]))
    mspe.GALL <- c(mspe.GALL, mean(GALL.pred[,loc]))
    mspe.GALLRS <- c(mspe.GALLRS, mean(GALLRS.pred[,loc]))
  }
  
  rec2 <- cbind(mspe.ALL,mspe.GALL,mspe.GALLRS)
  
  ################################################
  mse.bt0.ALL <- mse.bt0.GALL<- mse.bt0.GALLRS <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.bt0.ALL <- c(mse.bt0.ALL, mean(ALL.bt0.dif[,loc]))
    mse.bt0.GALL <- c(mse.bt0.GALL, mean(GALL.bt0.dif[,loc]))
    mse.bt0.GALLRS <- c(mse.bt0.GALLRS, mean(GALLRS.bt0.dif[,loc]))
    
  }
  
  rec3 <- cbind(mse.bt0.ALL,mse.bt0.GALL,mse.bt0.GALLRS)
  
  
  ##################################################
  Vara.GALL<-  Vara.ALL <-Vara.GALLRS<- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    Vara.ALL <- c(Vara.ALL, mean(ALL.var_a[,loc]))
    Vara.GALL <- c(Vara.GALL, mean(GALL.var_a[,loc]))
    Vara.GALLRS <- c(Vara.GALLRS, mean(GALLRS.var_a[,loc]))
  }
  
  rec4 <- cbind(Vara.ALL,Vara.GALL,Vara.GALLRS)
  
  ##################################################
  Vare.GALL<-  Vare.ALL <-Vare.GALLRS<- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    Vare.ALL <- c(Vare.ALL, mean(ALL.var_e[,loc]))
    Vare.GALL <- c(Vare.GALL, mean(GALL.var_e[,loc]))
    Vare.GALLRS <- c(Vare.GALLRS, mean(GALLRS.var_e[,loc]))
  }
  
  rec5 <- cbind(Vare.ALL,Vare.GALL,Vare.GALLRS)
  
  
  
  ##################################################
  
  
  
  
  
  save(rec1, rec2, rec3,rec4,rec5, file = paste0("fixedslope_", dist_x,".Rdata"))
  
  return(list(rec1,rec2,rec3,rec4,rec5))
}

##################################
Comp_RS=function(N,p, R_all, Var.e, nloop, n, dist_x="case1", dist_a="N.ori",groupsize,obj.c=0.5){
  big_column_vector<-c()
  beta=rep(1, p)
  sigma=diag(0.5,p,p)+matrix(0.5,p,p)
  #sigma=diag(1,p,p)
  lrs=length(R_all)
  names=c("ALL.bt.mat","GALL.bt.mat", "GALLRS.bt.mat","ALL.var_a","GALL.var_a","GALLRS.var_a",
          "ALL.pred","GALL.pred", "GALLRS.pred","ALL.var_e","GALL.var_e","GALLRS.var_e",
          "ALL.bt0.dif","GALL.bt0.dif","GALLRS.bt0.dif"
  )
  mat_names=c( "ALL.bt", "GALL.bt", "GALLRS.bt")
  for(name in names) {
    assign(name, matrix(NA, 1, nloop*lrs), envir = .GlobalEnv)
  }
  for(name in mat_names) {
    assign(name, matrix(NA, p+1, nloop*lrs), envir = .GlobalEnv)
  }
  
  itr = 0
  
  
  #######
  for (j in 1:lrs) {
    time.CGOSS=0
    meanR=0
    for (k in 1:nloop) {
      R <- R_all[j]
      m <- N/(10*R)
      random_numbers <- generate_groups(R,m,N,groupsize)
      C.test <- round(random_numbers)
      C.train<- 3*C.test
      SC.test = c(0, cumsum(C.test))
      SC.train  = c(0, cumsum(C.train))
      
      
      if (k%/%100 == k/100) cat(k, "-")
      itr <- itr+1
      set.seed(k* 100000)
      if(dist_a == "N.ori") {Var.a = 2.25; Fa.test = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C.test)
      Var.a = 2.25; Fa.train = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C.train)
      }
      if(dist_a == "N.ML") { Var.a <- 0; Fa.train<-Fa.test <- 0 }
      if(dist_a=="T"){Var.a = 3;Fa.test = rep(rt(R,3), C.test)
      Fa.train = rep(rt(R,3), C.train)}
      
      Fe.train = rnorm(max(SC.train),mean = 0,sd = sqrt(Var.e))
      Fe.test = rnorm(max(SC.test),mean = 0,sd = sqrt(Var.e))
      FXX.train = matrix(0, nrow = max(SC.train), ncol = p)
      FXX.test = matrix(0, nrow = max(SC.test), ncol = p)
      
      
      index.knowGOSS <-index.CGOSS<- index.GOSS<- index.GIBOSS<- index.knowGIBOSS <- c()
      cpu_time_index_goss<-0
      
      
      
      sigma2 = matrix(0, p, p)
      for(a in 1:p){
        for(b in 1:p){
          sigma2[a, b] = 0.8^abs(a - b)
        }
      }
      
      
      ##############
      for (i in 1:R) {
        
        setseed =  k * 100000 + i * 100
        set.seed(setseed)
        if(dist_x=="case1") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=matrix(runif(C.train[i]*p, -1, 1),C.train[i],p)
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=matrix(runif(C.test[i]*p, -1, 1),C.test[i],p) }
        
        
        if(dist_x=="case2") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=mvrnorm(C.train[i], rep(0, p), sigma)
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=mvrnorm(C.test[i], rep(0, p), sigma)}
        
        
        if(dist_x=="case3") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=matrix(runif(C.train[i]*p, -1.55+i/20, 0.45+i/20),C.train[i],p)
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=matrix(runif(C.test[i]*p, -1.55+i/20, 0.45+i/20),C.test[i],p)}
        
        
        if(dist_x=="case4") {FXX.train[(SC.train[i] + 1):(SC.train[i+1]),]=mvrnorm(C.train[i], rep(-2+(i-1)/5, p), sigma) 
        FXX.test[(SC.test[i] + 1):(SC.test[i+1]),]=mvrnorm(C.test[i], rep(-2+(i-1)/5, p), sigma)}
        
        if(dist_x=="case5") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            matrix(rt(C.train[i] * p, df = 3), C.train[i], p)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            matrix(rt(C.test[i] * p, df = 3), C.test[i], p)
        }
        
        if(dist_x=="case6") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            matrix(rlnorm(C.train[i] * p, meanlog = 0, sdlog = 1), C.train[i], p)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            matrix(rlnorm(C.test[i] * p, meanlog = 0, sdlog = 1), C.test[i], p)
        }
        
        if(dist_x=="case7") {
          id.train = rbinom(C.train[i], 1, 0.5)
          id.test  = rbinom(C.test[i], 1, 0.5)
          
          X1.train = matrix(rnorm(C.train[i] * p, mean = -2, sd = 1), C.train[i], p)
          X2.train = matrix(rnorm(C.train[i] * p, mean =  2, sd = 1), C.train[i], p)
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            X1.train * (1 - id.train) + X2.train * id.train
          
          X1.test = matrix(rnorm(C.test[i] * p, mean = -2, sd = 1), C.test[i], p)
          X2.test = matrix(rnorm(C.test[i] * p, mean =  2, sd = 1), C.test[i], p)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            X1.test * (1 - id.test) + X2.test * id.test
        }
        
        if(dist_x=="case8") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            mvrnorm(C.train[i], rep(0, p), sigma)
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            mvrnorm(C.test[i], rep(1.5, p), sigma)
        }
        
        if(dist_x=="case9") {
          FXX.train[(SC.train[i] + 1):(SC.train[i+1]), ] =
            mvrnorm(C.train[i], rep(0, p), diag(p))
          FXX.test[(SC.test[i] + 1):(SC.test[i+1]), ] =
            mvrnorm(C.test[i], rep(0, p), 2*sigma2)
        }
        
      }
      
      
      
      for (col_idx in 1:p) {
        # 联合训练集和测试集，寻找当前维度的全局最值
        global_min <- min(FXX.train[, col_idx], FXX.test[, col_idx])
        global_max <- max(FXX.train[, col_idx], FXX.test[, col_idx])
        
        # 避免分母为零的极端情况
        if (global_max > global_min) {
          # 使用极大极小线性缩放公式
          FXX.train[, col_idx] <- 2 * (FXX.train[, col_idx] - global_min) / (global_max - global_min) - 1
          FXX.test[, col_idx]  <- 2 * (FXX.test[, col_idx] - global_min) / (global_max - global_min) - 1
        } else {
          # 如果该列所有值都一样，统一置为 0
          FXX.train[, col_idx] <- 0
          FXX.test[, col_idx]  <- 0
        }
      }
      
      
      
      
      sigma.a <- 2.25              # 随机截距标准差
      sigma.b <- 0.1            # 随机斜率标准差
      sigma.e <- 9              # 误差标准差
      
      Fa.train <- rep(NA, sum(C.train))
      Fa.test  <- rep(NA, sum(C.test))
      FY.train <- rep(NA, sum(C.train))
      FY.test  <- rep(NA, sum(C.test))
      
      for(i in 1:R){
        
        idx.train <- (SC.train[i] + 1):(SC.train[i+1])
        idx.test  <- (SC.test[i] + 1):(SC.test[i+1])
        
        # 第 i 组随机截距
        a.i.train <- rnorm(1, mean = 0, sd = sqrt(sigma.a))
        a.i.test <- rnorm(1, mean = 0, sd = sqrt(sigma.a))
        # 第 i 组随机斜率向量
        
        b.i <- rnorm(p, mean = 0, sd = sqrt(sigma.b))
        
        # 误差项
        e.train <- rnorm(C.train[i], mean = 0, sd = sqrt(sigma.e))
        e.test  <- rnorm(C.test[i], mean = 0, sd = sqrt(sigma.e))
        
        # 保存随机截距
        Fa.train[idx.train] <- a.i.train
        Fa.test[idx.test]   <- a.i.test
        
        # 生成响应
        FY.train[idx.train] <-1+ a.i.train +
          FXX.train[idx.train, ] %*% (beta + b.i) + e.train
        
        FY.test[idx.test] <-1+ a.i.test +
          FXX.test[idx.test, ] %*% (beta + b.i) + e.test
      }
      
      
      
      
      
      
      ####################################################################################################
      
      T.initial<-50
      Cn=2
      K_mat <- calculate_K_dynamic(FXX.test)
      time2.start<-Sys.time()
      ################### 标准 SA 初始化 (必须在循环外) ###################
      # 1. 计算初始状态 (Current State)
      cluster.curr <- mbky(setseed, FXX.train, FY.train, Cn)
      R_CGOSS.curr <- cluster.curr$R_CGOSS
      FXXXX.curr   <- cluster.curr$data_matrix_sorted
      FYYY.curr    <- cluster.curr$sorted_y
      C.curr       <- cluster.curr$cluster_sizes_vector
      centroids.curr <- cluster.curr$centroids
      # 计算初始目标函数值 (避免重复调用 C++ 函数)
      info_res_curr <- count_info_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)
      I.curr <- -sum(diag( solve( info_res_curr$Information) %*% K_mat ))  
      D.curr<- info_res_curr$D
      A.curr<- info_res_curr$A
      obj.curr <- 0.7*log(D.curr)/p+0.3*log(A.curr/p)
      
      # 2. 初始化全局最优记录 (Global Best)
      obj.best <- obj.curr
      FXX.best <- FXXXX.curr
      FY.best <- FYYY.curr  # 注意：这里你原本写的是 FY.bestM，但下面更新是 FY.best，建议统一命名为 FY.best
      C.best   <- C.curr
      R.best   <- R_CGOSS.curr
      Cn.best  <- Cn
      centroids.best <- centroids.curr
      # 3. SA 参数设置
      T.curr <- T.initial        
      alpha  <- 0.95            
      iter   <- 0
      max_iter <- 50           
      
      ################### SA 主循环 ###################
      repeat {
        iter <- iter + 1
        if (T.curr < 1e-4 || iter > max_iter) break 
        
        # 修正步长逻辑：允许加 1 或减 1 的随机游走
        step <- sample(c(-1,1, 2), 1) 
        Cn.candi <- Cn + step
        
        # 防止聚类数异常 (不能小于 1)
        if (Cn.candi < 1) Cn.candi <- 1
        # 如果因为小于 1 被拉回导致 Cn.candi 等于 Cn，强制让它向右走
        if (Cn.candi == Cn) { 
          Cn.candi <- Cn + 1 
        }
        
        cluster.candi <- mbky(setseed + iter, FXX.train, FY.train, Cn.candi)
        
        R.candi <- cluster.candi$R_CGOSS
        F.candi <- cluster.candi$data_matrix_sorted
        Y.candi <- cluster.candi$sorted_y
        C.candi <- cluster.candi$cluster_sizes_vector
        
        # 计算候选目标函数值 (同样提取出来避免重复调用)
        info_res_candi <- count_info_cpp(F.candi, Y.candi, C.candi, R.candi, p)
        I.candi <- -sum(diag( solve(info_res_candi$Information) %*% K_mat ))
        D.candi<- info_res_candi$D
        A.candi<- info_res_candi$A
        obj.candi <- 0.7*log(D.candi)/p+0.3*log(A.candi/p)
        
        
        delta <- obj.candi - obj.curr
        print(delta)
        accept <- FALSE
        
        if (delta > 0) {
          accept <- TRUE
        } else {
          prob <- exp(delta / T.curr)
          if (runif(1) < prob) {
            accept <- TRUE
          }
        }
        
        if (accept) {
          # 补全所有的状态同步，防止新旧状态错乱
          Cn           <- Cn.candi
          obj.curr     <- obj.candi
          FXXXX.curr   <- F.candi
          FYYY.curr    <- Y.candi      # 新增补全
          C.curr       <- C.candi      # 新增补全
          R_CGOSS.curr <- R.candi      # 新增补全
          centroids.curr <- cluster.candi$centroids
          if (obj.candi > obj.best) {
            obj.best <- obj.candi
            FXX.best <- F.candi
            FY.best  <- Y.candi      # 与上面的 FY.bestM 统一一下比较好
            C.best   <- C.candi
            R.best   <- R.candi
            Cn.best  <- Cn.candi
            centroids.best <- centroids.curr
          }
        }
        
        # -------------------------------------------------------
        # E. 降温 (Cooling)
        # -------------------------------------------------------
        T.curr <- T.curr * alpha
        
        # 可选：打印进度
      }
      
      meanR <- meanR + R.best
      time2.end <- Sys.time()
      time.CGOSS <- time.CGOSS + as.numeric(difftime(time2.end, time2.start, units = "secs"))
      
      print(time.CGOSS)
      
      
      
      
      
      ##############GALLL##############
      
      GALL.Est <- Est_hat_cpp(xx=FXX.best, yy=FY.best, 
                              beta, Var.a, Var.e, C.best, R.best, p)
      GALL.pred[,itr] <- MSPE_fn(FY.test, FXX.test, FXX.best, FY.best, 
                                 GALL.Est[[5]], GALL.Est[[6]], GALL.Est[[7]], 
                                 C.best, centroids.best)
      GALL.bt.mat[,itr] <- GALL.Est[[1]]
      GALL.bt0.dif[,itr] <- GALL.Est[[4]]
      GALL.bt[,itr] <- GALL.Est[[5]]
      GALL.var_a[,itr] <- GALL.Est[[6]]
      GALL.var_e[,itr] <- GALL.Est[[7]]
      
      
      
      ##############GALLLRS##############
      GALLRS.Est <- Est_hat_RS_cpp(xx=FXX.best, yy=FY.best, 
                                   beta, Var.a, Var.e, C.best, R.best, p)
      GALLRS.pred[,itr] <- MSPE_fn_RS(FY.test, FXX.test, FXX.best, FY.best, 
                                      GALLRS.Est$beta2, GALLRS.Est$G.hat, GALLRS.Est$Var.e, 
                                      C.best, centroids.best)
      GALLRS.bt.mat[,itr] <- GALLRS.Est[[1]]
      GALLRS.bt0.dif[,itr] <- GALLRS.Est[[4]]
      GALLRS.bt[,itr] <- GALLRS.Est[[5]]
      GALLRS.var_a[,itr] <- GALLRS.Est[[6]]
      GALLRS.var_e[,itr] <- GALLRS.Est[[7]]
      
      
      
      ##############ALLL##############
      ALL.Est <- Est_hat_RS_cpp(xx=FXX.train, yy=FY.train, 
                                beta, Var.a, Var.e, C.train, R, p)
      ALL.pred[,itr] <- MSPE_tru_RS(FY.test, FXX.test, FXX.train, FY.train, 
                                    ALL.Est$beta2, ALL.Est$G.hat, ALL.Est$Var.e, 
                                    C.train, C.test, R)
      ALL.bt.mat[,itr] <- ALL.Est[[1]]
      ALL.bt0.dif[,itr] <- ALL.Est[[4]]
      ALL.bt[,itr] <- ALL.Est[[5]]
      ALL.var_a[,itr] <- ALL.Est[[6]]
      ALL.var_e[,itr] <- ALL.Est[[7]]
      
      
      
      cat(j,"-",k,"\n")
    }
    
    
    cat("mean time is",time.CGOSS/nloop,"\n")
    cat("mean CGOSS R is",meanR/nloop,"\n")
    
    
    cat("\n\n")
  }
  
  
  ##########################################################
  mse.ALL<-mse.GALL<-mse.GALLRS<-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.ALL <- c(mse.ALL, mean(ALL.bt.mat[,loc]))
    mse.GALL <- c(mse.GALL, mean(GALL.bt.mat[,loc]))
    mse.GALLRS <- c(mse.GALLRS, mean(GALLRS.bt.mat[,loc]))
    
  }
  
  rec1<-cbind(mse.ALL,mse.GALL,mse.GALLRS)
  
  
  
  ##################################################
  mspe.GALL<-  mspe.ALL<-  mspe.GALLRS <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mspe.ALL <- c(mspe.ALL, mean(ALL.pred[,loc]))
    mspe.GALL <- c(mspe.GALL, mean(GALL.pred[,loc]))
    mspe.GALLRS <- c(mspe.GALLRS, mean(GALLRS.pred[,loc]))
  }
  
  rec2 <- cbind(mspe.ALL,mspe.GALL,mspe.GALLRS)
  
  ################################################
  mse.bt0.ALL <- mse.bt0.GALL<- mse.bt0.GALLRS <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.bt0.ALL <- c(mse.bt0.ALL, mean(ALL.bt0.dif[,loc]))
    mse.bt0.GALL <- c(mse.bt0.GALL, mean(GALL.bt0.dif[,loc]))
    mse.bt0.GALLRS <- c(mse.bt0.GALLRS, mean(GALLRS.bt0.dif[,loc]))
    
  }
  
  rec3 <- cbind(mse.bt0.ALL,mse.bt0.GALL,mse.bt0.GALLRS)
  
  
  ##################################################
  Vara.GALL<-  Vara.ALL <-Vara.GALLRS<- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    Vara.ALL <- c(Vara.ALL, mean(ALL.var_a[,loc]))
    Vara.GALL <- c(Vara.GALL, mean(GALL.var_a[,loc]))
    Vara.GALLRS <- c(Vara.GALLRS, mean(GALLRS.var_a[,loc]))
  }
  
  rec4 <- cbind(Vara.ALL,Vara.GALL,Vara.GALLRS)
  
  ##################################################
  Vare.GALL<-  Vare.ALL <-Vare.GALLRS<- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    Vare.ALL <- c(Vare.ALL, mean(ALL.var_e[,loc]))
    Vare.GALL <- c(Vare.GALL, mean(GALL.var_e[,loc]))
    Vare.GALLRS <- c(Vare.GALLRS, mean(GALLRS.var_e[,loc]))
  }
  
  rec5 <- cbind(Vare.ALL,Vare.GALL,Vare.GALLRS)
  
  
  
  ##################################################
  
  
  
  
  
  save(rec1, rec2, rec3,rec4,rec5, file = paste0("randomslope_", dist_x,".Rdata"))
  
  return(list(rec1,rec2,rec3,rec4,rec5))
}

#########################



R_all=c(10,20,50)
N=2500
modeltype="N.ori"
result = Comp(N,p=50,R_all,Var.e=9,nloop=50,n=100,dist_x =filename, dist_a=modeltype,groupsize="large",obj.c=0.1)

result_RS = Comp_RS(N,p=50,R_all,Var.e=9,nloop=50,n=100,dist_x =filename, dist_a=modeltype,groupsize="large",obj.c=0.1)

result

result_RS

















