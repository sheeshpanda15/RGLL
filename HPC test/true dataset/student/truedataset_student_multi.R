
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
  if (N <= R * m) stop("N must be greater than R * m")
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
      idx <- sample(1:R, 1); result[idx] <- result[idx] + 1; difference <- difference - 1
    } else {
      idx <- sample(1:R, 1)
      if (result[idx] - 1 > m) { result[idx] <- result[idx] - 1; difference <- difference + 1 }
    }
  }
  return(result)
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
    current_Cn <- current_Cn - 1; setseed <- setseed + 1
  }
  R_CGOSS <- length(cluster_sizes)
  sort_idx <- order(batchs)
  return(list(R_CGOSS = R_CGOSS,
              data_matrix_sorted = FXX[sort_idx, , drop = FALSE],
              sorted_y = y[sort_idx],
              cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
              sorted_indices = (1:nrow(FXX))[sort_idx],
              centroids = mini_batch_kmeans$centroids,
              final_Cn = current_Cn))
}

# ================================================================
# supervised_kmy：残差增广聚类 + 保存组参数供软分配使用
# ================================================================
supervised_kmy <- function(setseed, FXX, FY, Cn,
                            lambda    = 1.0,
                            threshold = 5L * (ncol(FXX) + 1L)) {
  p          <- ncol(FXX)
  current_Cn <- Cn
  km_centers_X <- NULL; batchs <- NULL; cluster_sizes <- NULL
  ols_resid    <- as.vector(residuals(lm(FY ~ FXX)))
  resid_scaled <- scale(ols_resid) * lambda
  FXX_aug      <- cbind(FXX, resid_scaled)
  repeat {
    if (current_Cn < 1L) current_Cn <- 1L
    if (current_Cn == 1L) {
      batchs <- rep(1L, nrow(FXX)); cluster_sizes <- table(batchs)
      km_centers_X <- matrix(colMeans(FXX), nrow = 1L); break
    }
    set.seed(setseed)
    km <- stats::kmeans(FXX_aug, centers = current_Cn, nstart = 25, iter.max = 300)
    batchs <- km$cluster; cluster_sizes <- table(batchs)
    km_centers_X <- km$centers[, seq_len(p), drop = FALSE]
    if (min(cluster_sizes) >= threshold) break
    current_Cn <- current_Cn - 1L; setseed <- setseed + 1L
  }
  R_now <- max(batchs); sort_idx <- order(batchs)
  group_params <- vector("list", R_now)
  for (g in seq_len(R_now)) {
    idx_g <- which(batchs == g); Xg <- FXX[idx_g, , drop = FALSE]
    mu_g  <- colMeans(Xg)
    Sigma_g <- if (nrow(Xg) > p + 1L) cov(Xg) + diag(1e-6, p) else diag(1.0, p)
    group_params[[g]] <- list(mu = mu_g, Sigma = Sigma_g, prior = nrow(Xg) / nrow(FXX))
  }
  return(list(R_CGOSS = R_now,
              data_matrix_sorted = FXX[sort_idx, , drop = FALSE],
              sorted_y = FY[sort_idx],
              cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
              sorted_indices = (1:nrow(FXX))[sort_idx],
              centroids = km_centers_X,
              group_params = group_params,
              final_Cn = current_Cn))
}

# ================================================================
# soft_assign / MSPE_fn_soft / MSPE_fn_RS_soft
# ================================================================
soft_assign <- function(x_new, group_params) {
  log_probs <- sapply(group_params, function(par) {
    log(par$prior + 1e-12) + mvtnorm::dmvnorm(x_new, mean = par$mu, sigma = par$Sigma, log = TRUE)
  })
  log_probs <- log_probs - max(log_probs)
  w <- exp(log_probs); w / sum(w)
}
MSPE_fn_soft <- function(fy_test, fx_test, sx_train, sy_train,
                          beta_hat, Var.a, Var.e, nc_train, group_params) {
  R <- length(nc_train); mv_hat <- numeric(R); index <- 1L
  fixed_pred_train <- cbind(1, sx_train) %*% beta_hat
  for (i in seq_len(R)) {
    idx <- index:(index + nc_train[i] - 1L)
    mv_hat[i] <- (Var.a / (Var.e + nc_train[i] * Var.a)) * sum(sy_train[idx] - fixed_pred_train[idx])
    index <- index + nc_train[i]
  }
  n_test <- nrow(fx_test); fixed_pred_test <- cbind(1, fx_test) %*% beta_hat; y_hat <- numeric(n_test)
  for (j in seq_len(n_test)) {
    w <- soft_assign(as.numeric(fx_test[j, ]), group_params)
    y_hat[j] <- fixed_pred_test[j] + sum(w * mv_hat)
  }
  mean((fy_test - y_hat)^2)
}
MSPE_fn_RS_soft <- function(fy_test, fx_test, sx_train, sy_train,
                             beta_hat, G_hat, Var.e, nc_train, group_params) {
  R <- length(nc_train); q <- ncol(fx_test) + 1L
  u_vecs <- matrix(0, nrow = R, ncol = q)
  X_train_full <- cbind(1, sx_train); G_inv <- solve(G_hat + diag(1e-7, q)); index <- 1L
  for (i in seq_len(R)) {
    idx <- index:(index + nc_train[i] - 1L); Xg <- X_train_full[idx, , drop = FALSE]; yg <- sy_train[idx]
    M <- solve(Var.e * G_inv + crossprod(Xg))
    u_vecs[i, ] <- M %*% (t(Xg) %*% (yg - Xg %*% beta_hat)); index <- index + nc_train[i]
  }
  n_test <- nrow(fx_test); X_test_full <- cbind(1, fx_test); y_hat <- numeric(n_test)
  for (j in seq_len(n_test)) {
    w <- soft_assign(as.numeric(fx_test[j, ]), group_params)
    y_hat[j] <- sum(X_test_full[j, ] * (beta_hat + colSums(w * u_vecs)))
  }
  mean((fy_test - y_hat)^2)
}
kmy <- function(setseed, FXX, y, Cn, threshold = 5 * (ncol(FXX) + 1)) {
  current_Cn <- Cn; km_centers <- NULL
  repeat {
    if (current_Cn < 1) current_Cn <- 1
    if (current_Cn == 1) {
      batchs <- rep(1, nrow(FXX)); cluster_sizes <- table(batchs)
      km_centers <- matrix(colMeans(FXX), nrow = 1); break
    }
    set.seed(setseed)
    km <- stats::kmeans(FXX, centers = current_Cn, nstart = 25, iter.max = 300)
    batchs <- km$cluster; cluster_sizes <- table(batchs); km_centers <- km$centers
    if (min(cluster_sizes) >= threshold) break
    current_Cn <- current_Cn - 1; setseed <- setseed + 1
  }
  sort_idx <- order(batchs)
  return(list(R_CGOSS = length(cluster_sizes),
              data_matrix_sorted = FXX[sort_idx, , drop = FALSE],
              sorted_y = y[sort_idx],
              cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
              sorted_indices = (1:nrow(FXX))[sort_idx],
              centroids = km_centers, final_Cn = current_Cn))
}
findsubforCGOSS<-function(n,R){
  if (n %% R != 0) { me=floor(n/R); loss=n-me*R; mCGOSS=c(rep(me+1,loss),rep(me,R-loss)) }
  else { mCGOSS=c(rep(n/R,R)) }
  return(mCGOSS)
}
GOSS<-function(setseed,FXX,FY,n,Cn,p){
  cluster=mbky(setseed,FXX,FY,n,Cn); R_CGOSS=cluster$R_CGOSS
  FXXXX <- cluster$data_matrix_sorted; FYYY<- as.matrix(cluster$sorted_y)
  SCC <- c(0, cumsum(cluster$cluster_sizes_vector)); mcgoss<-findsubforCGOSS(n,R_CGOSS)
  index.CGOSS <- integer(0)
  for (i in 1:(length(SCC) - 1)) {
    current_indices <- (SCC[i] + 1):SCC[i + 1]
    index.CGOSS <- c(index.CGOSS, OAJ2_cpp(apply(FXXXX[current_indices, ], 2, scalex), mcgoss[i], tPow=2) + SCC[i])
  }
  index_CGOSS_interation <- cluster$sorted_indices[index.CGOSS]; ncCGOSS <- mcgoss
  D.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)$D
  A.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)$A
  return(list(index=index_CGOSS_interation,D=D.after,A=A.after,R=R_CGOSS,nc=ncCGOSS,C=cluster$cluster_sizes_vector,FX=FXXXX,FY=FYYY))
}
MSE_LM<-function(xx,yy,beta){
  beta0 <- as.matrix(lm(yy ~ xx)$coefficients)
  mse<-sum((beta0[-1]-beta)^2); bt0<-(beta0[1]-1)^2; return(list(bt0,mse,beta0))
}
MSPE_LM<-function(xx,yy,beta){ y.est<-cbind(1,xx)%*%beta; mspe<- mean((yy-y.est)^2) }

Comp=function(dataset,obj.c=0.5){
  FXX.train=dataset$FXX.train; FXX.test=dataset$FXX.test
  FY.test<-dataset$FY.test; FY.train<-dataset$FY.train
  R<-dataset$R; p<-dataset$p; C.train<-dataset$nc.train; C.test<-dataset$nc.test
  setseed <- 42
  names=c("ALL.bt.mat","GALL.bt.mat","GALLRS.bt.mat","ALL.var_a","GALL.var_a","GALLRS.var_a",
          "ALL.pred","GALL.pred","GALLRS.pred","ALL.var_e","GALL.var_e","GALLRS.var_e",
          "ALL.bt0.dif","GALL.bt0.dif","GALLRS.bt0.dif")
  mat_names=c("ALL.bt","GALL.bt","GALLRS.bt")
  for(name in names) assign(name, matrix(NA,1,1), envir=.GlobalEnv)
  for(name in mat_names) assign(name, matrix(NA,p+1,1), envir=.GlobalEnv)
  time.CGOSS=0; meanR=0
  index.knowGOSS<-index.CGOSS<-index.GOSS<-index.GIBOSS<-index.knowGIBOSS<-c()
  K_mat <- calculate_K_dynamic(FXX.test); Cn=2; time2.start<-Sys.time()
  cluster.curr <- supervised_kmy(setseed, FXX.train, FY.train, Cn)
  R_CGOSS.curr<-cluster.curr$R_CGOSS; FXXXX.curr<-cluster.curr$data_matrix_sorted
  FYYY.curr<-cluster.curr$sorted_y 
  C.curr<-cluster.curr$cluster_sizes_vector
  centroids.curr<-cluster.curr$centroids 
  group_params.curr<-cluster.curr$group_params
  info_res_curr <- count_info_rs_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)
  I.curr <- -sum(diag( solve( info_res_curr$Information) %*% K_mat ))  
  #D.curr<- info_res_curr$D
  #A.curr<- info_res_curr$A
  #obj.curr <- 0.7*log(D.curr)/p+0.3*log(A.curr/p)
  obj.curr<-I.curr
  obj.best<-obj.curr; FXX.best<-FXXXX.curr; FY.best<-FYYY.curr; C.best<-C.curr
  R.best<-R_CGOSS.curr; Cn.best<-cluster.curr$final_Cn
  centroids.best<-centroids.curr; group_params.best<-group_params.curr
  T.curr<-500; alpha<-0.97; iter<-0; max_iter<-100
  repeat {
    iter <- iter + 1
    if (T.curr < 1e-4 || iter > max_iter) break
    set.seed(setseed + iter * 7)
    step <- sample(c(-3,-2,-1,1,2,3), 1); Cn.candi <- max(2, Cn + step)
    cluster.candi <- supervised_kmy(setseed + iter, FXX.train, FY.train, Cn.candi)
    R.candi<-cluster.candi$R_CGOSS; F.candi<-cluster.candi$data_matrix_sorted
    Y.candi<-cluster.candi$sorted_y; C.candi<-cluster.candi$cluster_sizes_vector
    min_per_group <- 5L * (p + 1L)
    if (min(C.candi) < min_per_group) {
      T.curr <- T.curr * alpha
      cat(sprintf("Iter: %d, T: %.4f, Cn: %d, REJECTED (min_size=%d < %d), Best: %.4f, Cn.best: %d\n",
                  iter, T.curr, Cn.candi, min(C.candi), min_per_group, obj.best, Cn.best)); next
    }
    info_res_candi <- count_info_rs_cpp(F.candi, Y.candi, C.candi, R.candi, p)
    
    #D.candi <- info_res_candi$D
    # A.candi <- info_res_candi$A
    #obj.candi <- 0.7 * log(D.candi) / p + 0.3 * log(A.candi / p)
    obj.candi<-  -sum(diag( solve( info_res_candi$Information) %*% K_mat )) 
    delta <- obj.candi - obj.curr; print(delta); accept <- FALSE
    if (delta > 0) { accept <- TRUE } else { if (runif(1) < exp(delta/T.curr)) accept <- TRUE }
    if (accept) {
      Cn<-Cn.candi; obj.curr<-obj.candi; FXXXX.curr<-F.candi; FYYY.curr<-Y.candi
      C.curr<-C.candi; R_CGOSS.curr<-R.candi
      centroids.curr<-cluster.candi$centroids; group_params.curr<-cluster.candi$group_params
      if (obj.candi > obj.best) {
        obj.best<-obj.candi; FXX.best<-F.candi; FY.best<-Y.candi; C.best<-C.candi; R.best<-R.candi
        Cn.best<-cluster.candi$final_Cn; centroids.best<-centroids.curr; group_params.best<-group_params.curr
      }
    }
    T.curr <- T.curr * alpha
    cat(sprintf("Iter: %d, T: %.4f, Cn: %d, Obj: %.4f, Best: %.4f, Cn.best: %d\n",
                iter, T.curr, Cn, obj.curr, obj.best, Cn.best))
  }
  meanR <- meanR + R.best; time2.end <- Sys.time()
  time.CGOSS <- time.CGOSS + as.numeric(difftime(time2.end, time2.start, units="secs"))
  GALL.Est  <- Est_hat_cpp(xx=FXX.best, yy=FY.best, C.best, R.best, p)
  GALL.pred <- MSPE_fn_soft(FY.test, FXX.test, FXX.best, FY.best,
                            GALL.Est[[5]], GALL.Est[[6]], GALL.Est[[7]], C.best, group_params.best)
  GALL.bt.mat<-GALL.Est[[1]]; GALL.bt0.dif<-GALL.Est[[4]]
  GALL.bt<-GALL.Est[[5]]; GALL.var_a<-GALL.Est[[6]]; GALL.var_e<-GALL.Est[[7]]
  GALLRS.Est  <- Est_hat_RS_cpp(xx=FXX.best, yy=FY.best, C.best, R.best, p)
  GALLRS.pred <- MSPE_fn_RS_soft(FY.test, FXX.test, FXX.best, FY.best,
                                 GALLRS.Est$beta2, GALLRS.Est$G.hat, GALLRS.Est$Var.e,
                                 C.best, group_params.best)
  GALLRS.bt.mat<-GALLRS.Est[[1]]; GALLRS.bt0.dif<-GALLRS.Est[[4]]
  GALLRS.bt<-GALLRS.Est[[5]]; GALLRS.var_a<-GALLRS.Est[[6]]; GALLRS.var_e<-GALLRS.Est[[7]]
  ALL.Est <- Est_hat_RS_cpp(xx=FXX.train, yy=FY.train, C.train, R, p)
  ALL.pred <- MSPE_tru_RS(FY.test, FXX.test, FXX.train, FY.train,
                          ALL.Est$beta2, ALL.Est$G.hat, ALL.Est$Var.e, C.train, C.test, R)
  ALL.bt.mat<-ALL.Est[[1]]; ALL.bt0.dif<-ALL.Est[[4]]
  ALL.bt<-ALL.Est[[5]]; ALL.var_a<-ALL.Est[[6]]; ALL.var_e<-ALL.Est[[7]]
  rec1<-cbind(ALL.bt.mat, GALL.bt.mat, GALLRS.bt.mat)
  rec2<-cbind(ALL.pred,   GALL.pred,   GALLRS.pred)
  rec3<-cbind(ALL.bt0.dif, GALL.bt0.dif, GALLRS.bt0.dif)
  rec4<-cbind(ALL.var_a,  GALL.var_a,  GALLRS.var_a)
  rec5<-cbind(ALL.var_e,  GALL.var_e,  GALLRS.var_e)
  return(list(rec1,rec2,rec3,rec4,rec5))
}
##################################


# ── 读取数据 ───────────────────────────────────────────────
por <- read.csv("student-por.csv", header = TRUE, sep = ";")
mat <- read.csv("student-mat.csv", header = TRUE, sep = ";")
por$subject <- 1L; mat$subject <- 2L
dat <- rbind(por, mat)
cat("合并后样本量:", nrow(dat), "\n")

dat$parent_edu <- as.integer(dat$Medu + dat$Fedu)
dat$parent_edu[dat$parent_edu < 2] <- 2L
cat("父母教育分组分布:", table(dat$parent_edu), "\n")

cat_cols <- c("school","sex","address","famsize","Pstatus",
              "Mjob","Fjob","reason","guardian",
              "schoolsup","famsup","paid","activities",
              "nursery","higher","internet","romantic")
for (col in cat_cols) dat[[col]] <- as.integer(as.factor(dat[[col]]))

# ── 提取 y 和 X ────────────────────────────────────────────
y <- as.numeric(dat$G3)
y <- (y - mean(y)) / sd(y)

exclude_cols <- c("G3", "G1", "G2", "subject", "parent_edu")
X_cols <- setdiff(colnames(dat), exclude_cols)
X <- as.matrix(dat[, X_cols])
cat("X 列数 p:", ncol(X), "\n")
cat("X 列名:", paste(X_cols, collapse = ", "), "\n\n")

for (j in 1:ncol(X)) {
  lo <- min(X[, j]); hi <- max(X[, j])
  X[, j] <- if (hi > lo) 2 * (X[, j] - lo) / (hi - lo) - 1 else 0
}

# ── 固定全局 70/30 划分 ────────────────────────────────────
set.seed(42)
n <- nrow(dat)
train_global <- sort(sample(1:n, round(n * 0.7)))
test_global  <- setdiff(1:n, train_global)

# ── 过滤：保证所有分组方案下 train 每组 >= 2，test 每组 >= 1 ──
all_group_cols <- c("subject","sex","school","address","Pstatus",
                    "studytime","reason","goout","parent_edu")
for (col in all_group_cols) {
  repeat {
    valid_train <- names(table(dat[[col]][train_global])[table(dat[[col]][train_global]) >= 2])
    valid_test  <- names(table(dat[[col]][test_global]) [table(dat[[col]][test_global])  >= 1])
    valid_both  <- intersect(valid_train, valid_test)
    keep_train  <- dat[[col]][train_global] %in% valid_both
    keep_test   <- dat[[col]][test_global]  %in% valid_both
    if (all(keep_train) && all(keep_test)) break
    train_global <- train_global[keep_train]; test_global <- test_global[keep_test]
  }
}
cat("过滤后 Train N:", length(train_global), "| Test N:", length(test_global), "\n\n")

X_train <- X[train_global, ]; X_test  <- X[test_global,  ]
y_train <- y[train_global];   y_test  <- y[test_global]

# ================================================================
# make_group_data_multi：支持单列或多列组合分组
# 会动态排除参与分组的列（如果它们在 X 中）
# ================================================================
make_group_data_multi <- function(group_cols, min_train = 2L, min_test = 1L) {

  # ── 1. 构造组合标签 ──────────────────────────────────────
  if (length(group_cols) == 1L) {
    raw_label_train <- as.character(dat[[group_cols]][train_global])
    raw_label_test  <- as.character(dat[[group_cols]][test_global])
  } else {
    raw_label_train <- do.call(paste, c(
      lapply(group_cols, function(col) as.character(dat[[col]][train_global])), sep="_"))
    raw_label_test  <- do.call(paste, c(
      lapply(group_cols, function(col) as.character(dat[[col]][test_global])),  sep="_"))
  }

  # ── 2. 迭代过滤 ──────────────────────────────────────────
  keep_train <- rep(TRUE, length(train_global))
  keep_test  <- rep(TRUE, length(test_global))
  repeat {
    lbl_tr <- raw_label_train[keep_train]; lbl_te <- raw_label_test[keep_test]
    valid_both <- intersect(names(which(table(lbl_tr) >= min_train)),
                            names(which(table(lbl_te) >= min_test)))
    nk_tr <- keep_train; nk_tr[keep_train] <- lbl_tr %in% valid_both
    nk_te <- keep_test;  nk_te[keep_test]  <- lbl_te %in% valid_both
    if (identical(nk_tr, keep_train) && identical(nk_te, keep_test)) break
    keep_train <- nk_tr; keep_test <- nk_te
  }
  lbl_tr <- raw_label_train[keep_train]; lbl_te <- raw_label_test[keep_test]

  # ── 3. 动态排除分组列 ────────────────────────────────────
  drop_idx <- which(colnames(X_train) %in% group_cols)
  if (length(drop_idx) > 0) {
    X_tr <- X[train_global[keep_train], -drop_idx, drop = FALSE]
    X_te <- X[test_global[keep_test],   -drop_idx, drop = FALSE]
  } else {
    X_tr <- X[train_global[keep_train], , drop = FALSE]
    X_te <- X[test_global[keep_test],   , drop = FALSE]
  }
  p_use <- ncol(X_tr)

  # ── 4. 按组标签排序 ──────────────────────────────────────
  g_train <- as.integer(as.factor(lbl_tr)); ord_train <- order(g_train); g_train_s <- g_train[ord_train]
  g_test  <- as.integer(as.factor(lbl_te)); ord_test  <- order(g_test);  g_test_s  <- g_test[ord_test]

  list(
    FXX.train    = X_tr[ord_train, , drop = FALSE],
    FY.train     = as.matrix(y_train[keep_train][ord_train]),
    nc.train     = as.integer(table(g_train_s)),
    FXX.test     = X_te[ord_test,  , drop = FALSE],
    FY.test      = as.matrix(y_test[keep_test][ord_test]),
    nc.test      = as.integer(table(g_test_s)),
    R            = length(unique(g_train_s)),
    p            = p_use,
    group_cols   = group_cols,
    group_labels = sort(unique(lbl_tr))
  )
}

# ── 打印分组信息 ────────────────────────────────────────────
print_group_info <- function(name, d) {
  cat(sprintf("%-30s R=%-4d p=%-3d Train=%-6d Test=%-6d train组大小[%d-%d] test组大小[%d-%d]\n",
              name, d$R, d$p, sum(d$nc.train), sum(d$nc.test),
              min(d$nc.train), max(d$nc.train), min(d$nc.test), max(d$nc.test)))
}

# ── 分组方案 ────────────────────────────────────────────────
# 单列（原有方案）
d_subject    <- make_group_data_multi("subject")     # 科目,        2 组
d_sex        <- make_group_data_multi("sex")         # 性别,        2 组
d_school     <- make_group_data_multi("school")      # 学校,        2 组
d_address    <- make_group_data_multi("address")     # 城乡,        2 组
d_Pstatus    <- make_group_data_multi("Pstatus")     # 父母状态,    2 组
d_studytime  <- make_group_data_multi("studytime")   # 每周学习,    4 组
d_reason     <- make_group_data_multi("reason")      # 选校原因,    4 组
d_goout      <- make_group_data_multi("goout")       # 外出频率,    5 组
d_parent_edu <- make_group_data_multi("parent_edu")  # 父母教育,    7 组

# 多列组合（新增方案）
d_subject_sex        <- make_group_data_multi(c("subject", "sex"))         # 科目×性别,       最多 2×2=4  组
d_subject_school     <- make_group_data_multi(c("subject", "school"))      # 科目×学校,       最多 2×2=4  组
d_subject_address    <- make_group_data_multi(c("subject", "address"))     # 科目×城乡,       最多 2×2=4  组
d_subject_parent_edu <- make_group_data_multi(c("subject", "parent_edu"))  # 科目×父母教育,   最多 2×7=14 组
d_sex_parent_edu     <- make_group_data_multi(c("sex", "parent_edu"))      # 性别×父母教育,   最多 2×7=14 组
d_reason_studytime   <- make_group_data_multi(c("reason", "studytime"))    # 选校原因×学习时间, 最多 4×4=16 组
d_goout_studytime    <- make_group_data_multi(c("goout", "studytime"))     # 外出×学习时间,   最多 5×4=20 组
d_subject_sex_school <- make_group_data_multi(c("subject","sex","school")) # 三列组合, 最多 2×2×2=8 组

cat("=== 分组信息 ===\n")
all_datasets <- list(
  # 单列
  d_subject    = d_subject,    d_sex       = d_sex,
  d_school     = d_school,     d_address   = d_address,
  d_Pstatus    = d_Pstatus,    d_studytime = d_studytime,
  d_reason     = d_reason,     d_goout     = d_goout,
  d_parent_edu = d_parent_edu,
  # 多列组合
  d_subject_sex        = d_subject_sex,
  d_subject_school     = d_subject_school,
  d_subject_address    = d_subject_address,
  d_subject_parent_edu = d_subject_parent_edu,
  d_sex_parent_edu     = d_sex_parent_edu,
  d_reason_studytime   = d_reason_studytime,
  d_goout_studytime    = d_goout_studytime,
  d_subject_sex_school = d_subject_sex_school
)
for (nm in names(all_datasets)) print_group_info(nm, all_datasets[[nm]])

# ── 批量运行并汇总结果 ──────────────────────────────────────
results <- list()
for (nm in names(all_datasets)) {
  cat("\n>>> 正在运行:", nm, "\n")
  results[[nm]] <- tryCatch(
    Comp(all_datasets[[nm]], obj.c = 0.1),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
}

mspe_summary <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]; if (is.null(r)) return(NULL)
  rec2 <- r[[2]]
  data.frame(dataset=nm, R=all_datasets[[nm]]$R, p=all_datasets[[nm]]$p,
             mspe.ALL=as.numeric(rec2[1]), mspe.GALL=as.numeric(rec2[2]),
             mspe.GALLRS=as.numeric(rec2[3]))
}))

print(mspe_summary)
