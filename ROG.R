#This R code shows the comparison between different subsample methods on simulated full dataset.
#People can choose 4 types of generated dataset from CASE1-CASE4 by "CASE" in line 23
#People can also choose 3 different setting of random effect: N.ori, N.ML, T by "modeltype" to generate the dataset.
#People can use "large" or "small" to control the differences between the size of generated groups.
#People can use "setted_cluster" to control the cluster number we choose in the algorithm.


#The function Comp is the main function that do the comparison between 5 subsample methods. 
#People can change the value of "n" to have different size of subsample set.


#We use CPP to help speed up the code. People must put 'myoss.cpp','lmm_fast.cpp' to the environment.



#######################
rm(list=ls())
library(devtools)
library(ClusterR)
library(MASS)
Rcpp::sourceCpp('myoss.cpp')
Rcpp::sourceCpp("lmm_fast.cpp")

filename<-CASE<-"case1"

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
MSPE_fn = function(fy, fx, sx, sy, beta, Var.a, Var.e, nc, C, R){
  # nc: 训练集/子样本的组大小
  # C:  在此逻辑中被忽略 (被 nc 的比例投影替代)
  
  index <- 1
  mv_hat <- c()
  
  # 1. 基于训练结构计算随机效应
  for (i in 1:R) {
    if(i <= length(nc)){
      current_indices <- index:(index + nc[i] - 1)
      if(max(current_indices) <= length(sy)){
        term1 <- Var.a / (Var.e + nc[i] * Var.a)
        term2 <- sum((sy - cbind(1, sx) %*% beta)[current_indices])
        mv_hat[i] <- term1 * term2
        index <- index + nc[i]
      } else {
        mv_hat[i] <- 0
      }
    } else {
      mv_hat[i] <- 0
    }
  }
  
  # 2. 比例投影：按比例放大 nc 以适应测试集大小
  # 我们根据 nc 的比例构造一个新的 C 向量 (C_projected)
  
  N_test <- length(fy)          # 测试集的总样本量
  valid_nc <- nc[1:R]           # 确保我们只取对应于 R 个估计效应的组大小
  
  # 计算比例
  props <- valid_nc / sum(valid_nc)
  
  # 将这些比例投影到测试集大小
  C_projected <- floor(props * N_test)
  
  # 3. 修正舍入误差
  # 由于使用了 floor()，总和可能略小于 N_test。
  # 我们将余数分配给前几个组，以确保长度完全匹配。
  remainder <- N_test - sum(C_projected)
  if(remainder > 0){
    C_projected[1:remainder] <- C_projected[1:remainder] + 1
  }
  
  # 4. 预测
  # 现在 sum(C_projected) == length(fy)，所以 rep() 可以完美运行
  y_hat <- cbind(1, fx)%*%beta + rep(mv_hat, times=C_projected)
  
  mspe <- mean((fy - y_hat)^2)
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
mbky <- function(setseed, FXX, y, Cn) {
  set.seed(setseed)
  
  mini_batch_kmeans <- ClusterR::MiniBatchKmeans(FXX, clusters = Cn, batch_size = 500, 
                                                 num_init = 3, max_iters = 5, 
                                                 initializer = 'kmeans++')
  
  # 2. 【关键修改】直接使用 C++ 接口预测簇，替代了原来的 assign_clusters 循环
  # 这一步是秒出的，不会卡顿
  batchs <- ClusterR::predict_KMeans(FXX, mini_batch_kmeans$centroids)
  
  
  
  cluster_sizes <- table(batchs)
  R_CGOSS <- length(cluster_sizes)
  sort_idx <- order(batchs)
  
  # 利用索引直接重排矩阵和向量，速度更快，内存更省
  data_matrix_sorted <- FXX[sort_idx, , drop = FALSE]
  sorted_y <- y[sort_idx]
  sorted_indices <- (1:nrow(FXX))[sort_idx]
  
  # 重新计算排序后的 cluster sizes (其实和上面 table(batchs) 是一样的，但这保证顺序对应)
  # 注意：table 默认按因子水平排序，这里为了保险起见，按出现的顺序或数值统计
  cluster_sizes_vector <- as.vector(table(batchs[sort_idx]))
  
  return(list(R_CGOSS = R_CGOSS, 
              data_matrix_sorted = data_matrix_sorted, 
              sorted_y = sorted_y, 
              cluster_sizes_vector = cluster_sizes_vector, 
              sorted_indices = sorted_indices))
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
  D.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)[1]
  A.after=count_info_cpp(FXX[index_CGOSS_interation,],FY[index_CGOSS_interation,],ncCGOSS,R_CGOSS,p)[2]
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

Comp=function(N_all,p, R, Var.e, nloop, n, dist_x="case1", dist_a="N.ori",groupsize,setted_cluster=1,obj.c=0.5){
  big_column_vector<-c()
  beta=rep(1, p)
  sigma=diag(0.5,p,p)+matrix(0.5,p,p)
  #sigma=diag(1,p,p)
  lrs=length(N_all)
  names=c("ALL.bt.mat","GALL.bt.mat", 
          "ALL.pred","GALL.pred", 
          "ALL.bt0.dif","GALL.bt0.dif"
  )
  mat_names=c( "ALL.bt", "GALL.bt")
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
      m<-N/(10*R)
      N<-N_all[j]
      random_numbers <- generate_groups(R,m,N,groupsize)
      C.test <- round(random_numbers)
      C.train<- 3*C.test
      SC.test = c(0, cumsum(C.test))
      SC.train  = c(0, cumsum(C.train))
      
      
      if (k%/%100 == k/100) cat(k, "-")
      itr <- itr+1
      set.seed(k* 100000)
      if(dist_a == "N.ori") {Var.a = 0.5; Fa.test = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C.test)
      Var.a = 0.5; Fa.train = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C.train)
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
        
      }
      
      FY.test <- 1 + FXX.test%*%beta + Fa.test + Fe.test
      FY.train <- 1 + FXX.train%*%beta + Fa.train + Fe.train
      
      ####################################################################################################
      
      T.initial<-50
      Cn=1
      time2.start<-Sys.time()
      ################### 标准 SA 初始化 (必须在循环外) ###################
      # 1. 计算初始状态 (Current State)
      cluster.curr <- mbky(setseed, FXX.train, FY.train, Cn)
      R_CGOSS.curr <- cluster.curr$R_CGOSS
      FXXXX.curr   <- cluster.curr$data_matrix_sorted
      FYYY.curr    <- cluster.curr$sorted_y
      C.curr       <- cluster.curr$cluster_sizes_vector
      
      # 计算初始目标函数值
      D.curr <- count_info_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)[1]
      A.curr <- count_info_cpp(FXXXX.curr, FYYY.curr, C.curr, R_CGOSS.curr, p)[2]
      obj.curr <- (obj.c/p)*log(D.curr) - (1-obj.c)*(log(A.curr/p))
      
      # 2. 初始化全局最优记录 (Global Best)
      obj.best <- obj.curr
      FXX.best <- FXXXX.curr
      FY.bestM <- FYYY.curr
      C.best   <- C.curr
      R.best   <- R_CGOSS.curr
      Cn.best  <- Cn
      
      # 3. SA 参数设置
      T.curr <- T.initial       
      alpha  <- 0.95            
      iter   <- 0
      max_iter <- 100           
      
      ################### SA 主循环 ###################
      repeat {
        iter <- iter + 1
        if (T.curr < 1e-4 || iter > max_iter) break 
        
        step <- sample(c(0, 1), 1) 
        Cn.candi <- Cn + step
        
        if (Cn.candi < 1) Cn.candi <- 1
        if (Cn.candi == Cn) { 
          Cn.candi <- Cn + 1 
        }
        
        cluster.candi <- mbky(setseed, FXX.train, FY.train, Cn.candi)
        
        R.candi <- cluster.candi$R_CGOSS
        F.candi <- cluster.candi$data_matrix_sorted
        Y.candi <- cluster.candi$sorted_y
        C.candi <- cluster.candi$cluster_sizes_vector
        
        D.candi <- count_info_cpp(F.candi, Y.candi, C.candi, R.candi, p)[1]
        A.candi <- count_info_cpp(F.candi, Y.candi, C.candi, R.candi, p)[2]
        
        obj.candi <- (obj.c/p)*log(D.candi) - (1-obj.c)*(log(A.candi/p))
        
        delta <- obj.candi - obj.curr
        
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
          Cn <- Cn.candi
          obj.curr <- obj.candi
          FXXXX.curr <- F.candi
          
          
          if (obj.candi > obj.best) {
            obj.best <- obj.candi
            FXX.best <- F.candi
            FY.best  <- Y.candi
            C.best   <- C.candi
            R.best   <- R.candi
            Cn.best  <- Cn.candi
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
      time2.end<-Sys.time()
      time.CGOSS<-time.CGOSS+as.numeric(difftime(time2.end, time2.start, units = "secs"))
      
      print(time.CGOSS)
      
      
      
      
      
      ##############GALLL##############
      
      GALL.Est <- Est_hat_cpp(xx=FXX.best, yy=FY.best, 
                              beta, Var.a, Var.e, C.best, R.best, p)
      
      # 修改说明：交换了 C.test 和 C.best 的位置
      # nc 应该是 C.best (子样本组大小), C 应该是 C.test (测试集组大小)
      GALL.pred[,itr] <- MSPE_fn(FY.test, FXX.test, FXX.best, FY.best, 
                                 GALL.Est[[5]], GALL.Est[[6]], GALL.Est[[7]], 
                                 C.best, C.test, R.best)
      
      GALL.bt.mat[,itr] <- GALL.Est[[1]]
      GALL.bt0.dif[,itr] <- GALL.Est[[4]]
      GALL.bt[,itr] <- GALL.Est[[5]]
      
      
      
      ##############ALLL##############
      
      # 修改说明：此处原代码使用了未定义的变量 Train.FX 等，已修正为 FXX.train 等
      # 同时补充了缺失的 C 参数，这里假设全样本训练对应的组大小是 C.train
      ALL.Est <- Est_hat_cpp(xx=FXX.train, yy=FY.train, 
                             beta, Var.a, Var.e, C.train, R, p)
      
      ALL.pred[,itr] <- MSPE_fn(FY.test, FXX.test, FXX.train, FY.train, 
                                ALL.Est[[5]], ALL.Est[[6]], ALL.Est[[7]], 
                                C.train, C.test, R)
      
      ALL.bt.mat[,itr] <- ALL.Est[[1]]
      ALL.bt0.dif[,itr] <- ALL.Est[[4]]
      ALL.bt[,itr] <- ALL.Est[[5]]
      
      
      
      
      
      
      
      
      
      cat(j,"-",k,"\n")
    }
    
    
    cat("mean time is",time.CGOSS/nloop,"\n")
    cat("mean CGOSS R is",meanR/nloop,"\n")
    
    
    cat("\n\n")
  }
  
  
  ##########################################################
  mse.ALL<-mse.GALL<-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.ALL <- c(mse.ALL, mean(ALL.bt.mat[,loc]))
    mse.GALL <- c(mse.GALL, mean(GALL.bt.mat[,loc]))
    
    
  }
  
  rec1<-cbind(mse.ALL,mse.GALL)
  
  
  
  ##################################################
  mspe.GALL<-  mspe.ALL <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mspe.ALL <- c(mspe.ALL, mean(ALL.pred[,itr]))
    mspe.GALL <- c(mspe.GALL, mean(GALL.pred[,itr]))
    
  }
  
  rec2 <- cbind(mspe.ALL,mspe.GALL)
  
  ################################################
  mse.bt0.ALL <- mse.bt0.GALL <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.bt0.ALL <- c(mse.bt0.ALL, mean(ALL.bt0.dif[,itr]))
    mse.bt0.GALL <- c(mse.bt0.GALL, mean(GALL.bt0.dif[,itr]))
    
    
  }
  
  rec3 <- cbind(mse.bt0.ALL,mse.bt0.GALL)
  
  #save(rec1, rec2, rec3, file = paste0(dist_a,"_", dist_x,"NEW.Rdata"))
  
  return(list(rec1,rec2,rec3))
}

#########################



N=c(2500)
modeltype="N.ori"
result = Comp(N,p=50,R=20,Var.e=9,nloop=1,n=100,dist_x =filename, dist_a=modeltype,groupsize="large",setted_cluster=20,obj.c=0.1)
result

















