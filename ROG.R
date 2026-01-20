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

filename<-CASE<-"case4"

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
MSPE_fn=function(fy,fx, sx, sy, beta, Var.a, Var.e, nc,C, R){
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

mbky <- function(setseed, FXX, y, n, Cn) {
  set.seed(setseed)
  
  repeat {
    # 1. 运行 Mini-Batch K-means
    mini_batch_kmeans <- ClusterR::MiniBatchKmeans(FXX, clusters = Cn, batch_size = 4096, 
                                                   num_init = 3, max_iters = 5, 
                                                   initializer = 'kmeans++')
    
    # 2. 【关键修改】直接使用 C++ 接口预测簇，替代了原来的 assign_clusters 循环
    # 这一步是秒出的，不会卡顿
    batchs <- ClusterR::predict_KMeans(FXX, mini_batch_kmeans$centroids)
    
    # 3. 检查簇大小是否满足条件
    cluster_sizes <- table(batchs)
    threshold <- n / Cn
    
    if (any(cluster_sizes < threshold)) {
      Cn <- Cn - 1  
    } else {
      break  
    }
  }
  
  R_CGOSS <- length(cluster_sizes)
  
  # 4. 【优化排序】直接获取排序索引，避免创建巨大的 data.frame
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
  m=ceiling(n / R)
  sigma=diag(0.5,p,p)+matrix(0.5,p,p)
  #sigma=diag(1,p,p)
  lrs=length(N_all)
  names=c("CGOSS.bt.mat",  "IBOSS.bt.mat","ALL.bt.mat","GALL.bt.mat", "OSS.bt.mat", 
          "knowGOSS.bt.mat", "GOSS.bt.mat","GIBOSS.bt.mat","knowGIBOSS.bt.mat",
          "CGOSS.pred",  "IBOSS.pred","ALL.pred","GALL.pred", "OSS.pred", 
          "knowGOSS.pred", "GOSS.pred","GIBOSS.pred","knowGIBOSS.pred",
          "CGOSS.bt0.dif","IBOSS.bt0.dif","ALL.bt0.dif","GALL.bt0.dif","OSS.bt0.dif",
          "knowGOSS.bt0.dif","GOSS.bt0.dif","GIBOSS.bt0.dif","knowGIBOSS.bt0.dif",
          "CGOSS.Var.a","IBOSS.Var.a","ALL.Var.a","GALL.Var.a","OSS.Var.a",
          "knowGOSS.Var.a","GOSS.Var.a","GIBOSS.Var.a","knowGIBOSS.Var.a",
          "CGOSS.Var.e","IBOSS.Var.e","ALL.Var.e","GALL.Var.e","OSS.Var.e",
          "knowGOSS.Var.e","GOSS.Var.e","GIBOSS.Var.e","knowGIBOSS.Var.e")
  mat_names=c("CGOSS.bt",  "IBOSS.bt", "ALL.bt", "GALL.bt", "OSS.bt", 
              "knowGOSS.bt", "GOSS.bt","GIBOSS.bt","knowGIBOSS.bt")
  for(name in names) {
    assign(name, matrix(NA, 1, nloop*lrs), envir = .GlobalEnv)
  }
  for(name in mat_names) {
    assign(name, matrix(NA, p+1, nloop*lrs), envir = .GlobalEnv)
  }
  
  itr = 0
  
  
  #######
  for (j in 1:lrs) {
    D.CGOSS=0
    A.CGOSS=0
    A.OSS=0
    time.CGOSS=0
    meanR=0
    for (k in 1:nloop) {
      
      N<-N_all[j]
      random_numbers <- generate_groups(R,m,N,groupsize)
      C <- round(random_numbers)
      SC = c(0, cumsum(C))
      
      
      
      if (k%/%100 == k/100) cat(k, "-")
      itr <- itr+1
      set.seed(k* 100000)
      if(dist_a == "N.ori") {Var.a = 0.5; Fa = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C)}
      if(dist_a == "N.ML") { Var.a <- 0; Fa <- 0 }
      if(dist_a == "N.large") {Var.a = 100; Fa = rep(rnorm(R, mean = 0, sd = sqrt(Var.a)), C)}
      if(dist_a=="T"){Var.a = 3;Fa = rep(rt(R,3), C)}
      
      Fe = rnorm(max(SC),mean = 0,sd = sqrt(Var.e))
      FXX = matrix(0, nrow = max(SC), ncol = p)
      index.knowGOSS <-index.CGOSS<- index.GOSS<- index.GIBOSS<- index.knowGIBOSS <- c()
      cpu_time_index_goss<-0
      
      
      ##############
      for (i in 1:R) {
        
        setseed =  k * 100000 + i * 100
        set.seed(setseed)
        if(dist_x=="case1") {FXX[(SC[i] + 1):(SC[i+1]),]=matrix(runif(C[i]*p, -1, 1),C[i],p)}
        if(dist_x=="case2") {FXX[(SC[i] + 1):(SC[i+1]),]=mvrnorm(C[i], rep(0, p), sigma)}
        if(dist_x=="case3") {FXX[(SC[i] + 1):(SC[i+1]),]=matrix(runif(C[i]*p, -1.55+i/20, 0.45+i/20),C[i],p)}
        if(dist_x=="case4") {FXX[(SC[i] + 1):(SC[i+1]),]=mvrnorm(C[i], rep(-2+(i-1)/5, p), sigma) }
        index.knowGOSS <- c(index.knowGOSS, OAJ2_cpp(apply(FXX[(SC[i] + 1):(SC[i+1]),],2,scalex),m, tPow=2) + SC[i])
        index.knowGIBOSS <- c(index.knowGIBOSS, iboss(FXX[(SC[i] + 1):(SC[i+1]),],m) + SC[i])
        Fori <- FXX
      }
      
      
      FYori <- 1 + Fori%*%beta + Fa + Fe
      shuffled_indices <- sample(nrow(FXX))
      shuffled_df <- FXX[shuffled_indices, ]
      rownames(shuffled_df) <- rownames(FXX)[shuffled_indices]
      FXX <- shuffled_df
      
      CC=rep(N/R,R)
      SSC = c(0, cumsum(CC))
      for(i in 1:R){
        
        
        
        index.GOSS <- c(index.GOSS, OAJ2_cpp(apply(FXX[(SSC[i] + 1):(SSC[i+1]),],2,scalex),m, tPow=2) + SSC[i])
        
        index.GIBOSS <- c(index.GIBOSS, iboss(FXX[(SSC[i] + 1):(SSC[i+1]),],m) + SSC[i])
      }
      
      
      
      
      
      FY <- 1 + FXX%*%beta + Fa + Fe
      nc <- rep(m,R)
      
      
      ########################################## OSS
      nc2 <- c()
      index.OSS <- sort(OAJ2_cpp(apply(FXX,2,scalex),n, tPow=2))
      for (i in 1:R) {
        nc.OSS <- which(index.OSS >= (SSC[i] + 1) & index.OSS <= (SSC[i+1]))
        nc2[i] <- length(nc.OSS)
      }
      
      
      #########################CGOSS###########################################################################
      
      T.initial<-R
      time2.start<-Sys.time()
      repeat {
        Cn=1
        cluster= mbky(setseed,FXX,FY,n,Cn)
        R_CGOSS= cluster$R_CGOSS
        FXXXX <- cluster$data_matrix_sorted
        FYYY  <- as.matrix(cluster$sorted_y)
        C = cluster$cluster_sizes_vector
        
        D = count_info_cpp(FXXXX,FYYY,C,R_CGOSS,p)[1]
        A = count_info_cpp(FXXXX,FYYY,C,R_CGOSS,p)[2]
        
        obj.best<- obj.candi <- (obj.c/p)*log(D) - (1-obj.c)*(log(A/p))
        
        
        if (obj.candi >= obj.best) {
            obj.best  <- obj.candi
            obj.before<- obj.candi
            Cn <- Cn + 3
            next
        }
        alpha <- if (obj.before == obj.candi) 0.8 else if (Cn == informat$R) 0.95 else 0.85
        T.cool <- T.initial * alpha ^ Cn
        heatprob <- exp(-(obj.best - obj.candi)/T.cool)
        if (runif(1, min = 0.1, max = 0.9) < heatprob) {
          Cn <- Cn + 3
          obj.before <- obj.candi
          next
        } else {
          break
        }
      }
      
      meanR <- meanR + R_CGOSS
      time2.end<-Sys.time()
      time.CGOSS<-time.CGOSS+as.numeric(difftime(time2.end, time2.start, units = "secs"))
      
      print(time.CGOSS)
      
      ############################################## IBOSS
      nc3 <- c()
      index.IBOSS <- sort(iboss(FXX,n))
      for (i in 1:R) {
        nc.IBOSS <- which(index.IBOSS >= (SC[i] + 1) & index.IBOSS <= (SC[i+1]))
        nc3[i] <- length(nc.IBOSS)
      }
      
      
      ##########################################################  GOSS
    
      GOSS.Est <- Est_hat_cpp(xx=FXX[index.GOSS,], yy=FY[index.GOSS,], 
                              beta, Var.a, Var.e, nc, R, p)
      GOSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.GOSS,], FY[index.GOSS,], 
                                 GOSS.Est[[5]], GOSS.Est[[6]], GOSS.Est[[7]], nc,C, R)
      GOSS.bt.mat[,itr] <- GOSS.Est[[1]]
      GOSS.Var.a[,itr]<- GOSS.Est[[2]]
      GOSS.Var.e[,itr]<- GOSS.Est[[3]]
      GOSS.bt0.dif[,itr] <- GOSS.Est[[4]]
      GOSS.bt[,itr] <- GOSS.Est[[5]]
      
      ############################################################# estimate knowGOSS
      
      knowGOSS.Est <- Est_hat_cpp(xx=Fori[index.knowGOSS,], yy=FYori[index.knowGOSS,], 
                                  beta, Var.a, Var.e, nc, R, p)
      knowGOSS.pred[,itr] <- MSPE_fn(FYori, Fori, Fori[index.knowGOSS,], FYori[index.knowGOSS,], 
                                     knowGOSS.Est[[5]], knowGOSS.Est[[6]], knowGOSS.Est[[7]], nc,C, R)
      knowGOSS.bt.mat[,itr] <- knowGOSS.Est[[1]]
      knowGOSS.Var.a[,itr]<- knowGOSS.Est[[2]]
      knowGOSS.Var.e[,itr]<- knowGOSS.Est[[3]]
      knowGOSS.bt0.dif[,itr] <- knowGOSS.Est[[4]]
      knowGOSS.bt[,itr] <- knowGOSS.Est[[5]]
      
      ############################################################# estimate GIBOSS
      
      GIBOSS.Est <- Est_hat_cpp(xx=FXX[index.GIBOSS,], yy=FY[index.GIBOSS,], 
                                beta, Var.a, Var.e, nc, R, p)
      
      
      
      GIBOSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.GIBOSS,], FY[index.GIBOSS,], 
                                   GIBOSS.Est[[5]], GIBOSS.Est[[6]], GIBOSS.Est[[7]], nc,C, R)
      GIBOSS.bt.mat[,itr] <- GIBOSS.Est[[1]]
      GIBOSS.Var.a[,itr]<- GIBOSS.Est[[2]]
      GIBOSS.Var.e[,itr]<- GIBOSS.Est[[3]]
      GIBOSS.bt0.dif[,itr] <- GIBOSS.Est[[4]]
      GIBOSS.bt[,itr] <- GIBOSS.Est[[5]]
      
      ############################################################# estimate knowGIBOSS
      
      
      
      knowGIBOSS.Est <- Est_hat_cpp(xx=Fori[index.knowGIBOSS,], yy=FYori[index.knowGIBOSS,], 
                                    beta, Var.a, Var.e, nc, R, p)
      
      
      
      knowGIBOSS.pred[,itr] <- MSPE_fn(FYori, Fori, Fori[index.knowGIBOSS,], FYori[index.knowGIBOSS,], 
                                       knowGIBOSS.Est[[5]], knowGIBOSS.Est[[6]], knowGIBOSS.Est[[7]], nc,C, R)
      knowGIBOSS.bt.mat[,itr] <- knowGIBOSS.Est[[1]]
      knowGIBOSS.Var.a[,itr]<- knowGIBOSS.Est[[2]]
      knowGIBOSS.Var.e[,itr]<- knowGIBOSS.Est[[3]]
      knowGIBOSS.bt0.dif[,itr] <- knowGIBOSS.Est[[4]]
      knowGIBOSS.bt[,itr] <- knowGIBOSS.Est[[5]]
      
      
      
      ############################################################# estimate CGOSS
      print(length(FY.est))
      print(length(final_index_CGOSS))
      CGOSS.Est <- Est_hat_cpp(xx=FX.est[final_index_CGOSS,], yy=FY.est[final_index_CGOSS,], 
                               beta, Var.a, Var.e, ncCGOSS, R_CGOSS, p)
      CGOSS.pred[,itr]  <- MSPE_fn(FY.est, FX.est , FXX[final_index_CGOSS,], FY[final_index_CGOSS,], 
                                   CGOSS.Est[[5]], CGOSS.Est[[6]], CGOSS.Est[[7]], ncCGOSS,C.est, R_CGOSS)
      CGOSS.bt.mat[,itr] <- CGOSS.Est[[1]]
      CGOSS.Var.a[,itr]<- CGOSS.Est[[2]]
      CGOSS.Var.e[,itr]<- CGOSS.Est[[3]]
      CGOSS.bt0.dif[,itr] <- CGOSS.Est[[4]]
      CGOSS.bt[,itr] <- CGOSS.Est[[5]]
      
      
      
      ##############GALLL##############
      
      GALL.Est <- Est_hat_cpp(xx=FX.est, yy=FY.est, 
                              beta, Var.a, Var.e, C.est, R, p)
      #ALL.pred[,itr] <- MSPE_fn(FX.est, Fori, FXX[index.GOSS,], FY[index.GOSS,], 
       #                         ALL.Est[[5]], ALL.Est[[6]], ALL.Est[[7]], nc,C, R)
      GALL.pred[,itr]<- 0
      GALL.bt.mat[,itr] <- GALL.Est[[1]]
      GALL.Var.a[,itr]<- GALL.Est[[2]]
      GALL.Var.e[,itr]<- GALL.Est[[3]]
      GALL.bt0.dif[,itr] <- GALL.Est[[4]]
      GALL.bt[,itr] <- GALL.Est[[5]]
      
      
      
      ##############ALLL##############
      
      ALL.Est <- Est_hat_cpp(xx=Fori, yy=FYori, 
                             beta, Var.a, Var.e, C, R, p)
      #ALL.pred[,itr] <- MSPE_fn(FX.est, Fori, FXX[index.GOSS,], FY[index.GOSS,], 
      #                         ALL.Est[[5]], ALL.Est[[6]], ALL.Est[[7]], nc,C, R)
      ALL.pred[,itr]<- 0
      ALL.bt.mat[,itr] <- ALL.Est[[1]]
      ALL.Var.a[,itr]<- ALL.Est[[2]]
      ALL.Var.e[,itr]<- ALL.Est[[3]]
      ALL.bt0.dif[,itr] <- ALL.Est[[4]]
      ALL.bt[,itr] <- ALL.Est[[5]]
      
      
      
      
      
      ##########################################################  OSS
      # OSS.Est <- Est_hat_cpp(xx=FXX[index.OSS,], yy=FY[index.OSS,], 
      #                        beta, Var.a, Var.e, nc2, R, p)
      # OSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.OSS,], FY[index.OSS,], 
      #                           OSS.Est[[5]], OSS.Est[[6]], OSS.Est[[7]], nc2,(N/n)*nc2, R)
      EST_OSS_LM<-MSE_LM(FXX[index.OSS,],FY[index.OSS,],beta)
      OSS.pred[,itr]<-OSS_mspe<-MSPE_LM(FXX,FY,EST_OSS_LM[[3]])
      OSS.bt.mat[,itr] <- EST_OSS_LM[[2]]
      #OSS.Var.a[,itr]<- OSS.Est[[2]]
      #OSS.Var.e[,itr]<- OSS.Est[[3]]
      OSS.bt0.dif[,itr] <- EST_OSS_LM[[1]]
      OSS.bt[,itr] <- EST_OSS_LM[[3]]
      ############################################################# estimate IBOSS
      # IBOSS.Est <- Est_hat_cpp(xx=FXX[index.OSS,], yy=FY[index.OSS,], 
      #                        beta, Var.a, Var.e, nc2, R, p)
      # IBOSS.pred[,itr] <- MSPE_fn(FYori, Fori, FXX[index.OSS,], FY[index.OSS,], 
      #                           OSS.Est[[5]], OSS.Est[[6]], OSS.Est[[7]], nc2,(N/n)*nc2, R)
      EST_IBOSS_LM<-MSE_LM(FXX[index.IBOSS,],FY[index.IBOSS,],beta)
      IBOSS.pred[,itr]<-IBOSS_mspe<-MSPE_LM(FXX,FY,EST_IBOSS_LM[[3]])
      IBOSS.bt.mat[,itr] <- EST_IBOSS_LM[[2]]
      #OSS.Var.a[,itr]<- OSS.Est[[2]]
      #OSS.Var.e[,itr]<- OSS.Est[[3]]
      IBOSS.bt0.dif[,itr] <- EST_IBOSS_LM[[1]]
      IBOSS.bt[,itr] <- EST_IBOSS_LM[[3]]
      
      
      
      
      cat(j,"-",k,"\n")
    }
    
    
    cat("mean time is",time.CGOSS/nloop,"\n")
    cat("mean CGOSS R is",meanR/nloop,"\n")
    
    
    cat("\n\n")
  }
  
  
  ##########################################################
  mse.Goss<-mse.ALL<-mse.GALL<- mse.knowGOSS <- mse.CGOSS <- mse.oss <- mse.iboss <-mse.GIBOSS<-mse.knowGIBOSS<- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mse.knowGOSS <- c(mse.knowGOSS, mean(knowGOSS.bt.mat[,loc]))
    mse.ALL <- c(mse.ALL, mean(ALL.bt.mat[,loc]))
    mse.GALL <- c(mse.GALL, mean(GALL.bt.mat[,loc]))
    mse.Goss <- c(mse.Goss, mean(GOSS.bt.mat[,loc]))
    mse.CGOSS <- c(mse.CGOSS, mean(CGOSS.bt.mat[,loc]))
    mse.iboss <- c(mse.iboss, mean(IBOSS.bt.mat[,loc]))
    mse.oss <- c(mse.oss, mean(OSS.bt.mat[,loc]))
    mse.GIBOSS <- c(mse.GIBOSS, mean(GIBOSS.bt.mat[,loc]))
    mse.knowGIBOSS <- c(mse.knowGIBOSS, mean(knowGIBOSS.bt.mat[,loc]))
    
  }
  
  rec1<-cbind(mse.CGOSS,mse.ALL,mse.GALL, mse.iboss, mse.oss, mse.knowGOSS, mse.Goss,mse.GIBOSS,mse.knowGIBOSS)
  
  ###############################################
  mse.CGOSS.Var.a <- mse.oss.Var.a <- mse.iboss.Var.a  <- mse.Goss.Var.a <- mse.GIBOSS.Var.a<- mse.knowGIBOSS.Var.a<- mse.knowGOSS.Var.a <-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    mse.CGOSS.Var.a <- c(mse.CGOSS.Var.a, mean(CGOSS.Var.a[,loc]))
    mse.knowGOSS.Var.a <- c(mse.knowGOSS.Var.a, mean(knowGOSS.Var.a[,loc]))
    mse.Goss.Var.a <- c(mse.Goss.Var.a, mean(GOSS.Var.a[,loc]))
    mse.GIBOSS.Var.a <- c(mse.GIBOSS.Var.a, mean(GIBOSS.Var.a[,loc]))
    mse.knowGIBOSS.Var.a <- c(mse.knowGIBOSS.Var.a, mean(knowGIBOSS.Var.a[,loc]))
    
  }
  
  rec2<-cbind(mse.CGOSS.Var.a, mse.iboss.Var.a, mse.oss.Var.a, mse.knowGOSS.Var.a, mse.Goss.Var.a,mse.GIBOSS.Var.a,mse.knowGIBOSS.Var.a)
  
  ###############################################
  mse.CGOSS.Var.e <- mse.oss.Var.e <- mse.iboss.Var.e <- mse.Goss.Var.e<- mse.GIBOSS.Var.e<- mse.knowGIBOSS.Var.e <- mse.knowGOSS.Var.e <-c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    mse.CGOSS.Var.e <- c(mse.CGOSS.Var.e, mean(CGOSS.Var.e[,loc]))
    mse.knowGOSS.Var.e <- c(mse.knowGOSS.Var.e, mean(knowGOSS.Var.e[,loc]))
    mse.Goss.Var.e <- c(mse.Goss.Var.e, mean(GOSS.Var.e[,loc]))
    mse.GIBOSS.Var.e <- c(mse.GIBOSS.Var.e, mean(GIBOSS.Var.e[,loc]))
    mse.knowGIBOSS.Var.e <- c(mse.knowGIBOSS.Var.e, mean(knowGIBOSS.Var.e[,loc]))
    
  }
  
  rec3<-cbind(mse.CGOSS.Var.e, mse.iboss.Var.e, mse.oss.Var.e, mse.knowGOSS.Var.e, mse.Goss.Var.e,mse.GIBOSS.Var.e,mse.knowGIBOSS.Var.e)
  
  ##################################################
  mspe.Goss <- mspe.knowGOSS <- mspe.CGOSS<- mspe.oss<- mspe.GIBOSS<- mspe.knowGIBOSS<- mspe.iboss <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    mspe.knowGOSS <- c(mspe.knowGOSS, mean(knowGOSS.pred[,loc]))
    mspe.Goss <- c(mspe.Goss, mean(GOSS.pred[,loc]))
    mspe.CGOSS <- c(mspe.CGOSS, mean(CGOSS.pred[,loc]))
    mspe.iboss <- c(mspe.iboss, mean(IBOSS.pred[,loc]))
    mspe.oss <- c(mspe.oss, mean(OSS.pred[,loc]))
    mspe.GIBOSS <- c(mspe.GIBOSS, mean(GIBOSS.pred[,loc]))
    mspe.knowGIBOSS <- c(mspe.knowGIBOSS, mean(knowGIBOSS.pred[,loc]))
    
  }
  
  rec4 <- cbind(mspe.CGOSS, mspe.iboss, mspe.oss, mspe.knowGOSS, mspe.Goss,mspe.GIBOSS,mspe.knowGIBOSS)
  
  ################################################
  mse.bt0.Goss <- mse.bt0.knowGOSS <- mse.bt0.CGOSS<- mse.bt0.oss<- mse.bt0.GIBOSS<- mse.bt0.knowGIBOSS<- mse.bt0.iboss <- c()
  for (i in 1:lrs) {
    loc <- ((i-1)*nloop+1):(i*nloop)
    
    
    mse.bt0.knowGOSS <- c(mse.bt0.knowGOSS, mean(knowGOSS.bt0.dif[,loc]))
    mse.bt0.Goss <- c(mse.bt0.Goss, mean(GOSS.bt0.dif[,loc]))
    mse.bt0.CGOSS <- c(mse.bt0.CGOSS, mean(CGOSS.bt0.dif[,loc]))
    mse.bt0.iboss <- c(mse.bt0.iboss, mean(IBOSS.bt0.dif[,loc]))
    mse.bt0.oss <- c(mse.bt0.oss, mean(OSS.bt0.dif[,loc]))
    mse.bt0.GIBOSS <- c(mse.bt0.GIBOSS, mean(GIBOSS.bt0.dif[,loc]))
    mse.bt0.knowGIBOSS <- c(mse.bt0.knowGIBOSS, mean(knowGIBOSS.bt0.dif[,loc]))
    
  }
  
  rec5 <- cbind(mse.bt0.CGOSS, mse.bt0.iboss, mse.bt0.oss, mse.bt0.knowGOSS, mse.bt0.Goss,mse.bt0.GIBOSS,mse.bt0.knowGIBOSS)
  
  save(rec1, rec2, rec3, rec4, rec5, file = paste0(dist_a,"_", dist_x,".Rdata"))
  
  return(list(rec1,rec2,rec3,rec4,rec5))
}

#########################



N=c(1e4)
modeltype="N.ori"
result = Comp(N,p=50,R=20,Var.e=9,nloop=20,n=1e3,dist_x =filename, dist_a=modeltype,groupsize="large",setted_cluster=20,obj.c=0.5)
result


png(paste0(filename,"remake_mse_",modeltype,".png"),width=800,height = 600)
plot(log10(N), log10(result[[1]][,1]), xlab=expression(log["10"](N)), ylab = expression(log["10"](MSE)), pch=1,lwd=2,
     ylim=c(min(log10(result[[1]])), max(log10(result[[1]]))),xlim=c(min(log10(N)), max(log10(N))), type="o",main = "")
pp <- dim(result[[1]])[2]
for(i in 2:pp){
  lines(log10(N), log10(result[[1]][,i]), type="o", pch=(i), lty=i, col=i,lwd=2)
}

legend("topright", lty=1:pp, pch=(1:pp), col=1:pp,
       legend=c("CGOSS","IBOSS","OSS","known-GOSS","GOSS","GIBOSS","known-GIBOSS"),cex=0.75)

dev.off()

###########################################
png(paste0(filename,"remake_mspe_",modeltype,".png"),width=800,height = 600)
plot(log10(N), log10(result[[4]][,1]), xlab=expression(log["10"](N)), ylab = expression(log["10"](MSPE)), pch=1,lwd=2,
     ylim=c(min(log10(result[[4]])), max(log10(result[[4]]))),xlim=c(min(log10(N)), max(log10(N))), type="o",main = "")
pp <- dim(result[[4]])[2]
for(i in 2:pp){
  lines(log10(N), log10(result[[4]][,i]), type="o", pch=(i), lty=i, col=i,lwd=2)
}

legend("topright", lty=1:pp, pch=(1:pp), col=1:pp,
       legend=c("CGOSS","IBOSS","OSS","known-GOSS","GOSS","GIBOSS","known-GIBOSS"),cex=0.75)


dev.off()

















