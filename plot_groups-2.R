# ============================================================
# 可视化：SA 重分组前后每组的拟合线对比
#
# 输出 2 张大图（fixed / rs 各一张）
# 每张：3 列(p=2,10,50) × 2 行(Before/After) = 6 个子图
# 每个子图：散点(半透明) + 每组一条拟合线(Est_hat_RS_cpp+BLUP)
# 横轴：PC1，纵轴：Y，颜色区分组
# ============================================================

rm(list = ls())

library(ClusterR)
library(MASS)
library(ggplot2)
library(gridExtra)
library(grid)

Rcpp::sourceCpp("lmm_fast.cpp")

# ============================================================
# 1. 辅助函数
# ============================================================

generate_groups <- function(R, m, N, V) {
  if (N <= R * m) stop("N must be greater than R * m")
  Vu <- if (V == "large") 300 * N / (5 * R) else N / (5 * R)
  adjusted_sum   <- N - R * m
  random_numbers <- abs(rnorm(R, mean = N / R, sd = Vu))
  random_numbers <- random_numbers / sum(random_numbers) * adjusted_sum
  result     <- ceiling(random_numbers + m)
  difference <- N - sum(result)
  while (difference != 0) {
    idx <- sample(1:R, 1)
    if (difference > 0) {
      result[idx] <- result[idx] + 1; difference <- difference - 1
    } else if (result[idx] - 1 > m) {
      result[idx] <- result[idx] - 1; difference <- difference + 1
    }
  }
  result
}

mbky <- function(setseed, FXX, y, Cn, threshold = 2) {
  current_Cn <- Cn
  repeat {
    if (current_Cn < 1) current_Cn <- 1
    set.seed(setseed)
    mbk <- ClusterR::MiniBatchKmeans(
      FXX, clusters = current_Cn, batch_size = 1024,
      num_init = 3, max_iters = 5, initializer = "kmeans++")
    batchs        <- ClusterR::predict_KMeans(FXX, mbk$centroids)
    cluster_sizes <- table(batchs)
    if (min(cluster_sizes) >= threshold || current_Cn == 1) break
    current_Cn <- current_Cn - 1
    setseed    <- setseed + 1
  }
  sort_idx <- order(batchs)
  list(
    R_CGOSS              = length(cluster_sizes),
    data_matrix_sorted   = FXX[sort_idx, , drop = FALSE],
    sorted_y             = y[sort_idx],
    cluster_sizes_vector = as.vector(table(batchs[sort_idx])),
    sorted_indices       = (1:nrow(FXX))[sort_idx],
    centroids            = mbk$centroids,
    final_Cn             = current_Cn
  )
}

calculate_K_dynamic <- function(FXX_test) {
  X <- cbind(1, FXX_test); (t(X) %*% X) / nrow(FXX_test)
}

# ============================================================
# 2. 数据生成
# ============================================================

generate_data <- function(seed = 1, N = 2500, p = 50, R = 5,
                          Var.e = 9, dist_x = "case3",
                          groupsize = "large", mode = "fixed") {
  set.seed(seed * 100000)
  beta  <- rep(1, p)
  sigma <- diag(0.5, p, p) + matrix(0.5, p, p)

  m        <- N / (10 * R)
  C.test   <- round(generate_groups(R, m, N, groupsize))
  C.train  <- 3 * C.test
  SC.test  <- c(0, cumsum(C.test))
  SC.train <- c(0, cumsum(C.train))

  FXX.train <- matrix(0, max(SC.train), p)
  FXX.test  <- matrix(0, max(SC.test),  p)

  sigma2 <- matrix(0, p, p)
  for (a in 1:p) for (b in 1:p) sigma2[a, b] <- 0.8^abs(a - b)

  for (i in 1:R) {
    ss     <- seed * 100000 + i * 100
    set.seed(ss)
    idx_tr <- (SC.train[i] + 1):SC.train[i + 1]
    idx_te <- (SC.test[i]  + 1):SC.test[i  + 1]
    if (dist_x == "case4") {
      mu <- rep(-2 + (i-1)/5, p)
      FXX.train[idx_tr,] <- mvrnorm(C.train[i], mu, sigma)
      FXX.test[idx_te, ] <- mvrnorm(C.test[i],  mu, sigma)
    } else if (dist_x == "case1") {
      FXX.train[idx_tr,] <- matrix(runif(C.train[i]*p,-1,1),C.train[i],p)
      FXX.test[idx_te, ] <- matrix(runif(C.test[i]*p, -1,1),C.test[i], p)
    } else if (dist_x == "case3") {
      lo <- -1.55+i/20; hi <- 0.45+i/20
      FXX.train[idx_tr,] <- matrix(runif(C.train[i]*p,lo,hi),C.train[i],p)
      FXX.test[idx_te, ] <- matrix(runif(C.test[i]*p, lo,hi),C.test[i], p)
    }
  }

  for (col in 1:p) {
    gmin <- min(FXX.train[,col], FXX.test[,col])
    gmax <- max(FXX.train[,col], FXX.test[,col])
    if (gmax > gmin) {
      FXX.train[,col] <- 2*(FXX.train[,col]-gmin)/(gmax-gmin)-1
      FXX.test[,col]  <- 2*(FXX.test[,col] -gmin)/(gmax-gmin)-1
    } else {
      FXX.train[,col] <- FXX.test[,col] <- 0
    }
  }

  if (mode == "fixed") {
    Var.a    <- 2.25
    Fa.train <- rep(rnorm(R, 0, sqrt(Var.a)), C.train)
    Fa.test  <- rep(rnorm(R, 0, sqrt(Var.a)), C.test)
    Fe.train <- rnorm(max(SC.train), 0, sqrt(Var.e))
    Fe.test  <- rnorm(max(SC.test),  0, sqrt(Var.e))
    FY.train <- as.matrix(1 + FXX.train %*% beta + Fa.train + Fe.train)
    FY.test  <- as.matrix(1 + FXX.test  %*% beta + Fa.test  + Fe.test)
  } else {
    sigma.a <- 2.25; sigma.b <- 1.0
    FY.train <- numeric(max(SC.train))
    FY.test  <- numeric(max(SC.test))
    for (i in 1:R) {
      idx_tr <- (SC.train[i]+1):SC.train[i+1]
      idx_te <- (SC.test[i] +1):SC.test[i +1]
      b.i    <- rnorm(p, 0, sqrt(sigma.b))
      FY.train[idx_tr] <- 1 + rnorm(1,0,sqrt(sigma.a)) +
                          FXX.train[idx_tr,] %*% (beta+b.i) +
                          rnorm(C.train[i], 0, sqrt(Var.e))
      FY.test[idx_te]  <- 1 + rnorm(1,0,sqrt(sigma.a)) +
                          FXX.test[idx_te,]  %*% (beta+b.i) +
                          rnorm(C.test[i],  0, sqrt(Var.e))
    }
    FY.train <- as.matrix(FY.train)
    FY.test  <- as.matrix(FY.test)
  }

  list(FXX.train=FXX.train, FXX.test=FXX.test,
       FY.train=FY.train,   FY.test=FY.test,
       C.train=C.train,     C.test=C.test,
       SC.train=SC.train,   SC.test=SC.test,
       R=R, p=p, mode=mode, beta=beta)
}

# ============================================================
# 3. SA 聚类
# ============================================================

run_SA <- function(dat, seed=1, Cn_init=2,
                   max_iter=50, T_init=50, alpha=0.95) {
  FXX.train <- dat$FXX.train; FY.train <- dat$FY.train
  FXX.test  <- dat$FXX.test;  p <- dat$p; mode <- dat$mode
  setseed   <- seed * 100000
  K_mat     <- calculate_K_dynamic(FXX.test)

  calc_obj <- function(FXX_s, FY_s, C_s, R_s) {
    if (mode == "fixed") {
      info <- count_info_cpp(FXX_s, FY_s, C_s, R_s, p)
      0.7*log(info$D)/p + 0.3*log(info$A/p)
    } else {
      info <- count_info_rs_cpp(FXX_s, FY_s, C_s, R_s, p)
      -sum(diag(solve(info$Information) %*% K_mat))
    }
  }

  Cn  <- Cn_init
  cl  <- mbky(setseed, FXX.train, FY.train, Cn)
  obj <- calc_obj(cl$data_matrix_sorted, cl$sorted_y,
                  cl$cluster_sizes_vector, cl$R_CGOSS)
  obj.best <- obj; best <- cl; T.curr <- T_init

  for (iter in seq_len(max_iter)) {
    if (T.curr < 1e-4) break
    step     <- sample(c(-1,1,2), 1)
    Cn.candi <- max(1, Cn + step)
    if (Cn.candi == Cn) Cn.candi <- Cn + 1
    cl.c  <- mbky(setseed+iter, FXX.train, FY.train, Cn.candi)
    obj.c <- calc_obj(cl.c$data_matrix_sorted, cl.c$sorted_y,
                      cl.c$cluster_sizes_vector, cl.c$R_CGOSS)
    delta  <- obj.c - obj
    accept <- (delta > 0) || (runif(1) < exp(delta/T.curr))
    if (accept) {
      Cn <- Cn.candi; obj <- obj.c; cl <- cl.c
      if (obj.c > obj.best) { obj.best <- obj.c; best <- cl.c }
    }
    T.curr <- T.curr * alpha
  }

  group_after <- integer(nrow(FXX.train))
  group_after[best$sorted_indices] <- rep(seq_len(best$R_CGOSS),
                                          best$cluster_sizes_vector)
  list(group_after=group_after, R_after=best$R_CGOSS)
}

# ============================================================
# 4. 计算每组的拟合线（Est_hat_RS_cpp + BLUP，投影到 PC1）
#    返回 data.frame：Group, x_start, x_end, y_start, y_end
# ============================================================

compute_lines <- function(groups, FXX.train, FY.train,
                          pc1, Y, pca_rotation, beta_true, Var.e, p) {
  v1       <- pca_rotation[, 1]
  ord      <- order(groups)
  FXX_sort <- FXX.train[ord, , drop=FALSE]
  FY_sort  <- FY.train[ord]
  C_grp    <- as.integer(table(groups))
  R_grp    <- length(C_grp)

  est       <- Est_hat_RS_cpp(FXX_sort, FY_sort,
                              beta_true, 2.25, Var.e, C_grp, R_grp, p)
  beta_hat  <- est$beta2
  G.hat     <- est$G.hat
  Var.e_hat <- est$Var.e
  G_inv     <- solve(G.hat + diag(1e-7, p+1))
  X_full    <- cbind(1, FXX_sort)
  SC        <- c(0, cumsum(C_grp))

  do.call(rbind, lapply(seq_len(R_grp), function(i) {
    g       <- sort(unique(groups))[i]
    idx_grp <- (SC[i]+1):SC[i+1]
    idx_ori <- which(groups == g)
    Xg  <- X_full[idx_grp, , drop=FALSE]
    yg  <- FY_sort[idx_grp]
    M   <- solve(Var.e_hat * G_inv + t(Xg) %*% Xg)
    u_i <- as.numeric(M %*% (t(Xg) %*% (yg - Xg %*% beta_hat)))
    slope     <- as.numeric(v1 %*% (beta_hat[-1] + u_i[-1]))
    intercept <- mean(Y[idx_ori]) - slope * mean(pc1[idx_ori])
    x0 <- min(pc1[idx_ori]); x1 <- max(pc1[idx_ori])
    data.frame(Group=factor(g),
               x_start=x0, x_end=x1,
               y_start=intercept + slope*x0,
               y_end  =intercept + slope*x1)
  }))
}

# ============================================================
# 5. 单个子图：散点 + 拟合线
# ============================================================

make_panel <- function(pc1, Y, groups, lines_df, subtitle, y_range) {
  R_num <- length(unique(groups))
  pal   <- hcl.colors(R_num, palette = "Dark 2")
  df    <- data.frame(PC1=pc1, Y=Y, Group=factor(groups))
  lines_df$Group <- factor(lines_df$Group, levels=levels(df$Group))

  ggplot(df, aes(x=PC1, y=Y, color=Group)) +
    geom_point(alpha=0.2, size=0.5) +
    geom_segment(
      data=lines_df,
      aes(x=x_start, xend=x_end, y=y_start, yend=y_end, color=Group),
      linewidth=1.3, inherit.aes=FALSE
    ) +
    scale_color_manual(values=pal, name="Group") +
    scale_y_continuous(limits=y_range) +
    labs(subtitle=subtitle, x="PC1", y="Y") +
    theme_bw(base_size=10) +
    theme(plot.subtitle    =element_text(hjust=0.5, size=10, face="bold"),
          legend.position  ="right",
          legend.key.size  =unit(0.4,"cm"),
          legend.text      =element_text(size=8),
          legend.title     =element_text(size=9, face="bold"),
          panel.grid.minor =element_blank(),
          axis.title       =element_text(size=9)) +
    guides(color=guide_legend(
      override.aes=list(linewidth=1.2, size=2.5, alpha=1)))
}

# ============================================================
# 6. 主流程
# ============================================================

main <- function(seed=1, N=2500, R=5, Var.e=9,
                 dist_x="case3", groupsize="large",
                 p_vals=c(2, 10, 50), out_dir=".") {

  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

  for (mode in c("fixed", "rs")) {
    lbl <- if (mode=="fixed") "No Random Slope" else "Random Slope"
    cat(sprintf("\n===== %s =====\n", lbl))

    all_panels <- list()

    for (p in p_vals) {
      cat(sprintf("  p = %d ... ", p))
      beta_true <- rep(1, p)

      # 生成数据
      dat <- generate_data(seed=seed, N=N, p=p, R=R, Var.e=Var.e,
                           dist_x=dist_x, groupsize=groupsize, mode=mode)

      # PCA
      pca_res      <- prcomp(dat$FXX.train, center=TRUE, scale.=FALSE)
      pc1          <- pca_res$x[,1]
      pca_rotation <- pca_res$rotation
      Y            <- as.numeric(dat$FY.train)

      # 分组标签
      grp_before <- rep(seq_len(R), dat$C.train)
      sa         <- run_SA(dat, seed=seed)
      grp_after  <- sa$group_after
      R_after    <- sa$R_after
      cat(sprintf("R.best=%d\n", R_after))

      # 共享 Y 轴范围（before/after 用同一范围方便对比）
      y_range <- range(Y) + c(-1,1) * diff(range(Y)) * 0.05

      # 计算拟合线
      lines_bef <- compute_lines(grp_before, dat$FXX.train, dat$FY.train,
                                 pc1, Y, pca_rotation, beta_true, Var.e, p)
      lines_aft <- compute_lines(grp_after,  dat$FXX.train, dat$FY.train,
                                 pc1, Y, pca_rotation, beta_true, Var.e, p)

      # 两个子图
      panel_bef <- make_panel(pc1, Y, grp_before, lines_bef,
                              sprintf("p=%d  |  Before Regouping (%d groups)", p, R),
                              y_range)
      panel_aft <- make_panel(pc1, Y, grp_after,  lines_aft,
                              sprintf("p=%d  |  After Regouping (%d groups)", p, R_after),
                              y_range)

      all_panels[[paste0(p,"_bef")]] <- panel_bef
      all_panels[[paste0(p,"_aft")]] <- panel_aft
    }

    # 排列：行1=Before, 行2=After; 列=p=2,10,50
    grob_list <- c(
      lapply(p_vals, function(p) all_panels[[paste0(p,"_bef")]]),
      lapply(p_vals, function(p) all_panels[[paste0(p,"_aft")]])
    )

    combined <- arrangeGrob(
      grobs = grob_list,
      nrow  = 2,
      ncol  = length(p_vals),
      top   = textGrob(
        sprintf("Effect of re-clustering on fitted lines per group  [%s]", lbl),
        gp = gpar(fontsize=14, fontface="bold"))
    )

    out_path <- file.path(out_dir, paste0("plot_lines_", mode, ".png"))
    png(out_path, width=3600, height=2200, res=150)
    grid.draw(combined)
    dev.off()
    cat("  Saved:", out_path, "\n")
  }

  cat("\n=== Done! 2 PNG files saved to:", out_dir, "===\n")
}

# ---- 运行入口 ----
main(
  seed      = 1,
  N         = 2500,
  R         = 5,
  Var.e     = 9,
  dist_x    = "case3",
  groupsize = "large",
  p_vals    = c(2, 10, 50),
  out_dir   = "."
)
