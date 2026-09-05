# ROG_cases10_11_observation.R
# -----------------------------------------------------------------------------
# Case 10 and Case 11 rebuilt so that EVERY grouping method acts on observations.
# There is no unit-level clustering object in either scenario.
# -----------------------------------------------------------------------------

if (!exists("run_observation_comparison", mode = "function")) {
  if (!file.exists("ROG_observation_comparison.R")) stop("ROG_observation_comparison.R not found.")
  source("ROG_observation_comparison.R")
}

ar_cov_obs <- function(p, rho = 0.2, variance = 0.6^2) {
  idx <- seq_len(p)
  variance * outer(idx, idx, function(i, j) rho^abs(i - j))
}

draw_mvn_obs <- function(n, mu, Sigma) {
  Z <- matrix(rnorm(n * length(mu)), n, length(mu))
  ev <- eigen((Sigma + t(Sigma))/2, symmetric = TRUE)
  ev$values[ev$values < 1e-10] <- 1e-10
  L <- ev$vectors %*% diag(sqrt(ev$values), length(mu))
  sweep(Z %*% t(L), 2, mu, "+")
}

balanced_labels_obs <- function(n, K, seed) {
  lab <- rep(seq_len(K), length.out = n)
  set.seed(seed); sample(lab, length(lab), replace = FALSE)
}

continuous_group_r2_obs <- function(z, labels) {
  if (is.null(z) || is.null(labels)) return(NA_real_)
  labels <- as.integer(factor(labels)); z <- as.numeric(z)
  tot <- sum((z - mean(z))^2)
  if (!is.finite(tot) || tot < 1e-12) return(NA_real_)
  zh <- ave(z, labels, FUN = mean)
  1 - sum((z - zh)^2) / tot
}

pseudo_oracle_fit_obs <- function(beta, var_a, var_b, var_e, model_type) {
  list(beta = c(1, beta), Var.a = var_a,
       Var.b = if (model_type == "RS") var_b else NA_real_,
       Var.e = var_e,
       G = if (model_type == "RS") diag(c(var_a, rep(var_b, length(beta)))) else NULL,
       u_hat = if (model_type == "RS") matrix(0, 0, length(beta) + 1L) else numeric(0),
       converged = TRUE, iterations = 0L)
}

summarize_obs_structured <- function(raw) {
  num <- c("MSPE", "beta_MSE", "beta_SSE", "intercept_MSE",
           "Var_a_hat", "Var_b_hat", "Var_e_hat", "Var_a_MSE",
           "Var_b_MSE", "Var_e_MSE", "G_MSE", "selected_K", "runtime",
           "iterations", "group_ARI", "z_R2")
  parts <- lapply(split(raw, raw$method), function(d) {
    out <- data.frame(scenario = d$scenario[1], method = d$method[1],
                      n = nrow(d), convergence_rate = mean(d$converged))
    for (nm in num) {
      x <- d[[nm]]
      out[[paste0(nm, "_mean")]] <- if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
      out[[paste0(nm, "_MCSE")]] <- if (sum(is.finite(x)) <= 1L) NA_real_ else stats::sd(x, na.rm = TRUE)/sqrt(sum(is.finite(x)))
    }
    out
  })
  ans <- do.call(rbind, parts); rownames(ans) <- NULL; ans
}

# ------------------------ common observation evaluator -----------------------

evaluate_observation_dataset <- function(dat,
                                         methods = c("ORACLE","OBS","LM","KM","GMM","CPF","MIXREG","CIRG"),
                                         tau = 5 * (ncol(dat$X_train) + 1L),
                                         K_grid = 2:16,
                                         initial_Cn = 4L,
                                         sa_max_iter = 80,
                                         em_tol = 1e-5,
                                         em_max_iter = 500,
                                         seed = 12345L,
                                         cpf_graph_k = 8L,
                                         verbose = FALSE) {
  methods <- unique(toupper(methods)); rows <- list(); extra <- list()
  Xtr <- dat$X_train; Xrf <- dat$X_ref; Xte <- dat$X_test
  ytr <- dat$y_train; yte <- dat$y_test; beta <- dat$beta
  model_type <- dat$model_type; p <- ncol(Xtr)
  va <- dat$var_a_metric; vb <- dat$var_b_metric; ve <- dat$var_e
  true_te <- dat$true_partition_test
  G_true <- if (!is.null(dat$G_true)) dat$G_true else NULL
  min_group <- if (model_type == "RS") p + 1L else 2L
  K_grid <- sort(unique(as.integer(K_grid)))
  K_grid <- K_grid[K_grid >= 2L & K_grid <= floor(nrow(Xtr)/min_group)]
  add <- function(d, labels = NULL) {
    d$z_R2 <- continuous_group_r2_obs(dat$z_test, labels)
    rows[[length(rows)+1L]] <<- d
  }

  if ("ORACLE" %in% methods) {
    t0 <- proc.time()[3]
    if (!is.null(dat$oracle_pred)) {
      fit <- pseudo_oracle_fit_obs(beta, va, vb, ve, model_type)
      d <- metric_row("ORACLE", fit, dat$oracle_pred, yte, beta, va, vb, ve,
                      G_true, dat$oracle_K, proc.time()[3]-t0, NA_real_)
      add(d, NULL)
    } else {
      fp <- fit_predict_labeled(Xtr, ytr, Xte, dat$oracle_train, dat$oracle_test,
                                model_type, em_tol, em_max_iter)
      ari <- if (!is.null(true_te)) adjusted_rand_index(fp$labels_test, true_te) else NA_real_
      d <- metric_row("ORACLE", fp$fit, fp$pred, yte, beta, va, vb, ve,
                      G_true, fp$K, proc.time()[3]-t0, ari)
      add(d, fp$labels_test)
    }
  }

  if ("OBS" %in% methods) {
    t0 <- proc.time()[3]
    fp <- fit_predict_labeled(Xtr, ytr, Xte, dat$observed_train, dat$observed_test,
                              model_type, em_tol, em_max_iter)
    ari <- if (!is.null(true_te)) adjusted_rand_index(fp$labels_test, true_te) else NA_real_
    d <- metric_row("OBS", fp$fit, fp$pred, yte, beta, va, vb, ve,
                    G_true, fp$K, proc.time()[3]-t0, ari)
    add(d, fp$labels_test)
  }

  if ("LM" %in% methods) {
    t0 <- proc.time()[3]; fit <- fit_pooled_lm(Xtr, ytr); pred <- predict_pooled_lm(fit, Xte)
    d <- metric_row("LM", fit, pred, yte, beta, va, vb, ve, G_true, 1L, proc.time()[3]-t0, NA_real_)
    add(d, rep(1L, nrow(Xte)))
  }

  if ("KM" %in% methods) {
    t0 <- proc.time()[3]; km <- fit_kmeans_bic(Xtr, K_grid, seed+401L, tau=tau)
    trlab <- km$cluster$labels; telab <- assign_kmeans(Xte, km$cluster$centroids)
    fp <- fit_predict_labeled(Xtr,ytr,Xte,trlab,telab,model_type,em_tol,em_max_iter)
    ari <- if (!is.null(true_te)) adjusted_rand_index(telab,true_te) else NA_real_
    d <- metric_row("KM",fp$fit,fp$pred,yte,beta,va,vb,ve,G_true,fp$K,proc.time()[3]-t0,ari)
    add(d,telab); extra$KM <- km
  }

  if ("GMM" %in% methods && requireNamespace("mclust", quietly=TRUE)) {
    t0 <- proc.time()[3]; gm <- fit_gmm_bic(Xtr,K_grid,seed+503L)
    trlab <- gm$labels; telab <- predict_gmm_labels(gm,Xte)
    if (model_type == "RI" || min(tabulate(trlab,max(trlab))) >= p+1L) {
      fp <- fit_predict_labeled(Xtr,ytr,Xte,trlab,telab,model_type,em_tol,em_max_iter)
      ari <- if (!is.null(true_te)) adjusted_rand_index(telab,true_te) else NA_real_
      d <- metric_row("GMM",fp$fit,fp$pred,yte,beta,va,vb,ve,G_true,fp$K,proc.time()[3]-t0,ari)
      add(d,telab)
    }
    extra$GMM <- gm
  }

  if ("CPF" %in% methods) {
    t0 <- proc.time()[3]
    cpf <- fit_cpf_observation(Xtr,ytr,min_group_size=min_group,graph_k=cpf_graph_k)
    clf <- fit_x_label_classifier(Xtr,cpf$labels,seed+613L)
    telab <- predict_x_label_classifier(clf,Xte)
    fp <- fit_predict_labeled(Xtr,ytr,Xte,cpf$labels,telab,model_type,em_tol,em_max_iter)
    ari <- if (!is.null(true_te)) adjusted_rand_index(telab,true_te) else NA_real_
    d <- metric_row("CPF",fp$fit,fp$pred,yte,beta,va,vb,ve,G_true,fp$K,proc.time()[3]-t0,ari)
    add(d,telab); extra$CPF <- list(fit=cpf,classifier=clf)
  }

  if ("MIXREG" %in% methods) {
    t0 <- proc.time()[3]
    kk <- K_grid[K_grid <= min(10L,max(K_grid))]
    mr <- fit_mixreg_bic(Xtr,ytr,kk,min_group_size=min_group,seed=seed+641L)
    clf <- fit_x_label_classifier(Xtr,mr$labels,seed+647L)
    telab <- predict_x_label_classifier(clf,Xte)
    fp <- fit_predict_labeled(Xtr,ytr,Xte,mr$labels,telab,model_type,em_tol,em_max_iter)
    ari <- if (!is.null(true_te)) adjusted_rand_index(telab,true_te) else NA_real_
    d <- metric_row("MIXREG",fp$fit,fp$pred,yte,beta,va,vb,ve,G_true,fp$K,proc.time()[3]-t0,ari)
    add(d,telab); extra$MIXREG <- list(fit=mr,classifier=clf)
  }

  if ("CIRG" %in% methods) {
    t0 <- proc.time()[3]
    sr <- cirg_search(Xtr,ytr,Xrf,initial_Cn=initial_Cn,tau=tau,
                      model_type=model_type,max_iter=sa_max_iter,seed=seed+701L,
                      em_tol=em_tol,em_max_iter=em_max_iter,verbose=verbose)
    best <- sr$best; fit <- best$imspe_fit
    pred <- if (model_type=="RI") predict_ri_soft(fit,Xte,best$cluster$params) else predict_rs_soft(fit,Xte,best$cluster$params)
    telab <- predict_soft_labels(Xte,best$cluster$params)
    ari <- if (!is.null(true_te)) adjusted_rand_index(telab,true_te) else NA_real_
    d <- metric_row("CIRG",fit,pred,yte,beta,va,vb,ve,G_true,best$K,proc.time()[3]-t0,ari,best$objective)
    add(d,telab); extra$CIRG <- sr
  }

  raw <- do.call(rbind,rows); raw$scenario <- dat$scenario
  list(metrics=raw,extra=extra)
}

# ------------------------------ CASE 10 --------------------------------------
# Observation-level latent predictive groups. There are K_type true random-
# effect groups but K_type*modes_per_type geometric X modes. OBS uses the mode
# labels, so it is deliberately over-refined. No unit ID is used anywhere.

generate_case10_observation <- function(p=50,N_test=2500,K_type=6,
                                        modes_per_type=2,train_multiplier=3L,
                                        Var.a=2.25,Var.b=0,Var.e=9,
                                        type_sep=1.35,mode_sep=1.55,
                                        within_sd=0.65,seed=12345L) {
  set.seed(seed); beta <- rep(1,p); model_type <- if (Var.b>0) "RS" else "RI"
  K_mode <- K_type*modes_per_type
  # Fixed random directions define type signal and within-type nuisance modes.
  set.seed(seed+100L)
  Tdir <- matrix(rnorm(K_type*p),K_type,p)
  Tdir <- Tdir/sqrt(rowSums(Tdir^2))
  Mdir <- matrix(rnorm(K_type*p),K_type,p)
  for(k in seq_len(K_type)) {
    Mdir[k,] <- Mdir[k,] - sum(Mdir[k,]*Tdir[k,])*Tdir[k,]
    Mdir[k,] <- Mdir[k,]/sqrt(sum(Mdir[k,]^2))
  }
  offs <- seq(-mode_sep,mode_sep,length.out=modes_per_type)
  centers <- matrix(0,K_mode,p)
  for(g in seq_len(K_type)) for(s in seq_len(modes_per_type)) {
    m <- (g-1L)*modes_per_type+s
    centers[m,] <- type_sep*Tdir[g,] + offs[s]*Mdir[g,]
  }
  Sigma <- ar_cov_obs(p,0.15,within_sd^2)

  make_split <- function(n,offset) {
    type <- balanced_labels_obs(n,K_type,seed+offset)
    subtype <- integer(n)
    for(g in seq_len(K_type)) {
      ii <- which(type==g)
      subtype[ii] <- balanced_labels_obs(length(ii),modes_per_type,seed+offset+1000L*g)
    }
    mode <- (type-1L)*modes_per_type+subtype
    X <- matrix(0,n,p)
    for(m in seq_len(K_mode)) {
      ii <- which(mode==m); if(length(ii)) {
        set.seed(seed+offset+10007L*m); X[ii,] <- draw_mvn_obs(length(ii),centers[m,],Sigma)
      }
    }
    list(X=X,type=type,mode=mode)
  }
  tr <- make_split(as.integer(train_multiplier*N_test),10000L)
  rf <- make_split(N_test,20000L); te <- make_split(N_test,30000L)

  set.seed(seed+40000L); a <- rnorm(K_type,sd=sqrt(Var.a))
  B <- if(Var.b>0) matrix(rnorm(K_type*p,sd=sqrt(Var.b)),K_type,p) else matrix(0,K_type,p)
  make_y <- function(obj,offset) {
    set.seed(seed+offset)
    re <- a[obj$type] + rowSums(obj$X*B[obj$type,,drop=FALSE])
    1 + drop(obj$X%*%beta) + re + rnorm(nrow(obj$X),sd=sqrt(Var.e))
  }
  ytr <- make_y(tr,50000L); yte <- make_y(te,60000L)
  list(scenario="case10_observation_predictive_vs_geometric",
       X_train=tr$X,X_ref=rf$X,X_test=te$X,y_train=ytr,y_test=yte,beta=beta,
       model_type=model_type,var_a_metric=Var.a,var_b_metric=Var.b,var_e=Var.e,
       G_true=if(Var.b>0) diag(c(Var.a,rep(Var.b,p))) else NULL,
       oracle_train=tr$type,oracle_test=te$type,oracle_K=K_type,
       true_partition_test=te$type,
       observed_train=tr$mode,observed_test=te$mode,
       z_test=NULL,K_type=K_type,K_mode=K_mode,centers=centers,a_true=a,B_true=B)
}

run_case10_observation <- function(nloop=50,p=50,N_test=2500,K_type=6,
                                   modes_per_type=2,Var.a=2.25,Var.b=0,Var.e=9,
                                   tau=5*(p+1L),K_grid=2:16,initial_Cn=4L,
                                   sa_max_iter=80,em_tol=1e-5,em_max_iter=500,
                                   seed=12345L,
                                   methods=c("ORACLE","OBS","LM","KM","GMM","CPF","MIXREG","CIRG"),
                                   keep_extra=FALSE,verbose=FALSE) {
  raws <- list(); extras <- if(keep_extra) list() else NULL
  for(r in seq_len(nloop)) {
    s <- as.integer(seed+1000003L*r)
    dat <- generate_case10_observation(p,N_test,K_type,modes_per_type,3L,Var.a,Var.b,Var.e,seed=s)
    ans <- evaluate_observation_dataset(dat,methods,tau,K_grid,initial_Cn,sa_max_iter,em_tol,em_max_iter,s+77L,verbose=verbose)
    d <- ans$metrics; d$replication <- r; raws[[r]] <- d; if(keep_extra) extras[[r]] <- ans$extra
  }
  raw <- do.call(rbind,raws); list(raw=raw,summary=summarize_obs_structured(raw),extra=extras)
}

# ------------------------------ CASE 11 --------------------------------------
# Continuous observation-level heterogeneity. Every observation has a latent z;
# there is no true discrete K. Grouped LMMs approximate the smooth latent shift.

generate_case11_observation <- function(p=50,N_test=2500,train_multiplier=3L,
                                        obs_bins=4L,proxy_sd=0.35,
                                        Var.a=2.25,Var.b=0,Var.e=9,
                                        within_sd=0.70,seed=12345L) {
  set.seed(seed); beta <- rep(1,p); model_type <- if(Var.b>0) "RS" else "RI"
  # Scale smooth a(z) to the requested variance under z~U[-1,1].
  grid <- seq(-1,1,length.out=20001)
  f0 <- 1.15*grid + 0.75*sin(pi*grid)
  scale_a <- sqrt(Var.a/stats::var(f0))
  a_fun <- function(z) scale_a*(1.15*z+0.75*sin(pi*z))

  # Smooth slope heterogeneity. Scale average coordinate variance to Var.b.
  b_raw <- function(z) {
    if(Var.b<=0) return(matrix(0,length(z),p))
    J <- seq_len(p)
    outer(z,J,function(zz,j) 0.65*zz*cos(2*pi*j/p)+0.35*sin(pi*zz+2*pi*j/p))
  }
  if(Var.b>0) {
    Bg <- b_raw(grid[seq(1,length(grid),by=20)])
    sb <- sqrt(Var.b/mean(apply(Bg,2,var)))
  } else sb <- 0

  mu_fun <- function(z) {
    n <- length(z); M <- matrix(0,n,p)
    q1 <- seq_len(min(10,p)); M[,q1] <- outer(z,seq(1.3,0.5,length.out=length(q1)))
    if(p>10) {q2 <- 11:min(20,p); M[,q2] <- outer(z^2-1/3,seq(1.0,0.4,length.out=length(q2)))}
    if(p>20) {q3 <- 21:min(30,p); M[,q3] <- outer(sin(pi*z),seq(0.9,0.35,length.out=length(q3)))}
    M
  }
  Sigma <- ar_cov_obs(p,0.15,within_sd^2)
  make_split <- function(n,offset) {
    set.seed(seed+offset); z <- runif(n,-1,1); Mu <- mu_fun(z)
    X <- matrix(0,n,p)
    # Same covariance for all rows; draw centered noise once for efficiency.
    noise <- draw_mvn_obs(n,rep(0,p),Sigma); X <- Mu+noise
    list(X=X,z=z)
  }
  tr <- make_split(as.integer(train_multiplier*N_test),10000L)
  rf <- make_split(N_test,20000L); te <- make_split(N_test,30000L)
  effect <- function(obj) {
    a <- a_fun(obj$z); B <- if(Var.b>0) sb*b_raw(obj$z) else matrix(0,length(obj$z),p)
    a + rowSums(obj$X*B)
  }
  re_tr <- effect(tr); re_te <- effect(te)
  set.seed(seed+40000L); ytr <- 1+drop(tr$X%*%beta)+re_tr+rnorm(nrow(tr$X),sd=sqrt(Var.e))
  set.seed(seed+50000L); yte <- 1+drop(te$X%*%beta)+re_te+rnorm(nrow(te$X),sd=sqrt(Var.e))

  # Observation-level coarse labels based on a noisy observable proxy for z.
  set.seed(seed+60000L); ptr <- tr$z+rnorm(length(tr$z),sd=proxy_sd)
  set.seed(seed+70000L); pte <- te$z+rnorm(length(te$z),sd=proxy_sd)
  cuts <- as.numeric(stats::quantile(ptr,probs=seq(0,1,length.out=obs_bins+1L),type=8))
  cuts[1] <- -Inf; cuts[length(cuts)] <- Inf
  obs_tr <- as.integer(cut(ptr,breaks=cuts,include.lowest=TRUE,labels=FALSE))
  obs_te <- as.integer(cut(pte,breaks=cuts,include.lowest=TRUE,labels=FALSE))
  va_emp <- stats::var(a_fun(c(tr$z,te$z)))
  vb_emp <- if(Var.b>0) {
    Bt <- sb*b_raw(c(tr$z,te$z)); mean(apply(Bt,2,var))
  } else 0
  oracle_pred <- 1+drop(te$X%*%beta)+re_te
  list(scenario="case11_observation_continuous_heterogeneity",
       X_train=tr$X,X_ref=rf$X,X_test=te$X,y_train=ytr,y_test=yte,beta=beta,
       model_type=model_type,var_a_metric=va_emp,var_b_metric=vb_emp,var_e=Var.e,
       G_true=NULL,oracle_pred=oracle_pred,oracle_K=NA_integer_,
       true_partition_test=NULL,observed_train=obs_tr,observed_test=obs_te,
       z_train=tr$z,z_test=te$z)
}

run_case11_observation <- function(nloop=50,p=50,N_test=2500,obs_bins=4L,
                                   proxy_sd=0.35,Var.a=2.25,Var.b=0,Var.e=9,
                                   tau=5*(p+1L),K_grid=2:15,initial_Cn=4L,
                                   sa_max_iter=80,em_tol=1e-5,em_max_iter=500,
                                   seed=12345L,
                                   methods=c("ORACLE","OBS","LM","KM","GMM","CPF","MIXREG","CIRG"),
                                   keep_extra=FALSE,verbose=FALSE) {
  raws <- list(); extras <- if(keep_extra) list() else NULL
  for(r in seq_len(nloop)) {
    s <- as.integer(seed+1000003L*r)
    dat <- generate_case11_observation(p,N_test,3L,obs_bins,proxy_sd,Var.a,Var.b,Var.e,seed=s)
    ans <- evaluate_observation_dataset(dat,methods,tau,K_grid,initial_Cn,sa_max_iter,em_tol,em_max_iter,s+91L,verbose=verbose)
    d <- ans$metrics; d$replication <- r; raws[[r]] <- d; if(keep_extra) extras[[r]] <- ans$extra
  }
  raw <- do.call(rbind,raws); list(raw=raw,summary=summarize_obs_structured(raw),extra=extras)
}

# No automatic execution.
