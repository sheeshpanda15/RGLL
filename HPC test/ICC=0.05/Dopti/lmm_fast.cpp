// lmm_fast.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;
using arma::vec;
using arma::uword;

// --------- helper: scalex each column to [-1, 1] -------------
static inline mat scalex_mat(const mat& X) {
  mat Z = X;
  const uword n = Z.n_rows, p = Z.n_cols;
  for (uword j = 0; j < p; ++j) {
    double mn = Z.col(j).min();
    double mx = Z.col(j).max();
    double den = mx - mn;
    if (den == 0.0) {
      Z.col(j).zeros();               // consistent safe handling
    } else {
      Z.col(j) = 2.0 * (Z.col(j) - mn) / den - 1.0;
    }
  }
  return Z;
}

// --------- helper: OLS beta for y ~ 1 + X -------------------
static inline vec ols_beta(const mat& X, const vec& y) {
  // X is n x p (NO intercept). Build design with intercept.
  mat D(X.n_rows, X.n_cols + 1, arma::fill::ones);
  D.cols(1, X.n_cols) = X;
  // Solve least squares
  return arma::solve(D, y);
}

// [[Rcpp::export]]
Rcpp::List find_sigma_cpp(const arma::mat& xx,
                          const arma::vec& yy,
                          const arma::vec& beta,
                          Rcpp::IntegerVector nc,
                          int R) {
  // eta = y - [1, X] beta
  vec eta = yy - (beta[0] + xx * beta.subvec(1, beta.n_elem - 1));

  const double N  = static_cast<double>(eta.n_elem);
  const double s1 = arma::accu(eta);
  const double s2 = arma::dot(eta, eta);

  // Ue = N*sum(eta^2) - (sum eta)^2
  const double Ue = N * s2 - s1 * s1;

  // drop nc==0
  std::vector<int> ncz;
  ncz.reserve(nc.size());
  for (int v : nc) if (v != 0) ncz.push_back(v);

  const int R1 = static_cast<int>(ncz.size());
  double Ua = 0.0;
  int offset = 0;
  for (int g = 0; g < R1; ++g) {
    const int m = ncz[g];
    vec eg = eta.subvec(offset, offset + m - 1);
    const double sg1 = arma::accu(eg);
    const double sg2 = arma::dot(eg, eg);
    // Ua group contribution: sum(eg^2) - (sum eg)^2 / m
    Ua += sg2 - (sg1 * sg1) / static_cast<double>(m);
    offset += m;
  }

  const double NF = static_cast<double>(offset); // sum of nonzero nc
  double NS = 0.0;
  for (int m : ncz) NS += static_cast<double>(m) * static_cast<double>(m);

  // Var.e
  const double denom_e = (NF - static_cast<double>(R));
  double Var_e = Ua / denom_e;

  // Var.a
  const double denom_a = (NF * NF - NS);
  double Var_a = 0.0;
  if (denom_a != 0.0) {
    Var_a = Ue / denom_a - Ua * (NF * NF - NF) / (denom_a * denom_e);
  }
  if (Var_a < 0.0) Var_a = 0.0;

  return Rcpp::List::create(
    Rcpp::Named("Var.a") = Var_a,
    Rcpp::Named("Var.e") = Var_e
  );
}

// [[Rcpp::export]]
arma::vec find_beta_cpp(const arma::mat& xx,
                        const arma::vec& yy,
                        double Var_a,
                        double Var_e,
                        Rcpp::IntegerVector nc,
                        int R,
                        int p) {
  // drop nc==0
  std::vector<int> ncz;
  ncz.reserve(nc.size());
  for (int v : nc) if (v != 0) ncz.push_back(v);

  const int R1 = static_cast<int>(ncz.size());

  // XVX, XVY
  mat XVX(p + 1, p + 1, arma::fill::zeros);
  vec XVY(p + 1, arma::fill::zeros);

  int offset = 0;
  for (int g = 0; g < R1; ++g) {
    const int m = ncz[g];
    mat Xg(m, p + 1, arma::fill::ones);
    Xg.cols(1, p) = xx.rows(offset, offset + m - 1);
    vec yg = yy.subvec(offset, offset + m - 1);

    // gamma = (m*Var_a)/(Var_e + m*Var_a)
    const double denom = (Var_e + static_cast<double>(m) * Var_a);
    const double gamma = (denom == 0.0) ? 0.0 : (static_cast<double>(m) * Var_a) / denom;

    // Apply invV = (1/Var_e) * (I - (gamma/m) 11^T) WITHOUT building it:
    // invV * X = (1/Var_e) * (X - (gamma/m) * 1 * (1^T X))
    // invV * y = (1/Var_e) * (y - (gamma/m) * 1 * (1^T y))

    const double invVe = 1.0 / Var_e;
    vec ones(m, arma::fill::ones);

    arma::rowvec sX = arma::sum(Xg, 0);     // 1 x (p+1)
    const double sy = arma::accu(yg);       // scalar

    mat invV_X = invVe * (Xg - (gamma / static_cast<double>(m)) * (ones * sX));
    vec invV_y = invVe * (yg - (gamma / static_cast<double>(m)) * (ones * sy));

    XVX += Xg.t() * invV_X;
    XVY += Xg.t() * invV_y;

    offset += m;
  }

  // Solve XVX * beta = XVY
  // (If singular, arma::solve will throw; that is usually desirable to surface numerical issues.)
  vec bt = arma::solve(XVX, XVY);
  return bt;
}



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;
using arma::vec;

// --------- 内部辅助：GLS 更新 beta (全随机系数版) -------------
static inline vec find_beta_full_rs_internal(const mat& xx, const vec& yy, 
                                            const mat& G, double Ve, 
                                            Rcpp::IntegerVector nc, int p) {
    int q = p + 1;
    mat XVX(q, q, arma::fill::zeros);
    vec XVY(q, arma::fill::zeros);
    int offset = 0;
    mat G_inv = arma::inv(G);

    for (int g = 0; g < nc.size(); ++g) {
        int m = nc[g];
        if (m <= 0) continue;
        mat Xg(m, q, arma::fill::ones);
        Xg.cols(1, p) = xx.rows(offset, offset + m - 1);
        vec yg = yy.subvec(offset, offset + m - 1);

        mat M_inv = arma::inv(Ve * G_inv + Xg.t() * Xg);
        mat invV_X = (1.0 / Ve) * (Xg - Xg * (M_inv * (Xg.t() * Xg)));
        vec invV_y = (1.0 / Ve) * (yg - Xg * (M_inv * (Xg.t() * yg)));

        XVX += Xg.t() * invV_X;
        XVY += Xg.t() * invV_y;
        offset += m;
    }
    return arma::solve(XVX, XVY);
}

// --------- 内部辅助：方差分量更新 (全随机系数版) -------------
static inline Rcpp::List find_sigma_full_rs_internal(const mat& xx, const vec& yy, 
                                                   const vec& beta, const mat& G_old, 
                                                   double Ve_old, Rcpp::IntegerVector nc, int p) {
    int q = p + 1;
    mat G_new(q, q, arma::fill::zeros);
    double Ve_new = 0.0;
    int offset = 0;
    mat G_inv_old = arma::inv(G_old);

    for (int g = 0; g < nc.size(); ++g) {
        int m = nc[g];
        if (m <= 0) continue;
        mat Xg(m, q, arma::fill::ones);
        Xg.cols(1, p) = xx.rows(offset, offset + m - 1);
        vec yg = yy.subvec(offset, offset + m - 1);
        vec res_g = yg - Xg * beta;

        mat M_inv = arma::inv(Ve_old * G_inv_old + Xg.t() * Xg);
        vec Vg_inv_res = (1.0 / Ve_old) * (res_g - Xg * (M_inv * (Xg.t() * res_g)));
        vec ug = G_old * Xg.t() * Vg_inv_res;
        
        G_new += ug * ug.t(); 
        Ve_new += arma::dot(res_g - Xg * ug, res_g - Xg * ug);
        offset += m;
    }
    return Rcpp::List::create(Rcpp::Named("G") = G_new / static_cast<double>(nc.size()),
                             Rcpp::Named("Var.e") = Ve_new / static_cast<double>(yy.n_elem));
}


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using arma::mat;
using arma::vec;

// ======== helper 1: OLS beta for y ~ 1 + X ========
static inline arma::vec rs_count_ols_beta(const arma::mat& X, const arma::vec& y) {
    arma::mat D(X.n_rows, X.n_cols + 1, arma::fill::ones);
    D.cols(1, X.n_cols) = X;
    return arma::solve(D, y);
}

// ======== helper 2: GLS update beta under full random slopes ========
// Model: y_g = X_g beta + X_g u_g + e_g
// u_g ~ N(0, G), e_g ~ N(0, Ve I)
// Here X_g includes intercept column.
static inline arma::vec rs_count_find_beta_internal(const arma::mat& xx,
                                                    const arma::vec& yy,
                                                    const arma::mat& G,
                                                    double Ve,
                                                    Rcpp::IntegerVector nc,
                                                    int p) {
    const int q = p + 1;

    arma::mat XVX(q, q, arma::fill::zeros);
    arma::vec XVY(q, arma::fill::zeros);

    arma::mat G_inv;
    try {
        G_inv = arma::inv_sympd(G);
    } catch (...) {
        G_inv = arma::inv(G);
    }

    int offset = 0;
    for (int g = 0; g < nc.size(); ++g) {
        const int m = nc[g];
        if (m <= 0) continue;

        arma::mat Xg(m, q, arma::fill::ones);
        Xg.cols(1, p) = xx.rows(offset, offset + m - 1);
        arma::vec yg = yy.subvec(offset, offset + m - 1);

        arma::mat XtX = Xg.t() * Xg;

        // Woodbury identity:
        // V_g^{-1} = (1/Ve) [ I - Xg (Ve G^{-1} + Xg'Xg)^{-1} Xg' ]
        arma::mat M = Ve * G_inv + XtX;
        arma::mat M_inv;
        try {
            M_inv = arma::inv_sympd(M);
        } catch (...) {
            M_inv = arma::inv(M);
        }

        arma::mat invV_X = (1.0 / Ve) * (Xg - Xg * (M_inv * XtX));
        arma::vec invV_y = (1.0 / Ve) * (yg - Xg * (M_inv * (Xg.t() * yg)));

        XVX += Xg.t() * invV_X;
        XVY += Xg.t() * invV_y;

        offset += m;
    }

    return arma::solve(XVX, XVY);
}

// 严谨的 EM 方差更新函数 (请用此替换原有的同名函数)
static inline Rcpp::List rs_count_find_sigma_internal(const arma::mat& xx,
                                                      const arma::vec& yy,
                                                      const arma::vec& beta,
                                                      const arma::mat& G_old,
                                                      double Ve_old,
                                                      Rcpp::IntegerVector nc,
                                                      int p) {
    const int q = p + 1;
    arma::mat G_new(q, q, arma::fill::zeros);
    double Ve_new = 0.0;

    arma::mat G_inv_old;
    try {
        G_inv_old = arma::inv_sympd(G_old);
    } catch (...) {
        G_inv_old = arma::inv(G_old);
    }

    int offset = 0;
    int n_groups_nonzero = 0;
    int n_total = 0;

    for (int g = 0; g < nc.size(); ++g) {
        const int m = nc[g];
        if (m <= 0) continue;

        ++n_groups_nonzero;
        n_total += m;

        arma::mat Xg(m, q, arma::fill::ones);
        Xg.cols(1, p) = xx.rows(offset, offset + m - 1);
        arma::vec yg = yy.subvec(offset, offset + m - 1);

        arma::vec res_g = yg - Xg * beta;
        arma::mat XtX = Xg.t() * Xg;

        arma::mat M = Ve_old * G_inv_old + XtX;
        arma::mat M_inv;
        try {
            M_inv = arma::inv_sympd(M);
        } catch (...) {
            M_inv = arma::inv(M);
        }

        // 计算条件方差 Sigma_u 
        arma::mat Sigma_u = Ve_old * M_inv;

        arma::vec Vg_inv_res = (1.0 / Ve_old) * (res_g - Xg * (M_inv * (Xg.t() * res_g)));
        arma::vec ug = G_old * Xg.t() * Vg_inv_res;

        // 补全期望的后验方差项 (解决低估问题)
        G_new += ug * ug.t() + Sigma_u;

        arma::vec err = res_g - Xg * ug;
        
        // 补全残差方差的迹项
        Ve_new += arma::dot(err, err) + arma::trace(XtX * Sigma_u);

        offset += m;
    }

    if (n_groups_nonzero <= 0) n_groups_nonzero = 1;
    if (n_total <= 0) n_total = 1;

    G_new /= static_cast<double>(n_groups_nonzero);
    
    // 为了和无随机斜率版的自由度保持完全一致，分母修改为 N 减 R
    double denom_e = static_cast<double>(n_total - n_groups_nonzero);
    if (denom_e <= 0.0) denom_e = 1.0;
    Ve_new /= denom_e;

    // 稳定性保护
    if (Ve_new <= 1e-8) Ve_new = 1e-8;
    G_new.diag() += 1e-7;

    return Rcpp::List::create(
        Rcpp::Named("G") = G_new,
        Rcpp::Named("Var.e") = Ve_new
    );
}

// [[Rcpp::export]]
Rcpp::List count_info_rs_cpp(const arma::mat& xx_in,
                             const arma::vec& yy,
                             Rcpp::IntegerVector nc,
                             int R,
                             int p) {
    arma::mat xx = xx_in;

    double D = 0.0;
    double A = 0.0;
    arma::mat infor;

    // ===== case 1: single group -> ordinary linear model =====
    if (nc.size() == 1) {
        arma::mat X(xx.n_rows, p + 1, arma::fill::ones);
        X.cols(1, p) = xx;

        arma::vec b = arma::solve(X, yy);
        arma::vec resid = yy - X * b;

        double RSS = arma::dot(resid, resid);
        double Var_lm = RSS / (static_cast<double>(X.n_rows) - p - 1.0);
        if (Var_lm <= 1e-8) Var_lm = 1e-8;

        infor = (X.t() * X) / Var_lm;

        D = arma::det(infor);
        if (!arma::is_finite(D) || D <= 0.0) D = 1e-50;

        try {
            A = arma::trace(arma::inv_sympd(infor));
        } catch (...) {
            try {
                A = arma::trace(arma::inv(infor));
            } catch (...) {
                A = NA_REAL;
            }
        }
    }
    // ===== case 2: mixed model with random slopes =====
    else {
        const int q = p + 1;

        // ----- initialize -----
        arma::vec beta_curr = rs_count_ols_beta(xx, yy);
        arma::mat G_curr = arma::eye(q, q) * 0.1;
        double Ve_curr = 1.0;

        // ----- iterate -----
        for (int it = 0; it < 3; ++it) {
            Rcpp::List sigs = rs_count_find_sigma_internal(xx, yy, beta_curr, G_curr, Ve_curr, nc, p);
            G_curr = Rcpp::as<arma::mat>(sigs["G"]);
            Ve_curr = Rcpp::as<double>(sigs["Var.e"]);

            G_curr.diag() += 1e-7;
            if (Ve_curr <= 1e-8) Ve_curr = 1e-8;

            beta_curr = rs_count_find_beta_internal(xx, yy, G_curr, Ve_curr, nc, p);
        }

        // ----- build information matrix -----
        infor = arma::mat(q, q, arma::fill::zeros);

        arma::mat G_inv;
        try {
            G_inv = arma::inv_sympd(G_curr);
        } catch (...) {
            G_inv = arma::inv(G_curr);
        }

        int offset = 0;
        for (int g = 0; g < nc.size(); ++g) {
            const int m = nc[g];
            if (m <= 0) continue;

            arma::mat Xg(m, q, arma::fill::ones);
            Xg.cols(1, p) = xx.rows(offset, offset + m - 1);

            arma::mat XtX = Xg.t() * Xg;
            arma::mat M = Ve_curr * G_inv + XtX;
            arma::mat M_inv;

            try {
                M_inv = arma::inv_sympd(M);
            } catch (...) {
                M_inv = arma::inv(M);
            }

            // V_g^{-1} X_g
            arma::mat invV_X = (1.0 / Ve_curr) * (Xg - Xg * (M_inv * XtX));

            infor += Xg.t() * invV_X;

            offset += m;
        }

        D = arma::det(infor);
        if (!arma::is_finite(D) || D <= 0.0) D = 1e-50;

        try {
            A = arma::trace(arma::inv_sympd(infor));
        } catch (...) {
            try {
                A = arma::trace(arma::inv(infor));
            } catch (...) {
                A = NA_REAL;
            }
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("D") = D,
        Rcpp::Named("A") = A,
        Rcpp::Named("Information") = infor
    );
}


// 内部辅助：随机截距版 EM 方差更新
static inline Rcpp::List find_sigma_ri_internal(const arma::vec& yy,
                                                const arma::mat& X,
                                                const arma::vec& beta,
                                                double Var_a_old,
                                                double Var_e_old,
                                                Rcpp::IntegerVector nc) {
    double Var_a_new = 0.0;
    double Var_e_new = 0.0;
    
    int offset = 0;
    int n_groups_nonzero = 0;
    int n_total = 0;
    
    // 计算当前残差
    arma::vec resid = yy - X * beta;
    
    for (int g = 0; g < nc.size(); ++g) {
        int m = nc[g];
        if (m <= 0) continue;
        
        ++n_groups_nonzero;
        n_total += m;
        
        arma::vec res_g = resid.subvec(offset, offset + m - 1);
        double sum_res = arma::accu(res_g);
        
        // 计算 M 等于 Var_e 加上 m 乘以 Var_a
        double M = Var_e_old + static_cast<double>(m) * Var_a_old;
        
        // 组级别的条件方差 Sigma_u
        double Sigma_u = (Var_e_old * Var_a_old) / M;
        
        // 随机截距的经验贝叶斯估计 u_g
        double u_g = (Var_a_old / M) * sum_res;
        
        // 补全期望的后验方差项 (严格对齐 REML)
        Var_a_new += (u_g * u_g) + Sigma_u;
        
        // 补全残差方差项
        arma::vec err = res_g - u_g;
        Var_e_new += arma::dot(err, err) + static_cast<double>(m) * Sigma_u;
        
        offset += m;
    }
    
    if (n_groups_nonzero <= 0) n_groups_nonzero = 1;
    if (n_total <= 0) n_total = 1;
    
    Var_a_new /= static_cast<double>(n_groups_nonzero);
    
    // 自由度修正 N 减去 R
    double denom_e = static_cast<double>(n_total - n_groups_nonzero);
    if (denom_e <= 0.0) denom_e = 1.0;
    Var_e_new /= denom_e;
    
    // 稳定性保护
    if (Var_a_new <= 1e-8) Var_a_new = 1e-8;
    if (Var_e_new <= 1e-8) Var_e_new = 1e-8;
    
    return Rcpp::List::create(
        Rcpp::Named("Var.a") = Var_a_new,
        Rcpp::Named("Var.e") = Var_e_new
    );
}

// 替换原有的 Est_hat_cpp 函数
// [[Rcpp::export]]
Rcpp::List Est_hat_cpp(const arma::mat& xx_in,
                       const arma::vec& yy,
                       const arma::vec& beta_true,
                       double Var_a_true,
                       double Var_e_true,
                       Rcpp::IntegerVector nc,
                       int R,
                       int p) {
    arma::mat xx = xx_in;
    int q = p + 1;
    
    // 构建包含截距的设计矩阵
    arma::mat X(xx.n_rows, q, arma::fill::ones);
    X.cols(1, p) = xx;
    
    // 1. 初始值 (OLS 与固定初始方差，严格对齐另一套代码)
    arma::vec beta_curr = arma::solve(X, yy);
    double Var_a_curr = 0.1;
    double Var_e_curr = 1.0;
    
    // 2. 迭代 (固定 3 轮，匹配随机斜率逻辑深度)
    for(int i = 0; i < 3; ++i) {
        Rcpp::List sigs = find_sigma_ri_internal(yy, X, beta_curr, Var_a_curr, Var_e_curr, nc);
        Var_a_curr = Rcpp::as<double>(sigs["Var.a"]);
        Var_e_curr = Rcpp::as<double>(sigs["Var.e"]);
        
        beta_curr = find_beta_cpp(xx, yy, Var_a_curr, Var_e_curr, nc, R, p);
    }
    
    // 3. 计算 MSE
    arma::vec beta2_no_intercept = beta_curr.subvec(1, p);
    double bt_mse  = arma::accu(arma::square(beta2_no_intercept - beta_true));
    double bt0_mse = std::pow(beta_curr[0] - 1.0, 2.0);
    double va_mse  = std::pow(Var_a_curr - Var_a_true, 2.0);
    double ve_mse  = std::pow(Var_e_curr - Var_e_true, 2.0);
    
    return Rcpp::List::create(
        Rcpp::Named("bt.mse")  = bt_mse,
        Rcpp::Named("va.mse")  = va_mse,
        Rcpp::Named("ve.mse")  = ve_mse,
        Rcpp::Named("bt0.mse") = bt0_mse,
        Rcpp::Named("beta2")   = beta_curr,
        Rcpp::Named("Var.a")   = Var_a_curr,
        Rcpp::Named("Var.e")   = Var_e_curr
    );
}


// 替换原有的 count_info_cpp，请确保将其放置在 find_sigma_ri_internal 之后
// [[Rcpp::export]]
Rcpp::List count_info_cpp(const arma::mat& xx_in,
                          const arma::vec& yy,
                          Rcpp::IntegerVector nc,
                          int R,
                          int p) {
    arma::mat xx = xx_in;
    int q = p + 1;
    double D = 0.0;
    double A = 0.0;
    arma::mat infor;

    // 单组情况 (普通线性模型)
    if (nc.size() == 1) {
        arma::mat X(xx.n_rows, q, arma::fill::ones);
        X.cols(1, p) = xx;

        arma::vec b = arma::solve(X, yy);
        arma::vec resid = yy - X * b;

        double RSS = arma::dot(resid, resid);
        double Var_lm = RSS / (static_cast<double>(X.n_rows) - p - 1.0);
        
        if (Var_lm <= 1e-8) Var_lm = 1e-8;

        infor = (X.t() * X) / Var_lm;
        
        D = arma::det(infor);
        if (D <= 0.0) D = 1e-50; 
        
        try {
            A = arma::trace(arma::inv(infor));
        } catch (...) {
            A = NA_REAL; 
        }
    } 
    // 多组情况 (混合效应模型)
    else {
        arma::mat X(xx.n_rows, q, arma::fill::ones);
        X.cols(1, p) = xx;
        
        // 初始值与 Est_hat_cpp 绝对对齐
        arma::vec beta_curr = arma::solve(X, yy);
        double Var_a_curr = 0.1;
        double Var_e_curr = 1.0;

        // 严格套用三轮 EM 迭代
        for(int i = 0; i < 3; ++i) {
            Rcpp::List sigs = find_sigma_ri_internal(yy, X, beta_curr, Var_a_curr, Var_e_curr, nc);
            Var_a_curr = Rcpp::as<double>(sigs["Var.a"]);
            Var_e_curr = Rcpp::as<double>(sigs["Var.e"]);
            
            beta_curr = find_beta_cpp(xx, yy, Var_a_curr, Var_e_curr, nc, R, p);
        }

        if (Var_e_curr <= 1e-8) Var_e_curr = 1e-8;

        // 构建信息矩阵
        arma::mat XVX(q, q, arma::fill::zeros);
        int offset = 0;
        
        for (int g = 0; g < nc.size(); ++g) {
            int m = nc[g];
            if (m <= 0) continue;

            arma::mat Xg(m, q, arma::fill::ones);
            Xg.cols(1, p) = xx.rows(offset, offset + m - 1);

            double denom = Var_e_curr + static_cast<double>(m) * Var_a_curr;
            double weight = (denom <= 0) ? 0.0 : Var_a_curr / (Var_e_curr * denom);

            arma::mat term1 = (Xg.t() * Xg) / Var_e_curr;
            arma::mat sumX = arma::sum(Xg, 0).t();
            arma::mat term2 = weight * (sumX * sumX.t());

            XVX += (term1 - term2);
            offset += m;
        }

        infor = XVX;
        
        D = arma::det(infor);
        if (D <= 0.0) D = 1e-50; 
        
        try {
            A = arma::trace(arma::inv(infor));
        } catch (...) {
            A = NA_REAL;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("D") = D,
        Rcpp::Named("A") = A,
        Rcpp::Named("Information") = infor
    );
}

// [[Rcpp::export]]
Rcpp::List Est_hat_RS_cpp(const arma::mat& xx_in,
                          const arma::vec& yy,
                          const arma::vec& beta_true, 
                          double Var_a_true,          
                          double Var_e_true,
                          Rcpp::IntegerVector nc,
                          int R,
                          int p) {
    arma::mat xx = xx_in;
    int q = p + 1;

    // 1. 初始值 (OLS)
    arma::mat D(xx.n_rows, q, arma::fill::ones);
    D.cols(1, p) = xx;
    arma::vec beta_curr = arma::solve(D, yy);
    arma::mat G_curr = arma::eye(q, q) * 0.1;
    double Ve_curr = 1.0;

    // 2. 迭代 (固定 3 轮，匹配逻辑深度)
    for(int i = 0; i < 3; ++i) {
        // 【关键修复：已替换为带后验方差修正的新函数】
        Rcpp::List sigs = rs_count_find_sigma_internal(xx, yy, beta_curr, G_curr, Ve_curr, nc, p);
        G_curr = Rcpp::as<arma::mat>(sigs["G"]);
        Ve_curr = Rcpp::as<double>(sigs["Var.e"]);
        
        G_curr.diag() += 1e-7; 
        if (Ve_curr <= 1e-8) Ve_curr = 1e-8;
        
        // 【关键修复：已替换为对应的新版 Beta 更新函数】
        beta_curr = rs_count_find_beta_internal(xx, yy, G_curr, Ve_curr, nc, p);
    }

    // 3. 计算 MSE
    arma::vec beta2_no_intercept = beta_curr.subvec(1, p);
    double bt_mse  = arma::accu(arma::square(beta2_no_intercept - beta_true));
    double bt0_mse = std::pow(beta_curr[0] - 1.0, 2.0); 
    double va_mse  = std::pow(G_curr(0,0) - Var_a_true, 2.0); 
    double ve_mse  = std::pow(Ve_curr - Var_e_true, 2.0);

    return Rcpp::List::create(
        Rcpp::Named("bt.mse")  = bt_mse,
        Rcpp::Named("va.mse")  = va_mse, 
        Rcpp::Named("ve.mse")  = ve_mse,
        Rcpp::Named("bt0.mse") = bt0_mse,
        Rcpp::Named("beta2")   = beta_curr,
        Rcpp::Named("Var.a")   = G_curr(0,0), 
        Rcpp::Named("Var.e")   = Ve_curr,
        Rcpp::Named("G.hat")   = G_curr       
    );
}