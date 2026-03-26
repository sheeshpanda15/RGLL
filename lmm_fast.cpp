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

// [[Rcpp::export]]
Rcpp::List Est_hat_cpp(const arma::mat& xx_in,
                       const arma::vec& yy,
                       const arma::vec& beta_true, // length p (NO intercept)
                       double Var_a_true,
                       double Var_e_true,
                       Rcpp::IntegerVector nc,
                       int R,
                       int p) {
  // scale X as in R
  mat xx = xx_in;

  // beta0 from OLS
  vec beta0 = ols_beta(xx, yy);

  // sigma0 -> beta1 -> sigma1 -> beta2 -> sigma2
  Rcpp::List s0 = find_sigma_cpp(xx, yy, beta0, nc, R);
  vec beta1 = find_beta_cpp(xx, yy,
                            Rcpp::as<double>(s0["Var.a"]),
                            Rcpp::as<double>(s0["Var.e"]),
                            nc, R, p);

  Rcpp::List s1 = find_sigma_cpp(xx, yy, beta1, nc, R);
  vec beta2 = find_beta_cpp(xx, yy,
                            Rcpp::as<double>(s1["Var.a"]),
                            Rcpp::as<double>(s1["Var.e"]),
                            nc, R, p);

  Rcpp::List s2 = find_sigma_cpp(xx, yy, beta2, nc, R);
  const double Var_a_hat = Rcpp::as<double>(s2["Var.a"]);
  const double Var_e_hat = Rcpp::as<double>(s2["Var.e"]);

  // MSEs (match your R definition)
  vec beta2_no_intercept = beta2.subvec(1, p);
  const double bt_mse  = arma::accu(arma::square(beta2_no_intercept - beta_true));
  const double bt0_mse = std::pow(beta2[0] - 1.0, 2.0);
  const double va_mse  = std::pow(Var_a_hat - Var_a_true, 2.0);
  const double ve_mse  = std::pow(Var_e_hat - Var_e_true, 2.0);

  return Rcpp::List::create(
    Rcpp::Named("bt.mse")  = bt_mse,
    Rcpp::Named("va.mse")  = va_mse,
    Rcpp::Named("ve.mse")  = ve_mse,
    Rcpp::Named("bt0.mse") = bt0_mse,
    Rcpp::Named("beta2")   = beta2,
    Rcpp::Named("Var.a")   = Var_a_hat,
    Rcpp::Named("Var.e")   = Var_e_hat
  );
}

// ---------------- 替换从这里开始 ----------------

// [[Rcpp::export]]
Rcpp::List count_info_cpp(const arma::mat& xx_in,
                          const arma::vec& yy,
                          Rcpp::IntegerVector nc,
                          int R,
                          int p) {
    mat xx = xx_in; 

    double D = 0.0;
    double A = 0.0;
    mat infor;

    // 情况 1: Single group (Linear Model)
    if (nc.size() == 1) {
        mat X(xx.n_rows, p + 1, arma::fill::ones);
        X.cols(1, p) = xx;

        vec b = arma::solve(X, yy);
        vec resid = yy - X * b;

        double RSS = arma::dot(resid, resid);
        double Var_lm = RSS / (static_cast<double>(X.n_rows) - p - 1.0);
        
        // 增加安全保护，防止方差极小导致除零
        if (Var_lm <= 1e-8) Var_lm = 1e-8;

        // 信息矩阵
        infor = (X.t() * X) / Var_lm;
        
        // 计算准则
        D = arma::det(infor);
        if (D <= 0.0) D = 1e-50; // 强制矫正，防止 R 端 log() 报非数值错误
        
        try {
            A = arma::trace(arma::inv(infor));
        } catch (...) {
            A = NA_REAL; // 如果完全不可逆则返回 NA
        }
    } 
    // 情况 2: Mixed Model (LMM)
    else {
        // 迭代估计方差分量
        vec beta0 = ols_beta(xx, yy);

        Rcpp::List s0 = find_sigma_cpp(xx, yy, beta0, nc, R);
        vec beta1 = find_beta_cpp(xx, yy, 
                                  Rcpp::as<double>(s0["Var.a"]), 
                                  Rcpp::as<double>(s0["Var.e"]), 
                                  nc, R, p);

        Rcpp::List s1 = find_sigma_cpp(xx, yy, beta1, nc, R);
        vec beta2 = find_beta_cpp(xx, yy, 
                                  Rcpp::as<double>(s1["Var.a"]), 
                                  Rcpp::as<double>(s1["Var.e"]), 
                                  nc, R, p);

        Rcpp::List s2 = find_sigma_cpp(xx, yy, beta2, nc, R);
        double Var_a = Rcpp::as<double>(s2["Var.a"]);
        double Var_e = Rcpp::as<double>(s2["Var.e"]);

        // 增加安全保护：防止 Var_e 趋近于 0 导致矩阵出现 NaN
        if (Var_e <= 1e-8) Var_e = 1e-8;

        // 构建 XVX (信息矩阵)
        mat XVX(p + 1, p + 1, arma::fill::zeros);
        int offset = 0;
        
        for (int g = 0; g < nc.size(); ++g) {
            int m = nc[g];
            if (m <= 0) continue;

            mat Xg(m, p + 1, arma::fill::ones);
            Xg.cols(1, p) = xx.rows(offset, offset + m - 1);

            double denom = Var_e + static_cast<double>(m) * Var_a;
            double weight = (denom <= 0) ? 0.0 : Var_a / (Var_e * denom);

            mat term1 = (Xg.t() * Xg) / Var_e;
            mat sumX = arma::sum(Xg, 0).t();
            mat term2 = weight * (sumX * sumX.t());

            XVX += (term1 - term2);
            offset += m;
        }

        infor = XVX;
        
        D = arma::det(infor);
        if (D <= 0.0) D = 1e-50; // 强制矫正，防止 R 端 log() 报非数值错误
        
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

// ---------------- 替换到这里结束 ----------------

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

// [[Rcpp::export]]
Rcpp::List Est_hat_RS_cpp(const arma::mat& xx_in,
                                  const arma::vec& yy,
                                  const arma::vec& beta_true, // 长度 p (不含截距)
                                  double Var_a_true,          // 仅用于兼容接口，内部对比 G[0,0]
                                  double Var_e_true,
                                  Rcpp::IntegerVector nc,
                                  int R,
                                  int p) {
    mat xx = xx_in;
    int q = p + 1;

    // 1. 初始值 (OLS)
    mat D(xx.n_rows, q, arma::fill::ones);
    D.cols(1, p) = xx;
    vec beta_curr = arma::solve(D, yy);
    mat G_curr = arma::eye(q, q) * 0.1;
    double Ve_curr = 1.0;

    // 2. 迭代 (固定 3 轮，匹配原代码逻辑深度)
    for(int i = 0; i < 3; ++i) {
        Rcpp::List sigs = find_sigma_full_rs_internal(xx, yy, beta_curr, G_curr, Ve_curr, nc, p);
        G_curr = Rcpp::as<mat>(sigs["G"]);
        Ve_curr = Rcpp::as<double>(sigs["Var.e"]);
        G_curr.diag() += 1e-7; 
        beta_curr = find_beta_full_rs_internal(xx, yy, G_curr, Ve_curr, nc, p);
    }

    // 3. 计算 MSE (完全对齐原输出格式)
    vec beta2_no_intercept = beta_curr.subvec(1, p);
    double bt_mse  = arma::accu(arma::square(beta2_no_intercept - beta_true));
    double bt0_mse = std::pow(beta_curr[0] - 1.0, 2.0); // 假设真截距为 1.0
    double va_mse  = std::pow(G_curr(0,0) - Var_a_true, 2.0); // 对比截距方差
    double ve_mse  = std::pow(Ve_curr - Var_e_true, 2.0);

    return Rcpp::List::create(
        Rcpp::Named("bt.mse")  = bt_mse,
        Rcpp::Named("va.mse")  = va_mse, // 注意：此处对比的是 G 矩阵的第一个对角元素
        Rcpp::Named("ve.mse")  = ve_mse,
        Rcpp::Named("bt0.mse") = bt0_mse,
        Rcpp::Named("beta2")   = beta_curr,
        Rcpp::Named("Var.a")   = G_curr(0,0), // 返回截距的随机方差
        Rcpp::Named("Var.e")   = Ve_curr,
        Rcpp::Named("G.hat")   = G_curr       // 额外附带全矩阵供参考
    );
}