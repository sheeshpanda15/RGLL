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

// [[Rcpp::export]]
arma::vec count_info_cpp(const arma::mat& xx_in,
                         const arma::vec& yy,
                         Rcpp::IntegerVector nc,
                         int R,
                         int p) {
  // scale X as in R
  mat xx = scalex_mat(xx_in);

  // If single group -> LM case
  if (nc.size() == 1) {
    // Design matrix with intercept
    mat X(xx.n_rows, p + 1, arma::fill::ones);
    X.cols(1, p) = xx;

    vec b = arma::solve(X, yy);
    vec resid = yy - X * b;

    double RSS = arma::dot(resid, resid);
    double Var_lm = RSS / (static_cast<double>(X.n_rows) - p - 1.0);

    mat XtX = X.t() * X;
    double D = arma::det(XtX) / std::pow(Var_lm, p);
    double A = Var_lm * arma::trace(arma::inv(XtX));

    vec out(2);
    out[0] = D;
    out[1] = A;
    return out;
  }

  // LMM case: estimate Var.a, Var.e using same 0-1-2 iteration as your R code
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

  // Build XVX only (same as your count.info else-branch)
  // drop nc==0
  std::vector<int> ncz;
  ncz.reserve(nc.size());
  for (int v : nc) if (v != 0) ncz.push_back(v);

  const int R1 = static_cast<int>(ncz.size());
  mat XVX(p + 1, p + 1, arma::fill::zeros);

  int offset = 0;
  for (int g = 0; g < R1; ++g) {
    const int m = ncz[g];
    mat Xg(m, p + 1, arma::fill::ones);
    Xg.cols(1, p) = xx.rows(offset, offset + m - 1);

    const double denom = (Var_e + static_cast<double>(m) * Var_a);
    const double gamma = (denom == 0.0) ? 0.0 : (static_cast<double>(m) * Var_a) / denom;

    const double invVe = 1.0 / Var_e;
    vec ones(m, arma::fill::ones);
    arma::rowvec sX = arma::sum(Xg, 0);

    mat invV_X = invVe * (Xg - (gamma / static_cast<double>(m)) * (ones * sX));
    XVX += Xg.t() * invV_X;

    offset += m;
  }

  double D = arma::det(XVX);
  double A = arma::trace(arma::inv(XVX));

  vec out(2);
  out[0] = D;
  out[1] = A;
  return out;
}
