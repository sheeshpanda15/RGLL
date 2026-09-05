// lmm_fast_aligned.cpp
// Corrected implementation for the CIRG experiments.
// Main statistical choices:
//   1) convergence-based ML-EM estimation (not fixed three iterations),
//   2) RI and RS models have separate IMSPE objectives,
//   3) the RS covariance is structured as
//        G = diag(var_a, var_b, ..., var_b),
//      matching the simulation DGP and avoiding an unidentifiable full
//      (p+1)x(p+1) covariance with only a small number of groups,
//   4) full random-intercept + random-slope BLUP prediction.

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using arma::mat;
using arma::vec;
using arma::uword;

namespace {

constexpr double EPS_VAR = 1e-8;
constexpr double RIDGE0  = 1e-8;

mat add_intercept(const mat& X) {
  mat F(X.n_rows, X.n_cols + 1, arma::fill::ones);
  if (X.n_cols > 0) F.cols(1, X.n_cols) = X;
  return F;
}

mat safe_inv_spd(const mat& A) {
  mat B = 0.5 * (A + A.t());
  mat out;
  double ridge = RIDGE0;
  for (int k = 0; k < 8; ++k) {
    if (arma::inv_sympd(out, B)) return out;
    B.diag() += ridge;
    ridge *= 10.0;
  }
  return arma::pinv(B);
}

vec safe_solve(const mat& A, const vec& b) {
  vec out;
  bool ok = arma::solve(out, A, b, arma::solve_opts::fast);
  if (!ok || !out.is_finite()) out = arma::pinv(A) * b;
  return out;
}

double rel_vec_change(const vec& a, const vec& b) {
  return arma::norm(a - b, 2) / (1.0 + arma::norm(b, 2));
}

double rel_mat_change(const mat& a, const mat& b) {
  return arma::norm(a - b, "fro") / (1.0 + arma::norm(b, "fro"));
}

std::vector<int> positive_sizes(Rcpp::IntegerVector nc) {
  std::vector<int> out;
  out.reserve(nc.size());
  for (int v : nc) if (v > 0) out.push_back(v);
  return out;
}

struct RIFit {
  vec beta;
  double va;
  double ve;
  mat information;
  vec uhat;
  bool converged;
  int iterations;
};

vec ri_gls_beta(const mat& X, const vec& y, double va, double ve,
                const std::vector<int>& sizes) {
  const int q = static_cast<int>(X.n_cols);
  mat XVX(q, q, arma::fill::zeros);
  vec XVY(q, arma::fill::zeros);
  int offset = 0;
  for (int m : sizes) {
    mat Xg = X.rows(offset, offset + m - 1);
    vec yg = y.subvec(offset, offset + m - 1);
    double denom = ve + static_cast<double>(m) * va;
    double gamma = (denom <= 0.0) ? 0.0 : (static_cast<double>(m) * va / denom);
    vec one(m, arma::fill::ones);
    arma::rowvec sx = arma::sum(Xg, 0);
    double sy = arma::accu(yg);
    mat VinvX = (Xg - (gamma / static_cast<double>(m)) * one * sx) / ve;
    vec Vinvy = (yg - (gamma / static_cast<double>(m)) * one * sy) / ve;
    XVX += Xg.t() * VinvX;
    XVY += Xg.t() * Vinvy;
    offset += m;
  }
  return safe_solve(XVX, XVY);
}

mat ri_information(const mat& X, double va, double ve,
                   const std::vector<int>& sizes) {
  mat info(X.n_cols, X.n_cols, arma::fill::zeros);
  int offset = 0;
  for (int m : sizes) {
    mat Xg = X.rows(offset, offset + m - 1);
    double denom = ve + static_cast<double>(m) * va;
    double weight = (denom <= 0.0) ? 0.0 : va / (ve * denom);
    vec sx = arma::sum(Xg, 0).t();
    info += (Xg.t() * Xg) / ve - weight * (sx * sx.t());
    offset += m;
  }
  return 0.5 * (info + info.t());
}

RIFit fit_ri_internal(const mat& xx, const vec& y, Rcpp::IntegerVector nc,
                      double tol, int max_iter) {
  mat X = add_intercept(xx);
  std::vector<int> sizes = positive_sizes(nc);
  const int K = static_cast<int>(sizes.size());
  if (K < 1) Rcpp::stop("No non-empty groups supplied to RI fit.");

  vec beta = safe_solve(X.t() * X, X.t() * y);
  vec r0 = y - X * beta;
  double ve = std::max(EPS_VAR, arma::dot(r0, r0) / std::max(1.0, static_cast<double>(y.n_elem)));
  double va = std::max(0.1, 0.1 * ve);

  bool converged = false;
  int used = 0;
  for (int it = 1; it <= max_iter; ++it) {
    double va_new = 0.0;
    double ve_num = 0.0;
    int offset = 0;
    vec resid = y - X * beta;

    for (int m : sizes) {
      vec rg = resid.subvec(offset, offset + m - 1);
      double denom = ve + static_cast<double>(m) * va;
      denom = std::max(denom, EPS_VAR);
      double post_var = (ve * va) / denom;
      double ug = (va / denom) * arma::accu(rg);
      va_new += ug * ug + post_var;
      vec err = rg - ug;
      ve_num += arma::dot(err, err) + static_cast<double>(m) * post_var;
      offset += m;
    }

    va_new = std::max(EPS_VAR, va_new / static_cast<double>(K));
    // This is an ML-EM update. It is deliberately N in the denominator;
    // the legacy N-K denominator was not an exact EM/REML update.
    double ve_new = std::max(EPS_VAR, ve_num / static_cast<double>(y.n_elem));
    vec beta_new = ri_gls_beta(X, y, va_new, ve_new, sizes);

    double ch = std::max({
      std::abs(va_new - va) / (1.0 + std::abs(va)),
      std::abs(ve_new - ve) / (1.0 + std::abs(ve)),
      rel_vec_change(beta_new, beta)
    });

    va = va_new;
    ve = ve_new;
    beta = beta_new;
    used = it;
    if (ch < tol) {
      converged = true;
      break;
    }
  }

  vec uhat(K, arma::fill::zeros);
  vec resid = y - X * beta;
  int offset = 0;
  for (int k = 0; k < K; ++k) {
    int m = sizes[k];
    vec rg = resid.subvec(offset, offset + m - 1);
    double denom = std::max(EPS_VAR, ve + static_cast<double>(m) * va);
    uhat[k] = (va / denom) * arma::accu(rg);
    offset += m;
  }

  RIFit ans;
  ans.beta = beta;
  ans.va = va;
  ans.ve = ve;
  ans.information = ri_information(X, va, ve, sizes);
  ans.uhat = uhat;
  ans.converged = converged;
  ans.iterations = used;
  return ans;
}

struct RSFit {
  vec beta;
  mat G;
  double ve;
  mat information;
  mat uhat; // K x q, one row per group
  bool converged;
  int iterations;
};

double rs_var_b(const mat& G) {
  if (G.n_rows <= 1) return NA_REAL;
  double total = 0.0;
  for (uword j = 1; j < G.n_rows; ++j) total += G(j, j);
  return total / static_cast<double>(G.n_rows - 1);
}

vec rs_gls_beta(const mat& X, const vec& y, const mat& G, double ve,
                const std::vector<int>& sizes) {
  const int q = static_cast<int>(X.n_cols);
  mat Ginv = safe_inv_spd(G);
  mat XVX(q, q, arma::fill::zeros);
  vec XVY(q, arma::fill::zeros);
  int offset = 0;

  for (int m : sizes) {
    mat Xg = X.rows(offset, offset + m - 1);
    vec yg = y.subvec(offset, offset + m - 1);
    mat XtX = Xg.t() * Xg;
    mat Minv = safe_inv_spd(ve * Ginv + XtX);
    mat VinvX = (Xg - Xg * Minv * XtX) / ve;
    vec Vinvy = (yg - Xg * Minv * (Xg.t() * yg)) / ve;
    XVX += Xg.t() * VinvX;
    XVY += Xg.t() * Vinvy;
    offset += m;
  }
  return safe_solve(XVX, XVY);
}

mat rs_information(const mat& X, const mat& G, double ve,
                   const std::vector<int>& sizes) {
  mat Ginv = safe_inv_spd(G);
  mat info(X.n_cols, X.n_cols, arma::fill::zeros);
  int offset = 0;
  for (int m : sizes) {
    mat Xg = X.rows(offset, offset + m - 1);
    mat XtX = Xg.t() * Xg;
    mat Minv = safe_inv_spd(ve * Ginv + XtX);
    mat VinvX = (Xg - Xg * Minv * XtX) / ve;
    info += Xg.t() * VinvX;
    offset += m;
  }
  return 0.5 * (info + info.t());
}

RSFit fit_rs_internal(const mat& xx, const vec& y, Rcpp::IntegerVector nc,
                      double tol, int max_iter) {
  mat X = add_intercept(xx); // same design for fixed and random effects
  std::vector<int> sizes = positive_sizes(nc);
  const int K = static_cast<int>(sizes.size());
  const int q = static_cast<int>(X.n_cols);
  if (K < 1) Rcpp::stop("No non-empty groups supplied to RS fit.");
  if (q < 2) Rcpp::stop("RS fit requires at least one slope column.");

  vec beta = safe_solve(X.t() * X, X.t() * y);
  vec r0 = y - X * beta;
  double ve = std::max(EPS_VAR, arma::dot(r0, r0) /
                                  std::max(1.0, static_cast<double>(y.n_elem)));

  // Structured covariance matching the DGP:
  // G = diag(var_a, var_b, ..., var_b).
  double va = std::max(0.1, 0.10 * ve);
  double vb = std::max(0.01, 0.01 * ve);
  mat G(q, q, arma::fill::zeros);
  G(0, 0) = va;
  for (int j = 1; j < q; ++j) G(j, j) = vb;

  bool converged = false;
  int used = 0;

  for (int it = 1; it <= max_iter; ++it) {
    mat Ginv = safe_inv_spd(G);
    double va_num = 0.0;
    double vb_num = 0.0;
    double ve_num = 0.0;
    int offset = 0;
    vec resid = y - X * beta;

    for (int m : sizes) {
      mat Xg = X.rows(offset, offset + m - 1);
      vec rg = resid.subvec(offset, offset + m - 1);
      mat XtX = Xg.t() * Xg;
      mat Minv = safe_inv_spd(ve * Ginv + XtX);
      mat Sigma_u = ve * Minv;
      vec ug = Minv * (Xg.t() * rg);
      mat Euu = ug * ug.t() + Sigma_u;

      va_num += Euu(0, 0);
      for (int j = 1; j < q; ++j) vb_num += Euu(j, j);

      vec err = rg - Xg * ug;
      ve_num += arma::dot(err, err) + arma::trace(XtX * Sigma_u);
      offset += m;
    }

    double va_new = std::max(EPS_VAR, va_num / static_cast<double>(K));
    double vb_new = std::max(EPS_VAR,
      vb_num / static_cast<double>(K * (q - 1)));
    double ve_new = std::max(EPS_VAR,
      ve_num / static_cast<double>(y.n_elem));

    mat G_new(q, q, arma::fill::zeros);
    G_new(0, 0) = va_new;
    for (int j = 1; j < q; ++j) G_new(j, j) = vb_new;

    vec beta_new = rs_gls_beta(X, y, G_new, ve_new, sizes);

    double ch = std::max({
      std::abs(va_new - va) / (1.0 + std::abs(va)),
      std::abs(vb_new - vb) / (1.0 + std::abs(vb)),
      std::abs(ve_new - ve) / (1.0 + std::abs(ve)),
      rel_vec_change(beta_new, beta)
    });

    va = va_new;
    vb = vb_new;
    G = G_new;
    ve = ve_new;
    beta = beta_new;
    used = it;
    if (ch < tol) {
      converged = true;
      break;
    }
  }

  mat U(K, q, arma::fill::zeros);
  mat Ginv = safe_inv_spd(G);
  vec resid = y - X * beta;
  int offset = 0;
  for (int k = 0; k < K; ++k) {
    int m = sizes[k];
    mat Xg = X.rows(offset, offset + m - 1);
    vec rg = resid.subvec(offset, offset + m - 1);
    mat Minv = safe_inv_spd(ve * Ginv + Xg.t() * Xg);
    vec ug = Minv * (Xg.t() * rg);
    U.row(k) = ug.t();
    offset += m;
  }

  RSFit ans;
  ans.beta = beta;
  ans.G = G;
  ans.ve = ve;
  ans.information = rs_information(X, G, ve, sizes);
  ans.uhat = U;
  ans.converged = converged;
  ans.iterations = used;
  return ans;
}

Rcpp::List ri_to_list(const RIFit& fit) {
  return Rcpp::List::create(
    Rcpp::Named("beta") = fit.beta,
    Rcpp::Named("Var.a") = fit.va,
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("Information") = fit.information,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

Rcpp::List rs_to_list(const RSFit& fit) {
  return Rcpp::List::create(
    Rcpp::Named("beta") = fit.beta,
    Rcpp::Named("G") = fit.G,
    Rcpp::Named("Var.a") = fit.G(0, 0),
    Rcpp::Named("Var.b") = rs_var_b(fit.G),
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("Information") = fit.information,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List fit_ri_lmm_cpp(const arma::mat& xx,
                          const arma::vec& yy,
                          Rcpp::IntegerVector nc,
                          double tol = 1e-6,
                          int max_iter = 100) {
  return ri_to_list(fit_ri_internal(xx, yy, nc, tol, max_iter));
}

// [[Rcpp::export]]
Rcpp::List fit_rs_lmm_cpp(const arma::mat& xx,
                          const arma::vec& yy,
                          Rcpp::IntegerVector nc,
                          double tol = 1e-6,
                          int max_iter = 100) {
  return rs_to_list(fit_rs_internal(xx, yy, nc, tol, max_iter));
}

// Full empirical IMSPE for the random-intercept model.
// For a target point f(x) in region k, prediction is f(x)' beta + a_k.
// The empirical target summaries are:
//   W   = E[f f'],
//   m_k = E[f 1{x in region k}],
//   p_k = P(x in region k).
// [[Rcpp::export]]
Rcpp::List imspe_ri_cpp(const arma::mat& xx,
                        const arma::vec& yy,
                        Rcpp::IntegerVector nc,
                        const arma::mat& W,
                        Rcpp::List mk_list,
                        const arma::vec& pk,
                        double tol = 1e-6,
                        int max_iter = 100) {
  RIFit fit = fit_ri_internal(xx, yy, nc, tol, max_iter);
  mat X = add_intercept(xx);
  std::vector<int> sizes = positive_sizes(nc);
  const int K = static_cast<int>(sizes.size());
  const int q = static_cast<int>(X.n_cols);

  if (static_cast<int>(W.n_rows) != q || static_cast<int>(W.n_cols) != q)
    Rcpp::stop("W has incompatible dimension.");
  if (mk_list.size() != K)
    Rcpp::stop("mk_list length must equal the number of non-empty groups.");
  if (static_cast<int>(pk.n_elem) != K)
    Rcpp::stop("pk length must equal the number of non-empty groups.");

  mat A = X.t() * X / fit.ve;
  std::vector<vec> B(K);
  std::vector<double> Hinv(K);
  int offset = 0;
  for (int k = 0; k < K; ++k) {
    int m = sizes[k];
    mat Xg = X.rows(offset, offset + m - 1);
    B[k] = arma::sum(Xg, 0).t() / fit.ve;
    Hinv[k] = 1.0 / (1.0 / fit.va + static_cast<double>(m) / fit.ve);
    offset += m;
  }

  mat M = A;
  for (int k = 0; k < K; ++k)
    M -= Hinv[k] * (B[k] * B[k].t());
  M = 0.5 * (M + M.t());
  mat C11 = safe_inv_spd(M);

  double fixed_term = arma::trace(C11 * W);
  double cross_term = 0.0;
  double random_term = 0.0;

  for (int k = 0; k < K; ++k) {
    vec mk = Rcpp::as<vec>(mk_list[k]);
    if (static_cast<int>(mk.n_elem) != q)
      Rcpp::stop("A m_k vector has incompatible dimension.");

    vec C12k = -C11 * B[k] * Hinv[k];
    double C22kk = Hinv[k] +
      Hinv[k] * Hinv[k] * arma::as_scalar(B[k].t() * C11 * B[k]);
    cross_term += 2.0 * arma::dot(C12k, mk);
    random_term += C22kk * pk[k];
  }

  double imspe = fit.ve + fixed_term + cross_term + random_term;

  return Rcpp::List::create(
    Rcpp::Named("IMSPE") = imspe,
    Rcpp::Named("objective") = -imspe,
    Rcpp::Named("residual_term") = fit.ve,
    Rcpp::Named("fixed_term") = fixed_term,
    Rcpp::Named("cross_term") = cross_term,
    Rcpp::Named("random_term") = random_term,
    Rcpp::Named("beta") = fit.beta,
    Rcpp::Named("Var.a") = fit.va,
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("Information") = M,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

// Full empirical IMSPE corresponding to Eq. (10) in the paper for the
// random-intercept + random-slope model.  W_k must be scaled by n_ref, so
// sum_k W_k = W (up to empty numerical regions).
// [[Rcpp::export]]
Rcpp::List imspe_rs_cpp(const arma::mat& xx,
                        const arma::vec& yy,
                        Rcpp::IntegerVector nc,
                        const arma::mat& W,
                        Rcpp::List Wk_list,
                        double tol = 1e-6,
                        int max_iter = 100) {
  RSFit fit = fit_rs_internal(xx, yy, nc, tol, max_iter);
  mat X = add_intercept(xx);
  std::vector<int> sizes = positive_sizes(nc);
  const int K = static_cast<int>(sizes.size());
  const int q = static_cast<int>(X.n_cols);

  if (static_cast<int>(W.n_rows) != q || static_cast<int>(W.n_cols) != q)
    Rcpp::stop("W has incompatible dimension.");
  if (Wk_list.size() != K)
    Rcpp::stop("Wk_list length must equal the number of non-empty groups.");

  mat Ginv = safe_inv_spd(fit.G);
  std::vector<mat> B(K), Hinv(K);
  mat A(q, q, arma::fill::zeros);
  int offset = 0;
  for (int k = 0; k < K; ++k) {
    int m = sizes[k];
    mat Xg = X.rows(offset, offset + m - 1);
    mat XtX = Xg.t() * Xg;
    A += XtX / fit.ve;
    B[k] = XtX / fit.ve; // fixed and random designs are both Xg
    Hinv[k] = safe_inv_spd(Ginv + XtX / fit.ve);
    offset += m;
  }

  mat M = A;
  for (int k = 0; k < K; ++k)
    M -= B[k] * Hinv[k] * B[k].t();
  M = 0.5 * (M + M.t());
  mat C11 = safe_inv_spd(M);

  double fixed_term = arma::trace(C11 * W);
  double cross_term = 0.0;
  double random_term = 0.0;

  for (int k = 0; k < K; ++k) {
    mat Wk = Rcpp::as<mat>(Wk_list[k]);
    if (static_cast<int>(Wk.n_rows) != q || static_cast<int>(Wk.n_cols) != q)
      Rcpp::stop("A W_k matrix has incompatible dimension.");

    mat C12k = -C11 * B[k] * Hinv[k];
    mat C22kk = Hinv[k] + Hinv[k] * B[k].t() * C11 * B[k] * Hinv[k];
    cross_term += 2.0 * arma::trace(C12k * Wk);
    random_term += arma::trace(C22kk * Wk);
  }

  double classical_I = fixed_term;
  double imspe = fit.ve + fixed_term + cross_term + random_term;

  return Rcpp::List::create(
    Rcpp::Named("IMSPE") = imspe,
    Rcpp::Named("objective") = -imspe,
    Rcpp::Named("residual_term") = fit.ve,
    Rcpp::Named("fixed_term") = fixed_term,
    Rcpp::Named("cross_term") = cross_term,
    Rcpp::Named("random_term") = random_term,
    Rcpp::Named("classical_I") = classical_I,
    Rcpp::Named("beta") = fit.beta,
    Rcpp::Named("G") = fit.G,
    Rcpp::Named("Var.a") = fit.G(0,0),
    Rcpp::Named("Var.b") = rs_var_b(fit.G),
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("Information") = M,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

// -------------------------------------------------------------------------
// Backward-compatible wrappers.  They now use convergence-based fits.
// -------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List count_info_cpp(const arma::mat& xx,
                          const arma::vec& yy,
                          Rcpp::IntegerVector nc,
                          int R,
                          int p,
                          double tol = 1e-6,
                          int max_iter = 100) {
  RIFit fit = fit_ri_internal(xx, yy, nc, tol, max_iter);
  double detv = arma::det(fit.information);
  if (!arma::is_finite(detv) || detv <= 0.0) detv = 1e-50;
  double Aval = arma::trace(safe_inv_spd(fit.information));
  return Rcpp::List::create(
    Rcpp::Named("D") = detv,
    Rcpp::Named("A") = Aval,
    Rcpp::Named("Information") = fit.information,
    Rcpp::Named("beta") = fit.beta,
    Rcpp::Named("Var.a") = fit.va,
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

// [[Rcpp::export]]
Rcpp::List count_info_rs_cpp(const arma::mat& xx,
                             const arma::vec& yy,
                             Rcpp::IntegerVector nc,
                             int R,
                             int p,
                             double tol = 1e-6,
                             int max_iter = 100) {
  RSFit fit = fit_rs_internal(xx, yy, nc, tol, max_iter);
  double detv = arma::det(fit.information);
  if (!arma::is_finite(detv) || detv <= 0.0) detv = 1e-50;
  double Aval = arma::trace(safe_inv_spd(fit.information));
  return Rcpp::List::create(
    Rcpp::Named("D") = detv,
    Rcpp::Named("A") = Aval,
    Rcpp::Named("Information") = fit.information,
    Rcpp::Named("beta") = fit.beta,
    Rcpp::Named("G") = fit.G,
    Rcpp::Named("Var.a") = fit.G(0,0),
    Rcpp::Named("Var.b") = rs_var_b(fit.G),
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

// [[Rcpp::export]]
Rcpp::List Est_hat_cpp(const arma::mat& xx,
                       const arma::vec& yy,
                       const arma::vec& beta_true,
                       double Var_a_true,
                       double Var_e_true,
                       Rcpp::IntegerVector nc,
                       int R,
                       int p,
                       double tol = 1e-6,
                       int max_iter = 100) {
  RIFit fit = fit_ri_internal(xx, yy, nc, tol, max_iter);
  vec slope = fit.beta.subvec(1, p);
  double bt_mse = arma::accu(arma::square(slope - beta_true));
  double bt0_mse = std::pow(fit.beta[0] - 1.0, 2.0);
  double va_mse = std::pow(fit.va - Var_a_true, 2.0);
  double ve_mse = std::pow(fit.ve - Var_e_true, 2.0);
  return Rcpp::List::create(
    Rcpp::Named("bt.mse") = bt_mse,
    Rcpp::Named("va.mse") = va_mse,
    Rcpp::Named("ve.mse") = ve_mse,
    Rcpp::Named("bt0.mse") = bt0_mse,
    Rcpp::Named("beta2") = fit.beta,
    Rcpp::Named("Var.a") = fit.va,
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}

// [[Rcpp::export]]
Rcpp::List Est_hat_RS_cpp(const arma::mat& xx,
                          const arma::vec& yy,
                          const arma::vec& beta_true,
                          double Var_a_true,
                          double Var_e_true,
                          Rcpp::IntegerVector nc,
                          int R,
                          int p,
                          double tol = 1e-6,
                          int max_iter = 100) {
  RSFit fit = fit_rs_internal(xx, yy, nc, tol, max_iter);
  vec slope = fit.beta.subvec(1, p);
  double bt_mse = arma::accu(arma::square(slope - beta_true));
  double bt0_mse = std::pow(fit.beta[0] - 1.0, 2.0);
  double va_mse = std::pow(fit.G(0,0) - Var_a_true, 2.0);
  double ve_mse = std::pow(fit.ve - Var_e_true, 2.0);
  return Rcpp::List::create(
    Rcpp::Named("bt.mse") = bt_mse,
    Rcpp::Named("va.mse") = va_mse,
    Rcpp::Named("ve.mse") = ve_mse,
    Rcpp::Named("bt0.mse") = bt0_mse,
    Rcpp::Named("beta2") = fit.beta,
    Rcpp::Named("Var.a") = fit.G(0,0),
    Rcpp::Named("Var.b") = rs_var_b(fit.G),
    Rcpp::Named("Var.e") = fit.ve,
    Rcpp::Named("G.hat") = fit.G,
    Rcpp::Named("u_hat") = fit.uhat,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("iterations") = fit.iterations
  );
}
