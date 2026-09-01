# CIRG/RASC code alignment notes

Generated files:

- `lmm_fast_aligned.cpp`: convergence-based RI/RS mixed-model fitting, BLUPs, and full empirical IMSPE.
- `ROG_aligned.R`: RASC, SGA, independent reference sample, bidirectional SA, corrected simulation DGP, prediction and reporting.
- `smoke_test_aligned.R`: a small one-replication run for compilation/runtime validation on a machine with R.

The original `ROG.R`, `lmm_fast.cpp`, and `myoss.cpp` are left unchanged so old results remain reproducible.

## Method mapping

### Full IMSPE

The legacy search optimized only `-Tr(M^{-1}W)`. The aligned search calls `imspe_rs_cpp()`, which computes

`Var.e + Tr(C11 W) + sum_k {2 Tr(C12_k W_k) + Tr(C22_kk W_k)}`

and maximizes its negative.

### RASC

`rasc_cluster()` fits OLS, standardizes the residual, forms `[X, lambda*r_std]`, runs MiniBatch K-means, and decreases K until every group satisfies `tau` (default `5*(p+1)`). It also estimates `(mu_k, Sigma_k, pi_k)` from the original X coordinates.

### SGA

`soft_group_weights()` evaluates Gaussian posterior weights. `predict_ri_soft()` and `predict_rs_soft()` use these weights. The RS predictor evaluates the random-effect vector through the design vector, so the random contribution is `f(x)' u_soft` rather than a scalar intercept.

### SA

`cirg_search()` samples steps from `{-3,-2,-1,+1,+2,+3}` and retains a historical best state.

### Estimation

The C++ RI and RS fits iterate until a relative-change tolerance or `max_iter`. These are ML-EM updates, not exact REML. The old fixed three-iteration loops are not used in the aligned pipeline.

### Simulation fixes

- True train/test observations in the same original group share the same random intercept/slope.
- A separate reference X sample estimates W/W_k; the test X is untouched until final evaluation.
- Case 5 defaults to t(10), with df=3 available as a sensitivity setting.
- Case 6 uses `sdlog=sqrt(0.2)` for `log(X) ~ N(0,0.2)`.
- Random-slope variance defaults to 0.1 in `Comp_RS()`.
- Case 2 keeps the intended legacy equicorrelation matrix (diag=1, off-diagonal=0.5); the manuscript formula should be written explicitly to match this.

## One remaining manuscript-level issue

Equation (10) uses hard regions `chi_k`, while Algorithm 2 predicts with soft posterior weights. The aligned code uses MAP posterior regions to construct `W_k` for Eq. (10), then uses soft posterior weights for prediction. This is faithful to the current equation plus SGA algorithm, but the manuscript sentence claiming the optimized and reported quantities are literally identical should be weakened or the theory should be extended to a soft-assignment IMSPE containing cross-group `C22_{kl}` terms.
