# Prediction-focused simulation cases

Cases 12--16 evaluate predictive grouping. The `R` data-generating groups
create repeated observations and shared random effects, but recovering those
identities is not the target. Primary comparisons should use test MSPE,
selected K, runtime, and paired Monte Carlo differences in MSPE.

## Case 12: continuous predictive heterogeneity

A continuous latent coordinate controls nonlinear covariate centers, random
intercepts, and random slopes. There is no prespecified finite predictive
partition; a selected K is a finite approximation chosen for prediction.

## Case 13: density-prediction conflict

Two weak predictive states control the random effects, while six stronger and
independent nuisance states control much of the covariate density. This checks
whether a method selects groups for response prediction rather than for density
reconstruction alone.

## Case 14: overlapping predictive states

Two prediction states have overlapping Gaussian covariate distributions. This
is a boundary case for comparing CIRG soft assignment with hard K-means
assignment. A near tie is informative when the covariates do not support a
confident assignment.

## Case 15: rare high-loss state

Approximately 15 percent of test observations belong to a covariate-overlapping
state with a large random effect. The case checks whether predictive grouping
retains a small region that matters disproportionately for squared prediction
loss.

## Case 16: target-distribution shift

Training observations have latent-state weights `(0.70, 0.20, 0.10)`, whereas
the independent reference sample and test observations use
`(0.10, 0.20, 0.70)`. This evaluates target-distribution-aware IMSPE selection.

## Running the cases

Submit the complete experiment (cases 1--16) with:

```bash
bash submit_compare_core.sh
```

The primary specification uses `LAMBDA=0`, matching the existing experiment.
Run nonzero lambda values as separately named sensitivity analyses; do not pick
lambda after inspecting test MSPE.
