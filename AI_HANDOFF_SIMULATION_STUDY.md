# CIRG Simulation Study: AI Handoff Summary

## 1. Research objective

This project evaluates **predictive ability**, with test-set mean squared prediction error (MSPE) as the primary metric.

The proposed method is named `CIRG` in the code (the conversation occasionally used `CRIG`). Its goal is **not** to recover the true data-generating groups. The scientific claim should instead be framed as:

> CIRG searches for a grouping that is sufficiently good for prediction, even when that grouping differs from the latent data-generating partition.

Accordingly, adjusted Rand index or exact recovery of the true groups is secondary and should not be used as the main criterion for judging CIRG. The central question is whether the grouping selected by CIRG produces low out-of-sample MSPE.

## 2. Methods compared

The current experiment contains these feasible methods:

- `CIRG`: proposed RASC + SGA method, with the grouping selected by an IMSPE criterion.
- `KM`: covariate-only K-means followed by the common LMM prediction pipeline.
- `BLM`: a Bonhomme-Lamadon-Manresa-style two-step discretization adapted to the LMM setting.
- `CPF`: a concave pairwise-fusion method adapted to the LMM setting.
- `OBS`: LMM based on supplied proxy labels.
- `LM`: pooled linear model that ignores grouping.

There is also:

- `ORACLE`: an infeasible benchmark that uses the true data-generating group labels.

Under `LABEL_POLICY=unknown`, no feasible method receives the true group labels. `OBS`, `CPF`, and `BLM` use random pseudo-labels in the main comparison. Only `ORACLE` uses true labels, and every figure containing it explicitly identifies it as an infeasible benchmark.

## 3. Fixed simulation settings

The expanded simulation retained the original settings except for adding new data-generating cases:

- Sample size: `N = 2500`
- Number of covariates: `p = 50`
- Residual variance: `Var(e) = 9`
- Replications per setting: `20`
- Numbers of data-generating groups: `R = 10, 20, 50`
- Random-intercept variances: `Var(a) = 0.5, 2.25`
- Models:
  - `RI`: random intercept, `Var(b) = 0`
  - `RS`: random intercept and slope, `Var(b) = 0.1`
- Main label policy: `unknown`
- Misspecification type: `none`
- `lambda = 0`
- CIRG selection criterion: `IMSPE`
- Base seed: `12345`

There are `2 models x 16 cases x 3 R values x 2 Var(a) values = 192` settings. Each method has 20 Monte Carlo replications in each setting.

## 4. Sixteen data-generating cases

The figure titles and broad interpretations are:

1. **Independent uniform covariates**: independent `Uniform(-1, 1)` covariates.
2. **Correlated Gaussian covariates**: centered correlated Gaussian covariates.
3. **Group-shifted uniform covariates**: uniform covariates with group-dependent locations.
4. **Group-shifted Gaussian covariates**: correlated Gaussian covariates with group-dependent means.
5. **Heavy-tailed t covariates**: t-distributed covariates; the current default uses 10 degrees of freedom.
6. **Log-normal covariates**: log-normal covariates with `log(X)` variance 0.2.
7. **Two-component Gaussian mixture**: observations drawn from a two-component Gaussian mixture.
8. **Mean covariate shift**: train and target/test distributions differ in their means.
9. **Covariance shift**: train and target/test distributions differ in covariance.
10. **Type structure with nuisance modes**: group type is present in the covariates together with strong nuisance mixture modes.
11. **Structured latent predictive effects**: covariate structure and random effects are linked through a latent continuous score.
12. **Continuous predictive heterogeneity**: predictive effects vary continuously, so there is no prespecified finite target partition.
13. **Density-prediction conflict**: strong nuisance density modes conflict with a simpler prediction-relevant state.
14. **Overlapping predictive states**: prediction-relevant states overlap substantially in covariate space.
15. **Rare high-loss state**: a relatively rare state carries disproportionate predictive importance.
16. **Target-distribution shift**: the proportions of latent predictive states differ sharply between training and target/test data.

Cases 12-16 were designed specifically around **predictive grouping rather than true-group recovery**. Several true groups may share the same prediction-relevant state, and case 12 deliberately has continuous heterogeneity.

## 5. Main code changes

### Simulation and orchestration

- `ROG_method_comparison.R`
  - Contains comparison implementations and generators for cases 10-16.
  - Cases 12-16 explicitly target prediction-relevant grouping.
  - `ORACLE` is now permitted under the unknown-label experiment only as an explicitly infeasible benchmark.
  - Default method lists now include `ORACLE`.
- `run_compare.R`
  - Supports `CASE_LIST`, `R_LIST`, and `VAR_A_LIST`.
  - Uses the same fixed seed and settings across cases.
- `compare_job.sh`
  - Maps the Slurm array to model, case, and `R`.
- `submit_compare_core.sh`
  - Submits all 96 array tasks: `2 models x 16 cases x 3 R values`.
  - Each task runs both `Var(a)` values.

### Oracle backfill

- `add_oracle_benchmark.R`
  - Computes only `ORACLE` for the same 192 settings and the same 20 replication seeds.
  - It does not rerun or overwrite the six feasible methods.
  - Oracle results are saved separately under `result/results_oracle/`.

### Publication figures

- `plot_case_comparison.R`
  - Reads the 192 main result files and 192 Oracle-only result files.
  - Produces two complete figure versions: with and without Oracle.
  - Produces vector PDF output suitable for a paper.

## 6. Figure design

There is one figure for each of the 16 cases. Each figure contains four panels:

- RI with `Var(a) = 0.50`
- RI with `Var(a) = 2.25`
- RS with `Var(a) = 0.50`
- RS with `Var(a) = 2.25`

Within each panel:

- x-axis: `R = 10, 20, 50`
- y-axis: mean test MSPE
- error bars: mean `+/- 1.96 x Monte Carlo standard error`
- all methods are displayed simultaneously
- CIRG is emphasized with a thicker orange line

**Important:** each of the four panels has its own independent y-axis range. The earlier use of `facet_grid(..., scales = "free_y")` still shared scales within rows and compressed visible differences. It was replaced by `facet_wrap(..., scales = "free_y")` so all four panels are independently scaled.

## 7. Output locations

### Version without Oracle

Directory: `result/figures_without_oracle/`

- `case_01_mspe.pdf` through `case_16_mspe.pdf`
- `all_cases_mspe.pdf`: 16-page combined PDF
- `case_mspe_summary.csv`

This version is best for closely comparing the feasible methods because the Oracle line cannot stretch the y-axis.

### Version with Oracle

Directory: `result/figures_with_oracle/`

- `case_01_mspe.pdf` through `case_16_mspe.pdf`
- `all_cases_mspe.pdf`: 16-page combined PDF
- `case_mspe_summary.csv`

This version shows the gap to the infeasible true-group benchmark. Its caption states that Oracle uses true generating groups while all other methods do not.

### Raw results

- Main six-method results: `result/results_comparison/`
- Oracle-only results: `result/results_oracle/`

The Oracle directory has 192 RDS files, each containing 20 replications.

## 8. Reproduction commands

From the project root:

```bash
# Only needed if Oracle result files are missing.
N_CORES=8 Rscript add_oracle_benchmark.R

# Regenerate both PDF figure versions and their summary CSV files.
Rscript plot_case_comparison.R
```

On the Slurm system, the complete main simulation is submitted with:

```bash
bash submit_compare_core.sh
```

Do not rerun the main six methods merely to regenerate figures. The plotting script reads the existing RDS output.

## 9. Current result snapshot

Using the version without Oracle and selecting the numerically lowest mean MSPE in each of the 192 setting cells:

| Method | Number of cell wins | Mean rank |
|---|---:|---:|
| CIRG | 89 | 1.927 |
| KM | 56 | 2.125 |
| LM | 23 | 2.656 |
| BLM | 22 | 4.464 |
| OBS | 2 | 4.536 |
| CPF | 0 | 5.292 |

Useful case-level patterns:

- CIRG is strongest overall by both win count and mean rank.
- CIRG wins all 12 settings in case 12.
- CIRG wins 11 of 12 settings in cases 11 and 13.
- CIRG is also strong in cases 14-16, which were designed around prediction-relevant grouping.
- KM wins all 12 settings in case 10.
- LM is often competitive in simpler cases such as cases 1, 2, 5, and 6.

These are numerical winner counts, not formal significance claims. The Monte Carlo error bars should be considered when methods are close or intervals overlap.

## 10. Recommended paper framing

The main paper presentation should emphasize:

1. CIRG optimizes prediction-relevant grouping rather than latent partition recovery.
2. Cases 12-16 demonstrate settings where density clusters and prediction-relevant groups differ, overlap, are continuous, are rare, or undergo target shift.
3. The without-Oracle figures are the clearest primary comparison among feasible methods.
4. The with-Oracle figures can be supplementary or used to show the remaining gap to an infeasible benchmark.
5. Oracle must never be described as a competing implementable method.
6. Case 10 should be retained because it transparently shows a setting favoring KM; this helps demonstrate that the simulation suite is not constructed to make CIRG win universally.

## 11. Constraints for future edits

- Preserve all original simulation settings unless the user explicitly requests a change.
- Adding or changing cases should not silently change `N`, `p`, variances, replication count, tuning budgets, methods, seeds, or label policy.
- Keep Oracle separate from feasible methods in interpretation.
- Do not treat recovery of the true groups as the objective of CIRG.
- Do not overwrite the existing main RDS files when only an Oracle benchmark or new visualization is needed.
- Preserve the four independent y-axis scales in each case figure.
- The workspace may contain user-generated or cluster-generated changes; do not revert unrelated files.

## 12. Validation already performed

- Confirmed 192 main result settings.
- Confirmed 192 Oracle result files with 20 replications each.
- Confirmed each figure directory contains 16 individual PDFs plus one combined PDF.
- Confirmed both combined PDFs contain 16 pages.
- Visually inspected representative cases after enabling independent panel scales and after adding Oracle.
- R emitted harmless Fontconfig cache warnings during PDF creation, but all PDFs were generated successfully.
