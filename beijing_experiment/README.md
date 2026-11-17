# Beijing air-quality experiment

This folder reproduces the Beijing results in Table 2. It contains the original
station files, a snapshot of the fitting and preprocessing code, and the
reference summary from the runs used for the table.

## Run

Use R with a working C++ toolchain for `Rcpp::sourceCpp()`. The local
verification environment uses R 4.4.0, Rcpp 1.0.12, RcppArmadillo 0.12.8.4.0,
ClusterR 1.3.2, and MASS 7.3-60.0.2.

```r
install.packages(c("Rcpp", "RcppArmadillo", "ClusterR", "MASS"))
```

From this folder, run:

```sh
Rscript run_beijing.R
```

The defaults are 50 repetitions and four workers. Set `N_CORES` to change the
worker count. Windows uses one worker. Use `NLOOP=1` for a one-repetition
check, and `RESULT_DIR` to choose a separate output directory.

The `results/` directory contains the test scores for each split in `raw.csv`,
the means, Monte Carlo standard errors, and convergence counts in
`table2_summary.csv`, and the split indices, preprocessing parameters, and
session information in `beijing_results.rds`. Compare the full run to
`reference_table2_summary.csv` using the setting, model, and method columns.

## Design

- The outcome is `log1p` of daily mean PM2.5. Each station-day needs at least
  18 observed hours for the outcome and each weather covariate.
- Excl uses observed daily weather. Full adds annual-cycle and year terms.
  Both settings compare PGR, KM, BLM, CPF, OBS, and LM under RI and RS.
  Full also fits models grouped by station, season, or station and season.
- All methods share each split. The training, reference, and test sets use
  approximately 50%, 20%, and 30% of the 28-day calendar blocks. Two days at
  each end of a block are purged before splitting. Scaling uses training data.
- Within each setting, the table uses repetitions with finite predictions for
  all displayed methods in both models. Fits that reach the iteration limit
  with finite predictions are included and counted in the summary.

## Data attribution

Chen, S. (2017). *Beijing Multi-Site Air Quality* [Dataset]. UCI Machine
Learning Repository. [DOI: 10.24432/C5RK5G](https://doi.org/10.24432/C5RK5G).

UCI distributes the dataset under the
[Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
`data/` contains the 12 unmodified station CSVs, plus their source URL and
SHA-256 hashes in `source.json`. See the
[UCI dataset record](https://archive.ics.uci.edu/dataset/501/beijing+multi+site+air+quality+data).
