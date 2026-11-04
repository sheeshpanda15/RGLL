# Main Experiment

Only the main CPU study is supported. The separate CPFO experiment, its solver
and all CUDA/build/submission paths have been removed. CPF is a different,
group-level adaptation and remains in the main comparison.

## Upload

Finish or cancel old jobs before replacing sources in their working directory.
In particular, cancel the abandoned CPFO job if it is still running:

```bash
scancel 38111202
```

Copy these nine root-level runtime files to ~/RGSS (or a clean project folder):

```text
ROG_aligned_fixed.R
ROG_method_comparison.R
literature_baselines.R
lmm_fast_aligned.cpp
experiment_plan.R
run_compare.R
compare_job.sh
submit_compare_core.sh
collect_results.R
```

RUNNING.md is optional. No archive, result files, compiled libraries or caches
are needed. In an existing remote folder, the obsolete files can be removed:

```bash
rm -f submit_cpfo.sh submit_compare_gpu.sh build_cpf_gpu.sh cpf_original.cpp cpf_cuda.cu
```

Existing result directories are not deleted. New runs default to
result/results_main so that outputs from different source versions are not
silently mixed with result/results_main_no_checkpoint. The two previously
completed jobs are not automatically imported: a fresh run submits all 64 jobs.

## Submit

```bash
cd ~/RGSS
bash submit_compare_core.sh
```

This submits 64 independent ordinary Slurm jobs:

- RI and RS x CASE1-CASE16 x Var.a=0.5/2.25.
- Each job fixes model, CASE and Var.a, runs R=10/20/50, 50 repetitions each.
- Methods: LM, KM, OBS, CPF, BLMX, BLM and CIRG.
- Eight CPUs, 100 GB RAM and 48 hours per job.
- Each job has its own log and END/FAIL/TIME_LIMIT mail to k21181837@kcl.ac.uk.

BLMX is recorded as not_applicable because informative individual moments are
unavailable in the unknown-identity design. It is not a successful fitted method.
True labels are used for evaluation only, not as inputs to these seven methods.
Data generators, fitting settings and replication seed rules are unchanged.

There are no job arrays, dependency chains, test jobs, controllers, automatic
resubmissions or code-level queue concurrency caps. Slurm schedules accepted
jobs when resources and account policy permit. There are no checkpoints;
each job saves one compressed final RDS only after its whole workload finishes.
An interrupted job must be rerun in full. A 48-hour limit is not a runtime guarantee.

## Environment And Submission Errors

The worker loads r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6.
If loading fails, it retries the same module with the Lmod cache disabled; it
does not silently substitute a different R version. R_MODULE can override it.
Required R packages are Rcpp, RcppArmadillo, ClusterR and MASS.

The previous failure on erc-hpc-vm-zz4qbq occurred before R started. To avoid
that known problematic node while its module setup is investigated:

```bash
CPU_EXCLUDE=erc-hpc-vm-zz4qbq bash submit_compare_core.sh
```

Use either this command or the standard submission above, not both.
CPU_PARTITION, MAIL_USER and RESULT_DIR are optional overrides.

If sbatch reports a socket timeout, the job may still have been accepted.
The script stops without retrying; check squeue and sacct before continuing.
Job numbering remains Var.a first, then CASE, then model. For example:

```bash
bash submit_compare_core.sh 4
bash submit_compare_core.sh 3 64
```

The second command is appropriate only if jobs 1 and 2 in this same run were
accepted and no others were. Rerunning a task replaces its final file only on
success; end previous attempts first.

## Collect

After all 64 jobs have completed successfully:

```bash
cd ~/RGSS
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6
Rscript collect_results.R result/results_main
```

Collection does not load packages, compile C++ or fit models. It checks all 64
files, identical configurations, task metadata, method/setting/replication
coverage and duplicate rows before writing result/results_main/main_results.rds.
At 50 repetitions the output has 67200 raw rows, including 9600 BLMX
not_applicable rows, plus per-setting summaries and Monte Carlo standard errors.
File completeness does not imply all fits converged: inspect status and
convergence rates before interpreting the results.

An older complete and internally consistent set of 64 final BASE files can
still be collected by passing its directory. Do not mix old and new source
fingerprints in one run. Old per-setting checkpoint files are not imported.
