RGSS / CIRG METHOD-COMPARISON EXTENSION
=======================================

New files
---------
ROG_method_comparison.R        main comparison methods
run_compare.R                  generic runner controlled by environment variables
compare_job.sh                 generic KCL Slurm job
smoke_test_comparison.R        quick local/HPC smoke test
submit_compare_core.sh         108 short unknown-group jobs:
                                RI/RS x cases 1-9 x R 10/20/50 x Var.a 0.5/2.25
submit_compare_misspec_grid.sh disabled for the main unknown-group comparison

Methods
-------
LM      pooled linear model ignoring groups.
KM      X-only MiniBatch K-means + LMM; K selected by a BIC-like criterion.
OBS     LMM using random pseudo labels when LABEL_POLICY=unknown.
CPF     CPF-style pairwise fusion, initialized from random pseudo labels when
        LABEL_POLICY=unknown.
BLM     BLM-style two-step discretization, initialized from random pseudo labels
        when LABEL_POLICY=unknown.
CIRG    proposed RASC + SGA, RI/RS-matched IMSPE regrouping.
        Training clusters are built on [X, lambda * standardized residual].
        Test predictions use soft Gaussian posterior assignment from X.

Group-label policy
------------------
The main paper premise is that group information is unknown. Therefore the
default runner uses LABEL_POLICY=unknown. No method receives true group labels.
Methods that normally require observed labels are fed balanced random pseudo
labels, so they remain runnable without leaking hidden grouping information.

ORACLE requires true group labels and is not a valid main-comparison competitor
under LABEL_POLICY=unknown:

ORACLE  correctly specified true-group LMM. Infeasible benchmark only.

GMM was removed from the main runner to keep the code simple and avoid optional
package failures. Re-add it only as a separate diagnostic if needed.

Estimator note
--------------
The C++ mixed-model routines are convergence-based ML-EM estimators, not exact
REML optimizers. Manuscript text should say ML-EM/convergent EM unless a REML
implementation is added.

IMPORTANT METHOD-NAMING NOTE
----------------------------
CPF and BLM are LMM-compatible adaptations, not literal executions of the
original authors' model-specific programs. In a paper/table call them
"CPF-style" / "BLM-style" unless you separately reproduce their exact original
models. This code deliberately holds the final LMM estimator fixed so the
comparison mainly concerns how the grouping is constructed.

Observed-label misspecification
-------------------------------
Only used when LABEL_POLICY=observed.

MIS_TYPE=correct  true labels
MIS_TYPE=contam   permutes a fraction RHO of observation labels
MIS_TYPE=merge    merges adjacent true groups (factor 2 by default)
MIS_TYPE=split    splits each true group in two using a training X1 threshold

First smoke test
----------------
cd /users/k21181837/RGSS
Rscript smoke_test_comparison.R

One real job
------------
sbatch --export=ALL,LABEL_POLICY=unknown,MODEL=RI,CASE=2,R_LIST=20,VAR_A_LIST=2.25,MIS_TYPE=none,RHO=0,LAMBDA=0,NLOOP=20 compare_job.sh

Useful method subsets
---------------------
# Small main subset:
sbatch --export=ALL,LABEL_POLICY=unknown,MODEL=RI,CASE=2,R_LIST=20,VAR_A_LIST=2.25,MIS_TYPE=none,RHO=0,LAMBDA=0,NLOOP=20,METHODS=LM,KM,CIRG compare_job.sh

# Include pseudo-label versions of the label-dependent competitors:
sbatch --export=ALL,LABEL_POLICY=unknown,MODEL=RI,CASE=2,R_LIST=20,VAR_A_LIST=2.25,MIS_TYPE=none,RHO=0,LAMBDA=0,NLOOP=20,METHODS=LM,KM,OBS,CPF,BLM,CIRG compare_job.sh

Outputs
-------
run_compare.R saves one RDS per setting and one combined RDS. By default it
runs R_LIST=10,20,50 and VAR_A_LIST=0.5,2.25 unless R or VAR_A is explicitly
provided. The HPC submit script passes one R_LIST and one VAR_A_LIST per job so
each job stays short. Override with comma-separated CASE_LIST, R_LIST,
VAR_A_LIST, METHODS.
LAMBDA defaults to 0 in the unknown-group main comparison. Positive lambda
values are useful sensitivity checks, but they can make training clusters depend
on residual information that is unavailable for new test points. Set
CIRG_CRITERION=I to run the direct I-optimal search; the default is IMSPE.

Speed defaults
--------------
NLOOP defaults to 20, SA_MAX to 25, EM_MAX to 300, K_GRID_POINTS to 8,
KMEANS_NUM_INIT to 1, KMEANS_MAX_ITERS to 25, CPF_LAMBDA_POINTS to 6,
CPF_MAX_ITER to 300, and BLM_NSTART to 5. On Slurm, run_compare.R uses
SLURM_CPUS_PER_TASK as N_CORES by default and parallelizes replications on
Linux. Increase these values only after the fast run finishes cleanly.

Example:
results_comparison/RI_case2_R20_vara2.25_lambda0_none.rds
results_comparison/RI_cases2_R20_vara2.25_grid_lambda0_none.rds

Quick lambda diagnostic
-----------------------
Run lambda_diagnostic.R for a paired local comparison of CIRG with LAMBDA_LIST
values, defaulting to 0, 0.25, 0.5, and 1:

  Rscript lambda_diagnostic.R

Recommended sequence
--------------------
1. smoke_test_comparison.R
2. one RI case2 job at nloop=5
3. one RS case2 job at nloop=5
4. submit_compare_core.sh
5. only after inspecting those results, expand R_LIST, VAR_A_LIST, NLOOP, or
   CIRG_CRITERION for full label-free experiments.
