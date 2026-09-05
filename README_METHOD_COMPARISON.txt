RGSS / CIRG METHOD-COMPARISON EXTENSION
=======================================

New files
---------
ROG_method_comparison.R        main comparison methods
run_compare.R                  generic runner controlled by environment variables
compare_job.sh                 generic KCL Slurm job
smoke_test_comparison.R        quick local/HPC smoke test
install_comparison_packages.R  installs optional mclust package for GMM
submit_compare_core.sh         18 jobs: RI/RS x cases 1-9, 25% contamination
submit_compare_misspec_grid.sh 90-job sensitivity grid (run later)

Methods
-------
ORACLE  correctly specified true-group LMM. Upper benchmark only.
OBS     LMM using deliberately misspecified observed labels.
LM      pooled linear model ignoring groups.
KM      X-only MiniBatch K-means + LMM; K selected by a BIC-like criterion.
GMM     model-based Gaussian mixture clustering of X + LMM (mclust).
CPF     Ma-Huang-style MCP pairwise-fusion adaptation for repeated grouped data.
        Pairwise fusion acts on observed-group intercepts, tuning lambda by the
        modified-BIC form used in the pairwise-fusion literature. The final
        partition is evaluated with exactly the same RI/RS LMM as other methods.
BLM     Bonhomme-Lamadon-Manresa-style two-step discretization. First estimate
        informative group-level random-effect moments; second classify these
        moments by K-means; third fit the common RI/RS LMM.
CIRG    proposed corrected X-only, RI/RS-matched IMSPE regrouping.
        The obsolete residual-weight parameter lambda has been removed entirely.

IMPORTANT METHOD-NAMING NOTE
----------------------------
CPF and BLM are LMM-compatible adaptations, not literal executions of the
original authors' model-specific programs. In a paper/table call them
"CPF-style" / "BLM-style" unless you separately reproduce their exact original
models. This code deliberately holds the final LMM estimator fixed so the
comparison mainly concerns how the grouping is constructed.

Observed-label misspecification
-------------------------------
MIS_TYPE=correct  true labels
MIS_TYPE=contam   permutes a fraction RHO of observation labels
MIS_TYPE=merge    merges adjacent true groups (factor 2 by default)
MIS_TYPE=split    splits each true group in two using a training X1 threshold

First smoke test
----------------
cd /users/k21181837/RGSS
Rscript smoke_test_comparison.R

Install GMM dependency once if needed
-------------------------------------
Rscript install_comparison_packages.R

One real job
------------
sbatch --export=ALL,MODEL=RI,CASE=2,MIS_TYPE=contam,RHO=0.25 compare_job.sh

Useful method subsets
---------------------
# Exclude GMM if mclust is not installed:
sbatch --export=ALL,MODEL=RI,CASE=2,MIS_TYPE=contam,RHO=0.25,METHODS=ORACLE,OBS,LM,KM,CPF,BLM,CIRG compare_job.sh

# Only expensive/important competitors:
sbatch --export=ALL,MODEL=RI,CASE=2,MIS_TYPE=contam,RHO=0.25,METHODS=ORACLE,OBS,CPF,BLM,CIRG compare_job.sh

Outputs
-------
results_comparison/RI_case2_contam_rho25.rds
results_comparison/RS_case2_contam_rho25.rds
etc.

Recommended sequence
--------------------
1. smoke_test_comparison.R
2. one RI case2 job at nloop=5: add NLOOP=5 via --export
3. one RS case2 job at nloop=5
4. submit_compare_core.sh (25% contamination)
5. only after inspecting those results, submit the full misspecification grid.
