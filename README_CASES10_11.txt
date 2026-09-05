CASE 10 / CASE 11 EXTENSION
===========================

Files to place together in /users/k21181837/RGSS:
  ROG_aligned_fixed.R
  lmm_fast_aligned.cpp
  ROG_method_comparison.R
  ROG_cases10_11.R
  run_case10.R
  run_case11.R
  case10_job.sh
  case11_job.sh
  smoke_test_cases10_11.R
  submit_cases10_11.sh

Before long jobs:
  Rscript smoke_test_cases10_11.R

Case 10 design
--------------
- Default: 6 latent predictive types.
- 2 geometric X modes per predictive type -> 12 geometric modes.
- 5 observed administrative units per mode -> 60 observed units.
- Random effects are shared at the 6-type level.
- OBS uses the 60 over-split administrative units.
- ORACLE uses the 6 latent predictive types.
- KM/GMM optimize covariate clustering; CIRG chooses K via model-matched IMSPE
  after RASC clustering and uses SGA for test prediction.
- CPF and BLM use repeated-measure administrative-unit IDs.

This intentionally creates a prediction-vs-geometry distinction. A geometric
criterion can prefer 12 modes even when prediction benefits from pooling each
paired mode into one predictive type.

Case 11 design
--------------
- Default: 100 repeated-measure units.
- Each unit has continuous z in [-1,1]. There is no true discrete K.
- X distributions vary smoothly/nonlinearly with z.
- Random intercepts (and optional slopes) vary smoothly with z.
- OBS uses 4 coarse groups constructed from a noisy proxy for z.
- UNIT uses the original repeated-measure unit IDs and is a strong fine-grained
  benchmark, not a latent-type oracle.
- ORACLE uses the true continuous random effects only to provide an unattainable
  lower bound near the irreducible noise level.
- z_R2 reports how much continuous z heterogeneity is explained by a method's
  finite grouping; unit_consistency reports whether observations from one unit
  are assigned consistently.

Suggested first production runs
-------------------------------
Use nloop=10 before nloop=50:

  sbatch --export=ALL,NLOOP=10,VAR_B=0 case10_job.sh
  sbatch --export=ALL,NLOOP=10,VAR_B=0 case11_job.sh

Then RS:
  sbatch --export=ALL,NLOOP=10,VAR_B=0.1 case10_job.sh
  sbatch --export=ALL,NLOOP=10,VAR_B=0.1 case11_job.sh

GMM needs mclust. If unavailable:
  export METHODS=ORACLE,OBS,LM,KM,CPF,BLM,CIRG        # Case 10
  export METHODS=ORACLE,UNIT,OBS,LM,KM,CPF,BLM,CIRG   # Case 11

Important interpretation
------------------------
Case 10 has a true predictive partition, so group_ARI is meaningful.
Case 11 deliberately has no true discrete partition, so group_ARI is NA; use
MSPE, selected_K, z_R2 and unit_consistency instead.

Default tau
-----------
The current scripts use tau = 5*(p+1) by default, matching the main simulation
setup. Set TAU explicitly in the environment for sensitivity runs.
