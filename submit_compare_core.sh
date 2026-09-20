#!/bin/bash
# Main unknown-group comparison for cases 1--16. OBS/CPF/BLM receive random
# pseudo labels. LAMBDA=0 is explicit so stale shell values cannot leak in.
# Each task is one model/case/R setting and runs both Var.a values inside
# run_compare.R: 2 models x 16 cases x 3 R settings = 96 array tasks.
COMMON="LABEL_POLICY=unknown,MIS_TYPE=none,RHO=0,LAMBDA=0,NLOOP=20"
COMMON="$COMMON,SA_MAX=25,EM_MAX=300,K_GRID_POINTS=8"
COMMON="$COMMON,KMEANS_NUM_INIT=1,KMEANS_MAX_ITERS=25"
COMMON="$COMMON,CPF_LAMBDA_POINTS=6,CPF_MAX_ITER=300,BLM_NSTART=5"

sbatch --array=1-96%6 --export=ALL,$COMMON compare_job.sh
