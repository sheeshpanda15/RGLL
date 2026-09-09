#!/bin/bash
# Main unknown-group comparison. OBS/CPF/BLM receive random pseudo labels.
# LAMBDA=0 is explicit so stale shell environment values cannot leak in.
# Each job is one model/case/R/Var.a setting, keeping wall time short.
COMMON="LABEL_POLICY=unknown,MIS_TYPE=none,RHO=0,LAMBDA=0,NLOOP=20"
COMMON="$COMMON,SA_MAX=25,EM_MAX=300,K_GRID_POINTS=8"
COMMON="$COMMON,KMEANS_NUM_INIT=1,KMEANS_MAX_ITERS=25"
COMMON="$COMMON,CPF_LAMBDA_POINTS=6,CPF_MAX_ITER=300,BLM_NSTART=5"

for model in RI RS; do
  for case in {1..9}; do
    for r in 10 20 50; do
      for vara in 0.5 2.25; do
        sbatch --export=ALL,MODEL=$model,CASE=$case,R_LIST=$r,VAR_A_LIST=$vara,$COMMON compare_job.sh
      done
    done
  done
done
