#!/bin/bash
# Full sensitivity grid. Run only after the 25% core grid looks correct.
for model in RI RS; do
  for case in {1..9}; do
    for rho in 0.10 0.25 0.50; do
      sbatch --export=ALL,MODEL=$model,CASE=$case,MIS_TYPE=contam,RHO=$rho compare_job.sh
    done
    sbatch --export=ALL,MODEL=$model,CASE=$case,MIS_TYPE=merge,RHO=0 compare_job.sh
    sbatch --export=ALL,MODEL=$model,CASE=$case,MIS_TYPE=split,RHO=0 compare_job.sh
  done
done
