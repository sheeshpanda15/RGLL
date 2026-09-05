#!/bin/bash
# One central misspecification setting first: 25% label contamination.
for model in RI RS; do
  for case in {1..9}; do
    sbatch --export=ALL,MODEL=$model,CASE=$case,MIS_TYPE=contam,RHO=0.25 compare_job.sh
  done
done
