#!/bin/bash
set -euo pipefail
# RI
sbatch --export=ALL,VAR_B=0 case10_job.sh
sbatch --export=ALL,VAR_B=0 case11_job.sh
# RS
sbatch --export=ALL,VAR_B=0.1 case10_job.sh
sbatch --export=ALL,VAR_B=0.1 case11_job.sh
