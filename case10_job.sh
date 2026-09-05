#!/bin/bash
#SBATCH --job-name=case10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=CASE10_%j.txt

set -euo pipefail
module purge
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6
cd "${PROJECT_DIR:-/users/k21181837/RGSS}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
Rscript run_case10.R
