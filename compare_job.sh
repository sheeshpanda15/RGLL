#!/bin/bash -l
#SBATCH --output=/users/k21181837/RGSS/COMPARE_%A.txt
#SBATCH --job-name=RGCOMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=k21181837@kcl.ac.uk

source /etc/profile 2>/dev/null || true
source /etc/profile.d/modules.sh 2>/dev/null || true
source /usr/share/lmod/lmod/init/bash 2>/dev/null || true
module load r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6

R_DIR=/users/k21181837/RGSS
cd "$R_DIR"
export PROJECT_DIR="$R_DIR"

echo "Host: $HOSTNAME"
echo "MODEL=${MODEL:-RI} CASE=${CASE:-1} MIS_TYPE=${MIS_TYPE:-contam} RHO=${RHO:-0.25}"
Rscript "$R_DIR/run_compare.R"
