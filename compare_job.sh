#!/bin/bash -l
#SBATCH --output=/users/k21181837/RGSS/COMPARE_%A_%a.txt
#SBATCH --job-name=RGCOMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
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

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  task_id=$((SLURM_ARRAY_TASK_ID - 1))
  model_block=$((task_id / 48))
  within_model=$((task_id % 48))
  case_id=$((within_model / 3 + 1))
  r_index=$((within_model % 3))
  r_values=(10 20 50)

  if (( model_block == 0 )); then
    export MODEL=RI
  else
    export MODEL=RS
  fi
  export CASE=$case_id
  export R_LIST=${r_values[$r_index]}
  export VAR_A_LIST=0.5,2.25
fi

echo "MODEL=${MODEL:-RI} CASE=${CASE:-1} R_LIST=${R_LIST:-default} VAR_A_LIST=${VAR_A_LIST:-default}"
echo "LABEL_POLICY=${LABEL_POLICY:-unknown} MIS_TYPE=${MIS_TYPE:-none} RHO=${RHO:-0.25} LAMBDA=${LAMBDA:-0} CIRG_CRITERION=${CIRG_CRITERION:-IMSPE}"
echo "NLOOP=${NLOOP:-20} SA_MAX=${SA_MAX:-25} EM_MAX=${EM_MAX:-300} K_GRID_POINTS=${K_GRID_POINTS:-8} N_CORES=${N_CORES:-${SLURM_CPUS_PER_TASK:-1}}"
Rscript "$R_DIR/run_compare.R"
