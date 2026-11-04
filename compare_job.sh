#!/bin/bash
#SBATCH --job-name=RGCOMP
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=k21181837@kcl.ac.uk

set +u
set -Eeo pipefail
trap 'status=$?; printf "[RGLL] EXIT status=%s job=%s\n" "$status" "${SLURM_JOB_ID:-local}"' EXIT
printf '[RGLL] START job=%s mode=base host=%s\n' "${SLURM_JOB_ID:-local}" "${HOSTNAME:-unknown}"
[[ "${RUN_MODE:-base}" == base ]] || { echo "Only main-study jobs are supported." >&2; exit 1; }
cd "${PROJECT_DIR:-/users/k21181837/RGSS}"
export PROJECT_DIR="$PWD"

# Site module initializers can reference unset shell variables.
if ! command -v module >/dev/null 2>&1; then
  for init in "${RGLL_MODULE_INIT:-/etc/profile.d/modules.sh}" /usr/share/lmod/lmod/init/bash /etc/profile; do
    [[ -r "$init" ]] || continue
    source "$init" || true
    set +u
    command -v module >/dev/null 2>&1 && break
  done
fi
set +u
r_module=${R_MODULE:-r/4.3.0-gcc-13.2.0-withx-rmath-standalone-python-3.11.6}
printf '[RGLL] Loading R module: %s\n' "$r_module"
if ! module load "$r_module"; then
  echo "[RGLL] Retrying the same R module without the module cache."
  module --ignore_cache load "$r_module" || {
    echo "[RGLL] R module unavailable on this node; no experiment was started." >&2
    exit 1
  }
fi
Rscript --version
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1

printf '[RGLL] MODEL=%s CASE=%s R=%s Var.a=%s repetitions=%s\n' \
  "${MODEL:-RI}" "${CASE:-all}" "${R_LIST:-10,20,50}" "${VAR_A_LIST:-0.5,2.25}" "${NLOOP:-50}"
printf '[RGLL] Replication range=%s-%s; no checkpoints; one final result per job.\n' \
  "${REP_FIRST:-1}" "${REP_LAST:-${NLOOP:-50}}"
Rscript run_compare.R
