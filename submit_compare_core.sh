#!/bin/bash
# Queue the 64 independent main-study jobs; only final results are written.
set -euo pipefail
cd "$(dirname "$0")"
export PROJECT_DIR="$PWD"
[[ "${1:-}" != base ]] || shift
export RESULT_DIR=${RESULT_DIR:-result/results_main}
[[ "$RESULT_DIR" == /* ]] || export RESULT_DIR="$PROJECT_DIR/$RESULT_DIR"
[[ $# -le 2 ]] || { echo "Usage: bash submit_compare_core.sh [FIRST [LAST]]" >&2; exit 1; }
export NLOOP=${NLOOP:-50}
[[ "$NLOOP" =~ ^[1-9][0-9]*$ ]] || { echo "NLOOP must be positive." >&2; exit 1; }
total=64
first=${1:-1}; last=${2:-${1:-$total}}
[[ "$first" =~ ^[1-9][0-9]*$ && "$last" =~ ^[1-9][0-9]*$ ]] || { echo "FIRST/LAST must be positive integers." >&2; exit 1; }
(( first <= last && last <= total )) || { echo "Task range must be within 1-$total." >&2; exit 1; }

unset RGLL_RCPP_CACHE REPLICATION_IDS CASE_LIST R VAR_A METHODS RUN_NO_QUIT MODEL CASE
unset RUN_BUDGET_SECONDS JOB_DEADLINE_EPOCH SLURM_ARRAY_TASK_ID TASK_OFFSET
unset BASH_ENV ENV BASE_TASK_LAYOUT SUBMIT_MODE
export LABEL_POLICY=unknown MIS_TYPE=none RHO=0 LAMBDA=0
export N=2500 P=50 VAR_B=0.1 TAU=255
export SA_MAX=25 EM_MAX=300 K_GRID_POINTS=8 KMEANS_NUM_INIT=1 KMEANS_MAX_ITERS=25
export CPF_LAMBDA_POINTS=6 CPF_MAX_ITER=300
export BLM_NSTART=5 CIRG_CRITERION=IMSPE
mkdir -p "$RESULT_DIR/logs" "$RESULT_DIR/tasks"
accepted=(); task=0
submit_task() {
  task=$((task + 1))
  (( task >= first && task <= last )) || return 0
  export RUN_MODE=base MODEL=$1 CASE=$2 VAR_A_LIST=$3 R_LIST=10,20,50
  export REP_FIRST=1 REP_LAST=$NLOOP N_CORES=8
  resources=("--partition=${CPU_PARTITION:-cpu}" --cpus-per-task=8 --mem=100G)
  [[ -z "${CPU_EXCLUDE:-}" ]] || resources+=("--exclude=$CPU_EXCLUDE")
  var_tag=${VAR_A_LIST//./p}; var_tag=${var_tag//,/-}
  printf -v name 'RG_%s_C%02d_A%s_BASE' "$MODEL" "$CASE" "$var_tag"
  command=(sbatch --parsable --export=ALL --chdir="$PROJECT_DIR" "${resources[@]}"
    --time=48:00:00 --job-name="$name" --mail-type=END,FAIL,TIME_LIMIT
    "--mail-user=${MAIL_USER:-k21181837@kcl.ac.uk}"
    --output="$RESULT_DIR/logs/${name}_%j.txt" --error="$RESULT_DIR/logs/${name}_%j.txt" compare_job.sh)
  if [[ "${DRY_RUN:-0}" == 1 ]]; then
    printf '[task %s/%s] ' "$task" "$total"; printf '%q ' "${command[@]}"; printf '\n'; return 0
  fi
  if reply=$("${command[@]}"); then job=${reply%%;*}; else
    echo "Submission stopped at task $task. Confirmed job IDs: ${accepted[*]:-none}" >&2
    echo "A timed-out sbatch may still have accepted this task. Check squeue/sacct before retrying." >&2
    exit 1
  fi
  [[ "$job" =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch reply: $reply" >&2; exit 1; }
  accepted+=("$job")
  printf 'Task %s/%s: %s submitted as %s\n' "$task" "$total" "$name" "$job"
  printf '%s\t%s\t%s\t%s\n' "$task" "$job" "$name" "$RUN_MODE" >> "$RESULT_DIR/jobs.tsv"
}
echo "Submitting $first-$last of $total independent experiment jobs. No checkpoints or test jobs."
for model in RI RS; do
  for case_id in {1..16}; do
    for var_a in 0.5 2.25; do submit_task "$model" "$case_id" "$var_a"; done
  done
done
echo "Submission complete. Slurm schedules the jobs independently; no concurrency cap or dependencies were added."
