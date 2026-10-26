#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="$SCRIPT_DIR/job_nav15_rmsd.sbatch"
mkdir -p "$SCRIPT_DIR/logs"

previous_job=""
for wrapper in \
  run_wt_masked_v1new.sh \
  run_wt_masked_v2new.sh \
  run_qqq_masked_v1new.sh \
  run_qqq_masked_v2new.sh; do
  if [[ -n "$previous_job" ]]; then
    job_id=$(sbatch --parsable --dependency="afterok:${previous_job}" "$JOB_SCRIPT" "$wrapper")
  else
    job_id=$(sbatch --parsable "$JOB_SCRIPT" "$wrapper")
  fi
  echo "$wrapper: $job_id"
  previous_job="$job_id"
done

echo "Jobs were chained with afterok; only one RMSD run will execute at a time."
