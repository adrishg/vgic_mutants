#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JOB_SCRIPT="$SCRIPT_DIR/run_nav15_recycle_filter.sbatch"
mkdir -p "$SCRIPT_DIR/logs"

previous_job=""
for spec in \
  "wt masked_v1new wt_masked_v1new" \
  "wt masked_v2new wt_masked_v2new" \
  "qqq masked_v1new qqq_masked_v1new" \
  "qqq masked_v2new qqq_masked_v2new"; do
  read -r condition protocol output_name <<< "$spec"
  if [[ -n "$previous_job" ]]; then
    job_id=$(sbatch --parsable --dependency="afterok:${previous_job}" "$JOB_SCRIPT" "$condition" "$protocol" "$output_name")
  else
    job_id=$(sbatch --parsable "$JOB_SCRIPT" "$condition" "$protocol" "$output_name")
  fi
  echo "${output_name}: ${job_id}"
  previous_job="$job_id"
done

echo "Submitted four recycle-filtering jobs as an afterok chain."
