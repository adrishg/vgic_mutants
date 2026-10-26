#!/usr/bin/env bash
# Submit one of the four Nav1.5 RMSD convergence wrappers.
# Usage: sbatch job_nav15_rmsd.sbatch <wrapper-name>
#SBATCH --job-name=nav15-rmsd
#SBATCH --partition=high
#SBATCH --account=genome-center-grp
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: sbatch $0 <wrapper-name>" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER="$SCRIPT_DIR/$1"

if [[ ! -x "$WRAPPER" ]]; then
  echo "Executable wrapper not found: $WRAPPER" >&2
  exit 1
fi

mkdir -p "$SCRIPT_DIR/logs"
export NAV15_RMSD_WORKERS="${NAV15_RMSD_WORKERS:-4}"
exec "$WRAPPER"
