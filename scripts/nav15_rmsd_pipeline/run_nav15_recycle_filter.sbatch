#!/usr/bin/env bash
#SBATCH --job-name=nav15-filter
#SBATCH --account=genome-center-grp
#SBATCH --partition=high
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/nav15_filter_%A_%a.out
#SBATCH --error=logs/nav15_filter_%A_%a.err

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: sbatch $0 <condition> <protocol_dir> <output_name>" >&2
  exit 2
fi

CONDITION="$1"
PROTOCOL_DIR="$2"
OUTPUT_NAME="$3"
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
PROJECT_ROOT="${NAV15_PROJECT_ROOT:-/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5}"
FILTER_ROOT="${AF2_FILTER_ROOT:-/quobyte/yarovoygrp/ahgz/vgic_mutants/scripts/af2_recycle_rmsd_convergence_final}"
INPUT_DIR="$PROJECT_ROOT/$CONDITION/$PROTOCOL_DIR"
OUTPUT_DIR="$PROJECT_ROOT/Nav15_rmsd_convergence_${OUTPUT_NAME}"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory does not exist: $INPUT_DIR" >&2
  exit 1
fi
if [[ -e "$OUTPUT_DIR" ]] && find "$OUTPUT_DIR" -mindepth 1 -print -quit | grep -q .; then
  echo "Output already contains files: $OUTPUT_DIR" >&2
  exit 1
fi

eval "$(micromamba shell hook --shell bash)"
micromamba activate /quobyte/yarovoygrp/ahgz/mamba_envs/envs/bioadri
mkdir -p "$OUTPUT_DIR" "$SCRIPT_DIR/logs"

python "$FILTER_ROOT/af2_recycle_rmsd_convergence.py" \
  --input-dir "$INPUT_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --config "$SCRIPT_DIR/config_nav15_filter_3A.json" \
  --condition "$CONDITION" \
  --protocol masked \
  --workers "${NAV15_RMSD_WORKERS:-4}"
