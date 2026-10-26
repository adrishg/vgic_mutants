#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 <condition> <protocol_dir> <output_name>" >&2
  exit 2
fi

CONDITION="$1"
PROTOCOL_DIR="$2"
OUTPUT_NAME="$3"

PROJECT_ROOT="${NAV15_PROJECT_ROOT:-/quobyte/yarovoygrp/ahgz/vgic_mutants/Nav1.5}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="$PROJECT_ROOT/$CONDITION/$PROTOCOL_DIR"
OUTPUT_DIR="$PROJECT_ROOT/Nav15_rmsd_convergence_${OUTPUT_NAME}"
CONFIG_PATH="${NAV15_RMSD_CONFIG:-$SCRIPT_DIR/config_nav15.json}"
WORKERS="${NAV15_RMSD_WORKERS:-4}"
PROGRAM="${AF2_RECYCLE_RMSD_PROGRAM:-af2_recycle_rmsd_convergence.py}"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Input directory does not exist: $INPUT_DIR" >&2
  exit 1
fi
if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file does not exist: $CONFIG_PATH" >&2
  exit 1
fi
if ! command -v "$PROGRAM" >/dev/null 2>&1; then
  echo "Cannot find $PROGRAM on PATH; activate the established analysis environment first." >&2
  exit 1
fi
if [[ -e "$OUTPUT_DIR" && "${ALLOW_EXISTING_OUTPUT:-0}" != 1 ]]; then
  echo "Output already exists: $OUTPUT_DIR" >&2
  echo "Set ALLOW_EXISTING_OUTPUT=1 only after checking the previous run." >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
echo "Running Nav1.5 recycle-RMSD convergence"
echo "  condition: $CONDITION"
echo "  protocol:  $PROTOCOL_DIR"
echo "  input:     $INPUT_DIR"
echo "  output:    $OUTPUT_DIR"
echo "  workers:   $WORKERS"

exec "$PROGRAM" \
  --input-dir "$INPUT_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --config "$CONFIG_PATH" \
  --workers "$WORKERS"
