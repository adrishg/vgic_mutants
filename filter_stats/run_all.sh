#!/usr/bin/env bash
# Run data preparation (step 1), the topology screen (when the config has [screen]), the regions (step 2)
# and the statistics for one dataset config.
# Usage, from this folder:  PYTHON=~/anaconda3/bin/python3 ./run_all.sh datasets/<dataset>/<dataset>.toml
set -euo pipefail
CONFIG_ARG=${1:?usage: ./run_all.sh datasets/<dataset>/<dataset>.toml}
CONFIG="$(cd "$(dirname "$CONFIG_ARG")" && pwd)/$(basename "$CONFIG_ARG")"
PY=${PYTHON:-python3}
cd "$(dirname "$0")"
"$PY" -B pipeline/build_dv.py --config "$CONFIG"
if grep -q '^\[screen\]' "$CONFIG"; then
  "$PY" -B pipeline/screen.py --config "$CONFIG"
fi
"$PY" -B pipeline/region_dv.py --config "$CONFIG"
"$PY" -B pipeline/stats.py --config "$CONFIG"
