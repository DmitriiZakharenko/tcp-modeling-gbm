#!/usr/bin/env bash
# Execute notebooks top-to-bottom; write log to reports/notebook_run_log.txt
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG="$ROOT/reports/notebook_run_log.txt"
TIMEOUT=900
PYTHON="${PYTHON:-python3}"
JUPYTER="$PYTHON -m jupyter nbconvert"

NOTEBOOKS=(
  "$ROOT/notebooks/01_cohort_overview.ipynb"
  "$ROOT/notebooks/03_tcp_models.ipynb"
  "$ROOT/notebooks/04_parameter_estimation.ipynb"
  "$ROOT/notebooks/05_survival_analysis.ipynb"
  "$ROOT/notebooks/06_rano_multivariable_40gy.ipynb"
)

# Include 02 when local RTDOSE/GTV NIfTI is present (~52 GB), or when forced.
RAW="$ROOT/data/raw"
if [[ "${RUN_NOTEBOOK_02:-0}" == "1" ]]; then
  NOTEBOOKS=("$ROOT/notebooks/02_feature_extraction.ipynb" "${NOTEBOOKS[@]}")
elif [[ -d "$RAW" ]] && compgen -G "$RAW/**/RTDOSE*.nii.gz" >/dev/null 2>&1; then
  echo "data/raw RTDOSE found — including notebook 02_feature_extraction.ipynb" | tee -a "$LOG"
  NOTEBOOKS=("$ROOT/notebooks/02_feature_extraction.ipynb" "${NOTEBOOKS[@]}")
fi

: > "$LOG"
echo "Notebook run started $(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$LOG"
FAIL=0

for nb in "${NOTEBOOKS[@]}"; do
  name="$(basename "$nb")"
  echo "=== Executing $name ===" | tee -a "$LOG"
  if $JUPYTER --execute --to notebook --inplace "$nb" \
      --ExecutePreprocessor.timeout="$TIMEOUT" >>"$LOG" 2>&1; then
    echo "OK  $name" | tee -a "$LOG"
  else
    echo "FAIL $name (see log)" | tee -a "$LOG"
    FAIL=1
  fi
done

echo "Finished $(date -u +%Y-%m-%dT%H:%M:%SZ) exit=$FAIL" | tee -a "$LOG"
exit "$FAIL"
