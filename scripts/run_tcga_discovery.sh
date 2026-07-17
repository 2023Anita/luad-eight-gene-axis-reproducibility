#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"

cd "${ROOT}"

if [[ "${1:-}" != "--skip-download" ]]; then
  "${PYTHON_BIN}" scripts/download_tcga_xena_only.py
fi

"${PYTHON_BIN}" scripts/explore_luad_tme_clinical.py
"${PYTHON_BIN}" scripts/km_logrank_screen.py
"${PYTHON_BIN}" scripts/axis_score_screen.py
"${PYTHON_BIN}" scripts/prepare_formal_modeling_inputs.py
"${RSCRIPT_BIN}" scripts/run_survival_analysis.R
"${PYTHON_BIN}" scripts/build_tcga_candidate_decisions.py
"${PYTHON_BIN}" scripts/verify_tcga_reproduction.py
