#!/usr/bin/env bash
# repro/run_all.sh — End-to-end reproduction driver for the paper headlines.
#
# Runs (in order):
#   1. Non-linear LHCb 2025 refit             → reports/wo016c_nonlinear_refit.{md,csv}
#   2. Non-linear refit across 5 datasets     → reports/wo016d_nonlinear_xdataset.{md,csv}
#   3. Akaike-weight stacking across 5 fits   → reports/wo016a_akaike_stack.{md,csv}
#   4. Pure-geometry variant-selection check  → reports/wo016b_variant_geometry.{md,csv}
#   5. Paper figures F1, F2, F3               → paper/figures/fig_F{1,2,3}_*.{pdf,png}
#
# Idempotent. Persistent on-disk caches (data/processed/flavio_cache.json,
# graph cache inside wo010_universality) make subsequent runs near-instant.
#
# Usage:    bash repro/run_all.sh
# Requires: pip install -e ".[dev,plotting]" run once first.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

export PYTHONPATH="src${PYTHONPATH:+:$PYTHONPATH}"

step() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }

step "1/5  Non-linear LHCb 2025 refit (~30s)"
python3 scripts/wo016c_nonlinear_refit.py

step "2/5  Non-linear refit across 5 datasets (~2min)"
python3 scripts/wo016d_nonlinear_xdataset.py

step "3/5  Akaike-weight stacking"
python3 scripts/wo016a_akaike_stack.py

step "4/5  Pure-geometry variant-selection check"
python3 scripts/wo016b_variant_geometry.py

step "5/5  Paper figures (F1, F2, F3)"
python3 scripts/wo017_paper_figures.py

printf '\n\033[1;32mDone.\033[0m  See reports/ and paper/figures/ for outputs.\n'
printf 'To recompile the paper PDF: tectonic -X compile paper/main.tex\n'
