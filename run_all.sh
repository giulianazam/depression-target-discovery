#!/usr/bin/env bash
# Rerun the whole analysis. Needs the GSE144136 files in data/raw (see README).
set -euo pipefail
cd "$(dirname "$0")"
python3 src/01_load_and_qc.py
python3 src/02_pseudobulk_de.py
python3 src/03_rank_targets.py
python3 src/04_figures.py
python3 src/05_weight_sensitivity.py
python3 src/06_gpr37l1_robustness.py data/raw results/tables
