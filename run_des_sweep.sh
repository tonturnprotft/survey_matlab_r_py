#!/usr/bin/env bash
set -euo pipefail
ROOT="${ROOT:-/Users/tree/Simu}"
MATLAB_PATH="${MATLAB_PATH:-$ROOT/matlab}"
cd "$ROOT"

LAM="${LAM:-0.9}"
MU="${MU:-1.0}"
T_LIST=(${T_LIST:-1e4 3e4 1e5 3e5 1e6})

# new seeds 6..10
SEEDS_EXTRA=(${SEEDS_EXTRA:-6 7 8 9 10})

for T in "${T_LIST[@]}"; do
  for S in "${SEEDS_EXTRA[@]}"; do
    python run_with_time.py --out results/results_des_mm1_python_v2.csv -- \
      python python/des_mm1.py --lam "$LAM" --mu "$MU" --horizon "$T" --seed "$S" --out results/results_des_mm1_python_v2.csv
  done
done

for T in "${T_LIST[@]}"; do
  for S in "${SEEDS_EXTRA[@]}"; do
    python run_with_time.py --out results/results_des_mm1_R_v2.csv -- \
      Rscript R/des_mm1.R "$LAM" "$MU" "$T" "$S" results/results_des_mm1_R_v2.csv
  done
done

for T in "${T_LIST[@]}"; do
  for S in "${SEEDS_EXTRA[@]}"; do
    python run_with_time.py --out results/results_des_mm1_matlab_v2.csv -- \
      matlab -batch "addpath('$MATLAB_PATH'); des_mm1($LAM,$MU,$T,$S,'/Users/tree/Simu/results/results_des_mm1_matlab_v2.csv')"
  done
done
echo "[OK] DES extra seeds complete."