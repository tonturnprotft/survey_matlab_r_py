#!/usr/bin/env bash
set -euo pipefail
ROOT="${ROOT:-/Users/tree/Simu}"
MATLAB_PATH="${MATLAB_PATH:-$ROOT/matlab}"
cd "$ROOT"

RTOL="${RTOL:-1e-6}"
ATOL="${ATOL:-1e-8}"
T_LIST=(${T_LIST:-1e3 1e4 1e5})
REPS_MORE=${REPS_MORE:-5}
REF="${REF:-$ROOT/ref_robertson_t1e5.csv}"

for T in "${T_LIST[@]}"; do
  for rep in $(seq 1 "$REPS_MORE"); do
    python run_with_time.py --out results/results_ode_robertson_python.csv -- \
      python python/ode_robertson.py --t_end "$T" --rtol "$RTOL" --atol "$ATOL" --method BDF \
        --out results/results_ode_robertson_python.csv --ref "$REF"

    python run_with_time.py --out results/results_ode_robertson_R.csv -- \
      Rscript R/ode_robertson.R "$T" "$RTOL" "$ATOL" lsoda results/results_ode_robertson_R.csv "$REF"

    python run_with_time.py --out results/results_ode_robertson_matlab.csv -- \
      matlab -batch "addpath('$MATLAB_PATH'); ode_robertson($T,$RTOL,$ATOL,'ode15s','/Users/tree/Simu/results_ode_robertson_matlab.csv','$REF')"
  done
done
echo "[OK] ODE scalability extra reps complete."