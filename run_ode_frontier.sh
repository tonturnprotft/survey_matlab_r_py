#!/usr/bin/env bash
set -euo pipefail
ROOT="${ROOT:-/Users/tree/Simu}"
MATLAB_PATH="${MATLAB_PATH:-$ROOT/matlab}"
cd "$ROOT"

TEND="${TEND:-1e5}"
RTOLS=(${RTOLS:-1e-3 1e-6 1e-9})
REPS_MORE=${REPS_MORE:-5}
REF="${REF:-$ROOT/ref_robertson_t1e5.csv}"

for R in "${RTOLS[@]}"; do
  A=$(python - <<PY
R=${R}; print(R*1e-2)
PY
)
  for rep in $(seq 1 "$REPS_MORE"); do
    python run_with_time.py --out results/results_ode_robertson_python.csv -- \
      python python/ode_robertson.py --t_end "$TEND" --rtol "$R" --atol "$A" --method BDF \
        --out results/results_ode_robertson_python.csv --ref "$REF"

    python run_with_time.py --out results/results_ode_robertson_R.csv -- \
      Rscript R/ode_robertson.R "$TEND" "$R" "$A" lsoda results/results_ode_robertson_R.csv "$REF"

    python run_with_time.py --out results/results_ode_robertson_matlab.csv -- \
      matlab -batch "addpath('$MATLAB_PATH'); ode_robertson($TEND,$R,$A,'ode15s','/Users/tree/Simu/results_ode_robertson_matlab.csv','$REF')"
  done
done
echo "[OK] ODE frontier extra reps complete."