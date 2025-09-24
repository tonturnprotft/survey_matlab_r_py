# add 5 more reps per rtol level
REPS_MORE=${REPS_MORE:-5}
for r in 1e-3 1e-6 1e-9; do
  a=$(python - <<PY
r=$r; print(r*0.01)
PY
)
  for rep in $(seq 1 "$REPS_MORE"); do
    python run_with_time.py --out results_hybrid_bounce_python.csv -- \
      python python/hybrid_bounce.py --t_end 10 --rtol "$r" --atol "$a" \
        --out results_hybrid_bounce_python.csv --ref ref_bounce.csv

    python run_with_time.py --out results_hybrid_bounce_R.csv -- \
      Rscript R/hybrid_bounce.R 10 "$r" "$a" 0.9 10 0 results_hybrid_bounce_R.csv ref_bounce.csv

    python run_with_time.py --out results_hybrid_bounce_matlab.csv -- \
      matlab -batch "addpath('/Users/tree/Simu/matlab');hybrid_bounce(10,$r,$a,0.9,10,0,'results_hybrid_bounce_matlab.csv','/Users/tree/Simu/ref_bounce.csv')"
  done
done