# add 5 more repeats per n_rep level; vary seeds so trajectories differ
REPS_MORE=${REPS_MORE:-5}
for R in 1 3 10 30 100; do
  for rep in $(seq 1 "$REPS_MORE"); do
    SEED=$((42 + rep))

    python run_with_time.py --out results/results_ssa_sir_python.csv -- \
      python3 python/stoch_sir.py --n_rep "$R" --seed "$SEED" \
        --out results/results_ssa_sir_python.csv --ref ref_sir_mean.csv

    python run_with_time.py --out results/results_ssa_sir_R.csv -- \
      Rscript R/stoch_sir.R 0.0003 0.1 1000 10 160 "$R" "$SEED" results/results_ssa_sir_R.csv ref_sir_mean.csv

    python run_with_time.py --out results_ssa_sir_matlab.csv -- \
      matlab -batch "addpath('/Users/tree/Simu/matlab');stoch_sir(0.0003,0.1,1000,10,160,$R,$SEED,'/User/tree/Simu/results_ssa_sir_matlab.csv','/Users/tree/Simu/ref_sir_mean.csv')"
  done
done