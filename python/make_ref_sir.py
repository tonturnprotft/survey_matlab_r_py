#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd
from stoch_sir import sir_ssa_once, series_to_grid

ap = argparse.ArgumentParser()
ap.add_argument("--beta", type=float, default=0.3/1000.0)
ap.add_argument("--gamma", type=float, default=0.1)
ap.add_argument("--N", type=int, default=1000)
ap.add_argument("--I0", type=int, default=10)
ap.add_argument("--t_end", type=float, default=160.0)
ap.add_argument("--n_rep", type=int, default=5000)
ap.add_argument("--seed", type=int, default=99)
ap.add_argument("--out", default="ref_sir_mean.csv")
args = ap.parse_args()

S0 = args.N - args.I0
grid = np.linspace(0, args.t_end, int(args.t_end)+1)
rng = np.random.default_rng(args.seed)
acc = np.zeros((grid.size,3), float)
for _ in range(args.n_rep):
    T, X = sir_ssa_once(args.beta, args.gamma, S0, args.I0, 0, args.t_end, rng)
    acc += series_to_grid(T, X, grid)
mean = acc/args.n_rep
pd.DataFrame({"t":grid,"S_mean":mean[:,0],"I_mean":mean[:,1],"R_mean":mean[:,2]}).to_csv(args.out, index=False)
print(f"[OK] wrote {args.out} (n_rep={args.n_rep})")