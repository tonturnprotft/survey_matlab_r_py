#!/usr/bin/env python3
import argparse, json, math, time
import numpy as np, pandas as pd
from pathlib import Path

def sir_ssa_once(beta, gamma, S0, I0, R0, t_end, rng):
    S, I, R = int(S0), int(I0), int(R0)
    t = 0.0
    times = [0.0]; states = [(S,I,R)]
    while t < t_end and I > 0:
        a1 = beta*S*I
        a2 = gamma*I
        a0 = a1 + a2
        if a0 <= 0: break
        u1 = rng.random()
        tau = -math.log(u1)/a0
        t = t + tau
        u2 = rng.random()*a0
        if u2 < a1 and S>0:
            S -= 1; I += 1
        else:
            I -= 1; R += 1
        times.append(min(t, t_end))
        states.append((S,I,R))
    if times[-1] < t_end:
        times.append(t_end); states.append((S,I,R))
    return np.array(times), np.array(states, float)

def series_to_grid(times, states, grid):
    S = np.interp(grid, times, states[:,0], left=states[0,0], right=states[-1,0])
    I = np.interp(grid, times, states[:,1], left=states[0,1], right=states[-1,1])
    R = np.interp(grid, times, states[:,2], left=states[0,2], right=states[-1,2])
    return np.column_stack([S,I,R])

def rmse_vs_ref(mean_series, ref_df):
    ref = ref_df[["S_mean","I_mean","R_mean"]].to_numpy()
    diff = mean_series - ref
    rmse = float(np.sqrt((diff**2).mean()))
    linf = float(np.max(np.abs(diff)))
    return rmse, linf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float, default=0.3/1000.0)
    ap.add_argument("--gamma", type=float, default=0.1)
    ap.add_argument("--N", type=int, default=1000)
    ap.add_argument("--I0", type=int, default=10)
    ap.add_argument("--t_end", type=float, default=160.0)
    ap.add_argument("--n_rep", type=int, default=10)
    ap.add_argument("--seed", type=int, default=123)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref", default=None)  # CSV: t,S_mean,I_mean,R_mean
    args = ap.parse_args()

    S0 = args.N - args.I0
    grid = np.linspace(0, args.t_end, int(args.t_end)+1)
    rng = np.random.default_rng(args.seed)

    t0 = time.perf_counter()
    acc = np.zeros((grid.size, 3), float)
    steps_total = 0
    for _ in range(args.n_rep):
        times, states = sir_ssa_once(args.beta, args.gamma, S0, args.I0, 0, args.t_end, rng)
        steps_total += len(times)-1
        acc += series_to_grid(times, states, grid)
    mean_series = acc / args.n_rep
    wall = time.perf_counter() - t0

    rmse = linf = None
    if args.ref and Path(args.ref).exists():
        ref_df = pd.read_csv(args.ref)
        rmse, linf = rmse_vs_ref(mean_series, ref_df)

    row = {
        "platform":"python","family":"STOCH","benchmark_id":"ssa_sir",
        "tool":"event_loop","solver":"gillespie","status":"OK",
        "wall_time_s":wall,"rmse":rmse,"linf":linf,"n_steps":steps_total,
        "params_json":json.dumps({"beta":args.beta,"gamma":args.gamma,"N":args.N,
                                  "I0":args.I0,"t_end":args.t_end,"n_rep":args.n_rep,"seed":args.seed}),
    }
    p = Path(args.out); header = not p.exists()
    import csv
    with p.open("a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if header: w.writeheader()
        w.writerow(row)
    print(f"status=OK n_steps={steps_total} wall={wall:.3f}s rmse={rmse} linf={linf}")

if __name__ == "__main__":
    main()