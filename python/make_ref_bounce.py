#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd
from hybrid_bounce import integrate, to_grid
ap = argparse.ArgumentParser()
ap.add_argument("--t_end", type=float, default=10.0)
ap.add_argument("--rtol", type=float, default=1e-12)
ap.add_argument("--atol", type=float, default=1e-14)
ap.add_argument("--e", type=float, default=0.9)
ap.add_argument("--h0", type=float, default=10.0)
ap.add_argument("--v0", type=float, default=0.0)
ap.add_argument("--out", default="ref_bounce.csv")
args = ap.parse_args()
T, Y = integrate(np.array([args.h0,args.v0], float), args.t_end, args.e, args.rtol, args.atol)
grid = np.linspace(0,args.t_end,1001); series = to_grid(T,Y,grid)
pd.DataFrame({"t":grid,"h":series[:,0],"v":series[:,1]}).to_csv(args.out,index=False)
print(f"[OK] wrote {args.out}")