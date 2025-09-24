# ode_robertson_make_ref.py — Build a high-accuracy reference trajectory and save to CSV.
# Usage:
#   python3 ode_robertson_make_ref.py --t_end 1e5 --rtol 1e-12 --atol 1e-14 --method BDF --out ref_robertson_t1e5.csv
import argparse, numpy as np, csv
from scipy.integrate import solve_ivp

def robertson(t, y):
    y1, y2, y3 = y
    return [-0.04*y1 + 1.0e4*y2*y3, 0.04*y1 - 1.0e4*y2*y3 - 3.0e7*(y2**2), 3.0e7*(y2**2)]

def sample_times(t_end, n=2000):
    if t_end <= 1e-6:
        return np.linspace(0.0, t_end, n)
    else:
        return np.concatenate(([0.0], np.logspace(-6, np.log10(t_end), num=n-1)))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t_end", type=float, default=1e5)
    ap.add_argument("--rtol", type=float, default=1e-12)
    ap.add_argument("--atol", type=float, default=1e-14)
    ap.add_argument("--method", type=str, default="BDF")
    ap.add_argument("--out", type=str, default="ref_robertson_t1e5.csv")
    args = ap.parse_args()

    y0 = [1.0, 0.0, 0.0]
    t_eval = sample_times(args.t_end, n=2000)
    sol = solve_ivp(robertson, (0.0, args.t_end), y0, method=args.method, rtol=args.rtol, atol=args.atol, t_eval=t_eval)
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t","y1","y2","y3"])
        for ti, yi in zip(sol.t, sol.y.T):
            w.writerow([ti, yi[0], yi[1], yi[2]])
    print(f"Wrote reference to {args.out} with {len(sol.t)} samples")

if __name__ == "__main__":
    main()