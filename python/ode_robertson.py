# ode_robertson.py — Solve stiff Robertson ODE and emit one CSV row (results_schema-compatible).
# Usage:
#   python3 ode_robertson.py --t_end 1e5 --rtol 1e-6 --atol 1e-8 --method BDF \
#       --out results_ode_robertson_python.csv --ref ref_robertson_t1e5.csv
import argparse, time, csv, json, platform, os
import numpy as np
from math import log10
from pathlib import Path

from scipy import __version__ as scipy_version
from scipy.integrate import solve_ivp

def robertson(t, y):
    y1, y2, y3 = y
    dy1 = -0.04*y1 + 1.0e4*y2*y3
    dy2 =  0.04*y1 - 1.0e4*y2*y3 - 3.0e7*(y2**2)
    dy3 =  3.0e7*(y2**2)
    return [dy1, dy2, dy3]

def sample_times(t_end, n=500):
    # 0 plus logspace from 1e-6 to t_end with n-1 points
    if t_end <= 1e-6:
        grid = np.linspace(0.0, t_end, n)
    else:
        ls = np.logspace(-6, np.log10(t_end), num=n-1)
        grid = np.concatenate(([0.0], ls))
    return grid

def load_reference_csv(path):
    p = Path(path)
    if not p.exists():
        return None
    arr = np.genfromtxt(p, delimiter=",", names=True)
    # Expect columns: t,y1,y2,y3
    return np.vstack([arr['t'], arr['y1'], arr['y2'], arr['y3']]).T

def compute_errors(sol_times, sol_vals, ref):
    # ref: Nx4 [t,y1,y2,y3]; interpolate ref onto sol_times then compute RMSE/Linf across all states
    if ref is None:
        return "", ""
    t_ref = ref[:,0]
    ref_interp = np.vstack([
        np.interp(sol_times, t_ref, ref[:,1]),
        np.interp(sol_times, t_ref, ref[:,2]),
        np.interp(sol_times, t_ref, ref[:,3]),
    ]).T
    diff = sol_vals - ref_interp
    rmse = float(np.sqrt(np.mean(diff**2)))
    linf = float(np.max(np.abs(diff)))
    return rmse, linf

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t_end", type=float, default=1e5)
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-8)
    ap.add_argument("--method", type=str, default="BDF", choices=["BDF","LSODA","Radau","DOP853"])
    ap.add_argument("--out", type=str, default="results_ode_robertson_python.csv")
    ap.add_argument("--ref", type=str, default="")
    args = ap.parse_args()

    y0 = [1.0, 0.0, 0.0]
    t_eval = sample_times(args.t_end, n=500)

    t0_wall = time.perf_counter()
    t0_cpu = time.process_time()
    sol = solve_ivp(robertson, (0.0, args.t_end), y0, method=args.method, rtol=args.rtol, atol=args.atol, t_eval=t_eval, dense_output=False)
    wall = time.perf_counter() - t0_wall
    cpu  = time.process_time() - t0_cpu

    status = "OK" if sol.success else "FAIL"
    err_msg = "" if sol.success else str(sol.message)
    n_steps = int(len(sol.t) - 1)

    rmse, linf = "", ""
    ref = load_reference_csv(args.ref) if args.ref else None
    if ref is not None and status == "OK":
        vals = sol.y.T  # shape (N,3)
        rmse, linf = compute_errors(sol.t, vals, ref)

    params = {"t_end": args.t_end, "rtol": args.rtol, "atol": args.atol, "method": args.method}

    row = {
        "family": "ODE",
        "benchmark_id": "ode_robertson",
        "platform": "python",
        "tool": "scipy.integrate.solve_ivp",
        "version": scipy_version,
        "solver": args.method,
        "step_mode": "adaptive",
        "abs_tol": args.atol,
        "rel_tol": args.rtol,
        "dt_init": "",
        "dt_max": "",
        "seed": "",
        "params_json": json.dumps(params),
        "scale_param": "t_end",
        "scale_value": args.t_end,
        "repeat_idx": 0,
        "run_id": f"py_{int(time.time())}",
        "start_ts": "",
        "end_ts": "",
        "wall_time_s": wall,
        "cpu_time_s": cpu,
        "peak_mem_mb": "",
        "n_steps": n_steps,
        "n_events": "",
        "rmse": rmse,
        "linf": linf,
        "event_time_err_mean": "",
        "event_time_err_max": "",
        "throughput": "",
        "L_timeavg": "",
        "status": status,
        "error_message": err_msg,
        "host": platform.platform(),
    }

    field_order = ["family","benchmark_id","platform","tool","version","solver","step_mode",
                   "abs_tol","rel_tol","dt_init","dt_max","seed","params_json",
                   "scale_param","scale_value","repeat_idx","run_id","start_ts","end_ts",
                   "wall_time_s","cpu_time_s","peak_mem_mb","n_steps","n_events","rmse","linf",
                   "event_time_err_mean","event_time_err_max","throughput","L_timeavg","status","error_message","host"]
    exists = os.path.exists(args.out)
    with open(args.out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=field_order)
        if not exists:
            w.writeheader()
        w.writerow(row)
    print(f"status={status} n_steps={n_steps} wall={wall:.3f}s rmse={rmse} linf={linf}")

if __name__ == "__main__":
    main()