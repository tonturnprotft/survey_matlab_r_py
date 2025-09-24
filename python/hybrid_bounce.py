#!/usr/bin/env python3
import argparse, json, time, numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from pathlib import Path

G = 9.81
def bounce_event(t, y): return y[0]
bounce_event.terminal = True; bounce_event.direction = -1

def integrate(y0, t_end, e, rtol, atol):
    t0 = 0.0; T=[t0]; Y=[y0.copy()]
    while t0 < t_end:
        sol = solve_ivp(lambda t,y: [y[1], -G], [t0, t_end], y0, rtol=rtol, atol=atol, events=bounce_event)
        T += sol.t[1:].tolist(); Y += sol.y.T[1:].tolist()
        if sol.t_events[0].size == 0: break
        t0 = float(sol.t_events[0][-1])
        y0 = np.array([0.0, -e*sol.y_events[0][-1][1]], float)
    return np.array(T), np.array(Y)

def to_grid(T, Y, grid):
    h = np.interp(grid, T, Y[:,0]); v = np.interp(grid, T, Y[:,1]); return np.column_stack([h,v])


# --- Piecewise-aligned RMSE helpers ---
def _detect_ref_events(t_ref, h_ref, v_ref):
    """Detect impact times from a dense reference series using height zero-crossings
    with descending-then-ascending velocity pattern. Returns a sorted list of times.
    """
    t_ref = np.asarray(t_ref)
    h_ref = np.asarray(h_ref)
    v_ref = np.asarray(v_ref)
    events = []
    for i in range(len(t_ref) - 1):
        h0, h1 = h_ref[i], h_ref[i+1]
        v0, v1 = v_ref[i], v_ref[i+1]
        # zero-crossing by height near ground; prefer descending->ascending
        if h0 == h1:
            continue
        if (h0 > 0 and h1 <= 0) or (h0 >= 0 and h1 < 0):
            # linear interpolate to h=0
            alpha = h0 / (h0 - h1)
            tz = t_ref[i] + alpha * (t_ref[i+1] - t_ref[i])
            # optional velocity flip check; keep even if noisy
            events.append(float(tz))
    # de-dup and sort
    if not events:
        return []
    events = sorted(set(round(t, 12) for t in events))
    return events


def _piecewise_aligned_rmse(T, Y, t_ref, h_ref, v_ref, t_end):
    """Compute NRMSE/Linf by aligning segments between reference impact times.
    Normalization: RMSE(h)/range(h_ref) and RMSE(v)/max(|v_ref|), averaged.
    """
    events = [t for t in _detect_ref_events(t_ref, h_ref, v_ref) if 0.0 < t < t_end]
    bounds = [0.0] + events + [t_end]
    sse_h = 0.0
    sse_v = 0.0
    n = 0
    linf_norm = 0.0

    # reference scales
    h_scale = float(max(1e-12, np.max(h_ref) - np.min(h_ref)))
    v_scale = float(max(1e-12, np.max(np.abs(v_ref))))

    for a, b in zip(bounds[:-1], bounds[1:]):
        if b <= a:
            continue
        eps = min(1e-6, 1e-6 * t_end)
        aa = a + eps
        bb = b - eps
        if bb <= aa:
            aa, bb = a, b
        grid = np.linspace(aa, bb, 201)
        cand = to_grid(T, Y, grid)
        ref = np.column_stack([
            np.interp(grid, t_ref, h_ref),
            np.interp(grid, t_ref, v_ref),
        ])
        dh = cand[:,0] - ref[:,0]
        dv = cand[:,1] - ref[:,1]
        sse_h += float(np.sum(dh * dh))
        sse_v += float(np.sum(dv * dv))
        n += dh.size
        linf_norm = max(linf_norm, float(max(np.max(np.abs(dh))/h_scale, np.max(np.abs(dv))/v_scale)))

    if n == 0:
        return None, None
    rmse_h = (sse_h / n) ** 0.5
    rmse_v = (sse_v / n) ** 0.5
    nrmse = 0.5 * (rmse_h / h_scale + rmse_v / v_scale)
    return nrmse, linf_norm

ap = argparse.ArgumentParser()
ap.add_argument("--t_end", type=float, default=10.0)
ap.add_argument("--rtol", type=float, default=1e-6)
ap.add_argument("--atol", type=float, default=1e-8)
ap.add_argument("--e", type=float, default=0.9)
ap.add_argument("--h0", type=float, default=10.0)
ap.add_argument("--v0", type=float, default=0.0)
ap.add_argument("--out", required=True)
ap.add_argument("--ref", default=None)  # CSV: t,h,v
args = ap.parse_args()

grid = np.linspace(0, args.t_end, 1001)
t0 = time.perf_counter()
T, Y = integrate(np.array([args.h0, args.v0], float), args.t_end, args.e, args.rtol, args.atol)
wall = time.perf_counter() - t0
series = to_grid(T, Y, grid)

rmse = linf = None
if args.ref and Path(args.ref).exists():
    ref_df = pd.read_csv(args.ref)
    t_ref = ref_df["t"].to_numpy()
    h_ref = ref_df["h"].to_numpy()
    v_ref = ref_df["v"].to_numpy()
    rmse, linf = _piecewise_aligned_rmse(T, Y, t_ref, h_ref, v_ref, args.t_end)

row = {
    "platform":"python","family":"HYBRID","benchmark_id":"hybrid_bounce",
    "tool":"scipy.integrate.solve_ivp","solver":"events","status":"OK",
    "wall_time_s":wall,"rmse":rmse,"linf":linf,"n_steps":Y.shape[0],
    "params_json":json.dumps({"t_end":args.t_end,"rtol":args.rtol,"atol":args.atol,"e":args.e,"h0":args.h0,"v0":args.v0}),
}
import csv
header = not Path(args.out).exists()
with open(args.out,"a",newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(row.keys()))
    if header: w.writeheader()
    w.writerow(row)
print(f"status=OK steps={Y.shape[0]} wall={wall:.3f}s rmse={rmse} linf={linf}")