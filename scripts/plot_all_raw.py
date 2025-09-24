# plot_all_raw.py — Read results_*.csv and plot raw-point graphs (no medians/CI).
# Outputs:
#   figs/des_frontier_raw.png         — DES: % error of E[N] vs time
#   figs/ode_frontier_raw.png         — ODE: RMSE vs time
#   figs/ode_scalability_raw.png      — ODE: wall time vs t_end
#   figs/memory_vs_time_raw.png       — Peak memory vs time (all families)
#
# Usage:
#   python plot_all_raw.py --outdir figs [optional list of CSVs]
#
import argparse, json, glob, os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator

def load_results(csvs):
    if not csvs:
        csvs = sorted(glob.glob("results_*.csv"))
    frames = []
    for p in csvs:
        try:
            df = pd.read_csv(p, engine="python")
            df["__source__"] = os.path.basename(p)
            frames.append(df)
        except Exception as e:
            print(f"[WARN] Failed to read {p}: {e}")
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    # Parse params_json to populate common fields
    if "params_json" in out.columns:
        def parse_params(s):
            try: return json.loads(s)
            except: return {}
        P = out["params_json"].apply(parse_params)
        for k in ["lambda","mu","horizon","t_end","rtol","atol","method","solver"]:
            if k not in out.columns:
                out[k] = P.apply(lambda d: d.get(k))
    return out

# --- aggregation helpers: median + 95% bootstrap CI ---------------------------------
import numpy as _np
import pandas as _pd

def _group_bootstrap_ci(series, n_boot=1000, q=(2.5, 97.5), agg="median", seed=0):
    x = _pd.to_numeric(series, errors="coerce").dropna().to_numpy()
    if x.size == 0:
        return _np.nan, _np.nan, _np.nan, 0
    rng = _np.random.default_rng(seed)
    stats = []
    for _ in range(n_boot):
        samp = rng.choice(x, size=x.size, replace=True)
        stats.append(_np.median(samp) if agg == "median" else _np.mean(samp))
    low, high = _np.percentile(stats, q)
    center = _np.median(x) if agg == "median" else _np.mean(x)
    return float(center), float(low), float(high), int(x.size)

# Helper: draw a smooth curve in log–log space using a monotone cubic interpolator (PCHIP).
# Falls back to connecting sorted points if there are too few unique x's or any error.

def plot_loglog_curve(ax, x, y, **kwargs):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x = x[m]; y = y[m]
    if x.size < 2:
        return
    # sort by x
    order = np.argsort(x)
    x = x[order]; y = y[order]
    # if many duplicate x's, average y in log space
    try:
        lx = np.log10(x); ly = np.log10(y)
        uniq, inv = np.unique(lx, return_inverse=True)
        if uniq.size < 3:
            ax.plot(x, y, linestyle='-', **kwargs)
            return
        ly_avg = np.zeros_like(uniq)
        for i in range(uniq.size):
            ly_avg[i] = ly[inv == i].mean()
        p = PchipInterpolator(uniq, ly_avg)
        gx = np.linspace(uniq.min(), uniq.max(), 200)
        gy = p(gx)
        ax.plot(10**gx, 10**gy, **kwargs)
    except Exception:
        # graceful fallback
        ax.plot(x, y, linestyle='-', **kwargs)

def des_frontier(df, out):
    df = df[(df["family"]=="DES") & (df["benchmark_id"]=="des_mm1")].copy()
    if df.empty:
        print("[INFO] No DES rows for plotting."); return
    df["lambda"] = df["lambda"].astype(float)
    df["mu"] = df["mu"].astype(float)
    df["rho"] = df["lambda"] / df["mu"]
    df["L_theory"] = df["rho"] / (1.0 - df["rho"])
    df["L_err_pct"] = (df["L_timeavg"].astype(float) - df["L_theory"]).abs() / df["L_theory"] * 100.0

    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    for (plat, tool), grp in df.groupby(["platform","tool"], dropna=False):
        sc = ax.scatter(grp["wall_time_s"].astype(float), grp["L_err_pct"].astype(float), label=f"{plat}:{tool}")
        # use same color for the curve
        col = sc.get_facecolor()[0]
        plot_loglog_curve(ax,
                          grp["wall_time_s"].astype(float).values,
                          grp["L_err_pct"].astype(float).values,
                          color=col, alpha=0.9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)")
    ax.set_ylabel("Percent error of E[N] vs theory")
    ax.set_title("DES Speed–Accuracy (raw + smooth curve)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")

def ode_frontier(df, out):
    df = df[(df["family"]=="ODE") & (df["benchmark_id"]=="ode_robertson")].copy()
    if df.empty:
        print("[INFO] No ODE rows for plotting."); return
    df = df.dropna(subset=["wall_time_s","rmse"])

    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    for (plat, solver), grp in df.groupby(["platform","solver"], dropna=False):
        sc = ax.scatter(grp["wall_time_s"].astype(float), grp["rmse"].astype(float), label=f"{plat}:{solver}")
        col = sc.get_facecolor()[0]
        plot_loglog_curve(ax,
                          grp["wall_time_s"].astype(float).values,
                          grp["rmse"].astype(float).values,
                          color=col, alpha=0.9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)")
    ax.set_ylabel("RMSE vs reference")
    ax.set_title("ODE Speed–Accuracy Frontier (raw + smooth curve)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")

def ode_scalability(df, out):
    df = df[(df["family"]=="ODE") & (df["benchmark_id"]=="ode_robertson")].copy()
    if df.empty:
        print("[INFO] No ODE rows for plotting."); return
    df = df.dropna(subset=["wall_time_s","t_end"])

    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    for (plat, solver), grp in df.groupby(["platform","solver"], dropna=False):
        # draw raw points and a simple connecting line by t_end to show trend
        sg = grp.sort_values(by="t_end")
        ax.plot(sg["t_end"].astype(float), sg["wall_time_s"].astype(float), marker="o", linestyle="-", label=f"{plat}:{solver}")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("t_end")
    ax.set_ylabel("Wall time (s)")
    ax.set_title("ODE Scalability (raw points)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")

def mem_vs_time(df, out):
    df = df.dropna(subset=["wall_time_s","peak_mem_mb"])
    if df.empty:
        print("[INFO] No rows with memory for plotting."); return
    fig = plt.figure(figsize=(6,4))
    ax = fig.gca()
    for (fam, plat, tool), grp in df.groupby(["family","platform","tool"], dropna=False):
        sc = ax.scatter(grp["wall_time_s"].astype(float), grp["peak_mem_mb"].astype(float), label=f"{fam}:{plat}:{tool}")
        col = sc.get_facecolor()[0]
        plot_loglog_curve(ax,
                          grp["wall_time_s"].astype(float).values,
                          grp["peak_mem_mb"].astype(float).values,
                          color=col, alpha=0.9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)")
    ax.set_ylabel("Peak memory (MB)")
    ax.set_title("Peak memory vs Wall time (raw + smooth curve)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")

# =========================== Aggregated, cleaner plots ==============================

def des_frontier_agg(df, out):
    """DES: %error vs wall time — points are per (platform,tool,horizon),
    using median ±95% CI, label shows replicate count.
    """
    d = df[(df["family"] == "DES") & (df["benchmark_id"] == "des_mm1")].copy()
    if d.empty:
        print("[INFO] No DES rows for des_frontier_agg."); return
    d["lambda"] = _pd.to_numeric(d["lambda"], errors="coerce")
    d["mu"] = _pd.to_numeric(d["mu"], errors="coerce")
    d["rho"] = d["lambda"] / d["mu"]
    d["L_theory"] = d["rho"] / (1.0 - d["rho"])
    d["L_err_pct"] = (d["L_timeavg"].astype(float) - d["L_theory"]).abs() / d["L_theory"] * 100.0
    d["horizon"] = _pd.to_numeric(d["horizon"], errors="coerce")

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6,4)); ax = fig.gca()

    for (plat, tool), grp in d.groupby(["platform","tool"], dropna=False):
        rows = []
        # Aggregate primarily by HORIZON (include NaN so we can detect missing values)
        for h, g in grp.groupby("horizon", dropna=False):
            xm, xlo, xhi, nx = _group_bootstrap_ci(g["wall_time_s"])  # x := wall time
            ym, ylo, yhi, _ = _group_bootstrap_ci(g["L_err_pct"])    # y := % error
            rows.append((h, xm, xlo, xhi, ym, ylo, yhi, nx))
        # Fallback: if horizons are all NaN or rows empty, log-bin by wall time
        if (not rows) or all((r[0] != r[0]) for r in rows):  # NaN != NaN
            g = grp.copy()
            wt = _pd.to_numeric(g["wall_time_s"], errors="coerce").dropna()
            if wt.size >= 2:
                lo, hi = wt.min(), wt.max()
                import numpy as _np
                edges = _np.geomspace(lo, hi, num=5)  # 4 bins
                labels = _pd.cut(g["wall_time_s"], bins=edges, include_lowest=True)
                rows = []
                for b, gb in g.groupby(labels, dropna=False):
                    if gb.empty: continue
                    xm, xlo, xhi, nx = _group_bootstrap_ci(gb["wall_time_s"])  # x := wall time
                    ym, ylo, yhi, _ = _group_bootstrap_ci(gb["L_err_pct"])    # y := % error
                    rows.append((str(b), xm, xlo, xhi, ym, ylo, yhi, nx))
        # Sort: by horizon if numeric, else by median wall time
        def _key(t):
            h = t[0]
            return (1, t[1]) if (h != h) else (0, h)  # NaN last
        rows = sorted(rows, key=_key)
        if not rows:
            continue
        xs = [r[1] for r in rows]; ys = [r[4] for r in rows]
        xerr = [[r[1]-r[2] for r in rows],[r[3]-r[1] for r in rows]]
        yerr = [[r[4]-r[5] for r in rows],[r[6]-r[4] for r in rows]]
        N_total = len(grp)
        ax.errorbar(xs, ys, xerr=xerr, yerr=yerr, fmt="o-", capsize=2.5,
                    elinewidth=0.9, capthick=0.9, markersize=5, linewidth=1.2,
                    label=f"{plat}:{tool} (N={N_total})")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)"); ax.set_ylabel("Percent error of E[N] vs theory")
    ax.set_title("DES Speed–Accuracy (median ±95% CI)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.minorticks_on()
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")


def ode_frontier_agg(df, out):
    """ODE: RMSE vs wall time — points per (platform,solver,rtol),
    median ±95% CI, label shows replicate count.
    """
    d = df[(df["family"] == "ODE") & (df["benchmark_id"] == "ode_robertson")].copy()
    if d.empty:
        print("[INFO] No ODE rows for ode_frontier_agg."); return
    d["rtol"] = _pd.to_numeric(d["rtol"], errors="coerce")

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6,4)); ax = fig.gca()

    for (plat, solver), grp in d.groupby(["platform","solver"], dropna=False):
        rows = []
        for rtol, g in grp.groupby("rtol"):
            xm, xlo, xhi, nx = _group_bootstrap_ci(g["wall_time_s"])  # x := wall time
            ym, ylo, yhi, _ = _group_bootstrap_ci(g["rmse"])         # y := RMSE
            rows.append((rtol, xm, xlo, xhi, ym, ylo, yhi, nx))
        rows = sorted(rows, key=lambda t: (t[0] if _np.isfinite(t[0]) else _np.inf))
        if not rows:
            continue
        xs = [r[1] for r in rows]; ys = [r[4] for r in rows]
        xerr = [[r[1]-r[2] for r in rows],[r[3]-r[1] for r in rows]]
        yerr = [[r[4]-r[5] for r in rows],[r[6]-r[4] for r in rows]]
        N_total = len(grp)
        ax.errorbar(xs, ys, xerr=xerr, yerr=yerr, fmt="o-", capsize=2.5,
                    elinewidth=0.9, capthick=0.9, markersize=5, linewidth=1.2,
                    label=f"{plat}:{solver} (N={N_total})")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)"); ax.set_ylabel("RMSE vs reference")
    ax.set_title("ODE Speed–Accuracy Frontier (median ±95% CI)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.minorticks_on()
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")


def ode_scalability_agg(df, out):
    """ODE scalability: wall time vs t_end (points per platform:solver:t_end).
    Median ±95% CI on wall time; t_end on x.
    """
    d = df[(df["family"] == "ODE") & (df["benchmark_id"] == "ode_robertson")].copy()
    if d.empty:
        print("[INFO] No ODE rows for ode_scalability_agg."); return
    d["t_end"] = _pd.to_numeric(d["t_end"], errors="coerce")

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6,4)); ax = fig.gca()

    for (plat, solver), grp in d.groupby(["platform","solver"], dropna=False):
        rows = []
        for T, g in grp.groupby("t_end"):
            ym, ylo, yhi, n = _group_bootstrap_ci(g["wall_time_s"])  # y := wall time
            rows.append((T, ym, ylo, yhi, n))
        rows = sorted(rows, key=lambda t: (t[0] if _np.isfinite(t[0]) else _np.inf))
        if not rows:
            continue
        xs = [r[0] for r in rows]; ys = [r[1] for r in rows]
        yerr = [[r[1]-r[2] for r in rows],[r[3]-r[1] for r in rows]]
        N_total = len(grp)
        ax.errorbar(xs, ys, yerr=yerr, fmt="o-", capsize=2.5,
                    elinewidth=0.9, capthick=0.9, markersize=5, linewidth=1.2,
                    label=f"{plat}:{solver} (N={N_total})")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("t_end"); ax.set_ylabel("Wall time (s)")
    ax.set_title("ODE Scalability (median ±95% CI)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.minorticks_on()
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")


def mem_vs_time_agg(df, out):
    """Peak memory vs wall time — one point per (family,platform,tool,solver)
    using median ±95% CI on both axes.
    """
    d = df.dropna(subset=["wall_time_s","peak_mem_mb"]) .copy()
    if d.empty:
        print("[INFO] No rows for mem_vs_time_agg."); return
    d["solver"] = d["solver"].fillna("")

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(6,4)); ax = fig.gca()

    for (fam, plat, tool, solver), grp in d.groupby(["family","platform","tool","solver"], dropna=False):
        xm, xlo, xhi, nx = _group_bootstrap_ci(grp["wall_time_s"])  # x := wall time
        ym, ylo, yhi, _  = _group_bootstrap_ci(grp["peak_mem_mb"])   # y := memory
        N_total = len(grp)
        ax.errorbar([xm],[ym], xerr=[[xm-xlo],[xhi-xm]], yerr=[[ym-ylo],[yhi-ym]],
                    fmt="o", capsize=2.5, elinewidth=0.9, capthick=0.9, markersize=5,
                    label=f"{fam}:{plat}:{tool}{(':'+solver) if solver else ''} (N={N_total})")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)"); ax.set_ylabel("Peak memory (MB)")
    ax.set_title("Peak memory vs Wall time (median ±95% CI)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.minorticks_on()
    ax.legend(fontsize=7); fig.tight_layout(); fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")
# ===================================================================================
# Simpler, pattern-forward views (raw points, light jitter)

def _mul_jitter(x, scale=0.03, seed=0):
    import numpy as np
    x = np.asarray(x, float)
    if x.size == 0:
        return x
    rng = np.random.default_rng(seed)
    jitter = 1.0 + scale * (rng.random(x.size) - 0.5) * 2.0  # multiplicative ±scale
    return x * jitter


def des_vs_horizon(df, out):
    """DES: percent error vs HORIZON (raw points with light jitter)."""
    df = df[(df["family"] == "DES") & (df["benchmark_id"] == "des_mm1")].copy()
    if df.empty:
        print("[INFO] No DES rows for des_vs_horizon.")
        return
    df["lambda"] = df["lambda"].astype(float)
    df["mu"] = df["mu"].astype(float)
    df["rho"] = df["lambda"] / df["mu"]
    df["L_theory"] = df["rho"] / (1.0 - df["rho"])
    df["L_err_pct"] = (df["L_timeavg"].astype(float) - df["L_theory"]).abs() / df["L_theory"] * 100.0
    df["horizon"] = df["horizon"].astype(float)

    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    for (plat, tool), grp in df.groupby(["platform", "tool"], dropna=False):
        x = _mul_jitter(grp["horizon"].values, scale=0.04, seed=1)
        y = grp["L_err_pct"].astype(float).values
        ax.scatter(x, y, label=f"{plat}:{tool}", alpha=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Horizon T")
    ax.set_ylabel("Percent error of E[N] vs theory")
    ax.set_title("DES: Error vs Horizon (raw points)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")


def ode_vs_rtol(df, out):
    """ODE Robertson: RMSE vs RTOL (raw points with light jitter)."""
    df = df[(df["family"] == "ODE") & (df["benchmark_id"] == "ode_robertson")].copy()
    if df.empty:
        print("[INFO] No ODE rows for ode_vs_rtol.")
        return
    if "rtol" not in df.columns:
        print("[INFO] No rtol column found.")
        return
    df["rtol"] = df["rtol"].astype(float)

    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    for (plat, solver), grp in df.groupby(["platform", "solver"], dropna=False):
        x = _mul_jitter(grp["rtol"].values, scale=0.03, seed=2)
        y = grp["rmse"].astype(float).values
        ax.scatter(x, y, label=f"{plat}:{solver}", alpha=0.9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Relative tolerance (rtol)")
    ax.set_ylabel("RMSE vs reference")
    ax.set_title("ODE Robertson: RMSE vs rtol (raw points)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    print(f"[OK] wrote {out}")
# SSA frontier (RMSE vs wall time; groups = platform:tool by n_rep)
def ssa_frontier_agg(df, out):
    import json as _json, numpy as _np, pandas as _pd, matplotlib.pyplot as plt
    d = df[(df["family"]=="STOCH") & (df["benchmark_id"]=="ssa_sir")].copy()
    if d.empty: print("[INFO] no STOCH rows."); return
    if "n_rep" not in d.columns:
        def get_nr(s):
            try: return _json.loads(s).get("n_rep")
            except: return _np.nan
        d["n_rep"] = d["params_json"].apply(get_nr)
    fig = plt.figure(figsize=(6,4)); ax = fig.gca()
    for (plat, tool), grp in d.groupby(["platform","tool"], dropna=False):
        rows=[]
        for nrep,g in grp.groupby("n_rep"):
            xm,xlo,xhi,nx = _group_bootstrap_ci(g["wall_time_s"])
            ym,ylo,yhi,_  = _group_bootstrap_ci(g["rmse"])
            rows.append((nrep,xm,xlo,xhi,ym,ylo,yhi,nx))
        rows = sorted(rows, key=lambda t: (t[0] if _np.isfinite(t[0]) else _np.inf))
        if not rows: continue
        xs=[r[1] for r in rows]; ys=[r[4] for r in rows]
        xerr=[[r[1]-r[2] for r in rows],[r[3]-r[1] for r in rows]]
        yerr=[[r[4]-r[5] for r in rows],[r[6]-r[4] for r in rows]]
        N_total = len(grp)
        ax.errorbar(xs,ys,xerr=xerr,yerr=yerr,fmt="o-",capsize=2.5,elinewidth=0.9,capthick=0.9,markersize=5,linewidth=1.2,
                    label=f"{plat}:{tool} (N={N_total})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)"); ax.set_ylabel("RMSE vs reference")
    ax.set_title("SSA SIR: Speed–Accuracy (median ±95% CI)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6); ax.minorticks_on()
    fig.tight_layout(); fig.savefig(out, dpi=200); print(f"[OK] wrote {out}")

# Hybrid frontier (RMSE vs wall time; groups = platform:solver by rtol)
def hybrid_frontier_agg(df, out):
    import json as _json, numpy as _np, pandas as _pd, matplotlib.pyplot as plt
    d = df[(df["family"]=="HYBRID") & (df["benchmark_id"]=="hybrid_bounce")].copy()
    if d.empty: print("[INFO] no HYBRID rows."); return
    if "rtol" not in d.columns:
        def get_rt(s):
            try: return _json.loads(s).get("rtol")
            except: return _np.nan
        d["rtol"] = d["params_json"].apply(get_rt)
    fig = plt.figure(figsize=(6,4)); ax = fig.gca()
    for (plat, solver), grp in d.groupby(["platform","solver"], dropna=False):
        rows=[]
        for rtol,g in grp.groupby("rtol"):
            xm,xlo,xhi,nx = _group_bootstrap_ci(g["wall_time_s"])
            ym,ylo,yhi,_  = _group_bootstrap_ci(g["rmse"])
            rows.append((rtol,xm,xlo,xhi,ym,ylo,yhi,nx))
        rows = sorted(rows, key=lambda t: (t[0] if _np.isfinite(t[0]) else _np.inf))
        if not rows: continue
        xs=[r[1] for r in rows]; ys=[r[4] for r in rows]
        xerr=[[r[1]-r[2] for r in rows],[r[3]-r[1] for r in rows]]
        yerr=[[r[4]-r[5] for r in rows],[r[6]-r[4] for r in rows]]
        N_total = len(grp)
        ax.errorbar(xs,ys,xerr=xerr,yerr=yerr,fmt="o-",capsize=2.5,elinewidth=0.9,capthick=0.9,markersize=5,linewidth=1.2,
                    label=f"{plat}:{solver} (N={N_total})")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Wall time (s)"); ax.set_ylabel("NRMSE vs reference")
    ax.set_title("Hybrid (Bouncing Ball): Speed–Accuracy (median ±95% CI)")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6); ax.minorticks_on()
    fig.tight_layout(); fig.savefig(out, dpi=200); print(f"[OK] wrote {out}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="figs")
    ap.add_argument("csvs", nargs="*")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = load_results(args.csvs)
    if df.empty:
        print("[ERROR] No results_*.csv found or readable.")
        return

    des_frontier_agg(df, os.path.join(args.outdir, "des_frontier_agg.png"))
    ode_frontier_agg(df, os.path.join(args.outdir, "ode_frontier_agg.png"))
    ode_scalability_agg(df, os.path.join(args.outdir, "ode_scalability_agg.png"))
    mem_vs_time_agg(df, os.path.join(args.outdir, "memory_vs_time_agg.png"))
    ssa_frontier_agg(df, os.path.join(args.outdir, "ssa_frontier_agg.png"))
    hybrid_frontier_agg(df, os.path.join(args.outdir, "hybrid_frontier_agg.png"))

if __name__ == "__main__":
    main()