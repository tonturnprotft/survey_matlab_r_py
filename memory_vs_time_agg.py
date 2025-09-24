# memory_vs_time_agg.py (paste into a temp cell/script and run from /Users/tree/Simu)
import os, json, numpy as np, pandas as pd, matplotlib.pyplot as plt
from pathlib import Path

def boot_med_ci(x, n_boot=1000, seed=0):
    x = pd.to_numeric(pd.Series(x), errors="coerce").dropna().to_numpy()
    if x.size == 0: return np.nan, np.nan, np.nan, 0
    rng = np.random.default_rng(seed)
    meds = np.array([np.median(rng.choice(x, size=x.size, replace=True)) for _ in range(n_boot)])
    lo, hi = np.percentile(meds, [2.5, 97.5])
    return float(np.median(x)), float(lo), float(hi), int(x.size)

# gather all results
csvs = [
  "results_des_mm1_python_v2.csv", "results_des_mm1_R_v2.csv", "results_des_mm1_matlab_v2.csv",
  "results_ode_robertson_python.csv", "results_ode_robertson_R.csv", "results_ode_robertson_matlab.csv",
  "results_ssa_sir_python.csv", "results_ssa_sir_R.csv", "results_ssa_sir_matlab.csv",
  "results_hybrid_bounce_python.csv", "results_hybrid_bounce_R.csv", "results_hybrid_bounce_matlab.csv",
]
frames=[]
for c in csvs:
    p = Path(c)
    if p.exists():
        try:
            df = pd.read_csv(p)
            frames.append(df)
        except Exception as e:
            print(f"[WARN] {c}: {e}")
if not frames:
    raise SystemExit("No CSVs found (run from your Simu root).")
df = pd.concat(frames, ignore_index=True)

# basic cleanup
for col in ["wall_time_s","peak_mem_mb"]:
    if col in df.columns: df[col] = pd.to_numeric(df[col], errors="coerce")
df["platform"] = df["platform"].astype(str)
df["tool"] = df["tool"].astype(str)
df["family"] = df["family"].astype(str)

# aggregate per (FAMILY:platform:tool)
rows=[]
for (fam, plat, tool), g in df.groupby(["family","platform","tool"], dropna=False):
    x_med, x_lo, x_hi, nx = boot_med_ci(g["wall_time_s"])
    y_med, y_lo, y_hi, ny = boot_med_ci(g["peak_mem_mb"])
    if not np.isfinite(x_med) or not np.isfinite(y_med): continue
    rows.append({
        "label": f"{fam}:{plat}:{tool}",
        "fam": fam, "plat": plat, "tool": tool,
        "N": len(g),
        "x": x_med, "xlo": x_lo, "xhi": x_hi,
        "y": y_med, "ylo": y_lo, "yhi": y_hi
    })
agg = pd.DataFrame(rows)
# nice ordering
fam_order = ["DES","ODE","STOCH","HYBRID"]  # STOCH = SSA
agg["fam"] = pd.Categorical(agg["fam"], fam_order)
agg = agg.sort_values(["fam","plat","tool"])

# plot
fig, ax = plt.subplots(figsize=(7.2,4.2))
for (fam), g in agg.groupby("fam"):
    xs=g["x"].to_numpy(); ys=g["y"].to_numpy()
    xerr=[g["x"]-g["xlo"], g["xhi"]-g["x"]]; yerr=[g["y"]-g["ylo"], g["yhi"]-g["y"]]
    ax.errorbar(xs, ys, xerr=xerr, yerr=yerr, fmt="o", capsize=2.5, elinewidth=0.9, capthick=0.9, markersize=5,
                label=fam)

# add text labels (small), one per point
for _,r in agg.iterrows():
    ax.annotate(f"{r['plat']}:{r['tool']} (N={r['N']})", (r["x"], r["y"]),
                textcoords="offset points", xytext=(4,4), fontsize=8)

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Wall time (s)")
ax.set_ylabel("Peak memory (MB)")
ax.set_title("Peak memory vs Wall time (median ±95% CI)")
ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
ax.legend(title="Family", fontsize=9)
fig.tight_layout()
plt.show()
# If you want to save: fig.savefig("memory_vs_time_agg.png", dpi=200)