# feature_heatmap.py
# Usage: python feature_heatmap.py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Edit this small matrix in-code as needed.
# 1 = supported / good; 0 = not supported / awkward; 0.5 = partial.
rows = [
    # tool,                         DES,     SSA,        ODE_stiff, Events(Hybrid), CSV_IO, Repro_scripts, Ecosystem_pkgs
    ("R:simmer / deSolve",          1,       1,          1,         1,              1,      1,             1),
    ("MATLAB:SimEvents/ode15s",     1,       0.5,        1,         1,              1,      1,             1),
    ("Python:simpy/scipy",          1,       1,          1,         1,              1,      1,             1),
]
df = pd.DataFrame(rows, columns=["Tool","DES","SSA","ODE (stiff)","Hybrid events","CSV I/O","Repro","Ecosystem"]).set_index("Tool")

fig, ax = plt.subplots(figsize=(7.6, 2.8))
data = df.values.astype(float)
im = ax.imshow(data, aspect="auto", interpolation="nearest")
ax.set_xticks(range(df.shape[1])); ax.set_xticklabels(df.columns, rotation=25, ha="right")
ax.set_yticks(range(df.shape[0])); ax.set_yticklabels(df.index)
# put numeric labels
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        s = {0:"–", 0.5:"±", 1:"✓"}.get(data[i,j], f"{data[i,j]:.1f}")
        ax.text(j, i, s, ha="center", va="center", fontsize=10)
ax.set_title("Feature / Interop Heatmap")
fig.tight_layout()
plt.savefig("feature_heatmap.png", dpi=200)
plt.show()