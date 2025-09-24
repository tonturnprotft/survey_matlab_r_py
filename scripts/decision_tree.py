# decision_tree.py
# Usage: python decision_tree.py
import matplotlib.pyplot as plt

# Nodes (x,y) on a simple grid; adjust if you want more branches
nodes = {
    "start":      (0, 0),
    "stochastic": (4,  2),   "deterministic": (4, -2),
    "events":     (8,  1),   "no_events":     (8, -1),
    "cps_hybrid": (12, 1.5), "queueing_des":  (12, 0.5),
    "ode_stiff":  (8, -2.5), "ode_smooth":    (8, -3.5),
    # leaves (recommendations)
    "SSA":        (16, 1.8),  "DES":          (16, 0.2),
    "ODE_stiff":  (16,-2.3),  "ODE_regular":  (16,-3.7),
}
edges = [
    ("start","stochastic","uncertainty in dynamics?"),
    ("start","deterministic","no"),
    ("stochastic","events","discrete events / resources?"),
    ("stochastic","no_events","no"),
    ("events","cps_hybrid","hybrid continuous+events?"),
    ("events","queueing_des","no"),
    ("deterministic","ode_stiff","stiff / wide scale sep?"),
    ("deterministic","ode_smooth","no"),
    # to leaves
    ("cps_hybrid","SSA","if stochastic at micro-level → SSA (else Hybrid)"),
    ("queueing_des","DES","DES (SimPy / simmer / SimEvents)"),
    ("ode_stiff","ODE_stiff","BDF/LSODA/ode15s"),
    ("ode_smooth","ODE_regular","RK/DOP853/ode45"),
]

def draw_node(ax, name, xy, kind="q"):
    rx, ry = 1.4, 0.6
    x, y = xy
    if name in {"SSA","DES","ODE_stiff","ODE_regular"}:
        box = plt.Rectangle((x-rx, y-ry), 2*rx, 2*ry, fill=False, lw=1.5)
        ax.add_patch(box)
    else:
        circ = plt.Circle((x,y), 0.6, fill=False, lw=1.5)
        ax.add_patch(circ)
    ax.text(x, y, name, ha="center", va="center", fontsize=10)

def arrow(ax, a, b, label=None):
    (x1,y1),(x2,y2)=nodes[a],nodes[b]
    ax.annotate("", xy=(x2,y2), xytext=(x1,y1),
                arrowprops=dict(arrowstyle="->", lw=1.2))
    if label:
        xm, ym = (x1+x2)/2, (y1+y2)/2
        ax.text(xm, ym+0.4, label, ha="center", va="bottom", fontsize=9)

fig, ax = plt.subplots(figsize=(9, 4.8))
for n,xy in nodes.items(): draw_node(ax, n, xy)
for a,b,lbl in edges: arrow(ax, a, b, lbl)

ax.set_xlim(-1, 17); ax.set_ylim(-5, 3)
ax.axis("off")
ax.set_title("Practitioner Decision Tree")
fig.tight_layout()
plt.savefig("decision_tree.png", dpi=200)
plt.show()