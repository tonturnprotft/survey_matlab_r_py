
# Simulation Tools Comparison — Project Scaffolding

**Generated:** 2025-09-22T04:44:26.336721Z

This folder contains CSV templates for recording results and usability metrics, and a benchmark catalog to standardize parameterization across MATLAB (R2024b/2025a), R (4.4.x), and Python (3.11).

## Files

- `results_schema.csv`: Columns for speed/accuracy/scalability results across all benchmarks/tools.
- `usability_schema.csv`: Columns for proxy usability measures (time-to-first-plot, LoC, errors, etc.).
- `benchmark_catalog.csv`: Canonical list of benchmarks, baseline params, scaling axes, and reference types.

## Suggested directory layout (outside this folder)

```
project-root/
  matlab/
    des_mm1.m
    ode_robertson.m
    hybrid_thermostat.slx
    run_all.m
  R/
    des_mm1.R
    ode_robertson.R
    run_all.R
    renv.lock
  python/
    des_mm1.py
    ode_robertson.py
    run_all.py
    requirements.txt
  results/
    raw_runs/*.csv
    combined/combined_results.csv
  figs/
    speed_accuracy/
    scalability/
    heatmap/
    decision_tree/
```

## Timing & memory (macOS)

Use BSD time for consistent wall time and MaxRSS:

```bash
/usr/bin/time -l python3 python/des_mm1.py      # Max RSS in bytes
/usr/bin/time -l Rscript R/des_mm1.R
/usr/bin/time -l matlab -batch "run('matlab/run_all.m')"
```

Convert MaxRSS bytes → MB in post-processing.

## RNG & reproducibility

- Fix seeds per run; record in `results_schema.csv`.
- R: set.seed(); Python: numpy.random.default_rng(seed); MATLAB: rng(seed,'twister').
- Distributional comparisons across tools for stochastic models; pathwise only within a tool.

## Usability study (N=3)

Each participant uses a different tool on the same three tasks. Log proxy metrics to `usability_schema.csv`. No SUS; report as proxies with limitations.

