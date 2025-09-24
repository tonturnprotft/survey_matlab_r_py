# des_mm1.py — M/M/1 discrete-event simulation (SimPy) producing one CSV row
# UPDATED: saves L_timeavg to CSV
# Usage:
#   python3 des_mm1.py --lam 0.9 --mu 1.0 --horizon 100000 --seed 42 --out results_des_mm1_python_v2.csv
import argparse, time, csv, json, platform, os
import numpy as np
import simpy

def mm1_sim(lam, mu, horizon, seed=42):
    rng = np.random.default_rng(seed)
    env = simpy.Environment()
    server = simpy.Resource(env, capacity=1)
    n_system = 0
    last_t = 0.0
    area = 0.0
    arrivals = 0
    departures = 0
    events_processed = 0

    def update_area():
        nonlocal area, last_t
        now = env.now
        area += n_system * (now - last_t)
        last_t = now

    def customer():
        nonlocal n_system, departures, events_processed
        update_area()
        n_system += 1
        events_processed += 1  # arrival event
        with server.request() as req:
            yield req
            service = rng.exponential(1.0 / mu)
            yield env.timeout(service)
            update_area()
            n_system -= 1
            departures += 1
            events_processed += 1  # departure event

    def arrival_proc():
        nonlocal arrivals, events_processed
        t = 0.0
        while True:
            iat = rng.exponential(1.0 / lam)
            t += iat
            if t > horizon:
                break
            yield env.timeout(t - env.now)
            env.process(customer())
            arrivals += 1

    env.process(arrival_proc())
    env.run(until=horizon)
    if last_t < horizon:
        area += n_system * (horizon - last_t)

    throughput = departures / horizon
    l_time_avg = area / horizon
    return {
        "arrivals": arrivals,
        "departures": departures,
        "events_processed": events_processed,
        "L_timeavg": l_time_avg,
        "throughput": throughput,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", type=float, default=0.9)
    ap.add_argument("--mu", type=float, default=1.0)
    ap.add_argument("--horizon", type=float, default=100000.0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=str, default="results_des_mm1_python_v2.csv")
    args = ap.parse_args()

    params = {"lambda": args.lam, "mu": args.mu, "horizon": args.horizon}

    t0_wall = time.perf_counter()
    t0_cpu = time.process_time()
    res = mm1_sim(args.lam, args.mu, args.horizon, seed=args.seed)
    wall = time.perf_counter() - t0_wall
    cpu = time.process_time() - t0_cpu

    rho = args.lam / args.mu
    run_id = f"py_{int(time.time())}"

    row = {
        "family": "DES",
        "benchmark_id": "des_mm1",
        "platform": "python",
        "tool": "simpy",
        "version": getattr(simpy, "__version__", ""),
        "solver": "event_loop",
        "step_mode": "event",
        "abs_tol": "",
        "rel_tol": "",
        "dt_init": "",
        "dt_max": "",
        "seed": args.seed,
        "params_json": json.dumps(params),
        "scale_param": "rho",
        "scale_value": rho,
        "repeat_idx": 0,
        "run_id": run_id,
        "start_ts": "",
        "end_ts": "",
        "wall_time_s": wall,
        "cpu_time_s": cpu,
        "peak_mem_mb": "",             # to be filled by wrapper
        "n_steps": res["events_processed"],
        "n_events": res["departures"],
        "rmse": "",
        "linf": "",
        "event_time_err_mean": "",
        "event_time_err_max": "",
        "throughput": res["throughput"],
        "L_timeavg": res["L_timeavg"], # NEW
        "status": "OK",
        "error_message": "",
        "host": platform.platform(),
    }

    field_order = [
        "family","benchmark_id","platform","tool","version","solver","step_mode",
        "abs_tol","rel_tol","dt_init","dt_max","seed","params_json",
        "scale_param","scale_value","repeat_idx","run_id","start_ts","end_ts",
        "wall_time_s","cpu_time_s","peak_mem_mb","n_steps","n_events","rmse","linf",
        "event_time_err_mean","event_time_err_max","throughput","L_timeavg","status","error_message","host"
    ]

    exists = os.path.exists(args.out)
    # If file exists but header lacks L_timeavg, stop to avoid corrupt append.
    if exists:
        try:
            with open(args.out, "r") as f:
                header = f.readline().strip().split(",")
            if "L_timeavg" not in header:
                raise SystemExit(f"[ERROR] Output file '{args.out}' exists without 'L_timeavg' column. "
                                 f"Please use a NEW out file (e.g., results_des_mm1_python_v2.csv) or regenerate the header.")
        except Exception as e:
            pass

    with open(args.out, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=field_order)
        if not exists:
            w.writeheader()
        w.writerow(row)

    print(f"RUN_ID={run_id}")
    print(f"rho={rho:.3f}, departures={res['departures']}, L_timeavg≈{res['L_timeavg']:.3f}, thr≈{res['throughput']:.4f}, wall={wall:.3f}s")

if __name__ == "__main__":
    main()