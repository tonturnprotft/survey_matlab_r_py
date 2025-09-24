# run_with_time.py — Run a command via `/usr/bin/time -l`, parse MaxRSS, and fill peak_mem_mb in the last CSV row.
# Usage:
#   python3 run_with_time.py --out results_des_mm1_python_v2.csv -- python3 python/des_mm1.py --lam 0.9 --mu 1.0 --horizon 100000 --seed 42 --out results_des_mm1_python_v2.csv
import argparse, subprocess, sys, re, pandas as pd, math
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description="Run a command with BSD time and inject peak_mem_mb into CSV last row.")
    ap.add_argument("--out", required=True, help="CSV file to update (must be the same one the child writes to).")
    ap.add_argument("cmd", nargs=argparse.REMAINDER, help="Command to run after '--'")
    args = ap.parse_args()

    if not args.cmd or args.cmd[0] != "--":
        print("[ERROR] Separate command with '--'. Example: --out results.csv -- python3 script.py ...", file=sys.stderr)
        sys.exit(2)
    cmd = args.cmd[1:]

    # Run command under /usr/bin/time -l
    time_cmd = ["/usr/bin/time", "-l"] + cmd
    proc = subprocess.Popen(time_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = proc.communicate()
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)
    rc = proc.returncode
    if rc != 0:
        sys.exit(rc)

    # Parse MaxRSS or peak memory footprint
    maxrss_bytes = None
    m1 = re.search(r"^\s*(\d+)\s+maximum resident set size", stderr, re.MULTILINE|re.IGNORECASE)
    if m1:
        maxrss_bytes = int(m1.group(1))
    else:
        m2 = re.search(r"^\s*(\d+)\s+peak memory footprint", stderr, re.MULTILINE|re.IGNORECASE)
        if m2:
            maxrss_bytes = int(m2.group(1))

    if maxrss_bytes is None:
        print("[WARN] Could not parse MaxRSS from /usr/bin/time -l output; not updating CSV.", file=sys.stderr)
        sys.exit(0)

    peak_mb = maxrss_bytes / (1024*1024)

    csv_path = Path(args.out)
    if not csv_path.exists():
        print(f"[WARN] CSV '{csv_path}' not found; skipping update.", file=sys.stderr)
        sys.exit(0)

    # Update last row's peak_mem_mb
    df = pd.read_csv(csv_path)
    if "peak_mem_mb" not in df.columns:
        df["peak_mem_mb"] = float("nan")
    if len(df) == 0:
        print("[WARN] CSV is empty; nothing to update.", file=sys.stderr)
        sys.exit(0)
    df.loc[df.index[-1], "peak_mem_mb"] = round(peak_mb, 2)
    df.to_csv(csv_path, index=False)
    print(f"[INFO] Updated {csv_path.name} last row: peak_mem_mb={round(peak_mb,2)} MB")

if __name__ == "__main__":
    main()