# des_mm1.R — M/M/1 simulation with simmer, writes one CSV row.
# Run: Rscript des_mm1.R 0.9 1.0 100000 42 results_des_mm1_R.csv
suppressWarnings(suppressMessages({
  library(simmer)
  library(simmer.plot)
}))

args <- commandArgs(trailingOnly = TRUE)
lam     <- if (length(args) >= 1) as.numeric(args[1]) else 0.9
mu      <- if (length(args) >= 2) as.numeric(args[2]) else 1.0
horizon <- if (length(args) >= 3) as.numeric(args[3]) else 100000
seed    <- if (length(args) >= 4) as.integer(args[4]) else 42
out     <- if (length(args) >= 5) args[5] else "results_des_mm1_R.csv"

set.seed(seed)

traj <- trajectory("cust") %>%
  seize("server", 1) %>%
  timeout(function() rexp(1, mu)) %>%
  release("server", 1)

env <- simmer("mm1") %>%
  add_resource("server", capacity=1, queue_size=Inf) %>%
  add_generator("cust", traj, function() rexp(1, lam))

t0 <- proc.time()
env <- run(env, until=horizon)
elapsed <- (proc.time() - t0)[["elapsed"]]
cpu_time <- sum((proc.time() - t0)[c("user.self","sys.self")])

arr <- get_mon_arrivals(env)
finished <- sum(arr$finished)
throughput <- finished / horizon

events_processed <- finished*2L  # crude: arrivals+departures

# Compute time-average number in system E[N] via resource timeline
res <- get_mon_resources(env)
res_srv <- subset(res, resource == "server", select = c(time, server, queue))
res_srv <- res_srv[order(res_srv$time), ]
# Ensure timeline starts at 0 with known state
if (nrow(res_srv) == 0 || res_srv$time[1] > 0) {
  res_srv <- rbind(data.frame(time = 0, server = 0, queue = 0), res_srv)
}
# Close the interval at horizon
last_server <- tail(res_srv$server, 1)
last_queue  <- tail(res_srv$queue, 1)
res_srv <- rbind(res_srv, data.frame(time = horizon, server = last_server, queue = last_queue))
# Integrate N(t) = server + queue over time
n_vals <- head(res_srv$server + res_srv$queue, -1)
dt     <- diff(res_srv$time)
area   <- sum(n_vals * dt)
L_timeavg <- as.numeric(area / horizon)

rho <- lam / mu
params_json <- sprintf('{"lambda":%s,"mu":%s,"horizon":%s}', lam, mu, horizon)

# Append CSV row (schema compatible with results_schema.csv)
fields <- c("family","benchmark_id","platform","tool","version","solver","step_mode",
            "abs_tol","rel_tol","dt_init","dt_max","seed","params_json",
            "scale_param","scale_value","repeat_idx","run_id","start_ts","end_ts",
            "wall_time_s","cpu_time_s","peak_mem_mb","n_steps","n_events","rmse","linf",
            "event_time_err_mean","event_time_err_max","throughput","L_timeavg","status","error_message","host")

row <- list(
  family="DES",
  benchmark_id="des_mm1",
  platform="R",
  tool="simmer",
  version=as.character(packageVersion("simmer")),
  solver="event_loop",
  step_mode="event",
  abs_tol="",
  rel_tol="",
  dt_init="",
  dt_max="",
  seed=seed,
  params_json=params_json,
  scale_param="rho",
  scale_value=rho,
  repeat_idx=0,
  run_id=paste0("r_", as.integer(as.numeric(Sys.time()))),
  start_ts="",
  end_ts="",
  wall_time_s=elapsed,
  cpu_time_s=cpu_time,
  peak_mem_mb="",
  n_steps=events_processed,
  n_events=finished,
  rmse="",
  linf="",
  event_time_err_mean="",
  event_time_err_max="",
  throughput=throughput,
  L_timeavg=L_timeavg,
  status="OK",
  error_message="",
  host=R.version.string
)

df <- as.data.frame(row, stringsAsFactors = FALSE)
write_header <- !file.exists(out)
suppressWarnings(write.table(df, file=out, sep=",", row.names=FALSE, col.names=write_header, append=!write_header, qmethod="double"))

cat(sprintf("rho=%.3f finished=%d L_timeavg=%.3f thr=%.4f wall=%.3fs\n", rho, finished, L_timeavg, throughput, elapsed))