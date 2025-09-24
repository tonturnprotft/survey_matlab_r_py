# ode_robertson.R — Solve Robertson ODE with deSolve and emit one CSV row.
# Run:
#   Rscript ode_robertson.R 1e5 1e-6 1e-8 lsoda results_ode_robertson_R.csv ref_robertson_t1e5.csv
suppressWarnings(suppressMessages({
  library(deSolve)
}))

args <- commandArgs(trailingOnly = TRUE)
t_end <- if (length(args) >= 1) as.numeric(args[1]) else 1e5
rtol  <- if (length(args) >= 2) as.numeric(args[2]) else 1e-6
atol  <- if (length(args) >= 3) as.numeric(args[3]) else 1e-8
method <- if (length(args) >= 4) args[4] else "lsoda"
out     <- if (length(args) >= 5) args[5] else "results_ode_robertson_R.csv"
refcsv  <- if (length(args) >= 6) args[6] else ""

robertson <- function(t, y, parms) {
  y1 <- y[1]; y2 <- y[2]; y3 <- y[3]
  dy1 <- -0.04*y1 + 1.0e4*y2*y3
  dy2 <-  0.04*y1 - 1.0e4*y2*y3 - 3.0e7*(y2^2)
  dy3 <-  3.0e7*(y2^2)
  list(c(dy1,dy2,dy3))
}

sample_times <- function(t_end, n=500) {
  if (t_end <= 1e-6) return(seq(0, t_end, length.out=n))
  c(0, exp(seq(log(1e-6), log(t_end), length.out=n-1)))
}

times <- sample_times(t_end, n=500)
y0 <- c(y1=1, y2=0, y3=0)

t0 <- proc.time()
sol <- ode(y=y0, times=times, func=robertson, parms=NULL, method=method, rtol=rtol, atol=atol)
elapsed <- (proc.time() - t0)[["elapsed"]]
cpu_time <- sum((proc.time() - t0)[c("user.self","sys.self")])

status <- "OK"
if (any(!is.finite(sol[,"y1"])) || any(is.na(sol))) status <- "FAIL"

rmse <- ""; linf <- ""
if (refcsv != "" && file.exists(refcsv) && status == "OK") {
  ref <- read.csv(refcsv)
  # Interpolate ref onto our time grid
  y1r <- approx(ref$t, ref$y1, xout = sol[,"time"], rule = 2)$y
  y2r <- approx(ref$t, ref$y2, xout = sol[,"time"], rule = 2)$y
  y3r <- approx(ref$t, ref$y3, xout = sol[,"time"], rule = 2)$y
  diff1 <- sol[,"y1"] - y1r
  diff2 <- sol[,"y2"] - y2r
  diff3 <- sol[,"y3"] - y3r
  diffsq <- diff1^2 + diff2^2 + diff3^2
  rmse <- sqrt(mean(diffsq))
  linf <- max(abs(c(diff1,diff2,diff3)))
}

params_json <- sprintf('{"t_end":%s,"rtol":%s,"atol":%s,"method":"%s"}', t_end, rtol, atol, method)

# Prepare row (results_schema + L_timeavg column retained for compatibility)
fields <- c("family","benchmark_id","platform","tool","version","solver","step_mode",
            "abs_tol","rel_tol","dt_init","dt_max","seed","params_json",
            "scale_param","scale_value","repeat_idx","run_id","start_ts","end_ts",
            "wall_time_s","cpu_time_s","peak_mem_mb","n_steps","n_events","rmse","linf",
            "event_time_err_mean","event_time_err_max","throughput","L_timeavg","status","error_message","host")

row <- list(
  family="ODE",
  benchmark_id="ode_robertson",
  platform="R",
  tool="deSolve",
  version=as.character(packageVersion("deSolve")),
  solver=method,
  step_mode="adaptive",
  abs_tol=atol,
  rel_tol=rtol,
  dt_init="",
  dt_max="",
  seed="",
  params_json=params_json,
  scale_param="t_end",
  scale_value=t_end,
  repeat_idx=0,
  run_id=paste0("r_", as.integer(as.numeric(Sys.time()))),
  start_ts="",
  end_ts="",
  wall_time_s=elapsed,
  cpu_time_s=cpu_time,
  peak_mem_mb="",
  n_steps=length(times)-1,
  n_events="",
  rmse=rmse,
  linf=linf,
  event_time_err_mean="",
  event_time_err_max="",
  throughput="",
  L_timeavg="",
  status=status,
  error_message="",
  host=R.version.string
)

df <- as.data.frame(row, stringsAsFactors = FALSE)
write_header <- !file.exists(out)
suppressWarnings(write.table(df, file=out, sep=",", row.names=FALSE, col.names=write_header, append=!write_header, qmethod="double", quote=TRUE))

cat(sprintf("status=%s n_steps=%d wall=%.3fs rmse=%s linf=%s\n", status, length(times)-1, elapsed, as.character(rmse), as.character(linf)))