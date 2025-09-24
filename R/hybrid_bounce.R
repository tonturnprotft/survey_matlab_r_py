#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 7) stop("Usage: Rscript hybrid_bounce.R t_end rtol atol e h0 v0 out [ref_csv]")
t_end=as.numeric(args[1]); rtol=as.numeric(args[2]); atol=as.numeric(args[3])
e=as.numeric(args[4]); h0=as.numeric(args[5]); v0=as.numeric(args[6]); out=args[7]
ref <- if (length(args)>=8) args[8] else NA
library(deSolve)
# --- Helpers for piecewise-aligned RMSE (avoid phase error at impacts) ---
.detect_events <- function(t, h, v){
  ev <- c()
  n <- length(t)
  if (n < 2) return(ev)
  for (i in 1:(n-1)){
    h0 <- h[i]; h1 <- h[i+1]
    if ((h0 > 0 && h1 <= 0) || (h0 >= 0 && h1 < 0)){
      alpha <- h0 / (h0 - h1)
      tz <- t[i] + alpha * (t[i+1] - t[i])
      ev <- c(ev, tz)
    }
  }
  if (length(ev) == 0) return(ev)
  sort(unique(round(ev, 12)))
}

.piecewise_aligned_rmse <- function(t_cand, h_cand, v_cand, t_ref, h_ref, v_ref, t_end){
  ev <- .detect_events(t_ref, h_ref, v_ref)
  ev <- ev[ev > 0 & ev < t_end]
  bounds <- c(0.0, ev, t_end)
  sse_h <- 0.0; sse_v <- 0.0; n <- 0L; linf_norm <- 0.0
  if (length(bounds) < 2) return(list(rmse=NA_real_, linf=NA_real_))

  # reference scales
  h_scale <- max(1e-12, max(h_ref, na.rm=TRUE) - min(h_ref, na.rm=TRUE))
  v_scale <- max(1e-12, max(abs(v_ref), na.rm=TRUE))

  for (k in 1:(length(bounds)-1)){
    a <- bounds[k]; b <- bounds[k+1]
    if (b <= a) next
    eps <- min(1e-6, 1e-6*t_end)
    aa <- a + eps; bb <- b - eps
    if (bb <= aa){ aa <- a; bb <- b }
    grid <- seq(aa, bb, length.out=201)
    hc <- approx(t_cand, h_cand, xout=grid, rule=2)$y
    vc <- approx(t_cand, v_cand, xout=grid, rule=2)$y
    hrefg <- approx(t_ref, h_ref, xout=grid, rule=2)$y
    vrefg <- approx(t_ref, v_ref, xout=grid, rule=2)$y
    dh <- hc - hrefg; dv <- vc - vrefg
    sse_h <- sse_h + sum(dh^2)
    sse_v <- sse_v + sum(dv^2)
    n <- n + length(dh)
    linf_norm <- max(linf_norm, max(c(max(abs(dh))/h_scale, max(abs(dv))/v_scale), na.rm=TRUE))
  }
  if (n == 0L) return(list(rmse=NA_real_, linf=NA_real_))
  rmse_h <- sqrt(sse_h/n); rmse_v <- sqrt(sse_v/n)
  nrmse <- 0.5 * (rmse_h/h_scale + rmse_v/v_scale)
  list(rmse=nrmse, linf=linf_norm)
}
# -------------------------------------------------------------------------
G <- 9.81
derivs <- function(t,y,parms) list(c(y[2], -G))
rootfun <- function(t,y,parms) y[1]
eventfun <- function(t,y,parms){ y[1]<-0; y[2] <- -parms$e*y[2]; return(y) }
parms <- list(e=e); y0 <- c(h0,v0)
t_start <- proc.time()[["elapsed"]]
sol <- ode(y=y0, times=c(0,t_end), func=derivs, parms=parms, method="lsoda", rtol=rtol, atol=atol, rootfun=rootfun, events=list(func=eventfun, root=TRUE))
grid <- seq(0,t_end,length.out=1001)
h <- approx(sol[,"time"], sol[,"1"], xout=grid, rule=2)$y
v <- approx(sol[,"time"], sol[,"2"], xout=grid, rule=2)$y
rmse<-NA; linf<-NA
if(!is.na(ref) && file.exists(ref)){
  refdf <- read.csv(ref)
  pa <- .piecewise_aligned_rmse(sol[,"time"], sol[,"1"], sol[,"2"], refdf$t, refdf$h, refdf$v, t_end)
  rmse <- pa$rmse; linf <- pa$linf
}
wall <- proc.time()[["elapsed"]] - t_start
row <- data.frame(platform="R",family="HYBRID",benchmark_id="hybrid_bounce",tool="deSolve",solver="lsoda",status="OK",
  wall_time_s=wall, rmse=rmse, linf=linf, n_steps=nrow(sol),
  params_json=jsonlite::toJSON(list(t_end=t_end,rtol=rtol,atol=atol,e=e,h0=h0,v0=v0), auto_unbox=TRUE))
hdr <- !file.exists(out)
suppressWarnings(write.table(row, file=out, sep=",", row.names=FALSE, col.names=hdr, append=!hdr, qmethod="double"))
cat(sprintf("status=OK steps=%d wall=%.3fs rmse=%s linf=%s\n", nrow(sol), wall, as.character(rmse), as.character(linf)))
