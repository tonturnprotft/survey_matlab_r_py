#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 8) stop("Usage: Rscript stoch_sir.R beta gamma N I0 t_end n_rep seed out [ref_csv]")
beta=as.numeric(args[1]); gamma=as.numeric(args[2]); N=as.integer(args[3]); I0=as.integer(args[4])
t_end=as.numeric(args[5]); n_rep=as.integer(args[6]); seed=as.integer(args[7]); out=args[8]
ref <- if (length(args)>=9) args[9] else NA
set.seed(seed); S0 <- N - I0
sir_once <- function(beta,gamma,S,I,R,t_end){
  t<-0; times<-c(0); Sx<-c(S); Ix<-c(I); Rx<-c(R)
  while(t<t_end && I>0){
    a1<-beta*S*I; a2<-gamma*I; a0<-a1+a2; if(a0<=0) break
    tau <- -log(runif(1))/a0; t <- t + tau
    if (runif(1)*a0 < a1 && S>0) { S<-S-1; I<-I+1 } else { I<-I-1; R<-R+1 }
    times<-c(times, min(t,t_end)); Sx<-c(Sx,S); Ix<-c(Ix,I); Rx<-c(Rx,R)
  }
  if(tail(times,1)<t_end){ times<-c(times,t_end); Sx<-c(Sx,S); Ix<-c(Ix,I); Rx<-c(Rx,R) }
  list(times=times, states=cbind(Sx,Ix,Rx))
}
series_to_grid <- function(times, states, grid){
  cbind( approx(times,states[,1],xout=grid,rule=2)$y,
         approx(times,states[,2],xout=grid,rule=2)$y,
         approx(times,states[,3],xout=grid,rule=2)$y )
}
grid <- seq(0,t_end,by=1); acc <- matrix(0.0,nrow=length(grid),ncol=3); steps_total<-0L
t_start <- proc.time()[["elapsed"]]
for(r in 1:n_rep){ sim<-sir_once(beta,gamma,S0,I0,0L,t_end); steps_total <- steps_total + (length(sim$times)-1L); acc <- acc + series_to_grid(sim$times,sim$states,grid) }
mean_series <- acc/n_rep; rmse<-NA; linf<-NA
if(!is.na(ref) && file.exists(ref)){ refdf<-read.csv(ref); diff <- mean_series - as.matrix(refdf[,c("S_mean","I_mean","R_mean")]); rmse <- sqrt(mean(diff^2)); linf <- max(abs(diff)) }
wall <- proc.time()[["elapsed"]] - t_start
row <- data.frame(platform="R",family="STOCH",benchmark_id="ssa_sir",tool="event_loop",solver="gillespie",status="OK",
  wall_time_s=wall, rmse=rmse, linf=linf, n_steps=steps_total,
  params_json=jsonlite::toJSON(list(beta=beta,gamma=gamma,N=N,I0=I0,t_end=t_end,n_rep=n_rep,seed=seed), auto_unbox=TRUE))
hdr <- !file.exists(out)
suppressWarnings(write.table(row, file=out, sep=",", row.names=FALSE, col.names=hdr, append=!hdr, qmethod="double"))
cat(sprintf("status=OK n_steps=%d wall=%.3fs rmse=%s linf=%s\n", steps_total, wall, as.character(rmse), as.character(linf)))