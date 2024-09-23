### For piecewise-constant signals

library(devtools)
library(IDetect)
library(wbs)
library(breakfast)
library(changepoint)
library(changepoints)
library(changepoint.np)
library(not)
library(mscp)
library(mosum)

# The necessary mcode for seedBS can be downloaded from https://github.com/kovacssolt/SeedBinSeg

# The following code requires first running the R file Algorithms.R

options(expressions = 500000)

sim_thr <- function(x, const = 1.2, points, qmax_NOT, delta, verbose){
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric",det_time="numeric", time="numeric", mean="numeric"), 
           prototype(cpt=numeric(0), nocpt=0, det_time=numeric(0), time=numeric(0), mean=numeric(0)))
  
  if(verbose){print("DAIS")}
  DAIS <- new("cpt.est")
  start.time <- Sys.time()
  z <- DAIS(x, contrast = "mean", thr_const = const, points=points)
  end_time <- Sys.time()
  if(length(z) == 0){DAIS@cpt = 0
  DAIS@nocpt = 0}
  else{DAIS@cpt <- z
  DAIS@nocpt <- length(z)}
  DAIS@time <- as.numeric(end_time - start.time)
  DAIS@mean <- mean.from.cpt(x, DAIS@cpt)
  
  if(verbose){print("ID")}
  ID <- new("cpt.est")
  start.time <- Sys.time()
  z <- ID(x)
  end_time <- Sys.time()
  if(z$no_cpt == 0){ID@cpt = 0
  ID@nocpt = 0}
  else{ID@cpt <- z$cpt
  ID@nocpt <- z$no_cpt}
  ID@time <- as.numeric(end_time - start.time)
  ID@mean <- mean.from.cpt(x, ID@cpt)
  
  if(verbose){print("ID_th")}
  ID_th <- new("cpt.est")
  start.time <- Sys.time()
  z <- pcm_th(x, points = points, thr_const = 1.1)
  end_time <- Sys.time()
  if(length(z) == 0){ID_th@cpt = 0
  ID_th@nocpt = 0}
  else{ID_th@cpt <- z
  ID_th@nocpt <- length(z)}
  ID_th@time <- as.numeric(end_time - start.time)
  ID_th@mean <- mean.from.cpt(x, ID_th@cpt)
  
  if(verbose){print("ID_sic")}
  ID_sic <- new("cpt.est")
  start.time <- Sys.time()
  z <- pcm_ic(x)
  end_time <- Sys.time()
  if(any(is.na(z$cpt_ic$sic_pen))){ID_sic@cpt = 0
  ID_sic@nocpt = 0}
  else{ID_sic@cpt <- z$cpt_ic$sic_pen
  ID_sic@nocpt <- length(z$cpt_ic$sic_pen)}
  ID_sic@time <- as.numeric(end_time - start.time)
  ID_sic@mean <- mean.from.cpt(x, ID_sic@cpt)
  
  if(verbose){print("wbs")}
  WBS_th <- new("cpt.est")
  start.time <- Sys.time()
  z <- model.thresh(sol.wbs(x), th_const = 1)
  end_time <- Sys.time()
  WBS_th@cpt = z$cpts
  WBS_th@nocpt <- z$no.of.cpt
  WBS_th@time <- as.numeric(end_time - start.time)
  WBS_th@mean <- mean.from.cpt(x, WBS_th@cpt)
  
  wbssbic <- new("cpt.est")
  start.time <- Sys.time()
  z <- model.ic(sol.wbs(x))
  end_time <- Sys.time()
  wbssbic@cpt = z$cpts
  wbssbic@nocpt <- z$no.of.cpt
  wbssbic@time <- as.numeric(end_time - start.time)
  wbssbic@mean <- mean.from.cpt(x, wbssbic@cpt)
  
  if(verbose){print("wbs2")}
  wbs2 <- new("cpt.est")
  start.time <- Sys.time()
  z <-  model.sdll(sol.wbs2(x))
  end_time <- Sys.time()
  wbs2@nocpt <- z$no.of.cpt
  wbs2@cpt <- as.numeric(z$cpts)
  if(wbs2@nocpt == 0){wbs2@cpt <- 0}
  wbs2@mean <- mean.from.cpt(x, wbs2@cpt)
  wbs2@time <- as.numeric(end_time - start.time)
  
  if(verbose){print("pelt")}
  pelt <- new("cpt.est")
  start.time <- Sys.time()
  z <- cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT")
  end_time <- Sys.time()
  if (any(z@cpts[1:(length(z@cpts)-1)] == length(x))){pelt@cpt = 0
  pelt@nocpt = 0}else{
    pelt@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)])
    pelt@nocpt <- length(pelt@cpt)}
  pelt@mean <- mean.from.cpt(x, pelt@cpt)
  pelt@time <- as.numeric(end_time - start.time)
  
  if(verbose){print("not")}
  not <- new("cpt.est")
  start.time <- Sys.time()
  z <- not(x,method="not",contrast="pcwsConstMean")
  cpt.ic = features(z,q.max=qmax_NOT)
  end_time <- Sys.time()
  if(any(is.na(cpt.ic$cpt))) {not@cpt = 0
  not@nocpt <- 0}
  else {not@cpt= cpt.ic$cpt
  not@nocpt <- length(not@cpt)}
  not@mean <- mean.from.cpt(x, not@cpt)
  not@time <- as.numeric(end_time - start.time)
  
  if(verbose){print('mosum')}
  mosum <- new("cpt.est")
  start.time <- Sys.time()
  z <- multiscale.localPrune(x)$cpts
  end_time <- Sys.time()
  if(length(z) == 0) {mosum@cpt = 0
  mosum@nocpt <- 0}
  else {mosum@cpt= z
  mosum@nocpt <- length(z)}
  mosum@mean <- mean.from.cpt(x, z)
  mosum@time <- as.numeric(end_time - start.time)
  
  if(verbose){print("mscp")}
  mscp <- new("cpt.est")
  start.time <- Sys.time()
  z <- mscp(x, delta = delta)$cp
  end_time <- Sys.time()
  if(all(z == 0)) {mscp@cpt = 0
  mscp@nocpt <- 0}
  else {mscp@cpt= z
  mscp@nocpt <- length(z)}
  mscp@mean <- mean.from.cpt(x, z)
  mscp@time <- as.numeric(end_time - start.time)
  
  if(verbose){print("seedBS")}
  seedBS_th <- new("cpt.est")
  start.time <- Sys.time()
  z <- seedBS(x)
  cpt.z_1 = changepoints(z,th.const=1)
  end_time <- Sys.time()
  if(any(is.na(cpt.z_1$cpt.th[[1]]))){seedBS_th@cpt = 0}
  else{seedBS_th@cpt <- cpt.z_1$cpt.th[[1]]}
  seedBS_th@nocpt <- cpt.z_1$no.cpt.th
  seedBS_th@mean <- mean.from.cpt(x, seedBS_th@cpt)
  seedBS_th@time <- as.numeric(end_time - start.time)
  
  seedBS_ic <- new("cpt.est")
  if(any(is.na(cpt.z_1$cpt.ic[[1]]))) {seedBS_ic@cpt = 0}
  else {seedBS_ic@cpt= cpt.z_1$cpt.ic[[1]]}
  seedBS_ic@nocpt <- cpt.z_1$no.cpt.ic[[1]]
  seedBS_ic@mean <- mean.from.cpt(x, seedBS_ic@cpt) ## The computational time is the same as in seedBS_th
  
  if(verbose){print("DP_univar")}
  DP_univar <- new("cpt.est")
  start.time <- Sys.time()
  z <- DP.univar(x, gamma = 5, delta = 5)
  end_time <- Sys.time()
  if(length(z$cpt) == 0){DP_univar@cpt = 0}
  else{DP_univar@cpt <- z$cpt}
  DP_univar@nocpt <- length(z$cpt)
  DP_univar@mean <- mean.from.cpt(x, DP_univar@cpt)
  DP_univar@time <- as.numeric(end_time - start.time)
  
  list(DAIS = DAIS, ID = ID, ID_th = ID_th, ID_sic = ID_sic, WBS_th = WBS_th, 
       wbssbic = wbssbic, wbs2 = wbs2, pelt = pelt, not = not, mosum = mosum, mscp = mscp,
       seedBS_th = seedBS_th, seedBS_ic = seedBS_ic, DP_univar = DP_univar)
}

sim.study.thr <- function(signal, true.cpt=NULL,sigma, m = 100, seed = NULL, const = 1.2, points=3,
                          qmax_NOT = 25, delta = 20, verbose = FALSE) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric",
                                      mse="numeric", time= "numeric", det_delay="numeric"), 
           prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  ts <- list()
  DAIS <- new("est.eval")
  ID <- new("est.eval")
  ID_th <- new("est.eval")
  ID_sic <- new("est.eval")
  WBS_th <- new("est.eval")
  wbssbic <- new("est.eval")
  wbs2 <- new("est.eval")
  pelt <- new("est.eval")
  not <- new("est.eval")
  mosum <- new("est.eval")
  mscp <- new("est.eval")
  seedBS_th <- new("est.eval")
  seedBS_ic <- new("est.eval")
  DP_univar <- new("est.eval")
  
  no.of.cpt <- sum(abs(diff(signal)) > 0)
  n <- length(signal)
  ns <- max(c(diff(true.cpt), length(signal)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    ts[[i]] <- x
    
    est <- sim_thr(x, const = const, points=points, qmax_NOT = qmax_NOT, delta = delta, verbose = verbose)
    
    DAIS@dnc[i] <- est$DAIS@nocpt - no.of.cpt
    DAIS@mse[i] <- mean((est$DAIS@mean - signal)^2)
    DAIS@cpt[[i]] <- est$DAIS@cpt
    DAIS@diff <- abs(matrix(est$DAIS@cpt,nrow=no.of.cpt,ncol=length(est$DAIS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$DAIS@cpt),byr=F))
    DAIS@dh[i] <- max(apply(DAIS@diff,1,min),apply(DAIS@diff,2,min))/ns
    DAIS@time[i] <- est$DAIS@time
    
    ID@dnc[i] <- est$ID@nocpt - no.of.cpt
    ID@mse[i] <- mean((est$ID@mean - signal)^2)
    ID@cpt[[i]] <- est$ID@cpt
    ID@diff <- abs(matrix(est$ID@cpt,nrow=no.of.cpt,ncol=length(est$ID@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ID@cpt),byr=F))
    ID@dh[i] <- max(apply(ID@diff,1,min),apply(ID@diff,2,min))/ns
    ID@time[i] <- est$ID@time
    
    ID_th@dnc[i] <- est$ID_th@nocpt - no.of.cpt
    ID_th@mse[i] <- mean((est$ID_th@mean - signal)^2)
    ID_th@cpt[[i]] <- est$ID_th@cpt
    ID_th@diff <- abs(matrix(est$ID_th@cpt,nrow=no.of.cpt,ncol=length(est$ID_th@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ID_th@cpt),byr=F))
    ID_th@dh[i] <- max(apply(ID_th@diff,1,min),apply(ID_th@diff,2,min))/ns
    ID_th@time[i] <- est$ID_th@time
    
    ID_sic@dnc[i] <- est$ID_sic@nocpt - no.of.cpt
    ID_sic@mse[i] <- mean((est$ID_sic@mean - signal)^2)
    ID_sic@cpt[[i]] <- est$ID_sic@cpt
    ID_sic@diff <- abs(matrix(est$ID_sic@cpt,nrow=no.of.cpt,ncol=length(est$ID_sic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ID_sic@cpt),byr=F))
    ID_sic@dh[i] <- max(apply(ID_sic@diff,1,min),apply(ID_sic@diff,2,min))/ns
    ID_sic@time[i] <- est$ID_sic@time
    
    WBS_th@dnc[i] <- est$WBS_th@nocpt - no.of.cpt
    WBS_th@mse[i] <- mean((est$WBS_th@mean - signal)^2)
    WBS_th@cpt[[i]] <- est$WBS_th@cpt
    WBS_th@diff <- abs(matrix(est$WBS_th@cpt,nrow=no.of.cpt,ncol=length(est$WBS_th@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$WBS_th@cpt),byr=F))
    WBS_th@dh[i] <- max(apply(WBS_th@diff,1,min),apply(WBS_th@diff,2,min))/ns
    WBS_th@time[i] <- est$WBS_th@time
    
    wbssbic@dnc[i] <- est$wbssbic@nocpt - no.of.cpt
    wbssbic@mse[i] <- mean((est$wbssbic@mean - signal)^2)
    wbssbic@cpt[[i]] <- est$wbssbic@cpt
    wbssbic@diff <- abs(matrix(est$wbssbic@cpt,nrow=no.of.cpt,ncol=length(est$wbssbic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbssbic@cpt),byr=F))
    wbssbic@dh[i] <- max(apply(wbssbic@diff,1,min),apply(wbssbic@diff,2,min))/ns
    wbssbic@time[i] <- est$wbssbic@time
    
    wbs2@dnc[i] <- est$wbs2@nocpt - no.of.cpt
    wbs2@mse[i] <- mean((est$wbs2@mean - signal)^2)
    wbs2@cpt[[i]] <- est$wbs2@cpt
    wbs2@diff <- abs(matrix(est$wbs2@cpt,nrow=no.of.cpt,ncol=length(est$wbs2@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbs2@cpt),byr=F))
    wbs2@dh[i] <- max(apply(wbs2@diff,1,min),apply(wbs2@diff,2,min))/ns
    wbs2@time[i] <- est$wbs2@time
    
    pelt@dnc[i] <- est$pelt@nocpt - no.of.cpt
    pelt@mse[i] <- mean((est$pelt@mean - signal)^2)
    pelt@cpt[[i]] <- est$pelt@cpt
    pelt@diff <- abs(matrix(est$pelt@cpt,nrow=no.of.cpt,ncol=length(est$pelt@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$pelt@cpt),byr=F))
    pelt@dh[i] <- max(apply(pelt@diff,1,min),apply(pelt@diff,2,min))/ns
    pelt@time[i] <- est$pelt@time
    
    not@dnc[i] <- est$not@nocpt - no.of.cpt
    not@mse[i] <- mean((est$not@mean - signal)^2)
    not@cpt[[i]] <- est$not@cpt
    not@diff <- abs(matrix(est$not@cpt,nrow=no.of.cpt,ncol=length(est$not@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$not@cpt),byr=F))
    not@dh[i] <- max(apply(not@diff,1,min),apply(not@diff,2,min))/ns
    not@time[i] <- est$not@time
    
    mosum@dnc[i] <- est$mosum@nocpt - no.of.cpt
    mosum@mse[i] <- mean((est$mosum@mean - signal)^2)
    mosum@cpt[[i]] <- est$mosum@cpt
    mosum@diff <- abs(matrix(est$mosum@cpt,nrow=no.of.cpt,ncol=length(est$mosum@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$mosum@cpt),byr=F))
    mosum@dh[i] <- max(apply(mosum@diff,1,min),apply(mosum@diff,2,min))/ns
    mosum@time[i] <- est$mosum@time
    
    mscp@dnc[i] <- est$mscp@nocpt - no.of.cpt
    mscp@mse[i] <- mean((est$mscp@mean - signal)^2)
    mscp@cpt[[i]] <- est$mscp@cpt
    mscp@diff <- abs(matrix(est$mscp@cpt,nrow=no.of.cpt,ncol=length(est$mscp@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$mscp@cpt),byr=F))
    mscp@dh[i] <- max(apply(mscp@diff,1,min),apply(mscp@diff,2,min))/ns
    mscp@time[i] <- est$mscp@time
    
    seedBS_th@dnc[i] <- est$seedBS_th@nocpt - no.of.cpt
    seedBS_th@mse[i] <- mean((est$seedBS_th@mean - signal)^2)
    seedBS_th@cpt[[i]] <- est$seedBS_th@cpt
    seedBS_th@diff <- abs(matrix(est$seedBS_th@cpt,nrow=no.of.cpt,ncol=length(est$seedBS_th@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$seedBS_th@cpt),byr=F))
    seedBS_th@dh[i] <- max(apply(seedBS_th@diff,1,min),apply(seedBS_th@diff,2,min))/ns
    seedBS_th@time[i] <- est$seedBS_th@time
    
    seedBS_ic@dnc[i] <- est$seedBS_ic@nocpt - no.of.cpt
    seedBS_ic@mse[i] <- mean((est$seedBS_ic@mean - signal)^2)
    seedBS_ic@cpt[[i]] <- est$seedBS_ic@cpt
    seedBS_ic@diff <- abs(matrix(est$seedBS_ic@cpt,nrow=no.of.cpt,ncol=length(est$seedBS_ic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$seedBS_ic@cpt),byr=F))
    seedBS_ic@dh[i] <- max(apply(seedBS_ic@diff,1,min),apply(seedBS_ic@diff,2,min))/ns
    seedBS_ic@time[i] <- est$seedBS_th@time
    
    DP_univar@dnc[i] <- est$DP_univar@nocpt - no.of.cpt
    DP_univar@mse[i] <- mean((est$DP_univar@mean - signal)^2)
    DP_univar@cpt[[i]] <- est$DP_univar@cpt
    DP_univar@diff <- abs(matrix(est$DP_univar@cpt,nrow=no.of.cpt,ncol=length(est$DP_univar@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$DP_univar@cpt),byr=F))
    DP_univar@dh[i] <- max(apply(DP_univar@diff,1,min),apply(DP_univar@diff,2,min))/ns
    DP_univar@time[i] <- est$DP_univar@time
    
    gc()
  }
  
  list(ts = ts, DAIS = DAIS, ID = ID, ID_th = ID_th, ID_sic = ID_sic, WBS_th = WBS_th, 
       wbssbic = wbssbic, wbs2 = wbs2, pelt = pelt, not = not, mosum = mosum, mscp = mscp,
       seedBS_th = seedBS_th, seedBS_ic = seedBS_ic, DP_univar = DP_univar)
}

#S1
small_dist <- c(rep(0, 485), rep(1, 30), rep(0, 485))
SIMR13.small <- sim.study.thr(small_dist, true.cpt = c(485, 515), sigma = 1, delta = 15)

#S2
stairs10 = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),
             rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10))
SIMR5.small = sim.study.thr(stairs10, true.cpt = seq(11, 141, 10), sigma = 0.3, delta = 5)

#S3
mix2 = c(rep(7,11),rep(-7,10),rep(6,20),rep(-6,20),rep(5,30),rep(-5,30),rep(4,40),rep(-4,40),
         rep(3,50),rep(-3,50))
SIMR3.small2 = sim.study.thr(mix2,true.cpt=c(11, 21, 41, 61, 91, 121, 161, 201, 251), sigma = 4, delta = 5)

#S4
many_cpts_mix <- c(rep(0,5), rep(5,7), rep(0,5), rep(6,8), rep(0,6), rep(4,7), 
                   rep(0,6), rep(5,6), rep(0,6), rep(6,5), rep(0,6), rep(4,8))
many_cpts_mix_sim <- sim.study.thr(many_cpts_mix, true.cpt = which(abs(diff(many_cpts_long)) > 0), sigma = 1, delta = 3)

#S5
many_cpts <- rep(c(rep(0, 7), rep(4, 7)), 50)
SIMR.many <- sim.study.thr(many_cpts, true.cpt = which(diff(many_cpts) != 0), sigma = 1, qmax_NOT = 100, delta = 4)

#S6
many_cpts_long <- rep(c(rep(0, 5), rep(5, 5)), 60) #S6
many_cpts_long_sim <- sim.study.thr(many_cpts_long, true.cpt = which(abs(diff(many_cpts_long)) > 0), 
                                    sigma = 0.9, qmax_NOT = 120, delta = 3)

#S11
justnoise = rep(0,6000) 
NC.small = sim.study.thr(justnoise, sigma = 1, true.cpt = c(0))

#S12
long_signal <- c(rep(0,5500), rep(1.5,5500)) 
SIMR.large <- sim.study.thr(long_signal, true.cpt=5500, sigma=1)

# S13
small_dist3 <- c(rep(0,100), rep(1.5,30), rep(0, 355), rep(1, 30), rep(0, 355), 
                 rep(1.5,30), rep(0,100)) 
SIMR15.small <- sim.study.thr(small_dist3, true.cpt = c(100, 130, 485, 515, 870, 900), sigma = 1, delta = 15)

#S14
teeth10 = c(rep(0,11),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,20),
            rep(0,20),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,19)) 
SIMR4.small = sim.study.thr(teeth10, true.cpt = seq(11, 251, 20), sigma = 0.4, delta = 5)



### For piecewise-linear signals
library(earth)
library(stats)
library(IDetect)
library(cpop)
library(not)
#devtools::install_github("hadley/l1tf")
library("l1tf")
library(trendsegmentR)

options(expressions = 500000)

get.signal <- function(model){
  
  if(length(model$cpt) != 0){
    
    signal <- rep(0, model$n)
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    
    slope <- model$start[2]
    signal[segments[1,1]:segments[1,2]] <- model$start[1] + segments[1,1]:segments[1,2] * model$start[2]
    
    for(j in 2:nrow(segments)) {
      
      slope <- slope +  model$jump.size[j-1]
      
      for(k in segments[j,1]:segments[j,2]) signal[k] <- signal[k-1] + slope
    }
    
  }
  return(signal)
}


sim.model <- function(model, sigma=1){
  get.signal(model) + sigma * rnorm(model$n)
} 


## The following function is to get the estimated signal in the case of piecewise-linearity.
fit_lin_cont <- function(x, cpt) {
  lx = length(x)
  if (missing(cpt)) 
    cpt <- linear_find_changepoint_ic_fixedw3(x)
  if (!is.null(cpt)) 
    if (any(is.na(cpt))) 
      cpt <- cpt[!is.na(cpt)]
  cpt <- sort(unique(c(cpt, 0, lx)))
  fit <- rep(0, lx)
  cpt <- setdiff(cpt, c(0, lx))
  X <- bs(1:lx, knots = cpt, degree = 1, intercept = TRUE)
  fit <- lm.fit(X, x)$fitted.values
  return(fit)
}

## run trend-filtering input data
tf.run = function(z, sig = 1, lambdas = exp(seq(log(10), log(1000), length = 250))) {
  M = length(lambdas)
  n = length(z)
  SIC = rep(0, M)
  for (i in 1:M) {
    lambda <- lambdas[i]
    ans.tf <- l1tf(z, lambda) 
    SIC[i] = sum((z - ans.tf)^2)
    dans <- round(diff(ans.tf), digits = 5)
    CP <- c()
    for (ii in 1:(length(dans) - 1)) {
      if (dans[ii] != dans[ii + 1]) {
        CP <- c(CP, ii)
      }
    }
    SIC[i] = SIC[i]/sig^2 + log(n) * length(CP)
  }
  k = which.min(SIC)
  ans.tf <- l1tf(z, lambdas[k]) 
  dans <- round(diff(ans.tf), digits = 5)
  CP <- c()
  for (ii in 1:(length(dans) - 1)) {
    if (dans[ii] != dans[ii + 1]) {
      CP <- c(CP, ii)
    }
  }
  
  return(list(cpt = CP, fit = ans.tf, lam = lambdas[k]))
  
}

## Simulation functions

linear.rev.sim <- function(x, q_max_NOT = 25, FKS_knot = 10) {
  
  setClass("cpt.est", representation(cpt = "numeric", nocpt = "numeric", fit = "numeric", time = "numeric", dh = "numeric"), 
           prototype(cpt = numeric(0), nocpt = 0, fit = numeric(0), time = numeric(0), dh = numeric(0)))
  
  print("DAIS")
  DAIS <- new("cpt.est")
  z = DAIS(x, contrast = "slope")
  DAIS@cpt <- as.numeric(z)
  DAIS@nocpt <- length(z)
  DAIS@fit <- fit_lin_cont(x, DAIS@cpt)
  DAIS@time <- system.time(DAIS(x, contrast = "slope"))[[3]]
  
  print("ID_th")
  ID_th <- new("cpt.est")
  z <- cplm_th(x, thr_const = 1.1)
  if(length(z) == 0){ID_th@cpt = 0
  ID_th@nocpt = 0}
  else{ID_th@cpt <- z
  ID_th@nocpt <- length(z)}
  ID_th@fit <- fit_lin_cont(x, ID_th@cpt)
  ID_th@time <- system.time(cplm_th(x, thr_const = 1.1))[[3]]
  
  print("ID_sic")
  ID_sic <- new("cpt.est")
  z <- cplm_ic(x)
  if(any(is.na(z$cpt_ic$sic_pen))){ID_sic@cpt = 0
  ID_sic@nocpt = 0}
  else{ID_sic@cpt <- z$cpt_ic$sic_pen
  ID_sic@nocpt <- length(z$cpt_ic$sic_pen)}
  ID_sic@fit <- fit_lin_cont(x, ID_sic@cpt)
  ID_sic@time <- system.time(cplm_ic(x))[[3]]
  
  print("CPOP")
  cpop <- new("cpt.est")
  z <- cpop(x, beta = 2*log(length(x)), sd = rep(mad(diff(diff(x)))/sqrt(6), length(x)))
  if (length(z@changepoints) == 2) {
    cpop@cpt = 0
    cpop@nocpt = 0
  } else {
    cpop@cpt = z@changepoints[-c(1, length(z@changepoints))]
    cpop@nocpt <- length(cpop@cpt)
  }
  cpop@fit <- fit_lin_cont(x, cpop@cpt)
  cpop@time <- system.time(cpop(x, beta = 2*log(length(x)), sd = rep(mad(diff(diff(x)))/sqrt(6), length(x))))[[3]]
  
  print("not")
  not <- new("cpt.est")
  z <- not(x, method = "not", contrast = "pcwsLinContMean")
  cpt.ic = features(z, q.max = q_max_NOT)
  if (any(is.na(cpt.ic$cpt))) {
    not@cpt = 0
    not@nocpt = 0
  } else {
    not@cpt = cpt.ic$cpt
    not@nocpt <- length(not@cpt)
  }
  not@fit <- fit_lin_cont(x, not@cpt)
  not@time <- system.time(features(not(x, method = "not", contrast = "pcwsLinContMean"), q.max = q_max_NOT))[[3]]
  
  print("MARS")
  MARS <- new("cpt.est")
  z1 <- earth(1:length(x),y=x)
  z <- sort(unique(z1$cuts[z1$selected.terms,]))[-1]
  if (length(z) == 0) {
    MARS@cpt = 0
    MARS@nocpt = 0
  }
  else {
    MARS@cpt = z
    MARS@nocpt = length(z)
  }
  MARS@fit <- fit_lin_cont(x, MARS@cpt)
  MARS@time <- system.time(sort(unique(earth(1:length(x),y=x)$cuts[earth(1:length(x),y=x)$selected.terms,]))[-1])[[3]]
  
  print("t1f")
  t1f <- new("cpt.est")
  z <- tf.run(x/(mad(diff(diff(x)))/sqrt(6)))
  if (length(z$cpt) == 0) {
    t1f@cpt = 0
    t1f@nocpt = 0
  } else {
    t1f@cpt = z$cpt
    t1f@nocpt <- length(z$cpt)
  }
  t1f@fit <- fit_lin_cont(x, t1f@cpt)
  t1f@time <- system.time(tf.run(x/(mad(diff(diff(x)))/sqrt(6))))[[3]]
  
  print("trendsegment")
  trendsegment <- new("cpt.est")
  z = trendsegment(x, continuous = TRUE, indep = TRUE)$cpt
  if (is.null(z)) {
    trendsegment@cpt = 0
    trendsegment@nocpt = 0
  } else {
    trendsegment@cpt <- as.numeric(z)
    trendsegment@nocpt <- length(z)
  }
  trendsegment@fit <- fit_lin_cont(x, trendsegment@cpt)
  trendsegment@time <- system.time(trendsegment(x, continuous = TRUE, indep = TRUE))[[3]]
  
  list(DAIS = DAIS, ID_th = ID_th, ID_sic = ID_sic, cpop = cpop, not = not, 
       MARS = MARS, t1f = t1f, trendsegment = trendsegment)
}

linear_rev.sim.study <- function(model, sigma, m = 100, seed = NULL, gen_qmax = 25) {
  
  setClass("est.eval", representation(avg.signal = "numeric", fit = "list", cpt = "list", diff = "matrix", dh = "numeric", 
                                      cptall = "numeric", dnc = "numeric", mse = "numeric", time = "numeric"), 
           prototype(dnc = numeric(m), mse = numeric(m), time = numeric(m), dh = numeric(m)))
  
  DAIS <- new("est.eval")
  IsolateDetect <- new("est.eval")
  ID_th <- new("est.eval")
  ID_sic <- new("est.eval")
  not <- new("est.eval")
  cpop <- new("est.eval")
  t1f <- new("est.eval")
  MARS <- new("est.eval")
  trendsegment <- new("est.eval")
  
  signal = get.signal(model)
  
  no.of.cpt <- length(model$cpt)
  n <- length(signal)
  ns <- max(diff(c(0, model$cpt, n)))
  
  if (!is.null(seed)) 
    set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    
    est <- linear.rev.sim(x, q_max_NOT = gen_qmax, FKS_knot = no.of.cpt + 2)
    
    DAIS@dnc[i] <- est$DAIS@nocpt - no.of.cpt
    DAIS@cpt[[i]] <- est$DAIS@cpt
    DAIS@mse[i] <- mean((signal - est$DAIS@fit)^2)
    DAIS@diff <- abs(matrix(est$DAIS@cpt, nrow = no.of.cpt, ncol = length(est$DAIS@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                           ncol = length(est$DAIS@cpt), byr = F))
    DAIS@dh[i] <- max(apply(DAIS@diff, 1, min), apply(DAIS@diff, 2, min))/ns
    DAIS@time[i] <- est$DAIS@time
    
    ID_th@dnc[i] <- est$ID_th@nocpt - no.of.cpt
    ID_th@cpt[[i]] <- est$ID_th@cpt
    ID_th@mse[i] <- mean((signal - est$ID_th@fit)^2)
    ID_th@diff <- abs(matrix(est$ID_th@cpt, nrow = no.of.cpt, ncol = length(est$ID_th@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                              ncol = length(est$ID_th@cpt), byr = F))
    ID_th@dh[i] <- max(apply(ID_th@diff, 1, min), apply(ID_th@diff, 2, min))/ns
    ID_th@time[i] <- est$ID_th@time
    
    ID_sic@dnc[i] <- est$ID_sic@nocpt - no.of.cpt
    ID_sic@cpt[[i]] <- est$ID_sic@cpt
    ID_sic@mse[i] <- mean((signal - est$ID_sic@fit)^2)
    ID_sic@diff <- abs(matrix(est$ID_sic@cpt, nrow = no.of.cpt, ncol = length(est$ID_sic@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                                 ncol = length(est$ID_sic@cpt), byr = F))
    ID_sic@dh[i] <- max(apply(ID_sic@diff, 1, min), apply(ID_sic@diff, 2, min))/ns
    ID_sic@time[i] <- est$ID_sic@time
    
    cpop@dnc[i] <- est$cpop@nocpt - no.of.cpt
    cpop@cpt[[i]] <- est$cpop@cpt
    cpop@mse[i] <- mean((signal - est$cpop@fit)^2)
    cpop@diff <- abs(matrix(est$cpop@cpt, nrow = no.of.cpt, ncol = length(est$cpop@cpt), byr = T) - matrix(model$cpt,
                                                                                                           nrow = no.of.cpt, ncol = length(est$cpop@cpt), byr = F))
    cpop@dh[i] <- max(apply(cpop@diff, 1, min), apply(cpop@diff, 2, min))/ns
    cpop@time[i] <- est$cpop@time
    
    not@dnc[i] <- est$not@nocpt - no.of.cpt
    not@cpt[[i]] <- est$not@cpt
    not@mse[i] <- mean((signal - est$not@fit)^2)
    not@diff <- abs(matrix(est$not@cpt, nrow = no.of.cpt, ncol = length(est$not@cpt), byr = T) - matrix(model$cpt, 
                                                                                                              nrow = no.of.cpt, ncol = length(est$not@cpt), byr = F))
    not@dh[i] <- max(apply(not@diff, 1, min), apply(not@diff, 2, min))/ns
    not@time[i] <- est$not@time
    
    MARS@dnc[i] <- est$MARS@nocpt - no.of.cpt
    MARS@cpt[[i]] <- est$MARS@cpt
    MARS@mse[i] <- mean((signal - est$MARS@fit)^2)
    MARS@diff <- abs(matrix(est$MARS@cpt, nrow = no.of.cpt, ncol = length(est$MARS@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                           ncol = length(est$MARS@cpt), byr = F))
    MARS@dh[i] <- max(apply(MARS@diff, 1, min), apply(MARS@diff, 2, min))/ns
    MARS@time[i] <- est$MARS@time
    
    t1f@dnc[i] <- est$t1f@nocpt - no.of.cpt
    t1f@cpt[[i]] <- est$t1f@cpt
    t1f@mse[i] <- mean((signal - est$t1f@fit)^2)
    t1f@diff <- abs(matrix(est$t1f@cpt, nrow = no.of.cpt, ncol = length(est$t1f@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                        ncol = length(est$t1f@cpt), byr = F))
    t1f@dh[i] <- max(apply(t1f@diff, 1, min), apply(t1f@diff, 2, min))/ns
    t1f@time[i] <- est$t1f@time
    
    trendsegment@dnc[i] <- est$trendsegment@nocpt - no.of.cpt
    trendsegment@cpt[[i]] <- est$trendsegment@cpt
    trendsegment@mse[i] <- mean((signal - est$trendsegment@fit)^2)
    trendsegment@diff <- abs(matrix(est$trendsegment@cpt, nrow = no.of.cpt, ncol = length(est$trendsegment@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                                                   ncol = length(est$trendsegment@cpt), byr = F))
    trendsegment@dh[i] <- max(apply(trendsegment@diff, 1, min), apply(trendsegment@diff, 2, min))/ns
    trendsegment@time[i] <- est$trendsegment@time
    
    gc()
  }
  list(DAIS = DAIS, ID_th = ID_th, ID_sic = ID_sic, cpop = cpop, not = not, 
       MARS = MARS, t1f = t1f, trendsegment = trendsegment)
}

## S8
model.wave1 <- list(name = "wave1", cpt.type = "pcwsLinContMean", cpt = c(256, 512, 768, 1024, 1152, 1280, 1344), jump.size = (-1)^(1:7) * 
                      (1:7)/64, n = 1408, start = c(1, 1/256))
lin.SIMR1 = linear_rev.sim.study(model.wave1, 1, seed = 16)

## S9
model.wave2 <- list(name = "wave2", cpt.type = "pcwsLinContMean", cpt = (1:99) * 15, jump.size = (-1)^{1:100}, 
                    n = 15 * 100, start = c(-1/2, 1/40))
lin.SIMR2 = linear_rev.sim.study(model.wave2, 1, seed = 16, gen_qmax = 200)

## S10
model.wave3 <- list(name = "wave3", cpt.type = "pcwsLinContMean", cpt = (1:119) * 7, jump.size = (-1)^{1:120}, 
                    n = 840, start = c(-1/2, 1/32))
lin.SIMR3 = linear_rev.sim.study(model.wave3, m = 100, sigma = 0.3, seed = 16)

## S15
model.justnoise.wave <- list(name = "justnoise_wave", cpt.type = "pcwsLinContMean",
                             cpt = numeric(0), n = 1000, start = c(0,1))
lin.SIMR.justnoise = linear_rev.sim.study(model.justnoise.wave, m = 100, sigma = 1)

## S16
model.wave4 <- list(name = "wave4", cpt.type = "pcwsLinContMean", cpt = (1:9) * 20,
                    jump.size = c(1/6, 1/2,-3/4,-1/3, -2/3,1,1/4,3/4,-5/4), n = 200, start = c(1, 1/32))
lin.SIMR4 = linear_rev.sim.study(model.wave4, m = 100, sigma = 0.3, seed = 16)

## S15
model.wave5 <- list(name = "wave5", cpt.type = "pcwsLinContMean", cpt = (1:49) * 7,
                    jump.size = c(rep(c(-2.5,2.5),50)), n = 350, start = c(0, 1))
lin.SIMR5 = linear_rev.sim.study(model.wave5, m = 100, sigma = 1)

