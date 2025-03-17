### For piecewise-constant signals

library(devtools)
library(IDetect)
library(wbs)
library(breakfast)
library(changepoint)
library(changepoints)
library(mscp)
library(mosum)
library(stepR)

library(MsFPOP) # The necessary code for MsFPOP can be downloaded from https://github.com/aLiehrmann/MsFPOP
# The necessary code for seedBS can be downloaded from https://github.com/kovacssolt/SeedBinSeg

# The following code requires first running the R file Algorithms.R

options(expressions = 500000)

mean.from.cpt <- function(x, cpt) {
  n <- length(x)
  len.cpt <- length(cpt)
  if (len.cpt) cpt <- sort(cpt)
  beg <- endd <- rep(0, len.cpt+1)
  beg[1] <- 1
  endd[len.cpt+1] <- n
  if (len.cpt) {
    beg[2:(len.cpt+1)] <- cpt+1
    endd[1:len.cpt] <- cpt
  }
  means <- rep(0, len.cpt+1)
  for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])
  rep(means, endd-beg+1)
}

sim_thr <- function(x, const = 1.7/sqrt(2), points, q_max, delta, verbose){
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric",det_time="numeric", time="numeric", mean="numeric"), 
           prototype(cpt=numeric(0), nocpt=0, det_time=numeric(0), time=numeric(0), mean=numeric(0)))
  
  if(verbose){print("DAIS")}
  DAIS <- new("cpt.est")
  start_time <- Sys.time()
  z <- DAIS(x, contrast = "mean", thr_const = const, points=points)
  end_time <- Sys.time()
  if(length(z) == 0){DAIS@cpt = 0
  DAIS@nocpt = 0}
  else{DAIS@cpt <- z
  DAIS@nocpt <- length(z)}
  DAIS@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  DAIS@mean <- mean.from.cpt(x, DAIS@cpt)
  
  if(verbose){print("ID_th")}
  ID_th <- new("cpt.est")
  ID_th <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.thresh(sol.idetect_seq(x), th.const = 1.15)
  end_time <- Sys.time()
  ID_th@cpt = z$cpts
  ID_th@nocpt <- z$no.of.cpt
  ID_th@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  ID_th@mean <- mean.from.cpt(x, ID_th@cpt)
  
  if(verbose){print("ID_sic")}
  ID_sic <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.ic(sol.idetect(x), q.max = q_max)
  end_time <- Sys.time()
  ID_sic@cpt = z$cpts
  ID_sic@nocpt <- z$no.of.cpt
  ID_sic@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  ID_sic@mean <- mean.from.cpt(x, ID_sic@cpt)
  
  if(verbose){print("wbs")}
  WBS_th <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.thresh(sol.wbs(x), th.const = 1.15)
  end_time <- Sys.time()
  WBS_th@cpt = z$cpts
  WBS_th@nocpt <- z$no.of.cpt
  WBS_th@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  WBS_th@mean <- mean.from.cpt(x, WBS_th@cpt)
  
  wbssbic <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.ic(sol.wbs(x), q.max = q_max)
  end_time <- Sys.time()
  wbssbic@cpt = z$cpts
  wbssbic@nocpt <- z$no.of.cpt
  wbssbic@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  wbssbic@mean <- mean.from.cpt(x, wbssbic@cpt)
  
  if(verbose){print("wbs2")}
  wbs2 <- new("cpt.est")
  start_time <- Sys.time()
  z <-  model.sdll(sol.wbs2(x))
  end_time <- Sys.time()
  wbs2@nocpt <- z$no.of.cpt
  wbs2@cpt <- as.numeric(z$cpts)
  if(wbs2@nocpt == 0){wbs2@cpt <- 0}
  wbs2@mean <- mean.from.cpt(x, wbs2@cpt)
  wbs2@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("pelt")}
  pelt <- new("cpt.est")
  start_time <- Sys.time()
  z <- cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT")
  end_time <- Sys.time()
  if (any(z@cpts[1:(length(z@cpts)-1)] == length(x))){pelt@cpt = 0
  pelt@nocpt = 0}else{
    pelt@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)])
    pelt@nocpt <- length(pelt@cpt)}
  pelt@mean <- mean.from.cpt(x, pelt@cpt)
  pelt@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("not")}
  not <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.ic(sol.not(x), q.max = q_max)
  end_time <- Sys.time()
  not@nocpt <- z$no.of.cpt
  not@cpt <- as.numeric(z$cpts)
  not@mean <- mean.from.cpt(x, not@cpt)
  not@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print('mosum')}
  mosum <- new("cpt.est")
  start_time <- Sys.time()
  z <- multiscale.localPrune(x, G = bandwidths.default(length(x), d.min = 3, G.min = 10))$cpts
  end_time <- Sys.time()
  if(length(z) == 0) {mosum@cpt = 0
  mosum@nocpt <- 0}
  else {mosum@cpt= z
  mosum@nocpt <- length(z)}
  mosum@mean <- mean.from.cpt(x, z)
  mosum@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("mscp")}
  mscp <- new("cpt.est")
  start_time <- Sys.time()
  z <- mscp(x, delta = delta)$cp
  end_time <- Sys.time()
  if(all(z == 0)) {mscp@cpt = 0
  mscp@nocpt <- 0}
  else {mscp@cpt= z
  mscp@nocpt <- length(z)}
  mscp@mean <- mean.from.cpt(x, z)
  mscp@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("seedBS")}
  seedBS_th <- new("cpt.est")
  start_time <- Sys.time()
  z <- seedBS(x)
  cpt.z_1 = changepoints(z,th.const=1.15)
  end_time <- Sys.time()
  if(any(is.na(cpt.z_1$cpt.th[[1]]))){seedBS_th@cpt = 0}
  else{seedBS_th@cpt <- cpt.z_1$cpt.th[[1]]}
  seedBS_th@nocpt <- cpt.z_1$no.cpt.th
  seedBS_th@mean <- mean.from.cpt(x, seedBS_th@cpt)
  seedBS_th@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  seedBS_ic <- new("cpt.est")
  if(any(is.na(cpt.z_1$cpt.ic[[1]]))) {seedBS_ic@cpt = 0}
  else {seedBS_ic@cpt= cpt.z_1$cpt.ic[[1]]}
  seedBS_ic@nocpt <- cpt.z_1$no.cpt.ic[[1]]
  seedBS_ic@mean <- mean.from.cpt(x, seedBS_ic@cpt) 
  seedBS_ic@time <- seedBS_th@time ## The computational time is the same as in seedBS_th
  
  if(verbose){print("DP_univar")}
  DP_univar <- new("cpt.est")
  start_time <- Sys.time()
  z <- DP.univar(x, gamma = 5, delta = 5)
  end_time <- Sys.time()
  if(length(z$cpt) == 0){DP_univar@cpt = 0}
  else{DP_univar@cpt <- z$cpt}
  DP_univar@nocpt <- length(z$cpt)
  DP_univar@mean <- mean.from.cpt(x, DP_univar@cpt)
  DP_univar@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("SMUCE")}
  SMUCE <- new("cpt.est")
  start_time <- Sys.time()
  z <- stepFit(y = x, x = 1:length(x), alpha = 0.5, jumpint = FALSE, confband = FALSE)$rightIndex
  end_time <- Sys.time()
  if(length(z) == 1){SMUCE@cpt = 0
  SMUCE@nocpt = 0}
  else{SMUCE@cpt <- z[-length(z)]
  SMUCE@nocpt <- length(z)-1}
  SMUCE@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  SMUCE@mean <- mean.from.cpt(x, SMUCE@cpt)
  
  if(verbose){print("MsFPOP")}
  MsFPOP <- new("cpt.est")
  start_time <- Sys.time()
  z <- MsFPOP(y=x, beta=2.25, alpha=9+2.25*log(length(x)))$changepoint
  end_time <- Sys.time()
  if(length(z) == 1){MsFPOP@cpt = 0
  MsFPOP@nocpt = 0}
  else{MsFPOP@cpt <- z[-length(z)]
  MsFPOP@nocpt <- length(z)-1}
  MsFPOP@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  MsFPOP@mean <- mean.from.cpt(x, MsFPOP@cpt)
  
  list(DAIS = DAIS, ID_th = ID_th, ID_sic = ID_sic, WBS_th = WBS_th, 
       wbssbic = wbssbic, wbs2 = wbs2, pelt = pelt, not = not, mosum = mosum, mscp = mscp,
       seedBS_th = seedBS_th, seedBS_ic = seedBS_ic, DP_univar = DP_univar, SMUCE = SMUCE,
       MsFPOP = MsFPOP)
}

sim.study.thr <- function(signal, true.cpt=NULL, sigma, m = 100, seed = NULL, const = 1.7/sqrt(2), points=3,
                          delta = NULL, verbose = FALSE) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric",
                                      mse="numeric", time= "numeric", det_delay="numeric"), 
           prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  ts <- list()
  
  no.of.cpt <- sum(abs(diff(signal)) > 0)
  n <- length(signal)
  ns <- max(c(diff(true.cpt), length(signal)))
  d_T <- min(diff(c(0, true.cpt, n)))
  if(d_T<1){d_T = n}
  
  q_max = max(min(100, n/log(n)), ceiling(n/d_T))
  if(is.null(delta)){delta = min(ceiling(d_T/2), 20)}
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    ts[[i]] <- x
    
    ## perform the MC simulations before the time is recorded for SMUCE
    if(i == 1){stepFit(y = x, x = 1:length(x), alpha = 0.5, jumpint = FALSE, confband = FALSE)$rightIndex}
    
    est <- sim_thr(x, const = const, points=points, q_max = q_max, delta = delta, verbose = verbose)
    
    if(i == 1){
      for(j in names(est)){
        eval(parse(text=paste(j, " <- new('est.eval')",sep="")))
      }
    }
    
    for(j in names(est)){
      eval(parse(text=paste(j, "@dnc[", i, "] <- est$", j, "@nocpt - no.of.cpt",sep="")))
      eval(parse(text=paste(j, "@mse[", i, "] <- mean((est$", j, "@mean - signal)^2)",sep="")))
      eval(parse(text=paste(j, "@cpt[[", i, "]] <- est$", j, "@cpt",sep="")))
      eval(parse(text=paste(j, "@diff <- abs(matrix(est$", j, "@cpt,nrow=no.of.cpt,ncol=length(est$", j, "@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$", j, "@cpt),byr=F))",sep="")))
      eval(parse(text=paste(j, "@dh[i] <- max(apply(", j, "@diff,1,min),apply(", j, "@diff,2,min))/ns",sep="")))
      eval(parse(text=paste(j, "@time[i] <- est$", j, "@time",sep="")))
    }
    
    gc()
  }
  
  result <- list(ts = ts)
  for(i in names(est)){
    eval(parse(text=paste("result <- c(result, ", i, " = ", i, ")",sep="")))
  }
  return(result)
}

seed.temp = 15
#S1
small_dist <- c(rep(0, 485), rep(1, 30), rep(0, 485))
small_dist_sim <- sim.study.thr(small_dist, true.cpt = c(485, 515), sigma = 1, delta = 15, seed = seed.temp)

#S2
small_dist2 <- c(rep(0, 30), rep(2.3, 5), rep(8, 100))
small_dist2_sim <- sim.study.thr(small_dist, true.cpt = c(30, 35), sigma = 1, 
                                 delta = 15, seed = seed.temp)

#S3
stairs10 = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),
             rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10))
stairs10_sim = sim.study.thr(stairs10, true.cpt = seq(11, 141, 10), sigma = 0.3, 
                             delta = 5, seed = seed.temp)

#S4
mix2 = c(rep(7,11),rep(-7,10),rep(6,20),rep(-6,20),rep(5,30),rep(-5,30),rep(4,40),rep(-4,40),
         rep(3,50),rep(-3,50))
mix2_sim = sim.study.thr(mix2,true.cpt=c(11, 21, 41, 61, 91, 121, 161, 201, 251), 
                         sigma = 4, delta = 5, seed = seed.temp)

#S5
many_cpts_mix <- c(rep(0,5), rep(5,7), rep(0,5), rep(6,8), rep(0,6), rep(4,7), 
                   rep(0,6), rep(5,6), rep(0,6), rep(6,5), rep(0,6), rep(4,8))
many_cpts_mix_sim <- sim.study.thr(many_cpts_mix, true.cpt = which(abs(diff(many_cpts_long)) > 0), 
                                   sigma = 1, delta = 3, seed = seed.temp)

#S6
many_cpts <- rep(c(rep(0, 7), rep(4, 7)), 50)
many_cpts_sim <- sim.study.thr(many_cpts, true.cpt = which(diff(many_cpts) != 0), sigma = 1, 
                               delta = 4, seed = seed.temp)

#S7
many_cpts_long <- rep(c(rep(0, 5), rep(5, 5)), 60) #S6
many_cpts_long_sim <- sim.study.thr(many_cpts_long, true.cpt = which(abs(diff(many_cpts_long)) > 0), 
                                    sigma = 0.9, delta = 3, seed = seed.temp)

#S12
justnoise = rep(0,6000) 
justnoise_sim = sim.study.thr(justnoise, sigma = 1, true.cpt = c(0))

#S13
long_signal <- c(rep(0,5500), rep(1.5,5500)) 
long_signal_sim <- sim.study.thr(long_signal, true.cpt=5500, sigma=1, seed = seed.temp)

# S14
small_dist3 <- c(rep(0,100), rep(1.5,30), rep(0, 355), rep(1, 30), rep(0, 355), 
                 rep(1.5,30), rep(0,100)) 
small_dist3_sim <- sim.study.thr(small_dist3, true.cpt = c(100, 130, 485, 515, 870, 900), sigma = 1, 
                                 delta = 15, seed = seed.temp)

#S15
teeth10 = c(rep(0,11),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,20),
            rep(0,20),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,19)) 
teeth10_sim = sim.study.thr(teeth10, true.cpt = seq(11, 251, 20), sigma = 0.4, delta = 5, seed = seed.temp)



### For piecewise-linear signals
library(earth)
library(stats)
library(breakfast)
library(IDetect)
library(cpop)
library(not)
#devtools::install_github("hadley/l1tf")
library("l1tf")
library(trendsegmentR)
library(changepoints)

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

linear.rev.sim <- function(x, const = 2.1/sqrt(2), q_max, FKS_knot = 10, verbose = F, gamma_set) {
  
  setClass("cpt.est", representation(cpt = "numeric", nocpt = "numeric", fit = "numeric", time = "numeric", dh = "numeric"), 
           prototype(cpt = numeric(0), nocpt = 0, fit = numeric(0), time = numeric(0), dh = numeric(0)))
  
  if(verbose){print("DAIS")}
  DAIS <- new("cpt.est")
  start_time <- Sys.time()
  z = DAIS(x, contrast = "slope", thr_const = const)
  end_time <- Sys.time()
  DAIS@cpt <- as.numeric(z)
  DAIS@nocpt <- length(z)
  DAIS@fit <- fit_lin_cont(x, DAIS@cpt)
  DAIS@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("ID_th")}
  ID_th <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.thresh(sol.idetect_seq(x, type = 'lin.cont'))
  end_time <- Sys.time()
  ID_th@nocpt <- z$no.of.cpt
  ID_th@cpt <- as.numeric(z$cpts)
  ID_th@fit <- fit_lin_cont(x, ID_th@cpt)
  ID_th@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("ID_sic")}
  ID_sic <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.ic(sol.idetect(x, type = 'lin.cont'), q.max = q_max)
  end_time <- Sys.time()
  ID_sic@nocpt <- z$no.of.cpt
  ID_sic@cpt <- as.numeric(z$cpts)
  ID_sic@fit <- fit_lin_cont(x, ID_sic@cpt)
  ID_sic@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose) print("CPOP")
  cpop <- new("cpt.est")
  start_time <- Sys.time()
  z <- cpop(x, beta = 2*log(length(x)), sd = rep(mad(diff(diff(x)))/sqrt(6), length(x)))
  end_time <- Sys.time()
  if (length(z@changepoints) == 2) {
    cpop@cpt = 0
    cpop@nocpt = 0
  } else {
    cpop@cpt = z@changepoints[-c(1, length(z@changepoints))]
    cpop@nocpt <- length(cpop@cpt)
  }
  cpop@fit <- fit_lin_cont(x, cpop@cpt)
  cpop@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose) print("not")
  not <- new("cpt.est")
  start_time <- Sys.time()
  z <- model.ic(sol.not(x, type = 'lin.cont'), q.max = q_max)
  end_time <- Sys.time()
  not@nocpt <- z$no.of.cpt
  not@cpt <- as.numeric(z$cpts)
  not@fit <- fit_lin_cont(x, not@cpt)
  not@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose) print("MARS")
  MARS <- new("cpt.est")
  start_time <- Sys.time()
  z1 <- earth(1:length(x),y=x)
  z <- sort(unique(z1$cuts[z1$selected.terms,]))[-1]
  end_time <- Sys.time()
  if (length(z) == 0) {
    MARS@cpt = 0
    MARS@nocpt = 0
  }
  else {
    MARS@cpt = z
    MARS@nocpt = length(z)
  }
  MARS@fit <- fit_lin_cont(x, MARS@cpt)
  MARS@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose) print("t1f")
  t1f <- new("cpt.est")
  start_time <- Sys.time()
  z <- tf.run(x/(mad(diff(diff(x)))/sqrt(6)))
  end_time <- Sys.time()
  if (length(z$cpt) == 0) {
    t1f@cpt = 0
    t1f@nocpt = 0
  } else {
    t1f@cpt = z$cpt
    t1f@nocpt <- length(z$cpt)
  }
  t1f@fit <- fit_lin_cont(x, t1f@cpt)
  t1f@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose) print("trendsegment")
  trendsegment <- new("cpt.est")
  start_time <- Sys.time()
  z = trendsegment(x, continuous = TRUE, indep = TRUE)$cpt
  end_time <- Sys.time()
  if (is.null(z)) {
    trendsegment@cpt = 0
    trendsegment@nocpt = 0
  } else {
    trendsegment@cpt <- as.numeric(z)
    trendsegment@nocpt <- length(z)
  }
  trendsegment@fit <- fit_lin_cont(x, trendsegment@cpt)
  trendsegment@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  if(verbose){print("local_poly")}
  local_poly <- new("cpt.est")
  start_time <- Sys.time()
  DP_result = CV.search.DP.poly(x, r = 1, gamma_set, delta = 5)
  min_idx = which.min(DP_result$test_error)
  cpt_init = unlist(DP_result$cpt_hat[min_idx])
  if(length(cpt_init) != 0){
    z <- local.refine.poly(cpt_init, x, r = 1, delta_lr = 5)
  } else {
    z <- numeric(0)
  }
  end_time <- Sys.time()
  local_poly@cpt <- as.numeric(z)
  local_poly@nocpt <- length(z)
  local_poly@fit <- fit_lin_cont(x, local_poly@cpt)
  local_poly@time <- as.numeric(difftime(end_time,start_time, units = 'secs'))
  
  list(DAIS = DAIS, ID_th = ID_th, ID_sic = ID_sic, cpop = cpop, not = not, 
       MARS = MARS, t1f = t1f, trendsegment = trendsegment, local_poly = local_poly)
}

linear_rev.sim.study <- function(model, sigma, m = 100, seed = NULL, verbose = F, gamma_set = 3:7, const = 2.1/sqrt(2)) {
  
  setClass("est.eval", representation(avg.signal = "numeric", fit = "list", cpt = "list", diff = "matrix", dh = "numeric", 
                                      cptall = "numeric", dnc = "numeric", mse = "numeric", time = "numeric"), 
           prototype(dnc = numeric(m), mse = numeric(m), time = numeric(m), dh = numeric(m)))
  
  signal = get.signal(model)
  
  no.of.cpt <- length(model$cpt)
  n <- length(signal)
  ns <- max(diff(c(0, model$cpt, n)))
  d_T <- min(diff(c(0, model$cpt, n)))
  if(d_T<1){d_T = n}
  
  q_max = max(min(100, n/log(n)), ceiling(n/d_T))
  
  if (!is.null(seed)) 
    set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    
    est <- linear.rev.sim(x, q_max = q_max, FKS_knot = no.of.cpt + 2, verbose = verbose,
                          gamma_set = gamma_set, const = const)
    
    
    if(i == 1){
      for(j in names(est)){
        eval(parse(text=paste(j, " <- new('est.eval')",sep="")))
      }
    }
    
    for(j in names(est)){
      eval(parse(text=paste(j, "@dnc[", i, "] <- est$", j, "@nocpt - no.of.cpt",sep="")))
      eval(parse(text=paste(j, "@cpt[[", i, "]] <- est$", j, "@cpt",sep="")))
      eval(parse(text=paste(j, "@mse[", i, "] <- mean((est$", j, "@fit - signal)^2)",sep="")))
      eval(parse(text=paste(j, "@diff <- abs(matrix(est$", j, "@cpt,nrow=no.of.cpt,ncol=length(est$", j, "@cpt),byr=T)-matrix(model$cpt,nrow=no.of.cpt,ncol=length(est$", j, "@cpt),byr=F))",sep="")))
      eval(parse(text=paste(j, "@dh[i] <- max(apply(", j, "@diff,1,min),apply(", j, "@diff,2,min))/ns",sep="")))
      eval(parse(text=paste(j, "@time[i] <- est$", j, "@time",sep="")))
    }
    
    gc()
  }
  
  result <- list(ts = ts)
  for(i in names(est)){
    eval(parse(text=paste("result <- c(result, ", i, " = ", i, ")",sep="")))
  }
  return(result)
}

## S9
model.wave1 <- list(name = "wave1", cpt.type = "pcwsLinContMean", cpt = c(256, 512, 768, 1024, 1152, 1280, 1344), jump.size = (-1)^(1:7) * 
                      (1:7)/64, n = 1408, start = c(1, 1/256))
wave1_sim = linear_rev.sim.study(model.wave1, 1, seed = seed.temp)

## S10
model.wave2 <- list(name = "wave2", cpt.type = "pcwsLinContMean", cpt = (1:99) * 15, jump.size = (-1)^{1:100}, 
                    n = 15 * 100, start = c(-1/2, 1/40))
wave2_sim = linear_rev.sim.study(model.wave2, 1, seed = seed.temp)

## S11
model.wave3 <- list(name = "wave3", cpt.type = "pcwsLinContMean", cpt = (1:119) * 7, jump.size = (-1)^{1:120}, 
                    n = 840, start = c(-1/2, 1/32))
wave3_sim = linear_rev.sim.study(model.wave3, m = 100, sigma = 0.3, seed = seed.temp)

## S16
model.justnoise.wave <- list(name = "justnoise_wave", cpt.type = "pcwsLinContMean",
                             cpt = numeric(0), n = 1000, start = c(0,1))
justnoise_lin_sim = linear_rev.sim.study(model.justnoise.wave, sigma = 1, seed = seed.temp)

## S17
model.wave4 <- list(name = "wave4", cpt.type = "pcwsLinContMean", cpt = (1:9) * 20,
                    jump.size = c(1/6, 1/2,-3/4,-1/3, -2/3,1,1/4,3/4,-5/4), n = 200, start = c(1, 1/32))
wave4_sim = linear_rev.sim.study(model.wave4, sigma = 0.3, seed = seed.temp)

## S18
model.wave5 <- list(name = "wave5", cpt.type = "pcwsLinContMean", cpt = (1:49) * 7,
                    jump.size = c(rep(c(-2.5,2.5),50)), n = 350, start = c(0, 1))
wave5_sim = linear_rev.sim.study(model.wave5, sigma = seed.temp)

