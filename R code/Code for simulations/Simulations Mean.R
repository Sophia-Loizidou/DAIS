library(devtools)
library(IDetect)
library(wbs)
library(breakfast)
library(changepoint)
library(changepoint.np)
library(not)

## Functions for DAIS
#Finds largest difference in given sequence and returns the index
largest_diff <- function(x){
  sort.int(abs(diff(x)), decreasing=TRUE, index.return = TRUE)$ix[1]
}

#list of left and right expanding intervals around the point
endpoints <- function(l_diff, s, e, points = 3){
  intervals <- list()
  intervals[[1]] <- c(seq(l_diff, s, -points))
  if (intervals[[1]][length(intervals[[1]])] != s){intervals[[1]] = c(intervals[[1]], s)}
  intervals[[2]] <- seq(min(e, l_diff + points - 1), e, points)
  if (intervals[[2]][length(intervals[[2]])] != e){intervals[[2]] = c(intervals[[2]], e)}
  return(intervals)
}

#Calculate CUSUM
inner_prod_cumsum <- function(x, y = cumsum(x)) {
  if(length(x) == 1){
    res <- 0
  }
  else{
    n <- length(x)
    res <- sqrt( ( (n-1) : 1) / n / (1 : (n-1))) * y[1 : (n-1)] -
      sqrt( (1 : (n-1)) / n / ( (n-1):1)) * (y[n] - y[ 1 : (n-1)])
  }
  return(res)
}

DAIS_mean <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), thr_const = 1.2, 
                      thr_fin = sigma * thr_const * sqrt(2 * log(length(x))), s = 1, 
                      e = length(x), points = 3, k_l = 1, k_r = 1,
                      left_checked = numeric(), right_checked = numeric(), 
                      cpoints = numeric(), print = FALSE){
  points <- as.integer(points)
  l <- length(x)
  chp <- 0
  if (e - s <= 1) {
    cpt <- 0
  }
  else{
    cpoint <- largest_diff(x[s:(e-1)]) + s - 1 #find largest differences
    if(print == TRUE){cat('largest difference at', cpoint, 'when s is equal to', s, 
                          'and e is', e, '\n')}
    endpoints <- endpoints(l_diff=cpoint, s=s, e=e, points = points)
    left_points <- endpoints[[1]]
    right_points <- endpoints[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    k_l_temp <- 1
    k_r_temp <- 1
    if(any(cpoint %in% cpoints)){
      pos <- which(cpoints == cpoint)
      for(i in 1:length(pos)){
        if(left_checked[pos[i]] < left_points[lur]){
          k_l_temp <- lur + 1
          if(right_checked[pos[i]] > right_points[rur]){
            k_r_temp <- rur + 1
          } else {
            last_checked_r <- which(right_points == right_checked[pos[i]])
            k_r_temp <- max(last_checked_r, k_r)
            k_r_temp <- min(lur, k_r_temp)
          }
        } else {
          last_checked_l <- which(left_points == left_checked[pos[i]])
          k_l_temp <- max(last_checked_l, k_l)
          if(right_checked[pos[i]] > right_points[rur]){
            k_r_temp <- rur + 1
            k_l_temp <- min(k_l_temp, rur)
          } else {
            last_checked_r <- which(right_points == right_checked[pos[i]])
            k_r_temp <- max(last_checked_r, k_r)
          }
        }
      }
    }
    k_l_temp <- min(k_l_temp, lur)
    k_r_temp <- min(k_r_temp, rur)
    flag_l <- 0
    flag_r <- 0
    while((chp == 0) & ((k_l_temp <= lur) | (k_r_temp <= rur))){
      if(print == TRUE){cat(left_points[k_l_temp],right_points[k_r_temp], '\n')}
      x_temp <- x[left_points[k_l_temp]:right_points[k_r_temp]]
      ipc <- inner_prod_cumsum(x_temp)
      pos <- which.max(abs(ipc)) + left_points[k_l_temp] - 1
      CUSUM <- abs(ipc[pos - left_points[k_l_temp] + 1])
      if (CUSUM > thr_fin) {
        chp <- pos
        if(print == TRUE){cat('cpt', chp, "detected in the interval [", left_points[k_l_temp],
                              ',', right_points[k_r_temp], '] \n')}
      }
      else {
        flag_l <- 0
        flag_r <- 0
        if((k_l_temp < lur) & (k_r_temp < rur)){
          if(k_l_temp == k_r_temp){
            k_l_temp = k_l_temp + 1
            flag_l <- 1
          }
          else {
            k_r_temp = k_r_temp+1
            flag_r <- 1
          } 
        } else if(k_l_temp < lur){
          k_l_temp = k_l_temp + 1
          flag_l <- 1
        } else if(k_r_temp < rur){
          k_r_temp = k_r_temp + 1
          flag_r <- 1
        } else{
          k_l_temp = k_l_temp + 1
          k_r_temp = k_r_temp + 1
          flag_l <- 1
          flag_r <- 1
        }
      }
    }
    cpoints <- c(cpoints, cpoint)
    if(chp != 0){
      left_checked <- c(left_checked, left_points[k_l_temp])
      right_checked <- c(right_checked, right_points[k_r_temp])
    } else {
      if(flag_l == 1){
        left_checked <- c(left_checked, left_points[k_l_temp - 1])
      } else {
        left_checked <- c(left_checked, left_points[k_l_temp])
      } 
      if (flag_r == 1){
        right_checked <- c(right_checked, right_points[k_r_temp - 1])
      } else {
        right_checked <- c(right_checked, right_points[k_r_temp])
      }
    }
    if (chp != 0) {
      r_left <- DAIS_mean(x, s = s, e = chp, points = points, 
                          thr_fin = thr_fin, left_checked = left_checked, 
                          right_checked = right_checked, 
                          cpoints = cpoints, print = print)
      r_right <- DAIS_mean(x, s = chp+1, e = e, points = points, 
                           thr_fin = thr_fin, left_checked = left_checked, 
                           right_checked = right_checked, 
                           cpoints = cpoints, print = print)
      cpt <- c(chp, r_left, r_right)
    }
    else {
      cpt <- chp
    }
  }
  cpt <- cpt[cpt != 0]
  return(sort(cpt))
}

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

sim_thr <- function(x, const = 1.15, points, Kmax_wbs, qmax_NOT, max_two_sides = max_two_sides){
  setClass("cpt.est", representation(cpt="numeric", nocpt="numeric",det_time="numeric", time="numeric", mean="numeric"), 
           prototype(cpt=numeric(0), nocpt=0, det_time=numeric(0), time=numeric(0), mean=numeric(0)))
  
  print("ID")
  ID <- new("cpt.est")
  z <- ID(x)
  if(z$no_cpt == 0){ID@cpt = 0
  ID@nocpt = 0}
  else{ID@cpt <- z$cpt
  ID@nocpt <- z$no_cpt}
  ID@time <- system.time(ID(x))[[3]]
  ID@mean <- mean.from.cpt(x, ID@cpt)
  
  print("DAIS")
  DAIS <- new("cpt.est")
  z <- DAIS_mean(x, thr_const = const, points=points)
  if(length(z) == 0){DAIS@cpt = 0
  DAIS@nocpt = 0}
  else{DAIS@cpt <- z
  DAIS@nocpt <- length(z)}
  DAIS@time <- system.time(DAIS_mean(x, thr_const = const, points=points))[[3]]
  DAIS@mean <- mean.from.cpt(x, DAIS@cpt)
  
  print("wbs")
  z <- wbs(x)
  cpt.z_1 = changepoints(z,th.const=1,Kmax=Kmax_wbs)
  wbsth10 <- new("cpt.est")
  if(any(is.na(cpt.z_1$cpt.th[[1]]))){wbsth10@cpt = 0}
  else{wbsth10@cpt <- cpt.z_1$cpt.th[[1]]}
  wbsth10@nocpt <- cpt.z_1$no.cpt.th
  wbsth10@mean <- mean.from.cpt(x, wbsth10@cpt)
  wbsth10@time <- system.time(changepoints(wbs(x,5000),th.const=1,Kmax=Kmax_wbs))[[3]]
  
  wbssbic <- new("cpt.est")
  if(any(is.na(cpt.z_1$cpt.ic[[1]]))) {wbssbic@cpt = 0}
  else {wbssbic@cpt= cpt.z_1$cpt.ic[[1]]}
  wbssbic@nocpt <- cpt.z_1$no.cpt.ic[[1]]
  wbssbic@mean <- mean.from.cpt(x, wbssbic@cpt) ## The computational time is the same as in wbs10
  
  print("wbs2")
  wbs2 <- new("cpt.est")
  z <-  model.sdll(sol.wbs2(x))
  wbs2@nocpt <- z$no.of.cpt
  wbs2@cpt <- as.numeric(z$cpts)
  if(wbs2@nocpt == 0){wbs2@cpt <- 0}
  wbs2@mean <- mean.from.cpt(x, wbs2@cpt)
  wbs2@time <- system.time(model.sdll(sol.wbs2(x)))[[3]]
  
  print("pelt")
  z <- cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT")
  pelt <- new("cpt.est")
  if (any(z@cpts[1:(length(z@cpts)-1)] == length(x))){pelt@cpt = 0
  pelt@nocpt = 0}else{
    pelt@cpt <- as.numeric(z@cpts[1:(length(z@cpts)-1)])
    pelt@nocpt <- length(pelt@cpt)}
  pelt@mean <- mean.from.cpt(x, pelt@cpt)
  pelt@time <- system.time(cpt.mean(x/mad(diff(x)/sqrt(2)), method="PELT"))[[3]]
  
  print("not")
  z <- not(x,method="not",contrast="pcwsConstMean")
  cpt.ic = features(z,q.max=qmax_NOT)
  notic <- new("cpt.est")
  if(any(is.na(cpt.ic$cpt))) {notic@cpt = 0
  notic@nocpt <- 0}
  else {notic@cpt= cpt.ic$cpt
  notic@nocpt <- length(notic@cpt)}
  notic@mean <- mean.from.cpt(x, notic@cpt)
  notic@time <- system.time(features(not(x,method="not",contrast="pcwsConstMean"),q.max = qmax_NOT))[[3]]
  
  list(ID = ID, DAIS = DAIS, wbsth10 = wbsth10, 
       wbssbic = wbssbic, wbs2 = wbs2, pelt = pelt, notic = notic)
}

sim.study.thr <- function(signal, true.cpt=NULL,sigma, m = 100, seed = NULL, const = 1.15, points=3,
                          Kmax_wbs = 50, qmax_NOT = 25, max_two_sides = 1000) {
  
  setClass("est.eval", representation(avg.signal="numeric",mean="list",cpt="list",diff="matrix",dh="numeric",cptall="numeric",dnc="numeric",
                                      mse="numeric", time= "numeric", det_delay="numeric"), 
           prototype(dnc=numeric(m), mse=numeric(m),time=numeric(m),dh=numeric(m)))
  
  ID <- new("est.eval")
  DAIS <- new("est.eval")
  wbsth10 <- new("est.eval")
  wbssbic <- new("est.eval")
  wbs2 <- new("est.eval")
  pelt <- new("est.eval")
  notic <- new("est.eval")
  ts <- list()
  
  no.of.cpt <- sum(abs(diff(signal)) > 0)
  n <- length(signal)
  ns <- max(c(diff(true.cpt), length(signal)))
  
  if (!is.null(seed)) set.seed(seed)
  
  for (i in 1:m) {
    
    print(i)
    x <- signal + sigma * rnorm(n)
    ts[[i]] <- x
    
    est <- sim_thr(x, const = const, points=points, Kmax_wbs = Kmax_wbs, qmax_NOT = qmax_NOT, 
                   max_two_sides = max_two_sides)
    
    ID@dnc[i] <- est$ID@nocpt - no.of.cpt
    ID@mse[i] <- mean((est$ID@mean - signal)^2)
    ID@cpt[[i]] <- est$ID@cpt
    ID@diff <- abs(matrix(est$ID@cpt,nrow=no.of.cpt,ncol=length(est$ID@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$ID@cpt),byr=F))
    ID@dh[i] <- max(apply(ID@diff,1,min),apply(ID@diff,2,min))/ns
    ID@time[i] <- est$ID@time
    
    DAIS@dnc[i] <- est$DAIS@nocpt - no.of.cpt
    DAIS@mse[i] <- mean((est$DAIS@mean - signal)^2)
    DAIS@cpt[[i]] <- est$DAIS@cpt
    DAIS@diff <- abs(matrix(est$DAIS@cpt,nrow=no.of.cpt,ncol=length(est$DAIS@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$DAIS@cpt),byr=F))
    DAIS@dh[i] <- max(apply(DAIS@diff,1,min),apply(DAIS@diff,2,min))/ns
    DAIS@time[i] <- est$DAIS@time
    
    wbsth10@dnc[i] <- est$wbsth10@nocpt - no.of.cpt
    wbsth10@mse[i] <- mean((est$wbsth10@mean - signal)^2)
    wbsth10@cpt[[i]] <- est$wbsth10@cpt
    wbsth10@diff <- abs(matrix(est$wbsth10@cpt,nrow=no.of.cpt,ncol=length(est$wbsth10@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbsth10@cpt),byr=F))
    wbsth10@dh[i] <- max(apply(wbsth10@diff,1,min),apply(wbsth10@diff,2,min))/ns
    wbsth10@time[i] <- est$wbsth10@time
    
    wbssbic@dnc[i] <- est$wbssbic@nocpt - no.of.cpt
    wbssbic@mse[i] <- mean((est$wbssbic@mean - signal)^2)
    wbssbic@cpt[[i]] <- est$wbssbic@cpt
    wbssbic@diff <- abs(matrix(est$wbssbic@cpt,nrow=no.of.cpt,ncol=length(est$wbssbic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$wbssbic@cpt),byr=F))
    wbssbic@dh[i] <- max(apply(wbssbic@diff,1,min),apply(wbssbic@diff,2,min))/ns
    wbssbic@time[i] <- wbsth10@time[i]
    
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
    
    notic@dnc[i] <- est$notic@nocpt - no.of.cpt
    notic@mse[i] <- mean((est$notic@mean - signal)^2)
    notic@cpt[[i]] <- est$notic@cpt
    notic@diff <- abs(matrix(est$notic@cpt,nrow=no.of.cpt,ncol=length(est$notic@cpt),byr=T)-matrix(true.cpt,nrow=no.of.cpt,ncol=length(est$notic@cpt),byr=F))
    notic@dh[i] <- max(apply(notic@diff,1,min),apply(notic@diff,2,min))/ns
    notic@time[i] <- est$notic@time
    
    gc()
  }
  
  list(ts = ts, ID = ID, DAIS = DAIS, wbsth10 = wbsth10, 
       wbssbic = wbssbic, wbs2 = wbs2, pelt = pelt, notic = notic)
}

make_df <- function(x, decimals = 3){
  results <- data.frame('Method' = c('DAIS', 'ID', 'WBSTH', 'WBSIC', 'WBS2', 'PELT', 'NOT'),
                        'MSE' = rep(NA),
                        'd_h' = rep(NA),
                        'time' = rep(NA))
  
  results$MSE[1] <- signif(mean(x$DAIS@mse), decimals)
  results$MSE[2] <- signif(mean(x$ID@mse), decimals)
  results$MSE[3] <- signif(mean(x$wbsth10@mse), decimals)
  results$MSE[4] <- signif(mean(x$wbssbic@mse), decimals)
  results$MSE[5] <- signif(mean(x$wbs2@mse), decimals)
  results$MSE[6] <- signif(mean(x$pelt@mse), decimals)
  results$MSE[7] <- signif(mean(x$notic@mse), decimals)
  
  results$d_h[1] <- signif(mean(x$DAIS@dh), decimals)
  results$d_h[2] <- signif(mean(x$ID@dh), decimals)
  results$d_h[3] <- signif(mean(x$wbsth10@dh), decimals)
  results$d_h[4] <- signif(mean(x$wbssbic@dh), decimals)
  results$d_h[5] <- signif(mean(x$wbs2@dh), decimals)
  results$d_h[6] <- signif(mean(x$pelt@dh), decimals)
  results$d_h[7] <- signif(mean(x$notic@dh), decimals)
  
  results$time[1] <- signif(mean(x$DAIS@time), decimals)
  results$time[2] <- signif(mean(x$ID@time), decimals)
  results$time[3] <- signif(mean(x$wbsth10@time), decimals)
  results$time[4] <- signif(mean(x$wbssbic@time), decimals)
  results$time[5] <- signif(mean(x$wbs2@time), decimals)
  results$time[6] <- signif(mean(x$pelt@time), decimals)
  results$time[7] <- signif(mean(x$notic@time), decimals)
  
  return(results)
}

cpts_df <- function(x, breaks){
  x1 <- cut(x$DAIS@dnc, breaks = breaks)
  x2 <- cut(x$ID@dnc, breaks = breaks)
  x3 <- cut(x$wbsth10@dnc, breaks = breaks)
  x4 <- cut(x$wbssbic@dnc, breaks = breaks)
  x5 <- cut(x$wbs2@dnc, breaks = breaks)
  x6 <- cut(x$pelt@dnc, breaks = breaks)
  x7 <- cut(x$notic@dnc, breaks = breaks)
  cname <- c('Method', levels(x1))
  df <- data.frame(matrix(ncol = length(cname), nrow=7))
  colnames(df) <- cname
  df$Method <- c('DAIS', 'ID', 'WBSTH', 'WBSIC', 'WBS2', 'PELT', 'NOT')
  
  for(i in 2:(length(cname))){
    df[1,i] <- sum(x1 == cname[i])
    df[2,i] <- sum(x2 == cname[i])
    df[3,i] <- sum(x3 == cname[i])
    df[4,i] <- sum(x4 == cname[i])
    df[5,i] <- sum(x5 == cname[i])
    df[6,i] <- sum(x6 == cname[i])
    df[7,i] <- sum(x7 == cname[i])
  }
  return(df)
}

make_table <- function(x, breaks){
  df1 <- make_df(x)
  df2 <- cpts_df(x, breaks)
  df <- cbind(df2, df1[,-1])
  return(df)
}


seed.temp=15

justnoise = rep(0,6000)
NC.small = sim.study.thr(justnoise,sigma=1,true.cpt=c(0),seed=seed.temp, const = 1.2, max_two_sides = 500) ## GOOD RESULTS, with 1.2
make_table(NC.small, breaks = c(-1,0,1,2,50))

long_signal <- c(rep(0,5500), rep(1.5,5500))
SIMR.large <- sim.study.thr(long_signal,true.cpt=5500,sigma=1,seed=seed.temp, points=3, 
                            const = 1.2, max_two_sides = 100) ## 1.2
make_table(SIMR.large, breaks = c(-1,0,1,2,10))

small_dist <- c(rep(0, 485), rep(1, 30), rep(0, 485))
SIMR13.small <- sim.study.thr(small_dist,true.cpt=c(485, 515),sigma=1,seed=seed.temp,m=100, 
                              const = 1.2, max_two_sides =100)
make_table(SIMR13.small, breaks = c(-3,-2,-1,0,1,2,50))

small_dist2 <- c(rep(0, 485), rep(1, 30), rep(0, 385), rep(1.5,30), rep(0,70))
SIMR14.small <- sim.study.thr(small_dist2,true.cpt=c(485, 515, 900, 930),sigma=1,seed=seed.temp,m=100, 
                              const = 1.2, max_two_sides =100)
make_table(SIMR14.small, breaks = c(-10,-4,-3,-2,-1,0,1,2,50))

small_dist3 <- c(rep(0,100), rep(1.5,30), rep(0, 355), rep(1, 30), rep(0, 355), rep(1.5,30), rep(0,100))
SIMR15.small <- sim.study.thr(small_dist3,true.cpt=c(100, 130, 485, 515, 870, 900),sigma=1,seed=seed.temp,m=100, 
                              const = 1.2, max_two_sides =100)
make_table(SIMR15.small, breaks = c(-10,-3,-2,-1,0,1,2,50))

stairs10 = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10),rep(7,10),rep(8,10),rep(9,10),rep(10,10),rep(11,10),rep(12,10),rep(13,10),rep(14,10),rep(15,10))
SIMR5.small = sim.study.thr(stairs10,true.cpt=seq(11,141,10),sigma=0.3,seed=seed.temp, 
                            const = 1.2, max_two_sides = 100)
make_table(SIMR5.small, breaks = c(-10,-3,-2,-1,0,1,2,10))

blocks = c(rep(0,205),rep(14.64,62),rep(-3.66,41),rep(7.32,164),rep(-7.32,40),rep(10.98,308),
           rep(-4.39,82),rep(3.29,430),rep(19.03,225),rep(7.68,41),rep(15.37,61),rep(0,389))
SIMR1.small2 = sim.study.thr(blocks,true.cpt = which(diff(blocks) != 0),sigma=10,
                             seed=seed.temp, const = 1.2, m=100, max_two_sides = 100)
make_table(SIMR1.small2, breaks = c(-10,-5,-2,-1,0,1,4,10))

mix2 = c(rep(7,11),rep(-7,10),rep(6,20),rep(-6,20),rep(5,30),rep(-5,30),rep(4,40),rep(-4,40),
         rep(3,50),rep(-3,50))
SIMR3.small2 = sim.study.thr(mix2,true.cpt=c(11,21,41,61,91,121,161,201,251), sigma=4,seed=seed.temp, 
                             const = 1.2, max_two_sides = 100)
make_table(SIMR3.small2, breaks = c(-10,-3,-2,-1,0,1,2,10))

teeth10 = c(rep(0,11),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,20),
            rep(0,20),rep(1,20),rep(0,20),rep(1,20),rep(0,20),rep(1,19))
SIMR4.small = sim.study.thr(teeth10,true.cpt=seq(11,251,20),0.4,seed=seed.temp, 
                            const = 1.2, max_two_sides = 100)
make_table(SIMR4.small, breaks = c(-10,-3,-2,-1,0,1,2,3,10)) 

extr_short_2_cpt <- c(rep(0,15), rep(2,15), rep(-0.5,10))
SIMR8.small = sim.study.thr(extr_short_2_cpt,true.cpt=c(15,30),sigma=1,seed=seed.temp, 
                            const = 1.2, max_two_sides = 100)
make_table(SIMR8.small, breaks = c(-3,-2,-1,0,1,2,10,46))

one_spike_in_middle <- c(rep(0,999),rep(100, 2), rep(0,999))
SIMR11.small <- sim.study.thr(one_spike_in_middle,true.cpt=c(1000, 1001),sigma=1,seed=seed.temp, m=100)
make_table(SIMR11.small, breaks = c(-3,-2,-1,0,1,2,10,46))

