options(expressions = 500000)

library(earth)
library(stats)
library(IDetect)
library(not)
#devtools::install_github("hadley/l1tf")
library("l1tf")

get.signal <- function(model){
  if(model$cpt.type == "pcwsConstMean"){
    signal <- rep(0, model$n)
    segments <- cbind(c(1,model$cpt+1), c(model$cpt,model$n))
    signal[segments[1,1]:segments[1,2]] <- model$start[1]
    for(j in 2:nrow(segments)) signal[segments[j,1]:segments[j,2]] <- signal[segments[j,1]-1] + model$jump.size[j-1]
  }else if(model$cpt.type == "pcwsLinContMean"){
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

## Add noise to true signal
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
    ans.tf <- l1tf(z, lambda)  ##l1tf not there
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
  ans.tf <- l1tf(z, lambdas[k])  ##l1tf not there
  dans <- round(diff(ans.tf), digits = 5)
  CP <- c()
  for (ii in 1:(length(dans) - 1)) {
    if (dans[ii] != dans[ii + 1]) {
      CP <- c(CP, ii)
    }
  }
  return(list(cpt = CP, fit = ans.tf, lam = lambdas[k]))
}

## The SDLL related functions for ID
all.slopechanges.are.cpts <- function(x) {
  diff.x <- abs(diff(diff(x)))
  cpts <- which(diff.x > 0)
  no.of.cpt <- length(cpts)
  est <- x
  list(est=est, no.of.cpt=no.of.cpt, cpts=cpts)
}

idetect.th_linear <- function(x, sigma = stats::mad(diff(diff(x))) / sqrt(6), thr_const = 1.4,
                              thr_fin = sigma * thr_const * sqrt(2 * log(length(x))),
                              s = 1, e = length(x), points = 3, k_l = 1, k_r = 1) {
  l <- length(x)
  Res <- matrix(0, 1, 4)
  y <- c(0, cumsum(x))
  points <- as.integer(points)
  r_e_points <- seq(points, l, points)
  l_e_points <- seq(l - points + 1, 1, -points)
  chp <- 0
  if (e - s < 2) {
    Res_fin <- matrix(0, 1, 4)
    cpt <- 0
  } else {
    pos_r <- numeric()
    CUSUM_r <- numeric()
    pos_l <- numeric()
    CUSUM_l <- numeric()
    moving_points <- breakfast:::start_end_points(r_e_points, l_e_points, s, e)
    right_points <- moving_points[[1]]
    left_points <- moving_points[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    if (k_r < k_l) {
      while ( (chp == 0) & (k_r < min(k_l, rur))) {
        x_temp_r <- x[s:right_points[k_r]]
        ipcr <- IDetect:::cumsum_lin(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
          indic <- 0
        } else {
          k_r <- k_r + 1
        }
      }
    }
    if (k_l < k_r) {
      while ( (chp == 0) & (k_l < min(k_r, lur))) {
        x_temp_l <- x[left_points[k_l]:e]
        ipcl <- IDetect:::cumsum_lin(x_temp_l)
        pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
        CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
        Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
        if (CUSUM_l[k_l] > thr_fin) {
          chp <- pos_l[k_l]
          indic <- 1
        } else {
          k_l <- k_l + 1
        }
      }
    }
    if (chp == 0) {
      while ( (chp == 0) & (k_l <= lur) & (k_r <= rur)) {
        x_temp_r <- x[s:right_points[k_r]]
        ipcr <- IDetect:::cumsum_lin(x_temp_r)
        pos_r[k_r] <- which.max(abs(ipcr)) + s - 1
        CUSUM_r[k_r] <- abs(ipcr[pos_r[k_r] - s + 1])
        Res <- rbind(Res, c(s,right_points[k_r],pos_r[k_r],CUSUM_r[k_r]))
        if (CUSUM_r[k_r] > thr_fin) {
          chp <- pos_r[k_r]
          indic <- 0
        } else {
          x_temp_l <- x[left_points[k_l]:e]
          ipcl <- IDetect:::cumsum_lin(x_temp_l)
          pos_l[k_l] <- which.max(abs(ipcl)) + left_points[k_l] - 1
          CUSUM_l[k_l] <- abs(ipcl[pos_l[k_l] - left_points[k_l] + 1])
          Res <- rbind(Res, c(left_points[k_l], e, pos_l[k_l], CUSUM_l[k_l]))
          if (CUSUM_l[k_l] > thr_fin) {
            chp <- pos_l[k_l]
            indic <- 1
          } else {
            k_r <- k_r + 1
            k_l <- k_l + 1
          }
        }
      }
    }
    if (chp != 0) {
      if (indic == 1) {
        r <- idetect.th_linear(x, s = s, e = chp, points = points,
                               thr_fin = thr_fin, k_r = k_r, k_l = 1)
      } else {
        r <- idetect.th_linear(x, s = chp + 1, e = e, points = points,
                               thr_fin = thr_fin, k_r = 1, k_l = max(1, k_l - 1))
      }
      cpt <- c(chp, r[[1]])
      Res_fin <- rbind(Res, r[[2]])
    } else {
      cpt <- chp
      Res_fin <- Res
    }
  }
  cpt <- cpt[cpt != 0]
  Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
  return(list(changepoints = sort(cpt), full_information = Res_fin, y = y))
}
window.idetect.th_linear <- function(xd, sigma = stats::mad(diff(diff(xd))) / sqrt(6), thr_con = 1.4,
                                     c_win = 5000, w_points = 3) {
  lg <- length(xd)
  w_points <- as.integer(w_points)
  c_win <- min(lg, c_win)
  c_win <- as.integer(c_win)
  t <- sigma * thr_con * sqrt(2 * log(lg))
  if (lg <= c_win) {
    u <- idetect.th_linear(x = xd, thr_const = thr_con, points = w_points)
    return(u)
  } else {
    K <- ceiling(lg / c_win)
    tsm <- list()
    u <- list()
    ufin <- list()
    uaddition <- list()
    tsm[[1]] <- xd[1:c_win]
    ufin <- idetect.th_linear(tsm[[1]], thr_fin = t, points = w_points)
    uaddition[[1]] <- list()
    uaddition[[1]] <- idetect.th_linear(x = xd[(max(1, c_win - (10 * w_points) + 1)):min( (c_win + (10 * w_points)), lg)], thr_fin = t, points = 2)
    uaddition[[1]][[1]] <- uaddition[[1]][[1]] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,1] <- uaddition[[1]][[2]][,1] + c_win - (10 * w_points)
    uaddition[[1]][[2]][,2] <- min(uaddition[[1]][[2]][,2] + c_win - (10 * w_points),min( (c_win + (10 * w_points)), lg))
    uaddition[[1]][[2]][,3] <- uaddition[[1]][[2]][,3] + c_win - (10 * w_points)
    ufin[[1]] <- c(ufin[[1]], uaddition[[1]][[1]])
    i <- 2
    while (i < K) {
      tsm[[i]] <- xd[( (i - 1) * c_win + 1):(i * c_win)]
      u[[i]] <- list()
      u[[i]] <- idetect.th_linear(x = tsm[[i]], thr_fin = t, points = w_points)
      u[[i]][[1]] <- u[[i]][[1]] + (i - 1) * c_win
      u[[i]][[2]][,1] <- u[[i]][[2]][,1] + (i - 1) * c_win
      u[[i]][[2]][,2] <- min(u[[i]][[2]][,2] + (i - 1) * c_win, i * c_win)
      u[[i]][[2]][,3] <- u[[i]][[2]][,3] + (i - 1) * c_win
      uaddition[[i]] <- list()
      uaddition[[i]] <- idetect.th_linear(x = xd[(max(1, i * c_win - (10 * w_points) + 1)):(min(i * c_win + (10 * w_points), lg))], thr_fin = t, points = 2)
      uaddition[[i]][[1]] <- uaddition[[i]][[1]] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,1] <- uaddition[[i]][[2]][,1] + i * c_win - (10 * w_points)
      uaddition[[i]][[2]][,2] <- min(uaddition[[i]][[2]][,2] + i * c_win - (10 * w_points), min(i * c_win + (10 * w_points), lg))
      uaddition[[i]][[2]][,3] <- uaddition[[i]][[2]][,3] + i * c_win - (10 * w_points)
      ufin[[1]] <- c(ufin[[1]],u[[i]][[1]], uaddition[[i]][[1]])
      i <- i + 1
    }
    tsm[[K]] <- xd[( (K - 1) * c_win + 1):lg]
    u[[K]] <- list()
    u[[K]] <- idetect.th_linear(tsm[[K]], thr_fin = t, points = w_points)
    u[[K]][[1]] <- u[[K]][[1]]  + (K - 1) * c_win
    u[[K]][[2]][,1] <- u[[K]][[2]][,1]  + (K - 1) * c_win
    u[[K]][[2]][,2] <- min(u[[K]][[2]][,2]  + (K - 1) * c_win,lg)
    u[[K]][[2]][,3] <- u[[K]][[2]][,3]  + (K - 1) * c_win
    ufin_cpt <- c(ufin[[1]], u[[K]][[1]])
    Res_fin <- matrix(0, 1, 4)
    Res_fin <- rbind(Res_fin, ufin[[2]], uaddition[[1]][[2]])
    if (K > 2){
      for (i in 2:(K-1)){Res_fin <- rbind(Res_fin,u[[i]][[2]], uaddition[[i]][[2]])}}
    Res_fin <- rbind(Res_fin, u[[K]][[2]])
    Res_fin <- Res_fin[which(Res_fin[,3] != 0),]
    return(list(changepoints = sort(ufin_cpt), full_information = Res_fin,  y = c(0, cumsum(xd))))
  }
}
sol.idetect_linear <- function(x, thr_ic = 1.25, points = 3) {
  solutions.nested <- TRUE
  solution.set <- list()
  cands <- matrix(NA, 0, 4)
  lx <- length(x)
  if (lx < points) {solution.path <- integer()}
  else{
    points <- as.integer(points)
    step1 <- window.idetect.th_linear(x, thr_con = thr_ic, w_points = points)
    s1 <- as.matrix(step1$full_information)
    if (dim(s1)[2] == 1) {s1 <- t(s1)}
    ord <- order(s1[,4], decreasing = T)
    cands <- s1[ord, ,drop=FALSE]
    cpt_lower <- step1[[1]]
    lcpt_ic <- length(cpt_lower)
    seb_set <- c(unique(c(1, cpt_lower)), lx)
    lseb_set <- length(seb_set)
    min_C <- numeric()
    CS <- matrix(cpt_lower,1,lcpt_ic)
    while (lseb_set >= 3) {
      Rs <- IDetect:::linear_contr_one(x, seb_set[1:(lseb_set - 2)], seb_set[3:(lseb_set)],
                                       seb_set[2:(lseb_set - 1)])
      indic <- which.min(Rs)
      s1 <- setdiff(cpt_lower, seb_set[2:(lseb_set - 1)])
      d <- numeric(lcpt_ic)
      if(length(s1) == 0){d <- Rs}
      else{
        indic2 <- match(s1, cpt_lower)
        d[-indic2] <- Rs}
      CS <- rbind(CS,d)
      m_rs <- min(Rs)
      min_Rs <- seb_set[2:(lseb_set - 1)][indic]
      cands <- rbind(cands, c(seb_set[indic], seb_set[indic + 2], min_Rs, m_rs))
      min_C <- c(min_C, min_Rs)
      seb_set <- seb_set[-which(seb_set == min_Rs)]
      lseb_set <- lseb_set - 1
    }
    solution.path <- min_C[length(min_C):1]
    cusum_m <- apply(CS[-1,,drop = FALSE],2,max)
    indic3 <- match(cpt_lower, cands[,3])
    cands[indic3,4] <- cusum_m
    ord <- order(cands[,4], decreasing = T)
    cands <- cands[ord, ,drop=FALSE]#[-(length(solution.path)+1), ,drop = FALSE]
    cands <- cands[!duplicated(cands[,3]),,drop = FALSE]
    if(is.na(solution.path[1])){solution.path <- integer(0)}}
  ret = list(solutions.nested = solutions.nested, solution.path = solution.path, solution.set = solution.set, x = x, cands = cands, method = "idetect")
  
  class(ret) <- "cptpath"
  
  ret
}

model.sdll_linear <- function(cptpath.object, sigma = stats::mad(diff(diff(cptpath.object$x))) / sqrt(6), th.const = 1.25, th.const.min.mult = 0.3) {
  x <- cptpath.object$x
  n <- length(x)
  if (n <= 1) {
    est <- x
    no.of.cpt <- 0
    cpts <- integer(0)
  }
  else {
    if (sigma == 0) {
      s0 <- all.slopechanges.are.cpts(x)
      est <- s0$est
      no.of.cpt <- s0$no.of.cpt
      cpts <- s0$cpts
    } else {
      if (is.null(th.const)) {stop("th.const must be specified.")
      }
    }
    th.const.min <- th.const * th.const.min.mult
    th <- th.const * sqrt(2 * log(n)) * sigma
    th.min <- th.const.min * sqrt(2 * log(n)) * sigma
    if (cptpath.object$cands[1,4] < th) {
      no.of.cpt <- 0
      cpts <- integer(0)
    }
    else {
      indices <- which(cptpath.object$cands[,4] > th.min)
      if (length(indices) == 1) {
        cpts <- cptpath.object$cands[indices, 3]
        no.of.cpt <- 1
      }
      else {
        rc.sel <- cptpath.object$cands[indices,,drop=F]
        z <- cptpath.object$cands[indices,4]
        z.l <- length(z)
        dif <- -diff(log(z))
        dif.ord <- order(dif, decreasing=T)
        j <- 1
        while ((j < z.l) & (z[dif.ord[j]+1] > th)) j <- j+1
        if (j < z.l) no.of.cpt <- dif.ord[j] else no.of.cpt <- z.l
        cpts <- sort(cptpath.object$cands[1:no.of.cpt,3])			
      }
    }
    est <- IDetect:::est_signal(x, cpts, type = "slope")
  }
  ret <- list(est=est, no.of.cpt=no.of.cpt, cpts=cpts, model="sdll", solution.path=cptpath.object$method)
  class(ret) <- "cptmodel"
  ret
}


## Functions for DAIS
#Finds second largest difference in given sequence and returns the index
largest_diff_slope <- function(x){
  if(length(x) < 3){
    diff <- 1
  } else {
    diff <- sort.int(abs(diff(diff(x))), decreasing=TRUE, index.return = TRUE)$ix[1]
  }
  return(diff)
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

#Calculate contrast function
cumsum_lin <- function(x) {
  if (!(is.numeric(x))){
    stop("The input in `x' should be a numeric vector.")
  }
  res <- numeric()
  n <- length(x)
  if (n <= 2) {
    res <- 0
  } else {
    b <- 2:(n - 1)
    y1 <- cumsum(x * (1:n))
    y <- cumsum(x)
    a <- sqrt(6 / ( (n - 1) * n * (n + 1) * (2 - 2 * b ^ 2 + 2 * b * n - 1 + 2 * b - n)))
    be <- sqrt( ( (n - b + 1) * (n - b)) / ( (b - 1) * b))
    res[1] <- 0
    res[b] <- a * be * ( (2 * b + n - 1) * y1[b] - (n + 1) * b * y[b]) - (a / be) * ( ( 3 * n - 2 * b + 1) * (y1[n] - y1[b]) - (n + 1) * (2 * n - b) * (y[n] - y[b]))
  }
  return(res)
}

DAIS_slope <- function(x, sigma = stats::mad(diff(diff(x)))/sqrt(6), thr_const = 1.5, 
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
    cpoint <- largest_diff_slope(x[s:(e-1)]) + s - 1 #find largest differences
    if(print == TRUE){cat('largest difference at', cpoint, 'when s is equal to', s, 'and e is', e, '\n')}
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
      ipc <- cumsum_lin(x_temp)
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
      r_left <- DAIS_slope(x, s = s, e = chp, points = points, 
                           thr_fin = thr_fin, left_checked = left_checked, 
                           right_checked = right_checked, 
                           cpoints = cpoints, print = print)
      r_right <- DAIS_slope(x, s = chp+1, e = e, points = points, 
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

## Simulation functions

linear.rev.sim <- function(x, q_max_NOT = 25, FKS_knot = 10) {
  
  setClass("cpt.est", representation(cpt = "numeric", nocpt = "numeric", fit = "numeric", time = "numeric", dh = "numeric"), 
           prototype(cpt = numeric(0), nocpt = 0, fit = numeric(0), time = numeric(0), dh = numeric(0)))
  
  
  print("not")
  z <- not(x, method = "not", contrast = "pcwsLinContMean")
  cpt.ic = features(z, q.max = q_max_NOT)
  notic <- new("cpt.est")
  if (any(is.na(cpt.ic$cpt))) {
    notic@cpt = 0
    notic@nocpt = 0
  } else {
    notic@cpt = cpt.ic$cpt
    notic@nocpt <- length(notic@cpt)
  }
  notic@fit <- fit_lin_cont(x, notic@cpt)
  notic@time <- system.time(features(not(x, method = "not", contrast = "pcwsLinContMean"), q.max = q_max_NOT))[[3]]
  
  
  print("t1f")
  z <- tf.run(x/(mad(diff(diff(x)))/sqrt(6)))
  t1f <- new("cpt.est")
  if (length(z$cpt) == 0) {
    t1f@cpt = 0
    t1f@nocpt = 0
  } else {
    t1f@cpt = z$cpt
    t1f@nocpt <- length(z$cpt)
  }
  t1f@fit <- z$fit
  t1f@time <- system.time(tf.run(x/(mad(diff(diff(x)))/sqrt(6))))[[3]]
  
  print("MARS")
  z1 <- earth(1:length(x),y=x)
  z <- sort(unique(z1$cuts[z1$selected.terms,]))[-1]
  MARS <- new("cpt.est")
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
  
  print("ID")
  IsolateDetect <- new("cpt.est")
  z = ID_cplm(x)
  IsolateDetect@cpt <- as.numeric(z$cpt)
  IsolateDetect@nocpt <- z$no_cpt
  IsolateDetect@fit <- fit_lin_cont(x, IsolateDetect@cpt)
  IsolateDetect@time <- system.time(ID_cplm(x))[[3]]
  
  print("DAIS")
  DAIS <- new("cpt.est")
  z = DAIS_slope(x)
  DAIS@cpt <- as.numeric(z)
  DAIS@nocpt <- length(z)
  DAIS@fit <- fit_lin_cont(x, DAIS@cpt)
  DAIS@time <- system.time(DAIS_slope(x))[[3]]
  
  print("ID_SDLL")
  IsolateDetect_sdll <- new("cpt.est")
  z = model.sdll_linear(sol.idetect_linear(x))
  IsolateDetect_sdll@cpt <- as.numeric(z$cpts)
  IsolateDetect_sdll@nocpt <- z$no.of.cpt
  IsolateDetect_sdll@fit <- IDetect:::est_signal(x, IsolateDetect_sdll@cpt, type = "slope")
  IsolateDetect_sdll@time <- system.time(model.sdll_linear(sol.idetect_linear(x)))[[3]]
  
  list(notic = notic, MARS = MARS, t1f = t1f, IsolateDetect = IsolateDetect, 
       IsolateDetect_sdll = IsolateDetect_sdll, DAIS = DAIS)
}

linear_rev.sim.study <- function(model, sigma, m = 100, seed = NULL, gen_qmax = 25) {
  
  setClass("est.eval", representation(avg.signal = "numeric", fit = "list", cpt = "list", diff = "matrix", dh = "numeric", 
                                      cptall = "numeric", dnc = "numeric", mse = "numeric", time = "numeric"), 
           prototype(dnc = numeric(m), mse = numeric(m), time = numeric(m), dh = numeric(m)))
  
  DAIS <- new("est.eval")
  IsolateDetect <- new("est.eval")
  IsolateDetect_sdll <- new("est.eval")
  notic <- new("est.eval")
  t1f <- new("est.eval")
  MARS <- new("est.eval")
  
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
    
    notic@dnc[i] <- est$notic@nocpt - no.of.cpt
    notic@cpt[[i]] <- est$notic@cpt
    notic@mse[i] <- mean((signal - est$notic@fit)^2)
    notic@diff <- abs(matrix(est$notic@cpt, nrow = no.of.cpt, ncol = length(est$notic@cpt), byr = T) - matrix(model$cpt, 
                                                                                                              nrow = no.of.cpt, ncol = length(est$notic@cpt), byr = F))
    notic@dh[i] <- max(apply(notic@diff, 1, min), apply(notic@diff, 2, min))/ns
    notic@time[i] <- est$notic@time
    
    t1f@dnc[i] <- est$t1f@nocpt - no.of.cpt
    t1f@cpt[[i]] <- est$t1f@cpt
    t1f@mse[i] <- mean((signal - est$t1f@fit)^2)
    t1f@diff <- abs(matrix(est$t1f@cpt, nrow = no.of.cpt, ncol = length(est$t1f@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                        ncol = length(est$t1f@cpt), byr = F))
    t1f@dh[i] <- max(apply(t1f@diff, 1, min), apply(t1f@diff, 2, min))/ns
    t1f@time[i] <- est$t1f@time
    
    MARS@dnc[i] <- est$MARS@nocpt - no.of.cpt
    MARS@cpt[[i]] <- est$MARS@cpt
    MARS@mse[i] <- mean((signal - est$MARS@fit)^2)
    MARS@diff <- abs(matrix(est$MARS@cpt, nrow = no.of.cpt, ncol = length(est$MARS@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                           ncol = length(est$MARS@cpt), byr = F))
    MARS@dh[i] <- max(apply(MARS@diff, 1, min), apply(MARS@diff, 2, min))/ns
    MARS@time[i] <- est$MARS@time
    
    IsolateDetect@dnc[i] <- est$IsolateDetect@nocpt - no.of.cpt
    IsolateDetect@cpt[[i]] <- est$IsolateDetect@cpt
    IsolateDetect@mse[i] <- mean((signal - est$IsolateDetect@fit)^2)
    IsolateDetect@diff <- abs(matrix(est$IsolateDetect@cpt, nrow = no.of.cpt, ncol = length(est$IsolateDetect@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                                                      ncol = length(est$IsolateDetect@cpt), byr = F))
    IsolateDetect@dh[i] <- max(apply(IsolateDetect@diff, 1, min), apply(IsolateDetect@diff, 2, min))/ns
    IsolateDetect@time[i] <- est$IsolateDetect@time
    
    IsolateDetect_sdll@dnc[i] <- est$IsolateDetect_sdll@nocpt - no.of.cpt
    IsolateDetect_sdll@cpt[[i]] <- est$IsolateDetect_sdll@cpt
    IsolateDetect_sdll@mse[i] <- mean((signal - est$IsolateDetect_sdll@fit)^2)
    IsolateDetect_sdll@diff <- abs(matrix(est$IsolateDetect_sdll@cpt, nrow = no.of.cpt, ncol = length(est$IsolateDetect_sdll@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                                                                     ncol = length(est$IsolateDetect_sdll@cpt), byr = F))
    IsolateDetect_sdll@dh[i] <- max(apply(IsolateDetect_sdll@diff, 1, min), apply(IsolateDetect_sdll@diff, 2, min))/ns
    IsolateDetect_sdll@time[i] <- est$IsolateDetect_sdll@time
    
    DAIS@dnc[i] <- est$DAIS@nocpt - no.of.cpt
    DAIS@cpt[[i]] <- est$DAIS@cpt
    DAIS@mse[i] <- mean((signal - est$DAIS@fit)^2)
    DAIS@diff <- abs(matrix(est$DAIS@cpt, nrow = no.of.cpt, ncol = length(est$DAIS@cpt), byr = T) - matrix(model$cpt, nrow = no.of.cpt, 
                                                                                                           ncol = length(est$DAIS@cpt), byr = F))
    DAIS@dh[i] <- max(apply(DAIS@diff, 1, min), apply(DAIS@diff, 2, min))/ns
    DAIS@time[i] <- est$DAIS@time
    
  }
  list(notic = notic, MARS = MARS, t1f = t1f, IsolateDetect = IsolateDetect, 
       IsolateDetect_sdll = IsolateDetect_sdll, DAIS = DAIS)
}


## Present results
make_df <- function(x, decimals = 3){
  results <- data.frame('Method' = c('DAIS', 'ID', 'ID_sdll', 'MARS', 't1f', 'NOT'),
                        'MSE' = rep(NA),
                        'dh' = rep(NA),
                        'time' = rep(NA))
  
  results$MSE[1] <- signif(mean(x$DAIS@mse), decimals)
  results$MSE[2] <- signif(mean(x$IsolateDetect@mse), decimals)
  results$MSE[3] <- signif(mean(x$IsolateDetect_sdll@mse), decimals)
  results$MSE[4] <- signif(mean(x$MARS@mse), decimals)
  results$MSE[5] <- signif(mean(x$t1f@mse), decimals)
  results$MSE[6] <- signif(mean(x$notic@mse), decimals)
  
  results$dh[1] <- signif(mean(x$DAIS@dh), decimals)
  results$dh[2] <- signif(mean(x$IsolateDetect@dh), decimals)
  results$dh[3] <- signif(mean(x$IsolateDetect_sdll@dh), decimals)
  results$dh[4] <- signif(mean(x$MARS@dh), decimals)
  results$dh[5] <- signif(mean(x$t1f@dh), decimals)
  results$dh[6] <- signif(mean(x$notic@dh), decimals)
  
  results$time[1] <- signif(mean(x$DAIS@time), decimals)
  results$time[2] <- signif(mean(x$IsolateDetect@time), decimals)
  results$time[3] <- signif(mean(x$IsolateDetect_sdll@time), decimals)
  results$time[4] <- signif(mean(x$MARS@time), decimals)
  results$time[5] <- signif(mean(x$t1f@time), decimals)
  results$time[6] <- signif(mean(x$notic@time), decimals)
  
  return(results)
}

cpts_df <- function(x, breaks){
  x1 <- cut(x$DAIS@dnc, breaks = breaks)
  x2 <- cut(x$IsolateDetect@dnc, breaks = breaks)
  x3 <- cut(x$IsolateDetect_sdll@dnc, breaks = breaks)
  x4 <- cut(x$MARS@dnc, breaks = breaks)
  x5 <- cut(x$t1f@dnc, breaks = breaks)
  x6 <- cut(x$notic@dnc, breaks = breaks)
  cname <- c('Method', levels(x1))
  df <- data.frame(matrix(ncol = length(cname), nrow=6))
  colnames(df) <- cname
  df$Method <- c('DAIS', 'ID', 'ID_sdll', 'MARS', 't1f', 'NOT')
  
  for(i in 2:(length(cname))){
    df[1,i] <- sum(x1 == cname[i])
    df[2,i] <- sum(x2 == cname[i])
    df[3,i] <- sum(x3 == cname[i])
    df[4,i] <- sum(x4 == cname[i])
    df[5,i] <- sum(x5 == cname[i])
    df[6,i] <- sum(x6 == cname[i])
  }
  return(df)
}

make_table <- function(x, breaks){
  df1 <- make_df(x)
  df2 <- cpts_df(x, breaks)
  df <- cbind(df2, df1[,-1])
  return(df)
}


## Signals
model.wave1 <- list(name = "wave1", cpt.type = "pcwsLinContMean", cpt = c(256, 512, 768, 1024, 1152, 1280, 1344), jump.size = (-1)^(1:7) * 
                      (1:7)/64, n = 1408, start = c(1, 1/256))
lin.SIMR1 = linear_rev.sim.study(model.wave1, 1, seed = 16)
make_table(lin.SIMR1, breaks = c(-100,-15,-2,-1,0,1,14, 100))

model.wave2 <- list(name = "wave2", cpt.type = "pcwsLinContMean", cpt = (1:99) * 15, jump.size = (-1)^{1:100}, n = 15 * 100, start = c(-1/2, 1/40))
lin.SIMR2 = linear_rev.sim.study(model.wave2, 1, seed = 16)
make_table(lin.SIMR2, breaks = c(-100,-15,-2,-1,0,1,14, 115))

model.wave3 <- list(name = "wave3", cpt.type = "pcwsLinContMean", cpt = (1:119) * 7, jump.size = (-1)^{1:120}, n = 840, start = c(-1/2, 1/32))
lin.SIMR3 = linear_rev.sim.study(model.wave3, m = 100, sigma = 0.3, seed = 16)
make_table(lin.SIMR3, breaks = c(-500,-15,-2,-1,0,1,14, 115))

model.wave4 <- list(name = "wave4", cpt.type = "pcwsLinContMean", cpt = (1:9) * 20,
                    jump.size = c(1/6, 1/2,-3/4,-1/3, -2/3,1,1/4,3/4,-5/4), n = 200, start = c(1, 1/32))
lin.SIMR4 = linear_rev.sim.study(model.wave4, m = 100, sigma = 0.3, seed = 16)
make_table(lin.SIMR4, breaks = c(-100,-15,-2,-1,0,1,14, 115))

model.wave5 <- list(name = "wave5", cpt.type = "pcwsLinContMean", cpt = (1:19) * 50,
                    jump.size = c(-1/16,-5/16,-5/8,1, 5/16,15/32,-5/8,-7/32,-3/4,13/16,5/16,19/32,-1,-5/8,
                                  23/32,1/2,15/16,-25/16,-5/4), n = 1000, start = c(1, 1/32))
lin.SIMR5 = linear_rev.sim.study(model.wave5, m = 100, sigma = 0.6, seed = 5)
make_table(lin.SIMR5, breaks = c(-100,-15,-2,-1,0,1,14, 115))

