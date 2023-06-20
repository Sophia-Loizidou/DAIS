#Finds largest difference in given sequence and returns the index
largest_diff <- function(x){
  sort.int(abs(diff(x)), decreasing=TRUE, index.return = TRUE)$ix[1]
}

#list of left and right expanding intervals around the largest difference
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
    # create endpoints of intervals to be checked
    endpoints <- endpoints(l_diff=cpoint, s=s, e=e, points = points)
    left_points <- endpoints[[1]]
    right_points <- endpoints[[2]]
    lur <- length(left_points)
    rur <- length(right_points)
    k_l_temp <- 1
    k_r_temp <- 1
    # check if some intervals have already been checked
    if(any(cpoint %in% cpoints)){
      pos <- which(cpoints == cpoint)
      # check all end-points that were expanded around the same largest difference
      for(i in 1:length(pos)){
      #   if(left_checked[pos[i]] < left_points[lur]){
      #     k_l_temp <- lur + 1
      #   } else {
      #     last_checked_l <- which(left_points == left_checked[pos[i]])
      #     k_l_temp <- max(last_checked_l, k_l)
      #   }
      #   if(right_checked[pos[i]] > right_points[rur]){
      #     k_r_temp <- rur + 1
      #   } else {
      #     last_checked_r <- which(right_points == right_checked[pos[i]])
      #     k_r_temp <- max(last_checked_r, k_r)
      #   }
      # }
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

############################ Slope ####################################
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
#same as for changes in the mean
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
      #   if(left_checked[pos[i]] < left_points[lur]){
      #     k_l_temp <- lur + 1
      #   } else {
      #     last_checked_l <- which(left_points == left_checked[pos[i]])
      #     k_l_temp <- max(last_checked_l, k_l)
      #   }
      #   if(right_checked[pos[i]] > right_points[rur]){
      #     k_r_temp <- rur + 1
      #   } else {
      #     last_checked_r <- which(right_points == right_checked[pos[i]])
      #     k_r_temp <- max(last_checked_r, k_r)
      #   }
      #   k_l_temp <- min(k_l_temp, k_r_temp)
      #   k_r_temp <- min(k_l_temp, k_r_temp)
      # }
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

