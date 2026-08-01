expit <- function(x){
  to_ret <- exp(x) / (1+exp(x))
  if (is.na(to_ret)){return(1)}
  return(to_ret)
}

logit <- function(x){
  return(log(x/(1-x)))
}

#Reads in rcpp file
readCpp <- function(path) {
  tryCatch(
    {
      sourceCpp(file = path)
    },
    error = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of error
      NA
    },
    warning = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of warning
      NULL
    },
    finally = {
      message("Done")
    }
  )
}

#turns vector into dummy matrix
#useful for sociodemo covar
Vec2Mat <- function(vect){
  mat <- matrix(0,nrow = length(vect),ncol = max(vect))
  
  for(i in 1:length(vect)){
    mat[i,vect[i]] <- 1
  }
  return (mat)
}

#only used for 96 period len
#used in singleday
FirstDay2SingleDay <- function(first_day,target_day){
  
  day_to_keep_vec <- numeric(864)
  
  if (first_day == target_day){
    day_to_keep_vec[673:768] <- 1
  } else {
    if (first_day > target_day){
      day_ind <- 7 - (first_day-target_day)
    } else {
      day_ind <- target_day - first_day
    }
    first_day_ind <- (96 * (day_ind)) + 1
    last_day_ind <- first_day_ind + 95
    day_to_keep_vec[first_day_ind:last_day_ind] <- 1
  }
    
  return(day_to_keep_vec)
}

#determines week/weekend
FirstDay2WeekInd <- function(first_day,
                             period_len = DEFAULT_PERIODS_PER_DAY){
  
  if (period_len == DEFAULT_PERIODS_PER_DAY){
    weekday <- numeric(DEFAULT_PERIODS_PER_DAY)
    friday <- c(rep(0,68),rep(1,28))
    saturday <- numeric(DEFAULT_PERIODS_PER_DAY)+1
    sunday <- c(rep(1,68),rep(0,28))
  } else{
    weekday <- numeric(period_len)
    friday <- c(rep(0,period_len * 2 / 3),rep(1,period_len/3))
    saturday <- numeric(period_len)+1
    sunday <- c(rep(1,period_len * 2 / 3),rep(0,period_len/3))
  }

  if (first_day == 1){
    covar_vec <- c(sunday,rep(weekday,4),friday,saturday,sunday,weekday)
  } else if (first_day == 2) {
    covar_vec <- c(rep(weekday,4),friday,saturday,sunday,rep(weekday,2))
  } else if (first_day == 3) {
    covar_vec <- c(rep(weekday,3),friday,saturday,sunday,rep(weekday,3))
  } else if (first_day == 4) {
    covar_vec <- c(rep(weekday,2),friday,saturday,sunday,rep(weekday,4))
  } else if (first_day == 5) {
    covar_vec <- c(weekday,friday,saturday,sunday,rep(weekday,4),friday)
  } else if (first_day == 6) {
    covar_vec <- c(friday,saturday,sunday,rep(weekday,4),friday,saturday)
  } else if (first_day == 7) {
    covar_vec <- c(saturday,sunday,rep(weekday,4),friday,saturday,sunday)
  }
  
  return(covar_vec)
}

#loads on rcpp functionr
#works on both pc and cluster
