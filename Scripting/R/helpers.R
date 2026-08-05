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


make_em_reference_scales <- function(emit_act, emit_light){
  act_sd <- emit_act[, 2L, , , drop = FALSE]
  light_sd <- emit_light[, 2L, , , drop = FALSE]

  list(act_mean_sd = act_sd, light_mean_sd = light_sd)
}

pack_transition_probability_values <- function(
    params_tran_array,
    period_len){

  parameter_dimensions <- dim(params_tran_array)

  mix_num <- parameter_dimensions[1L]
  vcovar_num <- parameter_dimensions[3L]

  values <- vector("list",length = mix_num * vcovar_num)

  output_index <- 1L

  for (vcovar_ind in seq_len(vcovar_num)){

    params_tran <- params_tran_array[, , vcovar_ind]

    if (mix_num == 1L){params_tran <- matrix(params_tran,  nrow = 1L)}

    for (re_ind in seq_len(mix_num)){

      transition_matrix <- Params2TranVectorTresid(
        re_ind = re_ind,
        len = period_len,
        params_tran = params_tran,
        period_len = period_len
      )

      # Column 2 is p21 and column 3 is p12.
      # The diagonal probabilities are redundant.
      values[[output_index]] <- as.numeric(
        transition_matrix[, c(2L, 3L), drop = FALSE]
      )

      output_index <- output_index + 1L
    }
  }

  values <- unlist(values, recursive = FALSE, use.names = FALSE)

  values
}

pack_em_convergence_values <- function(
    init,
    params_tran_array,
    emit_act,
    emit_light,
    corr_mat,
    nu_mat,
    beta_vec,
    surv_coef,
    period_len,
    reference_scales){

  act_mean <- emit_act[, 1L, , , drop = FALSE]
  act_sd <- emit_act[, 2L, , , drop = FALSE]

  light_mean <- emit_light[, 1L, , , drop = FALSE]
  light_sd <- emit_light[, 2L, , , drop = FALSE]

  # Only one initial-state probability is free because
  # the two state probabilities sum to one.
  init_probability <- pmin(pmax(init[, 2L], 1e-10), 1 - 1e-10)

  transition_probability <-
    pack_transition_probability_values(
      params_tran_array = params_tran_array,
      period_len = period_len
    )

  # The first mixing class and first class-survival
  # coefficient are reference parameters.
  mixing_parameters <- if (ncol(nu_mat) > 1L){
    as.numeric(nu_mat[, -1L, drop = FALSE])
  } else {
    numeric(0)
  }

  class_survival_parameters <- if (length(beta_vec) > 1L){
    as.numeric(beta_vec[-1L])
  } else {
    numeric(0)
  }

  values <- c(
    qlogis(init_probability),
    transition_probability,
    as.numeric(act_mean / reference_scales$act_mean_sd),
    log(as.numeric(act_sd)),
    as.numeric(light_mean / reference_scales$light_mean_sd),
    log(as.numeric(light_sd)),
    atanh(as.numeric(corr_mat)),
    mixing_parameters,
    class_survival_parameters,
    as.numeric(unlist(surv_coef, use.names = FALSE))
  )

  values
}