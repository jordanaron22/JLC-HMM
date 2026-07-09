#calculates linear value of risk due to sociodemo covar, does not incl mixture
SurvCovarRiskVec <- function(surv_covar,surv_coef){
  surv_covar[[1]] <- matrix(surv_covar[[1]],ncol = 1)
  surv_coef[[1]] <- matrix(surv_coef[[1]],ncol = 1)
  surv_covar_risk_vec <- rowSums(mapply(function(x, y) x %*% y,  surv_covar,surv_coef, SIMPLIFY = T))
  return(surv_covar_risk_vec)
}

#calculates non-parametric baseline haz using breslow estimator
CalcBLHaz <- function(surv_coef,beta_vec, re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar){
  n <- length(surv_event)
  if (n == 0){
    return(list(numeric(0),numeric(0)))
  }

  if (!is.null(dim(re_prob)) && ncol(re_prob) != 1){
    risk_score <- drop(re_prob %*% exp(beta_vec)) * exp(surv_covar_risk_vec)
  } else {
    risk_score <- exp(surv_covar_risk_vec)
  }

  time_factor <- factor(surv_time,levels = sort(unique(surv_time)))
  risk_by_time <- as.numeric(tapply(risk_score,time_factor,sum))
  denom_by_time <- rev(cumsum(rev(risk_by_time)))
  denom_by_time[denom_by_time > .Machine$double.xmax] <-
    .Machine$double.xmax

  bline_vec <- unname(surv_event / denom_by_time[as.integer(time_factor)])
  bline_by_time <- as.numeric(tapply(bline_vec,time_factor,sum))
  cbline_by_time <- cumsum(bline_by_time)
  cbline_vec <- unname(cbline_by_time[as.integer(time_factor)])

  return(list(bline_vec,cbline_vec))
}

#Calculates likelihood only due to survival component
SurvLike <- function(beta_vec,surv_covar_risk_vec,surv_coef,survival_context){
  re_prob <- survival_context$re_prob
  surv_event <- survival_context$surv_event
  surv_time <- survival_context$surv_time
  surv_covar <- survival_context$surv_covar
  bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,
                        surv_event,surv_time,surv_covar)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]


  # loglike <- sum(log(bline_vec^surv_event) +
  #                  ((re_prob %*% beta_vec)+surv_covar_risk_vec) * surv_event -
  #                  cbline_vec*exp((re_prob %*% beta_vec)+surv_covar_risk_vec))

  #keeping above lines as found this bug after most of runs were done
  #exponentiating the probabilities in first, should only do that to betas in last line
  class_event_lp <- drop(re_prob %*% beta_vec)
  class_risk <- drop(re_prob %*% exp(beta_vec))

  loglike <- sum(
    log(bline_vec^surv_event) +
      (class_event_lp + surv_covar_risk_vec) * surv_event -
      cbline_vec * exp(surv_covar_risk_vec) * class_risk
  )


  return(-loglike)
}

#converts long vector of survival coef into beta and sociodemo list
OutofBetaSurvCoef <- function(beta_surv_coef,surv_coef_len,fit_mix_num){
  if (fit_mix_num > 1){
    beta_vec <- c(0,beta_surv_coef[2:fit_mix_num])
  } else {
    beta_vec <- 0
  }

  surv_coef_new <- list()
  surv_coef_new <- append(surv_coef_new,beta_surv_coef[[1]])


  surv_coef_len_alt <- cumsum(surv_coef_len-1)


  if (length(surv_coef_len) > 1){
    for (i in 1:(length(surv_coef_len)-1)) {
      coef_vec <- c(0,beta_surv_coef[(fit_mix_num+1+surv_coef_len_alt[i]):(fit_mix_num+surv_coef_len_alt[i+1])])
      surv_coef_new <- append(surv_coef_new,list(coef_vec))
    }
  }

  return(list(beta_vec,surv_coef_new))
}

#input is beta vector and sociodemo coef list, output is long combined vector
#age first
IntoBetaSurvCoef <- function(beta_vec,surv_coef,fit_mix_num){
  beta_surv_coef_len <- fit_mix_num + length(unlist(surv_coef)) - length(surv_coef)
  beta_surv_coef <- numeric(beta_surv_coef_len)

  beta_surv_coef[1] <- surv_coef[[1]]
  if(fit_mix_num > 1){
    beta_surv_coef[2:fit_mix_num] <- beta_vec[-1]
  }


  altered_surv_coef <- surv_coef[-1]
  if (length(altered_surv_coef) > 0){
    for (i in 1:length(altered_surv_coef)){
      altered_surv_coef[[i]] <- altered_surv_coef[[i]][-1]
    }

    beta_surv_coef[(fit_mix_num+1):beta_surv_coef_len] <- unlist(altered_surv_coef)
  }

  return(beta_surv_coef)
}

# Calculates the profiled Cox/Breslow survival score and positive information
# at fixed finite-dimensional survival parameters.
CalcSurvivalScoreInfo <- function(beta_surv_coef,survival_context,
                                  surv_coef_len,fit_mix_num,
                                  surv_covar_risk_vec = NULL){
  re_prob <- survival_context$re_prob
  surv_event <- survival_context$surv_event
  surv_time <- survival_context$surv_time
  surv_covar <- survival_context$surv_covar

  beta_surv_coef_list <- OutofBetaSurvCoef(beta_surv_coef,surv_coef_len,
                                           fit_mix_num)
  beta_vec <- beta_surv_coef_list[[1]]
  surv_coef <- beta_surv_coef_list[[2]]

  if (is.null(surv_covar_risk_vec)){
    surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
  }

  param_num <- fit_mix_num + length(unlist(surv_coef)) - length(surv_coef)
  grad <- numeric(param_num)
  hess <- matrix(0,param_num,param_num)

  if (is.null(dim(re_prob))){
    re_prob <- matrix(re_prob,ncol = fit_mix_num)
  }
  if (ncol(re_prob) != fit_mix_num){
    stop("survival_context$re_prob must have fit_mix_num columns")
  }

  build_event_z <- function(ind){
    z <- numeric(param_num)
    z[1] <- surv_covar[[1]][ind]

    if (fit_mix_num > 1){
      z[2:fit_mix_num] <- re_prob[ind,2:fit_mix_num]
    }

    if (length(surv_covar) > 1){
      pos <- fit_mix_num + 1
      for (surv_covar_ind in 2:length(surv_covar)){
        if (length(surv_coef[[surv_covar_ind]]) <= 1){
          next
        }
        covar_mat <- as.matrix(surv_covar[[surv_covar_ind]])
        non_ref_cols <- 2:length(surv_coef[[surv_covar_ind]])
        covar_values <- as.numeric(covar_mat[ind,non_ref_cols,drop = FALSE])
        idx <- pos:(pos + length(covar_values) - 1)
        z[idx] <- covar_values
        pos <- pos + length(covar_values)
      }
    }

    z
  }

  build_risk_set_z <- function(risk_set){
    risk_n <- sum(risk_set)
    subject_ind <- rep(seq_len(risk_n),each = fit_mix_num)
    class_ind <- rep(seq_len(fit_mix_num),times = risk_n)
    z <- matrix(0,nrow = risk_n * fit_mix_num,ncol = param_num)

    z[,1] <- surv_covar[[1]][risk_set][subject_ind]

    if (fit_mix_num > 1){
      for (re_ind in 2:fit_mix_num){
        z[,re_ind] <- as.numeric(class_ind == re_ind)
      }
    }

    if (length(surv_covar) > 1){
      pos <- fit_mix_num + 1
      for (surv_covar_ind in 2:length(surv_covar)){
        if (length(surv_coef[[surv_covar_ind]]) <= 1){
          next
        }
        covar_mat <- as.matrix(surv_covar[[surv_covar_ind]])
        non_ref_cols <- 2:length(surv_coef[[surv_covar_ind]])
        covar_block <- covar_mat[risk_set,non_ref_cols,drop = FALSE]
        idx <- pos:(pos + length(non_ref_cols) - 1)
        z[,idx] <- covar_block[subject_ind,,drop = FALSE]
        pos <- pos + length(non_ref_cols)
      }
    }

    z
  }

  for (ind in which(surv_event == 1)){
    risk_set <- surv_time >= surv_time[ind]
    linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]

    risk_weight_mat <- re_prob[risk_set,,drop = FALSE] *
      matrix(exp(linear_surv_covar_risk),nrow = sum(risk_set),
             ncol = fit_mix_num) *
      matrix(exp(beta_vec),nrow = sum(risk_set),ncol = fit_mix_num,
             byrow = TRUE)
    risk_n <- nrow(risk_weight_mat)
    subject_ind <- rep(seq_len(risk_n),each = fit_mix_num)
    class_ind <- rep(seq_len(fit_mix_num),times = risk_n)
    risk_weights <- risk_weight_mat[cbind(subject_ind,class_ind)]
    denom <- sum(risk_weights)

    if (!is.finite(denom) || denom <= 0){
      stop("Non-positive or non-finite survival risk-set denominator")
    }

    z <- build_risk_set_z(risk_set)
    risk_prob <- risk_weights / denom
    ez <- colSums(sweep(z,1,risk_prob,`*`))
    ezz <- crossprod(z,sweep(z,1,risk_prob,`*`))

    grad <- grad + build_event_z(ind) - ez
    hess <- hess + ezz - tcrossprod(ez)
  }

  hess <- 0.5 * (hess + t(hess))

  list(score = grad,
       info = hess,
       beta_vec = beta_vec,
       surv_coef = surv_coef,
       surv_covar_risk_vec = surv_covar_risk_vec)
}

#LM approach for calculating survival coefficients
CalcBetaManual <- function(beta_surv_coef,surv_covar_risk_vec,stop_crit,
                           survival_context,surv_coef_len,fit_mix_num){
  l2norm <- 101
  hess <- NULL
  while (l2norm > stop_crit){

    beta_surv_coef_list <- OutofBetaSurvCoef(
      beta_surv_coef,
      surv_coef_len,
      fit_mix_num
    )

    beta_vec <- beta_surv_coef_list[[1]]
    surv_coef <- beta_surv_coef_list[[2]]
    surv_covar_risk_vec <- SurvCovarRiskVec(
      survival_context$surv_covar,
      surv_coef
    )
    score_info <- CalcSurvivalScoreInfo(beta_surv_coef,
                                        survival_context,
                                        surv_coef_len,
                                        fit_mix_num,
                                        surv_covar_risk_vec)
    beta_vec <- score_info$beta_vec
    surv_coef <- score_info$surv_coef
    grad <- score_info$score
    hess <- score_info$info

    old_slike <- SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,
                          survival_context)
    slike_diff <- -1

    step_size <- .01
    step_mat <- matrix(0,dim(hess)[1],dim(hess)[2])
    nr_fact <- numeric(length(grad))-1
    max_val <- 6
    while(slike_diff < 0|all(nr_fact == -1) | max_val > 5){
      diag(step_mat) <- (diag(step_mat) + step_size - .01)
      nr_fact <- SolveCatch(hess + step_mat,-grad)
      step_size <- step_size * 10

      beta_surv_coef_new <- beta_surv_coef-nr_fact
      beta_vec_new <- OutofBetaSurvCoef(beta_surv_coef_new,
                                        surv_coef_len,fit_mix_num)[[1]]
      surv_coef_new <- OutofBetaSurvCoef(beta_surv_coef_new,
                                         surv_coef_len,fit_mix_num)[[2]]
      max_val <- abs(max(c(beta_vec_new,unlist(surv_coef_new))))

      #removed this and added following lines because was having an issue where old baselines were used
      #shouldnt change too much hopefully
      # new_slike <- SurvLike(beta_vec_new,surv_covar_risk_vec,surv_coef_new,
      #                       survival_context)

      surv_covar_risk_vec_new <- SurvCovarRiskVec(
        survival_context$surv_covar,
        surv_coef_new
      )

      new_slike <- SurvLike(
        beta_vec_new,
        surv_covar_risk_vec_new,
        surv_coef_new,
        survival_context
      )
      slike_diff <- old_slike - new_slike
    }

    l2norm <- sum(sqrt((beta_surv_coef_new - beta_surv_coef)^2))
    beta_surv_coef <- beta_surv_coef_new
  }

  list(beta_surv_coef,sqrt(diag(solve(hess))))
}

RemFirCol <- function(x){return(x[,-1])}

#Calculates beta, manual LM for JM, standard Cox for 2 stage
CalcBeta <- function(beta_surv_coef, combined_covar_mat,surv_covar_risk_vec,
                     incl_surv, survival_context, surv_coef_len, fit_mix_num,
                     stop_crit = .1){

  if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
    stop_crit <- JOINT_BETA_STOP_CRIT
    beta_surv_coef_new <- CalcBetaManual(beta_surv_coef, surv_covar_risk_vec,
                                         stop_crit, survival_context,
                                         surv_coef_len, fit_mix_num)
    return(beta_surv_coef_new)
  }

  surv_data <- data.frame(time = survival_context$surv_time,
                          status = survival_context$surv_event,
                          age = survival_context$surv_covar[[1]])

  surv_data <- cbind(surv_data,survival_context$re_prob,combined_covar_mat)
  colnames(surv_data)[4] <- "toRem"

  fit <- coxph(Surv(time, status) ~ .  - toRem, data = surv_data)
  beta_surv_coef_new <- fit$coefficients
  se <- sqrt(diag(fit$var))

  return(list(beta_surv_coef_new,se))


}

CalcS <- function(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec){
  surv_mat_ind <- matrix(NA,length(event_time),dim(re_prob)[1])
  surv_vec <- numeric(length(event_time))
  fitted_mix_num <- ncol(re_prob)

  for (ind in 1:dim(re_prob)[1]){
    for (t in 1:length(event_time)){
      rs <- 0


      for (re_ind in 1:fitted_mix_num){
        rs <- rs + exp(-cbline_vec_new[t]* exp(beta_vec[re_ind]+surv_covar_risk_vec[ind])) * re_prob[ind,re_ind]
      }

      surv_vec[t] <- rs
    }
    surv_mat_ind[,ind] <- surv_vec
  }

  return(surv_mat_ind)
}

SurvCovar2Coef <- function(covar_mat,max_val = .1){
  return(seq(0,max_val,length.out = dim(covar_mat)[2]))
}

SubsetSurvCovar <- function(surv_covar,leave_out_inds){
  for (i in 1:length(surv_covar)){
    if (is.null(dim(surv_covar[[i]]))){
      surv_covar[[i]] <- surv_covar[[i]][leave_out_inds]
    } else {
      surv_covar[[i]] <- surv_covar[[i]][leave_out_inds,]
    }
  }
  return(surv_covar)
}

