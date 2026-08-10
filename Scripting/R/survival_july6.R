#calculates linear value of risk due to sociodemo covar, does not incl mixture
SurvCovarRiskVec <- function(surv_covar,surv_coef){
  surv_covar[[1]] <- matrix(surv_covar[[1]],ncol = 1)
  surv_coef[[1]] <- matrix(surv_coef[[1]],ncol = 1)
  surv_covar_risk_vec <- rowSums(mapply(function(x, y) x %*% y,  surv_covar,surv_coef, SIMPLIFY = T))
  return(surv_covar_risk_vec)
}

# The two-stage Cox model treats posterior class probabilities as ordinary
# covariates.  Use the probabilities exactly as they were supplied to
# FitTwoStageCox(); silently renormalizing them here could make prediction
# differ from the fitted Cox design matrix.
CalcTwoStageLinearPredictor <- function(beta_vec,re_prob,
                                        surv_covar_risk_vec){
  re_prob <- as.matrix(re_prob)
  num_people <- nrow(re_prob)

  if (ncol(re_prob) != length(beta_vec)){
    stop("ncol(re_prob) must equal length(beta_vec)")
  }
  if (any(!is.finite(beta_vec))){
    stop("beta_vec contains nonfinite values")
  }
  if (length(surv_covar_risk_vec) != num_people){
    stop("re_prob and surv_covar_risk_vec must describe the same participants")
  }
  if (any(!is.finite(re_prob)) || any(re_prob < 0)){
    stop("re_prob must contain finite, nonnegative probabilities")
  }
  if (any(rowSums(re_prob) <= 0)){
    stop("Every row of re_prob must have positive probability mass")
  }
  if (any(!is.finite(surv_covar_risk_vec))){
    stop("surv_covar_risk_vec contains nonfinite values")
  }

  drop(re_prob %*% beta_vec) + surv_covar_risk_vec
}

# Generic Breslow baseline-hazard calculation for an already determined Cox
# linear predictor.  This is used only for the post-fit two-stage model; the
# established joint-model CalcBLHaz() path remains unchanged below.
CalcBLHazFromLinearPredictor <- function(linear_predictor,surv_event,
                                         surv_time,sweights_vec){
  num_people <- length(surv_event)

  if (length(surv_time) != num_people ||
      length(linear_predictor) != num_people ||
      length(sweights_vec) != num_people){
    stop("Survival vectors, linear predictor, and weights must have equal lengths")
  }
  if (num_people == 0L){
    return(list(bline_vec = numeric(0),cbline_vec = numeric(0)))
  }
  if (any(!surv_event %in% c(0,1))){
    stop("surv_event must contain only 0 and 1")
  }
  if (any(!is.finite(surv_time))){
    stop("surv_time contains nonfinite values")
  }
  if (any(!is.finite(linear_predictor))){
    stop("linear_predictor contains nonfinite values")
  }
  if (any(!is.finite(sweights_vec)) ||
      any(sweights_vec < 0) ||
      sum(sweights_vec) <= 0){
    stop("Invalid sample weights")
  }

  positive_weight <- sweights_vec > 0
  linear_predictor_offset <- max(linear_predictor[positive_weight])
  weighted_relative_risk <- sweights_vec *
    exp(linear_predictor - linear_predictor_offset)

  time_factor <- factor(surv_time,levels = sort(unique(surv_time)))
  risk_by_time <- as.numeric(
    tapply(weighted_relative_risk,time_factor,sum)
  )
  denominator_by_time <- rev(cumsum(rev(risk_by_time)))
  event_weight_by_time <- as.numeric(
    tapply(surv_event * sweights_vec,time_factor,sum)
  )

  if (any(!is.finite(denominator_by_time)) ||
      any(denominator_by_time <= 0)){
    stop("Two-stage baseline-hazard denominator is nonpositive or nonfinite")
  }

  log_denominator_by_time <-
    linear_predictor_offset + log(denominator_by_time)
  bline_by_time <- numeric(length(event_weight_by_time))
  event_times <- event_weight_by_time > 0
  bline_by_time[event_times] <- exp(
    log(event_weight_by_time[event_times]) -
      log_denominator_by_time[event_times]
  )
  cbline_by_time <- cumsum(bline_by_time)

  list(
    bline_vec = unname(bline_by_time[as.integer(time_factor)]),
    cbline_vec = unname(cbline_by_time[as.integer(time_factor)])
  )
}

CalcBLHazTwoStage <- function(beta_vec,re_prob,surv_covar_risk_vec,
                              surv_event,surv_time,sweights_vec){
  linear_predictor <- CalcTwoStageLinearPredictor(
    beta_vec = beta_vec,
    re_prob = re_prob,
    surv_covar_risk_vec = surv_covar_risk_vec
  )

  CalcBLHazFromLinearPredictor(
    linear_predictor = linear_predictor,
    surv_event = surv_event,
    surv_time = surv_time,
    sweights_vec = sweights_vec
  )
}

#calculates non-parametric baseline haz using breslow estimator
CalcBLHaz <- function(surv_coef,beta_vec, re_prob,surv_covar_risk_vec,surv_event,surv_time,surv_covar, sweights_vec){
  n <- length(surv_event)
  if (n == 0){
    return(list(numeric(0),numeric(0)))
  }

  if (!is.null(dim(re_prob)) && ncol(re_prob) != 1){
    risk_score <- sweights_vec * drop(re_prob %*% exp(beta_vec)) * exp(surv_covar_risk_vec)
  } else {
    risk_score <- sweights_vec * exp(surv_covar_risk_vec)
  }

  time_factor <- factor(surv_time,levels = sort(unique(surv_time)))
  risk_by_time <- as.numeric(tapply(risk_score,time_factor,sum))
  denom_by_time <- rev(cumsum(rev(risk_by_time)))
  denom_by_time[denom_by_time > .Machine$double.xmax] <-
    .Machine$double.xmax


  event_weight_by_time <- as.numeric(
    tapply(surv_event * sweights_vec, time_factor, sum)
  )

  bline_by_time <- event_weight_by_time / denom_by_time
  cbline_by_time <- cumsum(bline_by_time)

  bline_vec <- unname(
    bline_by_time[as.integer(time_factor)]
  )

  cbline_vec <- unname(
    cbline_by_time[as.integer(time_factor)]
  )

  return(list(bline_vec,cbline_vec))
}

#Calculates likelihood only due to survival component
SurvLike <- function(beta_vec,surv_covar_risk_vec,surv_coef,survival_context){
  re_prob <- survival_context$re_prob
  surv_event <- survival_context$surv_event
  surv_time <- survival_context$surv_time
  surv_covar <- survival_context$surv_covar
  sweights_vec <- survival_context$sweights_vec

  bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,surv_covar_risk_vec,
                        surv_event,surv_time,surv_covar,sweights_vec)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]

  # July 6 behavior: posterior-average the class log-risk first,
  # then exponentiate in the cumulative-hazard term.
  class_lp <- drop(re_prob %*% beta_vec)

  # Safe equivalent of surv_event * log(bline_vec): do not evaluate
  # log(0) for censored observations.
  event_log_bline <- numeric(length(surv_event))
  event_rows <- surv_event == 1
  event_log_bline[event_rows] <- log(bline_vec[event_rows])

  loglike_i <-
    event_log_bline +
    surv_event * (class_lp + surv_covar_risk_vec) -
    cbline_vec * exp(class_lp + surv_covar_risk_vec)

  loglike <- sum(sweights_vec * loglike_i)

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
  sweights_vec <- survival_context$sweights_vec

  beta_surv_coef_list <- OutofBetaSurvCoef(beta_surv_coef,surv_coef_len,
                                           fit_mix_num)
  beta_vec <- beta_surv_coef_list[[1]]
  surv_coef <- beta_surv_coef_list[[2]]

  if (is.null(surv_covar_risk_vec)){
    surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
  }

  if (is.null(dim(re_prob))){
    re_prob <- matrix(re_prob,ncol = fit_mix_num)
  }
  if (ncol(re_prob) != fit_mix_num){
    stop("survival_context$re_prob must have fit_mix_num columns")
  }

  param_num <- fit_mix_num + length(unlist(surv_coef)) - length(surv_coef)
  grad <- numeric(param_num)
  hess <- matrix(0,param_num,param_num)

  # July 6 score/information implementation, with current sample-weight plumbing.
  # When sweights_vec == 1 (all simulations), this is numerically the July 6 code.
  for (ind in which(surv_event == 1)){
    event_weight <- sweights_vec[ind]
    risk_set <- surv_time >= surv_time[ind]
    linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]
    ageadj_risk <- exp(linear_surv_covar_risk) %x% t(exp(beta_vec))

    weighted_re_prob <- sweep(re_prob[risk_set,,drop = FALSE],
                              1,sweights_vec[risk_set],"*")

    num_list <- list()

    num0 <- sum(weighted_re_prob * ageadj_risk *
                  surv_covar[[1]][risk_set] %x%
                  t(numeric(fit_mix_num) + 1))
    num02 <- sum(weighted_re_prob * ageadj_risk *
                   surv_covar[[1]][risk_set]^2 %x%
                   t(numeric(fit_mix_num) + 1))
    num_list[[1]] <- c(num0,num02)

    denom <- sum(weighted_re_prob * ageadj_risk)

    grad[1] <- grad[1] + event_weight *
      (surv_covar[[1]][ind] - num0 / denom)
    hess[1,1] <- hess[1,1] + event_weight *
      ((num02 / denom) - (num0 / denom)^2)

    if (fit_mix_num > 1){
      for (beta_ind in 2:fit_mix_num){
        num <- sum(weighted_re_prob[,beta_ind] *
                     exp(beta_vec[beta_ind]) *
                     exp(linear_surv_covar_risk))
        grad[beta_ind] <- grad[beta_ind] + event_weight *
          (re_prob[ind,beta_ind] - num / denom)

        hess[beta_ind,beta_ind] <- hess[beta_ind,beta_ind] + event_weight *
          ((num / denom) - (num / denom)^2)

        num_cross <- sum(weighted_re_prob[,beta_ind] *
                           exp(beta_vec[beta_ind]) *
                           exp(linear_surv_covar_risk) *
                           surv_covar[[1]][risk_set])
        hess[1,beta_ind] <- hess[1,beta_ind] + event_weight *
          (denom * num_cross - num * num0) / denom^2
        hess[beta_ind,1] <- hess[beta_ind,1] + event_weight *
          (denom * num_cross - num * num0) / denom^2
      }
    }

    if (length(surv_covar) > 1){
      for (surv_covar_ind in 2:length(surv_covar)){
        num_disc_covar_vec <- c()
        for (surv_covar_indicator in 2:length(surv_coef[[surv_covar_ind]])){
          num_disc_covar <- sum(weighted_re_prob * ageadj_risk *
                                  surv_covar[[surv_covar_ind]][risk_set,
                                                                surv_covar_indicator] %x%
                                  t(numeric(fit_mix_num) + 1))

          num_disc_covar_vec <- c(num_disc_covar_vec,num_disc_covar)
          num_list[[surv_covar_ind]] <- num_disc_covar_vec
        }
      }

      list_of_lens <- unlist(lapply(num_list,length))
      list_of_lens[1] <- 0
      list_of_cum_lens <- cumsum(list_of_lens)

      for (surv_covar_ind in 2:length(surv_covar)){
        starting_index <- fit_mix_num + list_of_cum_lens[surv_covar_ind - 1] + 1
        ending_index <- fit_mix_num + list_of_cum_lens[surv_covar_ind]
        curr_covar_inds <- starting_index:ending_index

        grad[curr_covar_inds] <- grad[curr_covar_inds] + event_weight *
          (surv_covar[[surv_covar_ind]][ind,-1] -
             num_list[[surv_covar_ind]] / denom)

        diag(hess)[curr_covar_inds] <- diag(hess)[curr_covar_inds] +
          event_weight *
          (num_list[[surv_covar_ind]] / denom -
             (num_list[[surv_covar_ind]] / denom)^2)

        for (ind_curr_covar_inds in seq_along(curr_covar_inds)){
          num_cross_age <- sum(weighted_re_prob * ageadj_risk *
                                 surv_covar[[1]][risk_set] *
                                 surv_covar[[surv_covar_ind]][risk_set,
                                                               ind_curr_covar_inds + 1] %x%
                                 t(numeric(fit_mix_num) + 1))
          hess[1,curr_covar_inds[ind_curr_covar_inds]] <-
            hess[1,curr_covar_inds[ind_curr_covar_inds]] + event_weight *
            (denom * num_cross_age -
               num_list[[surv_covar_ind]][[ind_curr_covar_inds]] * num0) /
            denom^2
          hess[curr_covar_inds[ind_curr_covar_inds],1] <-
            hess[curr_covar_inds[ind_curr_covar_inds],1] + event_weight *
            (denom * num_cross_age -
               num_list[[surv_covar_ind]][[ind_curr_covar_inds]] * num0) /
            denom^2
        }

        if (fit_mix_num > 1){
          for (beta_ind in 2:fit_mix_num){
            num <- sum(weighted_re_prob[,beta_ind] *
                         exp(beta_vec[beta_ind]) *
                         exp(linear_surv_covar_risk))

            for (ind_curr_covar_inds in seq_along(curr_covar_inds)){
              num_cross <- sum(weighted_re_prob[,beta_ind] *
                                 exp(beta_vec[beta_ind]) *
                                 exp(linear_surv_covar_risk) *
                                 surv_covar[[surv_covar_ind]][risk_set,
                                                               ind_curr_covar_inds + 1])
              hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] <-
                hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] + event_weight *
                (denom * num_cross -
                   num * num_list[[surv_covar_ind]][[ind_curr_covar_inds]]) /
                denom^2

              hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] <-
                hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] + event_weight *
                (denom * num_cross -
                   num * num_list[[surv_covar_ind]][[ind_curr_covar_inds]]) /
                denom^2
            }
          }
        }

        if (surv_covar_ind != length(surv_covar)){
          for (surv_covar_ind2 in (surv_covar_ind + 1):length(surv_covar)){
            starting_index2 <- fit_mix_num +
              list_of_cum_lens[surv_covar_ind2 - 1] + 1
            ending_index2 <- fit_mix_num + list_of_cum_lens[surv_covar_ind2]
            curr_covar_inds2 <- starting_index2:ending_index2

            for (ind_curr_covar_inds in seq_along(curr_covar_inds)){
              for (ind_curr_covar_inds2 in seq_along(curr_covar_inds2)){
                num_cross_covar <- sum(weighted_re_prob * ageadj_risk *
                                         surv_covar[[surv_covar_ind]][risk_set,
                                                                      ind_curr_covar_inds + 1] *
                                         surv_covar[[surv_covar_ind2]][risk_set,
                                                                       ind_curr_covar_inds2 + 1] %x%
                                         t(numeric(fit_mix_num) + 1))

                hess[curr_covar_inds2[ind_curr_covar_inds2],
                     curr_covar_inds[ind_curr_covar_inds]] <-
                  hess[curr_covar_inds2[ind_curr_covar_inds2],
                       curr_covar_inds[ind_curr_covar_inds]] + event_weight *
                  (denom * num_cross_covar -
                     num_list[[surv_covar_ind]][[ind_curr_covar_inds]] *
                     num_list[[surv_covar_ind2]][[ind_curr_covar_inds2]]) /
                  denom^2

                hess[curr_covar_inds[ind_curr_covar_inds],
                     curr_covar_inds2[ind_curr_covar_inds2]] <-
                  hess[curr_covar_inds[ind_curr_covar_inds],
                       curr_covar_inds2[ind_curr_covar_inds2]] + event_weight *
                  (denom * num_cross_covar -
                     num_list[[surv_covar_ind]][[ind_curr_covar_inds]] *
                     num_list[[surv_covar_ind2]][[ind_curr_covar_inds2]]) /
                  denom^2
              }
            }
          }
        }
      }
    }
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

    # July 6 behavior: use the surv_covar_risk_vec supplied to CalcBetaManual
    # throughout this internal optimization rather than recomputing it for
    # the current or proposed survival-covariate coefficients.
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

      new_slike <- SurvLike(beta_vec_new,surv_covar_risk_vec,surv_coef_new,
                            survival_context)
      slike_diff <- old_slike - new_slike
    }

    l2norm <- sum(sqrt((beta_surv_coef_new - beta_surv_coef)^2))
    beta_surv_coef <- beta_surv_coef_new
  }

  list(beta_surv_coef,sqrt(diag(solve(hess))))
}

RemFirCol <- function(x){return(x[,-1])}

FitTwoStageCox <- function(re_prob,survival_context,combined_covar_mat,
                           ties = "breslow",init = NULL,
                           expected_coef_names = NULL){
  if (!ties %in% c("breslow","efron","exact")){
    stop("ties must be one of breslow, efron, or exact")
  }

  re_prob <- as.matrix(re_prob)
  n <- length(survival_context$surv_time)
  if (!all(dim(re_prob) == c(n,ncol(re_prob)))){
    stop("re_prob must have one row per survival observation")
  }
  if (length(survival_context$surv_event) != n ||
      length(survival_context$surv_covar[[1]]) != n ||
      length(survival_context$sweights_vec) != n){
    stop("survival_context entries must have one row per survival observation")
  }

  surv_data <- data.frame(time = survival_context$surv_time,
                          status = survival_context$surv_event,
                          age = survival_context$surv_covar[[1]])

  surv_data <- cbind(surv_data,re_prob,combined_covar_mat)
  if (ncol(surv_data) < 4){
    stop("two-stage Cox data must include at least one class probability")
  }
  colnames(surv_data)[4] <- "toRem"

  if (nrow(surv_data) != n){
    stop("two-stage Cox data row count changed during construction")
  }

  if (!is.null(init)){
    init <- as.numeric(init)
  }

  cox_args <- list(
    formula = survival::Surv(time,status) ~ . - toRem,
    data = surv_data,
    weights = survival_context$sweights_vec,
    ties = ties,
    x = TRUE,
    y = TRUE,
    model = TRUE,
    robust = FALSE
  )
  if (!is.null(init)){
    cox_args$init <- init
  }

  fit <- do.call(survival::coxph,cox_args)

  if (!is.null(fit$na.action) && length(fit$na.action) > 0){
    stop("FitTwoStageCox omitted rows; missing-data handling is unsupported")
  }
  if (is.null(fit$model) || nrow(fit$model) != n){
    stop("FitTwoStageCox model frame does not match the input row count")
  }

  fit_coef <- stats::coef(fit)
  if (any(is.na(fit_coef)) || any(!is.finite(fit_coef))){
    stop("FitTwoStageCox returned aliased, missing, or nonfinite coefficients")
  }

  if (!is.null(expected_coef_names) &&
      !identical(names(fit_coef),expected_coef_names)){
    stop("FitTwoStageCox coefficient order does not match the expected order")
  }

  fit
}

#Calculates beta, manual LM for JM, standard Cox for 2 stage
CalcBeta <- function(beta_surv_coef, combined_covar_mat,surv_covar_risk_vec,
                     incl_surv, survival_context, surv_coef_len, fit_mix_num,
                     stop_crit = .1){

  if (incl_surv == MODEL_TYPE_CODES[["joint"]]) {


    if (incl_surv == MODEL_TYPE_CODES[["joint"]]) {
      beta_result <- CalcBetaManual(
        beta_surv_coef = beta_surv_coef,
        surv_covar_risk_vec = surv_covar_risk_vec,
        stop_crit = JOINT_BETA_STOP_CRIT,
        survival_context = survival_context,
        surv_coef_len = surv_coef_len,
        fit_mix_num = fit_mix_num
      )

      return(beta_result)
    }


    # #used to use calcbetamanual but this is equivalent  (maybe not so equivalent after all)
    # #safer to as it uses coxph
    # coxph_result <- CalcBetaCoxphJoint(
    #   beta_surv_coef = beta_surv_coef,
    #   survival_context = survival_context,
    #   surv_coef_len = surv_coef_len,
    #   fit_mix_num = fit_mix_num,
    #   one_step = FALSE
    # )

    # beta_surv_coef_new <-
    #   coxph_result$beta_surv_coef

    # conditional_se <-
    #   coxph_result$conditional_se

    # expected_parameter_num <-
    #   length(beta_surv_coef)

    # if (
    #   length(beta_surv_coef_new) !=
    #     expected_parameter_num
    # ) {
    #   stop(
    #     paste(
    #       "CalcBetaCoxphJoint returned",
    #       length(beta_surv_coef_new),
    #       "coefficients, but",
    #       expected_parameter_num,
    #       "were expected"
    #     )
    #   )
    # }

    # if (any(!is.finite(beta_surv_coef_new))) {
    #   stop(
    #     "CalcBetaCoxphJoint returned nonfinite coefficients"
    #   )
    # }

    # if (
    #   length(conditional_se) !=
    #     expected_parameter_num
    # ) {
    #   stop(
    #     paste(
    #       "CalcBetaCoxphJoint returned",
    #       length(conditional_se),
    #       "standard errors, but",
    #       expected_parameter_num,
    #       "were expected"
    #     )
    #   )
    # }

    # return(
    #   list(
    #     beta_surv_coef_new,
    #     conditional_se
    #   )
    # )
  }

  fit <- FitTwoStageCox(re_prob = survival_context$re_prob,
                        survival_context = survival_context,
                        combined_covar_mat = combined_covar_mat,
                        ties = "breslow")
  beta_surv_coef_new <- stats::coef(fit)
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

CalcSTwoStage <- function(event_time,cbline_vec_new,beta_vec,re_prob,
                          surv_covar_risk_vec){
  if (length(event_time) != length(cbline_vec_new)){
    stop("event_time and cbline_vec_new must have equal lengths")
  }
  if (any(!is.finite(cbline_vec_new)) || any(cbline_vec_new < 0)){
    stop("cbline_vec_new must contain finite, nonnegative values")
  }

  linear_predictor <- CalcTwoStageLinearPredictor(
    beta_vec = beta_vec,
    re_prob = re_prob,
    surv_covar_risk_vec = surv_covar_risk_vec
  )
  relative_risk <- exp(linear_predictor)

  if (any(!is.finite(relative_risk))){
    stop("Two-stage relative-risk calculation produced nonfinite values")
  }

  outer(
    cbline_vec_new,
    relative_risk,
    function(baseline_hazard,risk){
      exp(-baseline_hazard * risk)
    }
  )
}

CalcSurvivalProbabilities <- function(model_type,event_time,cbline_vec_new,
                                      beta_vec,re_prob,
                                      surv_covar_risk_vec){
  model_type <- match.arg(model_type,c("joint","two_stage"))

  if (model_type == "joint"){
    return(CalcS(event_time,cbline_vec_new,beta_vec,re_prob,
                 surv_covar_risk_vec))
  }

  CalcSTwoStage(event_time,cbline_vec_new,beta_vec,re_prob,
                surv_covar_risk_vec)
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

#need to validate this against calcbetamanual
#this uses coxph funciton and avoids doing a newton-raphson by hand (probably more accurate)
CalcBetaCoxphJoint <- function(
    beta_surv_coef,
    survival_context,
    surv_coef_len,
    fit_mix_num,
    one_step = FALSE) {

  re_prob <- as.matrix(survival_context$re_prob)
  surv_time <- survival_context$surv_time
  surv_event <- survival_context$surv_event
  sweights_vec <- survival_context$sweights_vec

  n <- length(surv_time)

  if (!all(dim(re_prob) == c(n, fit_mix_num))) {
    stop("re_prob dimensions do not match the survival data")
  }

  subject_design <- BuildSurvivalCovariateMatrix(
    survival_context$surv_covar
  )

  subject_index <- rep(seq_len(n), each = fit_mix_num)
  class_index <- rep(seq_len(fit_mix_num), times = n)

  expanded_data <- data.frame(
    time = surv_time[subject_index],
    status = surv_event[subject_index],
    latent_class = factor(
      class_index,
      levels = seq_len(fit_mix_num)
    ),
    subject_design[subject_index, , drop = FALSE],
    check.names = FALSE
  )

  # as.vector(t(re_prob)) gives:
  # subject 1, classes 1:G; subject 2, classes 1:G; ...
  expanded_weights <-
    sweights_vec[subject_index] *
    as.vector(t(re_prob))

  # Validate before passing weights to coxph().
  if (any(!is.finite(expanded_weights))) {
    stop("Expanded Cox weights contain nonfinite values")
  }

  if (any(expanded_weights < 0)) {
    stop("Expanded Cox weights contain negative values")
  }

  # coxph() requires strictly positive weights.
  #
  # Rows with zero posterior weight make exactly zero
  # contribution to the weighted partial likelihood, score,
  # and information, so dropping them is exact.
  positive_weight <-
    expanded_weights > 0

  if (!any(positive_weight)) {
    stop("All expanded Cox weights are zero")
  }

  expanded_data <-
    expanded_data[
      positive_weight,
      ,
      drop = FALSE
    ]

  expanded_weights <-
    expanded_weights[
      positive_weight
    ]
  



  # Order terms to match IntoBetaSurvCoef():
  # age, class 2:G, remaining survival covariates
  other_covariates <- setdiff(
    colnames(subject_design),
    "age"
  )

  formula_terms <- c(
    "age",
    if (fit_mix_num > 1L) {
      "latent_class"
    },
    other_covariates
  )

  survival_formula <- reformulate(
    formula_terms,
    response = "Surv(time, status)"
  )

  control <- if (one_step) {
    coxph.control(iter.max = 1)
  } else {
    coxph.control()
  }

  fit <- coxph(
    formula = survival_formula,
    data = expanded_data,
    weights = expanded_weights,
    ties = "breslow",
    init = beta_surv_coef,
    control = control,
    robust = FALSE,
    model = FALSE,
    x = FALSE,
    y = FALSE
  )

  list(
    beta_surv_coef = unname(coef(fit)),
    conditional_se = sqrt(diag(vcov(fit))),
    fit = fit
  )
}

BuildSurvivalCovariateMatrix <- function(surv_covar) {
  age <- matrix(
    as.numeric(surv_covar[[1]]),
    ncol = 1,
    dimnames = list(NULL, "age")
  )

  other_blocks <- list()

  if (length(surv_covar) > 1) {
    for (j in 2:length(surv_covar)) {
      covar_mat <- as.matrix(surv_covar[[j]])

      # First column is the reference category.
      if (ncol(covar_mat) > 1) {
        block <- covar_mat[, -1, drop = FALSE]

        block_name <- names(surv_covar)[j]
        if (is.null(block_name) || is.na(block_name) || block_name == "") {
          block_name <- paste0("covar", j)
        }

        colnames(block) <- paste0(
          make.names(block_name),
          "_",
          seq_len(ncol(block))
        )

        other_blocks[[length(other_blocks) + 1]] <- block
      }
    }
  }

  if (length(other_blocks) == 0) {
    return(age)
  }

  cbind(age, do.call(cbind, other_blocks))
}
