#calculates linear value of risk due to sociodemo covar, does not incl mixture
SurvCovarRiskVec <- function(surv_covar,surv_coef){
  surv_covar[[1]] <- matrix(surv_covar[[1]],ncol = 1)
  surv_coef[[1]] <- matrix(surv_coef[[1]],ncol = 1)
  surv_covar_risk_vec <- rowSums(mapply(function(x, y) x %*% y,  surv_covar,surv_coef, SIMPLIFY = T))
  return(surv_covar_risk_vec)
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


  # loglike <- sum(log(bline_vec^surv_event) +
  #                  ((re_prob %*% beta_vec)+surv_covar_risk_vec) * surv_event -
  #                  cbline_vec*exp((re_prob %*% beta_vec)+surv_covar_risk_vec))

  #keeping above lines as found this bug after most of runs were done
  #exponentiating the probabilities in first, should only do that to betas in last line
  class_event_lp <- drop(re_prob %*% beta_vec)
  class_risk <- drop(re_prob %*% exp(beta_vec))

  event_log_bline <- numeric(length(surv_event))
  event_rows <- surv_event == 1
  event_log_bline[event_rows] <- log(bline_vec[event_rows])

  loglike_i <-
    event_log_bline +
    surv_event * (
      class_event_lp +
      surv_covar_risk_vec
    ) -
    cbline_vec *
      exp(surv_covar_risk_vec) *
      class_risk

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
    event_weight <- sweights_vec[ind]

    risk_set <- surv_time >= surv_time[ind]
    linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]

    risk_weight_mat <- re_prob[risk_set,,drop = FALSE] *
      sweights_vec[risk_set]*
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

    grad <- grad + event_weight * (build_event_z(ind) - ez)
    hess <- hess + event_weight * (ezz - tcrossprod(ez))
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

  if (incl_surv == MODEL_TYPE_CODES[["joint"]]) {
    #used to use calcbetamanual but this is equivalent 
    #safer to as it uses coxph
    coxph_result <- CalcBetaCoxphJoint(
      beta_surv_coef = beta_surv_coef,
      survival_context = survival_context,
      surv_coef_len = surv_coef_len,
      fit_mix_num = fit_mix_num,
      one_step = FALSE
    )

    beta_surv_coef_new <-
      coxph_result$beta_surv_coef

    conditional_se <-
      coxph_result$conditional_se

    expected_parameter_num <-
      length(beta_surv_coef)

    if (
      length(beta_surv_coef_new) !=
        expected_parameter_num
    ) {
      stop(
        paste(
          "CalcBetaCoxphJoint returned",
          length(beta_surv_coef_new),
          "coefficients, but",
          expected_parameter_num,
          "were expected"
        )
      )
    }

    if (any(!is.finite(beta_surv_coef_new))) {
      stop(
        "CalcBetaCoxphJoint returned nonfinite coefficients"
      )
    }

    if (
      length(conditional_se) !=
        expected_parameter_num
    ) {
      stop(
        paste(
          "CalcBetaCoxphJoint returned",
          length(conditional_se),
          "standard errors, but",
          expected_parameter_num,
          "were expected"
        )
      )
    }

    return(
      list(
        beta_surv_coef_new,
        conditional_se
      )
    )
  }





  surv_data <- data.frame(time = survival_context$surv_time,
                          status = survival_context$surv_event,
                          age = survival_context$surv_covar[[1]])

  surv_data <- cbind(surv_data,survival_context$re_prob,combined_covar_mat)
  colnames(surv_data)[4] <- "toRem"

  fit <- coxph(Surv(time, status) ~ .  - toRem, data = surv_data, weights = survival_context$sweights_vec, ties = 'breslow')
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