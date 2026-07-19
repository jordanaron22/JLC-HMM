# Counts identifiable parameters that are updated by the current model.
CountModelParameters <- function(mix_num,vcovar_num,nu_covar_num,
                                 incl_act,incl_light,joint_model,surv_coef){
  state_num <- NUM_MARKOV_STATES
  marker_num <- as.integer(incl_act) + as.integer(incl_light)

  parameter_count <- c(
    initial_state = mix_num * (state_num - 1),
    transition = mix_num * vcovar_num * state_num *
      (state_num - 1) * 3,
    emission = marker_num * mix_num * vcovar_num * state_num * 2,
    correlation = as.integer(incl_act && incl_light) *
      mix_num * vcovar_num * state_num,
    class_membership = (mix_num - 1) * nu_covar_num,
    survival_class = 0,
    survival_covariates = 0
  )

  if (joint_model){
    if (length(surv_coef) == 0){
      stop("surv_coef is required when counting joint-model parameters")
    }
    parameter_count[["survival_class"]] <- mix_num - 1
    # The first survival covariate is continuous. Remaining coefficient
    # vectors use their first category as the reference level.
    parameter_count[["survival_covariates"]] <-
      1 + sum(pmax(lengths(surv_coef[-1]) - 1,0))
  }

  # The Breslow baseline hazard is profiled rather than included in this
  # regression-parameter count. Its dimension is constant across class counts
  # fitted to the same data, so it does not affect which class count minimizes BIC.
  c(parameter_count,total = sum(parameter_count))
}

# The likelihood is a sum over independent participants, so BIC uses the
# number of participants rather than the number of repeated measurements.
CalcBIC <- function(new_likelihood,num_of_people,parameter_count){
  if (num_of_people < 1){
    stop("num_of_people must be positive")
  }
  if (parameter_count < 1){
    stop("parameter_count must be positive")
  }
  parameter_count * log(num_of_people) - (2 * new_likelihood)
}

CalcAIC <- function(new_likelihood,parameter_count){
  if (parameter_count < 1){
    stop("parameter_count must be positive")
  }
  (2 * parameter_count) - (2 * new_likelihood)
}

#needed to calculate prob of censoring for brier score
Vec2StepPlot <- function(cens_dist,cens_dist_time,t_i){
  for (i in 1:(length(cens_dist_time)-1)){
    if ((t_i >= cens_dist_time[i]) & (t_i < cens_dist_time[i+1])){
      return(cens_dist[i])
    }
  }
  return(cens_dist[length(cens_dist_time)])

}

#calculates bries score
BrierScore <- function(bs_t,surv_event,surv_time,cens_dist,cens_dist_time,surv_mat_ind,event_time){
  bs_sum <- 0
  num_of_people <- length(surv_time)
  cens_dist[length(cens_dist)] <- cens_dist[length(cens_dist)-1]

  for (ind in 1:num_of_people){

    sprob <- Vec2StepPlot(surv_mat_ind[,ind],event_time,bs_t)

    num1 <- sprob^2 * (surv_time[ind] <= bs_t) * surv_event[ind]
    num2 <- (1-sprob)^2 * (surv_time[ind] > bs_t)


    denom1 <- Vec2StepPlot(cens_dist,cens_dist_time,surv_time[ind])
    denom2 <-  Vec2StepPlot(cens_dist,cens_dist_time,bs_t)

    bs_sum <- bs_sum + (num1/denom1) + (num2/denom2)
    if (is.na((num1/denom1) + (num2/denom2))){
      #debugging NA
      # print(ind)
      # print(bs_t)
      # break
    }
  }
  return(bs_sum/num_of_people)
}

#vectorizes bries score
vBrierScore <- Vectorize(BrierScore,vectorize.args = "bs_t")


#calculates c-index
#old
# CalcCindex <- function(surv_time,surv_event,beta_vec,re_prob,
#                        surv_covar_risk_vec){
#   # risk_score <- drop(re_prob %*% beta_vec) + surv_covar_risk_vec
#   risk_score <- log(drop(re_prob %*% exp(beta_vec))) + surv_covar_risk_vec
#   denom <- 0
#   num <- 0
#   for (j in which(surv_event == 1)){
#     comparable <- surv_time > surv_time[j]
#     denom <- denom + sum(comparable)
#     num <- num + sum(risk_score[j] > risk_score[comparable])
#   }

#   return(num/denom)
# }

CalcCindex <- function(
    surv_time,
    surv_event,
    beta_vec,
    re_prob,
    surv_covar_risk_vec,
    sweights_vec
) {

  num_people <- length(surv_time)

  if (is.null(dim(re_prob)) || length(dim(re_prob)) != 2){
    stop("re_prob must be a matrix")
  }

  if (length(surv_event) != num_people ||
      length(surv_covar_risk_vec) != num_people ||
      length(sweights_vec) != num_people ||
      nrow(re_prob) != num_people){
    stop(
      paste(
        "surv_time, surv_event, surv_covar_risk_vec,",
        "sweights_vec, and re_prob must describe the same participants"
      )
    )
  }

  if (length(beta_vec) != ncol(re_prob)){
    stop("length(beta_vec) must equal ncol(re_prob)")
  }

  if (any(!surv_event %in% c(0,1))){
    stop("surv_event must contain only 0 and 1")
  }

  if (any(!is.finite(sweights_vec)) ||
      any(sweights_vec < 0) ||
      sum(sweights_vec) <= 0){
    stop(
      paste(
        "sweights_vec must contain finite, nonnegative weights",
        "with a positive sum"
      )
    )
  }

  if (any(!is.finite(re_prob)) || any(re_prob < 0)){
    stop("re_prob must contain finite, nonnegative probabilities")
  }

  re_prob_sum <- rowSums(re_prob)
  if (any(!is.finite(re_prob_sum)) || any(re_prob_sum <= 0)){
    stop("Each row of re_prob must have positive probability mass")
  }

  re_prob <- re_prob/re_prob_sum

  risk_score <-
    log(drop(re_prob %*% exp(beta_vec))) +
    surv_covar_risk_vec

  if (any(!is.finite(risk_score))){
    stop("risk_score contains non-finite values")
  }

  concordance_data <- data.frame(
    surv_time = surv_time,
    surv_event = surv_event,
    risk_score = risk_score,
    sweights_vec = sweights_vec
  )

  result <- survival::concordance(
    survival::Surv(surv_time, surv_event) ~ risk_score,
    data = concordance_data,
    weights = sweights_vec,
    reverse = TRUE,
    timewt = "n"
  )

  unname(result$concordance)
}

CalcClassEntropy <- function(re_prob){
  if (is.null(dim(re_prob)) || length(dim(re_prob)) != 2){
    stop("re_prob must be a matrix")
  }
  if (any(!is.finite(re_prob)) || any(re_prob < 0)){
    stop("re_prob must contain finite, nonnegative probabilities")
  }

  probability_sums <- rowSums(re_prob)
  if (any(probability_sums <= 0)){
    stop("Each row of re_prob must have positive probability mass")
  }
  normalized_prob <- re_prob / probability_sums
  class_num <- ncol(normalized_prob)

  if (class_num == 1){
    return(list(
      normalized_mean = NA_real_,
      normalized_median = NA_real_,
      normalized_q25 = NA_real_,
      normalized_q75 = NA_real_,
      mean_max_posterior = 1
    ))
  }

  log_prob <- matrix(0,nrow(normalized_prob),ncol(normalized_prob))
  positive_prob <- normalized_prob > 0
  log_prob[positive_prob] <- log(normalized_prob[positive_prob])
  normalized_entropy <-
    -rowSums(normalized_prob * log_prob) / log(class_num)

  list(
    normalized_mean = mean(normalized_entropy),
    normalized_median = median(normalized_entropy),
    normalized_q25 = unname(quantile(normalized_entropy,0.25)),
    normalized_q75 = unname(quantile(normalized_entropy,0.75)),
    mean_max_posterior = mean(apply(normalized_prob,1,max))
  )
}

# Evaluates the fitted cumulative baseline hazard on a new time grid.
BaselineHazardAtTimes <- function(baseline_surv_time,cbline_vec,prediction_times){
  if (length(baseline_surv_time) != length(cbline_vec)){
    stop("baseline_surv_time and cbline_vec must have the same length")
  }

  baseline_times <- sort(unique(baseline_surv_time))
  baseline_values <- vapply(
    baseline_times,
    function(time) max(cbline_vec[baseline_surv_time == time]),
    numeric(1)
  )
  time_index <- findInterval(prediction_times,baseline_times)
  prediction_hazard <- numeric(length(prediction_times))
  prediction_hazard[time_index > 0] <-
    baseline_values[time_index[time_index > 0]]
  prediction_hazard
}

#official implementation with survex package
CalcIBS <- function(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,
                    surv_covar,re_prob,incl_surv,mix_assignment,
                    surv_covar_risk_vec,baseline_surv_time = surv_time){
  event_time <- unique(sort(surv_time))
  cbline_vec_new <- BaselineHazardAtTimes(
    baseline_surv_time,cbline_vec,event_time
  )

  surv_mat_ind <- CalcS(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec)



  IBS_times <- sort(unique(surv_time))
  # IBS_times <- sort(surv_time)
  # Calculate the IBS, again including the observed event times, censoring
  # variables, and the prediction timepoints corresponding to each row of the
  # survival prob matrix

  ibs <- integrated_brier_score(y_true = Surv(surv_time, surv_event),
                                surv = t(surv_mat_ind),
                                times = IBS_times)

  return(ibs)
}

#manual implementation
CalcIBS2 <- function(surv_time,surv_event,cbline_vec,beta_vec,re_prob,
                     surv_covar_risk_vec,baseline_surv_time = surv_time) {
  km_fit <- survfit(Surv(surv_time, 1 - surv_event) ~ 1)
  cens_dist <- c(1,summary(km_fit)$surv)
  cens_dist_time <- c(0,summary(km_fit)$time)
  G <- stepfun(km_fit$time, c(1, km_fit$surv))

  event_time <- unique(sort(surv_time))
  cbline_vec_new <- BaselineHazardAtTimes(
    baseline_surv_time,cbline_vec,event_time
  )

  surv_mat_ind <- CalcS(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec)

  ibs_score <- integrate(vBrierScore,lower = 0,upper = max(event_time),
                         surv_event = surv_event,surv_time = surv_time,cens_dist = cens_dist,
                         cens_dist_time = cens_dist_time,surv_mat_ind = surv_mat_ind,event_time = event_time,
                         rel.tol=.05)
  return(ibs_score$value/max(event_time))
}


########################

###############################################################################
# Summarize participant-specific log scores using NHANES sampling weights
###############################################################################

SummarizeWeightedLogScore <- function(individual_log_score,
                                      sweights_vec) {

  if (length(individual_log_score) != length(sweights_vec)) {
    stop(
      "individual_log_score and sweights_vec must have equal lengths"
    )
  }

  if (any(!is.finite(individual_log_score))) {
    stop("individual_log_score contains non-finite values")
  }

  if (any(!is.finite(sweights_vec)) ||
      any(sweights_vec < 0) ||
      sum(sweights_vec) <= 0) {
    stop(
      paste(
        "sweights_vec must contain finite, nonnegative weights",
        "with a positive sum"
      )
    )
  }

  weighted_sum <- sum(
    sweights_vec * individual_log_score
  )

  weight_sum <- sum(sweights_vec)

  list(
    weighted_mean = weighted_sum / weight_sum,
    weighted_sum = weighted_sum,
    weight_sum = weight_sum,
    individual_log_score = individual_log_score
  )
}


###############################################################################
# Survey-weighted held-out longitudinal log-likelihood
#
# alpha must come from a Forward() call that excludes survival.
###############################################################################

CalcCVLongitudinalLogLik <- function(alpha,
                                     pi_l,
                                     sweights_vec) {

  individual_loglik <- CalcLikelihoodIndVec(
    alpha = alpha,
    pi_l = pi_l
  )

  SummarizeWeightedLogScore(
    individual_log_score = individual_loglik,
    sweights_vec = sweights_vec
  )
}


###############################################################################
# Survey-weighted held-out conditional interval-survival log-likelihood
#
# Death contributions:
#   log{S(lower interval boundary) - S(upper interval boundary)}
#
# Censoring contributions:
#   log{S(observed censoring time)}
#
# cbline_vec and baseline_surv_time must both come from the training sample.
###############################################################################

CalcCVIntervalSurvivalLogLik <- function(
    surv_time,
    surv_event,
    cbline_vec,
    beta_vec,
    re_prob,
    surv_covar_risk_vec,
    sweights_vec,
    baseline_surv_time,
    interval_breaks = c(
      0, 12, 24, 36, 48, 60, 72, 84, 102
    )
) {

  num_people <- length(surv_time)

  ###########################################################################
  # Input checks
  ###########################################################################

  if (is.null(dim(re_prob)) || length(dim(re_prob)) != 2){
    stop("re_prob must be a matrix")
  }

  if (length(surv_event) != num_people ||
      length(surv_covar_risk_vec) != num_people ||
      length(sweights_vec) != num_people ||
      nrow(re_prob) != num_people) {
    stop(
      paste(
        "surv_time, surv_event, surv_covar_risk_vec,",
        "sweights_vec, and re_prob must describe the same participants"
      )
    )
  }

  if (length(beta_vec) != ncol(re_prob)){
    stop("length(beta_vec) must equal ncol(re_prob)")
  }

  if (any(!surv_event %in% c(0, 1))) {
    stop("surv_event must contain only 0 and 1")
  }

  if (any(!is.finite(sweights_vec)) ||
      any(sweights_vec < 0) ||
      sum(sweights_vec) <= 0) {
    stop(
      paste(
        "sweights_vec must contain finite, nonnegative weights",
        "with a positive sum"
      )
    )
  }

  if (length(baseline_surv_time) != length(cbline_vec)) {
    stop(
      "baseline_surv_time and cbline_vec must have equal lengths"
    )
  }

  interval_breaks <- as.numeric(interval_breaks)

  if (length(interval_breaks) < 2 ||
      any(!is.finite(interval_breaks)) ||
      any(diff(interval_breaks) <= 0)) {
    stop(
      "interval_breaks must be a strictly increasing finite vector"
    )
  }

  if (interval_breaks[1] != 0) {
    stop("The first interval boundary must be 0")
  }

  if (any(!is.finite(surv_time))){
    stop("surv_time must contain only finite values")
  }

  if (any(surv_time < interval_breaks[1])) {
    stop(
      paste0(
        "All survival times must be >= ",
        interval_breaks[1],
        "."
      )
    )
  }

  if (any(surv_event == 1 & surv_time > max(interval_breaks))) {
    stop(
      paste0(
        "All death times must fall within interval_breaks. ",
        "The last interval boundary is ",
        max(interval_breaks),
        ". Censored times may exceed this boundary."
      )
    )
  }

  if (any(!is.finite(re_prob)) || any(re_prob < 0)){
    stop("re_prob must contain finite, nonnegative probabilities")
  }

  re_prob_sum <- rowSums(re_prob)

  if (any(!is.finite(re_prob_sum)) ||
      any(re_prob_sum <= 0)) {
    stop("Each row of re_prob must have positive probability mass")
  }

  re_prob <- re_prob / re_prob_sum

  ###########################################################################
  # Determine the interval containing each observed death
  #
  # right = TRUE gives:
  #   [0,12], (12,24], ..., (84,102]
  ###########################################################################

  event_indicator <- surv_event == 1
  censored_indicator <- !event_indicator

  event_interval_index <- rep(
    NA_integer_,
    num_people
  )

  if (any(event_indicator)) {

    event_interval_index[event_indicator] <- cut(
      surv_time[event_indicator],
      breaks = interval_breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = FALSE
    )

    if (anyNA(event_interval_index[event_indicator])) {
      stop("At least one death was not assigned to an interval")
    }
  }

  event_lower_time <- rep(
    NA_real_,
    num_people
  )

  event_upper_time <- rep(
    NA_real_,
    num_people
  )

  if (any(event_indicator)) {

    current_interval <-
      event_interval_index[event_indicator]

    event_lower_time[event_indicator] <-
      interval_breaks[current_interval]

    event_upper_time[event_indicator] <-
      interval_breaks[current_interval + 1L]
  }

  ###########################################################################
  # Times at which survival probabilities are needed
  #
  # Events require interval boundaries.
  # Censored observations require their exact censoring time.
  ###########################################################################

  prediction_times <- sort(
    unique(
      c(
        interval_breaks,
        surv_time[censored_indicator]
      )
    )
  )

  cbline_prediction <- BaselineHazardAtTimes(
    baseline_surv_time = baseline_surv_time,
    cbline_vec = cbline_vec,
    prediction_times = prediction_times
  )

  ###########################################################################
  # Marginal survival predictions over latent-class uncertainty
  #
  # Rows are prediction times and columns are held-out participants.
  ###########################################################################

  survival_predictions <- CalcS(
    event_time = prediction_times,
    cbline_vec_new = cbline_prediction,
    beta_vec = beta_vec,
    re_prob = re_prob,
    surv_covar_risk_vec = surv_covar_risk_vec
  )

  if (!all(
    dim(survival_predictions) ==
      c(length(prediction_times), num_people)
  )) {
    stop(
      paste(
        "CalcS returned survival predictions with unexpected",
        "dimensions"
      )
    )
  }

  # Protect against negligible floating-point excursions.
  survival_predictions <- pmin(
    pmax(survival_predictions, 0),
    1
  )

  ###########################################################################
  # Probability assigned to each observed outcome
  ###########################################################################

  observed_probability <- numeric(num_people)

  if (any(event_indicator)) {

    event_people <- which(event_indicator)

    lower_row <- match(
      event_lower_time[event_people],
      prediction_times
    )

    upper_row <- match(
      event_upper_time[event_people],
      prediction_times
    )

    survival_lower <- survival_predictions[
      cbind(lower_row, event_people)
    ]

    survival_upper <- survival_predictions[
      cbind(upper_row, event_people)
    ]

    event_probability <-
      survival_lower -
      survival_upper

    # A substantial negative value would imply a nonmonotone predicted
    # survival curve and should not be silently corrected.
    if (any(event_probability < -1e-12)) {
      stop(
        "Predicted survival increased across at least one event interval"
      )
    }

    event_probability <- pmax(
      event_probability,
      0
    )

    observed_probability[event_people] <-
      event_probability
  }

  if (any(censored_indicator)) {

    censored_people <- which(censored_indicator)

    censor_row <- match(
      surv_time[censored_people],
      prediction_times
    )

    observed_probability[censored_people] <-
      survival_predictions[
        cbind(censor_row, censored_people)
      ]
  }

  ###########################################################################
  # Check for zero-probability observations
  ###########################################################################

  zero_probability <- which(
    !is.finite(observed_probability) |
      observed_probability <= 0
  )

  if (length(zero_probability) > 0) {

    zero_event <- zero_probability[
      surv_event[zero_probability] == 1
    ]

    if (length(zero_event) > 0) {

      zero_intervals <- unique(
        paste0(
          "(",
          event_lower_time[zero_event],
          ", ",
          event_upper_time[zero_event],
          "]"
        )
      )

      stop(
        paste0(
          "At least one held-out death received zero interval ",
          "probability. Problem interval(s): ",
          paste(zero_intervals, collapse = ", "),
          ". Check whether the training baseline cumulative hazard ",
          "increased within each interval."
        )
      )
    }

    stop(
      "At least one censored observation received zero survival probability"
    )
  }

  ###########################################################################
  # Participant-specific and weighted log scores
  ###########################################################################

  individual_loglik <- log(
    observed_probability
  )

  result <- SummarizeWeightedLogScore(
    individual_log_score = individual_loglik,
    sweights_vec = sweights_vec
  )

  result$observed_probability <-
    observed_probability

  result$surv_time <- surv_time
  result$surv_event <- surv_event
  result$event_interval_index <-
    event_interval_index

  result$event_lower_time <-
    event_lower_time

  result$event_upper_time <-
    event_upper_time

  result$interval_breaks <-
    interval_breaks

  result$prediction_times <-
    prediction_times

  result
}
