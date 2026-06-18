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
CalcCindex <- function(surv_time,surv_event,beta_vec,surv_coef,re_prob,surv_covar,surv_covar_risk_vec){
  denom <- 0
  num <- 0
  for (j in 1:length(surv_time)){


    for (i in 1:length(surv_time)){
      if (surv_event[j] == 1){
        if (i != j){

          risk_i <- (re_prob[i,] %*% beta_vec)+surv_covar_risk_vec[i]
          risk_j <- (re_prob[j,] %*% beta_vec)+surv_covar_risk_vec[j]

          if (surv_time[i] > surv_time[j]){
            denom <- denom + 1

            if (risk_j > risk_i){
              num <- num + 1
            }

          }


        }
      }

    }
  }

  return(num/denom)
}

#official implementation with survex package
CalcIBS <- function(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,surv_covar,re_prob,incl_surv,mix_assignment,surv_covar_risk_vec){
  event_time <- unique(sort(surv_time))
  cbline_vec_new <- unique(sort(cbline_vec))
  # cbline_vec_new <- basehaz(cox1,centered = F)[,1]

  stime_vec <- c()
  for (stime in event_time){
    stime_vec <- c(stime_vec,which(surv_time == stime)[1])
  }

  cbline_vec_new <- cbline_vec[stime_vec]

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
CalcIBS2 <- function(surv_time,surv_event,cbline_vec,beta_vec,re_prob,surv_covar_risk_vec) {
  km_fit <- survfit(Surv(surv_time, 1 - surv_event) ~ 1)
  cens_dist <- c(1,summary(km_fit)$surv)
  cens_dist_time <- c(0,summary(km_fit)$time)
  G <- stepfun(km_fit$time, c(1, km_fit$surv))

  event_time <- unique(sort(surv_time))
  cbline_vec_new <- unique(sort(cbline_vec))
  # cbline_vec_new <- basehaz(cox1,centered = F)[,1]

  stime_vec <- c()
  for (stime in event_time){
    stime_vec <- c(stime_vec,which(surv_time == stime)[1])
  }

  cbline_vec_new <- cbline_vec[stime_vec]

  surv_mat_ind <- CalcS(event_time,cbline_vec_new,beta_vec,re_prob,surv_covar_risk_vec)

  ibs_score <- integrate(vBrierScore,lower = 0,upper = max(event_time),
                         surv_event = surv_event,surv_time = surv_time,cens_dist = cens_dist,
                         cens_dist_time = cens_dist_time,surv_mat_ind = surv_mat_ind,event_time = event_time,
                         rel.tol=.05)
  return(ibs_score$value/max(event_time))
}
