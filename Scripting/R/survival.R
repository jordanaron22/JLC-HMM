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
  bline_vec <- numeric(n)
  cbline_vec <- numeric(n)

  for (time_ind in 1:n){
    risk_set <- surv_time >= surv_time[time_ind]
    linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]

    if (dim(re_prob)[2] != 1){
      denom <- sum((re_prob[risk_set,] %*% exp(beta_vec)) * exp(linear_surv_covar_risk))
    } else {
      denom <- sum(exp(linear_surv_covar_risk))
    }

    if (denom > .Machine$double.xmax){denom <- .Machine$double.xmax}
    bline_vec[time_ind] <- surv_event[time_ind]/denom
  }

  for(time_ind in 1:n){
    anti_risk_set <- surv_time <= surv_time[time_ind]
    cbline_vec[time_ind] <- sum(bline_vec[anti_risk_set])
  }

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


  loglike <- sum(log(bline_vec^surv_event) +
                   ((re_prob %*% beta_vec)+surv_covar_risk_vec) * surv_event -
                   cbline_vec*exp((re_prob %*% beta_vec)+surv_covar_risk_vec))

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


  for (i in 1:(length(surv_coef_len)-1)) {
    coef_vec <- c(0,beta_surv_coef[(fit_mix_num+1+surv_coef_len_alt[i]):(fit_mix_num+surv_coef_len_alt[i+1])])
    surv_coef_new <- append(surv_coef_new,list(coef_vec))
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
  for (i in 1:length(altered_surv_coef)){
    altered_surv_coef[[i]] <- altered_surv_coef[[i]][-1]
  }

  beta_surv_coef[(fit_mix_num+1):beta_surv_coef_len] <- unlist(altered_surv_coef)

  return(beta_surv_coef)
}

#LM approach for calculating survival coefficients
CalcBetaManual <- function(beta_surv_coef,surv_covar_risk_vec,stop_crit,
                           survival_context,surv_coef_len,fit_mix_num){
  re_prob <- survival_context$re_prob
  surv_event <- survival_context$surv_event
  surv_time <- survival_context$surv_time
  surv_covar <- survival_context$surv_covar
  l2norm <- 101
  while (l2norm > stop_crit){

    beta_vec <- OutofBetaSurvCoef(beta_surv_coef,surv_coef_len,fit_mix_num)[[1]]
    surv_coef <- OutofBetaSurvCoef(beta_surv_coef,surv_coef_len,fit_mix_num)[[2]]

    old_slike <- SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,
                          survival_context)
    slike_diff <- -1


    grad <- numeric(fit_mix_num + length(unlist(surv_coef)) - length(surv_coef))
    hess <- matrix(0,fit_mix_num + length(unlist(surv_coef)) - length(surv_coef),fit_mix_num + length(unlist(surv_coef)) - length(surv_coef))

    for (ind in which(surv_event == 1)){
      risk_set <- surv_time >= surv_time[ind]


      linear_surv_covar_risk <- surv_covar_risk_vec[risk_set]
      ageadj_risk <- exp(linear_surv_covar_risk) %x% t(exp(beta_vec))

      num_list <- list()

      #First survival covariate must be age
      #Age is only continuous surv covar so its treated differently
      num0 <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[1]][risk_set] %x% t(numeric(fit_mix_num)+1))
      num02 <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[1]][risk_set]^2 %x% t(numeric(fit_mix_num)+1))
      num_list[[1]] <- c(num0,num02)

      denom <- sum(re_prob[risk_set,] * ageadj_risk)

      #fix age to be first variable
      grad[1] <- grad[1] + (surv_covar[[1]][ind] - num0/denom)
      hess[1,1] <- hess[1,1] + (num02/denom) - (num0/denom)^2



      if (fit_mix_num > 1){
        for (beta_ind in 2:fit_mix_num){
          num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk))
          grad[beta_ind] <- grad[beta_ind] + (re_prob[ind,beta_ind] - num/denom)

          hess[beta_ind,beta_ind] <- hess[beta_ind,beta_ind] + (num/denom) - (num/denom)^2

          num_cross <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[1]][risk_set])
          hess[1,beta_ind] <- hess[1,beta_ind] + (denom*num_cross - num*num0)/denom^2
          hess[beta_ind,1] <- hess[beta_ind,1] + (denom*num_cross - num*num0)/denom^2
        }
      }


      if (length(surv_covar) > 1){
        for (surv_covar_ind in 2:length(surv_covar)){
          #need to only include people with same covariate status in numerator
          #num0 and num02 will be equal, need to adjust accordingly

          num_disc_covar_vec <- c()
          for (surv_covar_indicator in 2:length(surv_coef[[surv_covar_ind]])){

            num_disc_covar <- sum(re_prob[risk_set,] * ageadj_risk *
                                    surv_covar[[surv_covar_ind]][risk_set,surv_covar_indicator] %x% t(numeric(fit_mix_num)+1))

            num_disc_covar_vec <- c(num_disc_covar_vec,num_disc_covar)
            num_list[[surv_covar_ind]] <- num_disc_covar_vec
          }
        }
      }

      #Covariates
      list_of_lens <- unlist(lapply(num_list,length))
      list_of_lens[1] <- 0
      list_of_cum_lens <- cumsum(list_of_lens)
      for (surv_covar_ind in 2:length(surv_covar)){
        #Get indicies for covariates in gradient/hessian
        starting_index <- fit_mix_num + list_of_cum_lens[surv_covar_ind-1] +1
        ending_index <- fit_mix_num + list_of_cum_lens[surv_covar_ind]
        curr_covar_inds <- c(starting_index:ending_index)

        # grad[5] <- grad[5] + (surv_covar[[2]][ind,2] - num_list[[2]][2]/denom)
        grad[curr_covar_inds] <- grad[curr_covar_inds] + surv_covar[[surv_covar_ind]][ind,-1] - num_list[[surv_covar_ind]]/denom

        # hess[5,5] <- hess[5,5] + num_list[[2]][2]/denom - (num_list[[2]][2]/denom)^2
        diag(hess)[curr_covar_inds] <- diag(hess)[curr_covar_inds] + num_list[[surv_covar_ind]]/denom - (num_list[[surv_covar_ind]]/denom)^2

        #Cross terms for age
        for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
          # num_cross_age2 <- sum(re_prob[risk_set,] * ageadj_risk * age_vec[risk_set]*surv_covar[[2]][risk_set,2]%x% t(numeric(fit_mix_num)+1))
          # hess[1,5] <- hess[1,5] + (denom*num_cross_age2 - num_list[[2]][[2]]*num0)/denom^2
          num_cross_age <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[1]][risk_set]*surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1]%x% t(numeric(fit_mix_num)+1))
          hess[1,curr_covar_inds[ind_curr_covar_inds]] <- hess[1,curr_covar_inds[ind_curr_covar_inds]] + (denom*num_cross_age - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num0)/denom^2
          hess[curr_covar_inds[ind_curr_covar_inds],1] <- hess[curr_covar_inds[ind_curr_covar_inds],1] + (denom*num_cross_age - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num0)/denom^2
        }

        #Cross terms for betas
        if (fit_mix_num > 1){
          for (beta_ind in 2:fit_mix_num){

            num <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk))

            for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
              #sample code, vectorized this now
              # num_cross2 <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[2]][risk_set,2])
              # hess[5,beta_ind] <- hess[5,beta_ind] + (denom*num_cross2 - num*num_list[[2]][[2]])/denom^2

              num_cross <- sum(re_prob[risk_set,beta_ind] * exp(beta_vec[beta_ind])*exp(linear_surv_covar_risk)* surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1])
              hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] <- hess[curr_covar_inds[ind_curr_covar_inds],beta_ind] +
                (denom*num_cross - num*num_list[[surv_covar_ind]][[ind_curr_covar_inds]])/denom^2

              hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] <- hess[beta_ind,curr_covar_inds[ind_curr_covar_inds]] +
                (denom*num_cross - num*num_list[[surv_covar_ind]][[ind_curr_covar_inds]])/denom^2
            }
          }
        }


        #Cross terms for surv covariates
        if(surv_covar_ind != length(surv_covar)){
          for (surv_covar_ind2 in (surv_covar_ind+1):length(surv_covar)){

            starting_index2 <- fit_mix_num + list_of_cum_lens[surv_covar_ind2-1] +1
            ending_index2 <- fit_mix_num + list_of_cum_lens[surv_covar_ind2]
            curr_covar_inds2 <- c(starting_index2:ending_index2)

            for (ind_curr_covar_inds in 1:length(curr_covar_inds)){
              for (ind_curr_covar_inds2 in 1:length(curr_covar_inds2)){

                num_cross_covar <- sum(re_prob[risk_set,] * ageadj_risk * surv_covar[[surv_covar_ind]][risk_set,ind_curr_covar_inds+1]*surv_covar[[surv_covar_ind2]][risk_set,ind_curr_covar_inds2+1]%x% t(numeric(fit_mix_num)+1))

                hess[curr_covar_inds2[ind_curr_covar_inds2],curr_covar_inds[ind_curr_covar_inds]] <- hess[curr_covar_inds2[ind_curr_covar_inds2],curr_covar_inds[ind_curr_covar_inds]] +
                  (denom*num_cross_covar - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num_list[[surv_covar_ind2]][[ind_curr_covar_inds2]])/denom^2

                hess[curr_covar_inds[ind_curr_covar_inds],curr_covar_inds2[ind_curr_covar_inds2]] <- hess[curr_covar_inds[ind_curr_covar_inds],curr_covar_inds2[ind_curr_covar_inds2]] +
                  (denom*num_cross_covar - num_list[[surv_covar_ind]][[ind_curr_covar_inds]]*num_list[[surv_covar_ind2]][[ind_curr_covar_inds2]])/denom^2
              }
            }



          }
        }



      }


    }

    #LM aspect
    step_size <- .01
    step_mat <- matrix(0,dim(hess)[1],dim(hess)[2])
    nr_fact <- numeric(length(grad))-1
    max_val <- 6
    while(slike_diff < 0|all(nr_fact == -1) | max_val > 5){
      diag(step_mat) <- (diag(step_mat) + step_size - .01)
      nr_fact <- SolveCatch(hess + step_mat,-grad)
      step_size <- step_size * 10


      beta_surv_coef_new <- beta_surv_coef-nr_fact

      beta_vec_new <- OutofBetaSurvCoef(beta_surv_coef_new,surv_coef_len,fit_mix_num)[[1]]
      surv_coef_new <- OutofBetaSurvCoef(beta_surv_coef_new,surv_coef_len,fit_mix_num)[[2]]
      max_val <- abs(max(c(beta_vec_new,unlist(surv_coef_new))))


      new_slike <- SurvLike(beta_vec_new,surv_covar_risk_vec,surv_coef_new,
                            survival_context)
      slike_diff <- old_slike - new_slike
    }

    l2norm <- sum(sqrt((beta_surv_coef_new - beta_surv_coef)^2))



    beta_surv_coef <- beta_surv_coef_new

  }
  return(list(beta_surv_coef,sqrt(diag(solve(hess)))))
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
