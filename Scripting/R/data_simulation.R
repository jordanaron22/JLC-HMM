#Simulates hidden states
SimulateMC <- function(day_length,init,tran_list_ind,mixture_ind,vcovar_vec){
  hidden_states <- numeric(day_length)

  for (i in 1:day_length){

    tran <- tran_list_ind[[vcovar_vec[i]]][[(i-1)%%period_len+1]]

    if (i == 1) {
      hidden_states[1] <- rbinom(1,1,init[mixture_ind,2])
    } else {
      hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
    }
  }


  return(hidden_states)
}

#used to simulate survival/censor times
finv <- function(lam,time, randu, xb_ind){
  return((1-pexp(time,rate = lam))^(exp(xb_ind)) - randu)
}

#simulates survival by transforming runif var into rexponential
#uses: S(t) = 1-F(t)
#finds root of survival time: 1 - F(t|X) - u
SimSurvival <- function(mixture_vec,beta_vec,beta_age_true,beta_covar_sim,age_vec,surv_covar_sim,lam = 1/20){

  num_of_people <- length(mixture_vec)

  failure_times <- numeric(num_of_people)
  censor_times <- numeric(num_of_people)

  for(i in 1:num_of_people){
    xb <- beta_vec[mixture_vec[i]] + beta_age_true*age_vec[i] + beta_covar_sim[surv_covar_sim[i]]

    evnt<-runif(1)
    cens<-runif(1)
    failure_times[i] <- uniroot(finv, interval=c(0, 60), lam=lam, randu=evnt, xb_ind=xb, extendInt = "yes")$root
    censor_times[i] <- uniroot(finv, interval=c(0, 60), lam=lam, randu=cens, xb_ind=xb, extendInt = "yes")$root

  }

  time <- pmin(failure_times, censor_times)
  event <- as.integer(failure_times<censor_times)

  return(list(time = time,event = event))
}

#Simulates data
SimulateHMM <- function(day_length,num_of_people,init,params_tran_array,
                        emit_act,emit_light,corr_mat,
                        lod_act,lod_light, nu_mat,beta_vec_true,beta_age_true,beta_covar_sim,
                        missing_perc,lambda_act_mat,lambda_light_mat,
                        true_mix_num){

  mix_num <- dim(emit_act)[3]
  #this mix_num is confusing but is talking about simulation parameters, so we want that to match the truth
  if (mix_num != true_mix_num){
    stop("true_mix_num does not match the number of classes in the truth parameters")
  }

  #simulates age and single hypothetical categorical variable with 3 categories
  age_vec <- floor(runif(num_of_people,10,81))
  surv_covar_sim_wide <- t(rmultinom(num_of_people,1,c(1/3,1/3,1/3)))
  surv_covar_sim <- apply(surv_covar_sim_wide, 1, which.max)

  #simulates stationary and moderact activity and then calulates nu matrix
  statact_vec <- floor(runif(num_of_people,0,16))
  modact_vec <- rbinom(num_of_people,1,.5)
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,statact_vec,statact_vec^2)

  pi_l_true <- CalcPi(nu_mat,nu_covar_mat)

  #edge case for 1 mixture model
  if (dim(pi_l_true)[2] > 1){
    mixture_vec <- rMultinom(pi_l_true,1)
  } else {
    mixture_vec <- matrix(1,nrow = dim(pi_l_true)[1],ncol = 1)
  }

  #simulates week/weekend divide
  #only works for specific lengths
  if (day_length == 192){
    vcovar_vec <- c(rep(0,48),rep(1,period_len),rep(0,48))
  } else if (day_length == 96){
    vcovar_vec <- rep(0,96)
  } else if (day_length == 288){
    vcovar_vec <- c(rep(0,96*2),rep(1,96*1))
  } else if (day_length == 384){
    vcovar_vec <- c(rep(0,96*2),rep(1,96*2))
  } else if (day_length == 672){
    vcovar_vec <- c(rep(0,96*5),rep(1,96*2))
  } else {
    vcovar_vec <- c(rep(0,48),rep(1,96),rep(0,48))
    vcovar_vec <- rep(vcovar_vec,day_length/192)
  }

  vcovar_mat <- replicate(num_of_people, vcovar_vec)

  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                           period_len = period_len)


  #actually simulates data here
  for (ind in 1:num_of_people){
    activity <- numeric(day_length)
    light <- numeric(day_length)

    mixture_ind <- mixture_vec[ind]

    tran_vcovar_list <- tran_list[[mixture_ind]]
    vcovar_vec_ind <- vcovar_mat[,ind]+1

    hidden_states <- SimulateMC(day_length,init,tran_vcovar_list,mixture_ind,vcovar_vec_ind)

    for (i in 1:day_length){

      #sets normal parameters for current mixture and week status
      mu_act <- emit_act[hidden_states[i] + 1,1,mixture_ind,vcovar_vec_ind[i]]
      sig_act <- emit_act[hidden_states[i] + 1,2,mixture_ind,vcovar_vec_ind[i]]
      mu_light <- emit_light[hidden_states[i] + 1,1,mixture_ind,vcovar_vec_ind[i]]
      sig_light <- emit_light[hidden_states[i] + 1,2,mixture_ind,vcovar_vec_ind[i]]
      bivar_corr <- corr_mat[mixture_ind,hidden_states[i] + 1,vcovar_vec_ind[i]]

      sigma_mat <- matrix(c(sig_act^2,bivar_corr * sig_act* sig_light,
                            bivar_corr * sig_act* sig_light,sig_light^2),2,2,byrow = T)

      act_light <- mvrnorm(n = 1,
                           mu = c(mu_act,mu_light),
                           Sigma = sigma_mat)

      activity[i] <-act_light[1]
      light[i] <-act_light[2]

    }

    if (ind == 1){
      hidden_states_matrix <- hidden_states
      activity_matrix <- activity
      light_matrix <- light
    } else {
      hidden_states_matrix <- cbind(hidden_states_matrix,hidden_states)
      activity_matrix <- cbind(activity_matrix,activity)
      light_matrix <- cbind(light_matrix,light)
    }
  }

  #simulates survival data
  surv_list <- SimSurvival(mixture_vec,beta_vec_true,beta_age_true,beta_covar_sim,age_vec,surv_covar_sim)


  # Censor values below the limit of detection.
  light_matrix[light_matrix<lod_light] <- lod_light
  activity_matrix[activity_matrix<lod_act] <- lod_act

  #removes some data
  act_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                        ncol = num_of_people)
  light_missing <- matrix(rbinom(day_length * num_of_people,1,missing_perc),
                          ncol = num_of_people)

  activity_matrix[act_missing==1] <- NA
  light_matrix[light_missing==1] <- NA



  make_simulated_hmm_list(mc = hidden_states_matrix,
                          act = activity_matrix,
                          light = light_matrix,
                          mixture_mat = mixture_vec,
                          age_vec = age_vec,
                          nu_covar_mat = nu_covar_mat,
                          vcovar_mat = vcovar_mat,
                          survival = surv_list,
                          surv_covar_sim = surv_covar_sim)
}
