#turns F-B into simple prob of wake/sleep arrays
CondMarginalize <- function(alpha,beta,pi_l){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)

  #need to add prob of being in each mixture
  #not included directly in F-B
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[ind,re_ind])
    }
  }

  #individual likelihood of being in specific mixture
  ind_like_mat <- apply(alpha_beta,c(1,4),logSumExp)

  weight_array_wake <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  weight_array_sleep <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  for (ind in 1:dim(alpha_beta)[4]){
    for (t in 1:dim(alpha_beta)[1]){
      weight_array_wake[t,ind,] <- alpha_beta[t,1,,ind] - ind_like_mat[t,ind]
      weight_array_sleep[t,ind,] <- alpha_beta[t,2,,ind] - ind_like_mat[t,ind]
    }
  }

  return(list(weight_array_wake,weight_array_sleep))
}

#calculates initial probability
CalcInit <- function(alpha, beta,pi_l,sweights_vec){

  #setup
  num_obs <- dim(alpha[[1]][,,1])[1]
  time <- 1
  init_0_vec <- matrix(0,length(alpha),dim(pi_l)[2])
  init_1_vec <- matrix(0,length(alpha),dim(pi_l)[2])
  init_mat <- matrix(0,dim(pi_l)[2],2)

  #individual likelihood vector
  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)

  for(ind in 1:length(alpha)){
    ind_like <- ind_like_vec[ind]

    init_0_vec[ind,] <- alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(pi_l[ind,]) - ind_like + log(sweights_vec[ind])
    init_1_vec[ind,] <- alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(pi_l[ind,]) - ind_like + log(sweights_vec[ind])
  }

  #normalizes
  for (re_ind in 1:(dim(pi_l)[2])){
    init_0 <- logSumExp(init_0_vec[,re_ind])
    init_1 <- logSumExp(init_1_vec[,re_ind])
    init_vec <- exp(c(init_0,init_1) - logSumExp(c(init_0,init_1)))
    init_mat[re_ind,] <- init_vec
  }

  return(init_mat)
}

#calculates probability of an individual being in each mixture
CalcProbRE <- function(alpha,pi_l){

  len <- dim(alpha[[1]])[1]
  re_len <- dim(alpha[[1]])[3]
  re_weight_vec <- numeric(re_len)
  re_weights <- matrix(0,nrow = length(alpha),ncol = re_len)

  for (ind in 1:length(alpha)){
    for (re_ind in 1:re_len){
      #sums over latent states for last time and specific mixture
      re_weights[ind,re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(pi_l[ind,re_ind])
    }
    #normalizes
    re_weights[ind,] <- exp(re_weights[ind,] - logSumExp(c(re_weights[ind,])))
  }

  return(re_weights)
}

#calculates overall likelihood by summing individual likelihood
#accounts for sample weights vector
CalcLikelihood <- function(alpha, pi_l, sweigts_vec){
  individual_loglik <- CalcLikelihoodIndVec(alpha, pi_l)

  if (length(sweigts_vec) != length(individual_loglik)){
    stop("sample_weights must have one value per subject")
  }

  sum(sweigts_vec * individual_loglik)
}

#calculates likelihood vector, each individual as entry
CalcLikelihoodIndVec <- function(alpha,pi_l){
  num_obs <- dim(alpha[[1]][,,1])[1]
  like_vec <- c()
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(SumOverREIndTime(alpha,pi_l,i,num_obs)))
    like_vec <- c(like_vec,ind_like)
  }
  return(like_vec)
}

#adds pi to F-B output
SumOverREIndTime <- function(fb,pi_l,ind,time, add_re = T){

  fb_ind <- fb[[ind]]

  fb_sum <- numeric(2)
  if (add_re){
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,] + log(pi_l[ind,])))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,] + log(pi_l[ind,])))
  } else {
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,]))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,]))
  }

  return(fb_sum)
}

#Calculates likelihood of individual
IndLike <- function(alpha,pi_l,ind,len){
  likelihood <- logSumExp(SumOverREIndTime(alpha,pi_l,ind,len))
  return(likelihood)
}

#viterbi algorithm for global decoding
Viterbi <- function(act,light,vcovar_mat){
  decoded_array <- array(NA, dim = c(day_length,num_of_people,mix_num))
  for (ind in 1:num_of_people){
    for (clust_i in 1:mix_num){
      vit_ind_vec <- ViterbiIndHelper(ind,clust_i,act,light,vcovar_mat)
      decoded_array[,ind,clust_i] <- vit_ind_vec
    }
  }
  return(decoded_array)
}

#Viterbi algorithm on an individual given their parameters and conditioning on mixture
ViterbiIndHelper <- function(ind,clust_i,act,light,vcovar_mat){

  tran_list_clust <- tran_list[[clust_i]]

  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))

  vcovar_vec <- vcovar_mat[,ind]

  log_class_0_week <- logClassificationC( act[,ind], light[,ind],
                                          emit_act_week[1,1,clust_i],
                                          emit_act_week[1,2,clust_i],
                                          emit_light_week[1,1,clust_i],
                                          emit_light_week[1,2,clust_i],
                                          lod_act, lod_light, corr_mat[clust_i,1,1], lintegral_mat[clust_i,1,1],
                                          lambda_act_mat[clust_i,1,1],lambda_light_mat[clust_i,1,1],tobit)

  log_class_1_week <- logClassificationC( act[,ind], light[,ind],
                                          emit_act_week[2,1,clust_i],
                                          emit_act_week[2,2,clust_i],
                                          emit_light_week[2,1,clust_i],
                                          emit_light_week[2,2,clust_i],
                                          lod_act, lod_light, corr_mat[clust_i,2,1], lintegral_mat[clust_i,2,1],
                                          lambda_act_mat[clust_i,2,1],lambda_light_mat[clust_i,2,1],tobit)

  log_class_0_weekend <- logClassificationC( act[,ind], light[,ind],
                                          emit_act_weekend[1,1,clust_i],
                                          emit_act_weekend[1,2,clust_i],
                                          emit_light_weekend[1,1,clust_i],
                                          emit_light_weekend[1,2,clust_i],
                                          lod_act, lod_light, corr_mat[clust_i,1,2], lintegral_mat[clust_i,1,2],
                                          lambda_act_mat[clust_i,1,2],lambda_light_mat[clust_i,1,2],tobit)

  log_class_1_weekend <- logClassificationC( act[,ind], light[,ind],
                                          emit_act_weekend[2,1,clust_i],
                                          emit_act_weekend[2,2,clust_i],
                                          emit_light_weekend[2,1,clust_i],
                                          emit_light_weekend[2,2,clust_i],
                                          lod_act, lod_light, corr_mat[clust_i,2,2], lintegral_mat[clust_i,2,2],
                                          lambda_act_mat[clust_i,2,2],lambda_light_mat[clust_i,2,2],tobit)

  log_class_0 <- (log_class_0_week * (1-vcovar_vec)) + (log_class_0_weekend * vcovar_vec)
  log_class_1 <- (log_class_1_week * (1-vcovar_vec)) + (log_class_1_weekend * vcovar_vec)

  viterbi_mat <- matrix(NA,2,day_length)
  viterbi_mat[1,1] <- log(init[clust_i,1]) + log_class_0[1]
  viterbi_mat[2,1] <- log(init[clust_i,2]) + log_class_1[1]

  viterbi_ind_mat <- matrix(NA,2,day_length)

  for (time in 2:day_length){

    tran <- tran_list_clust[[vcovar_vec[time]+1]][[time]]

    viterbi_mat[1,time] <- log_class_0[[time]] +
      max(viterbi_mat[1,time-1] + log(tran[1,1]),
          viterbi_mat[2,time-1] + log(tran[2,1]))

    viterbi_mat[2,time] <- log_class_1[[time]] +
      max(viterbi_mat[1,time-1] + log(tran[1,2]),
          viterbi_mat[2,time-1] + log(tran[2,2]))

    viterbi_ind_mat[1,time] <-  which.max(c(viterbi_mat[1,time-1] + log(tran[1,1]),
                                            viterbi_mat[2,time-1] + log(tran[2,1])))

    viterbi_ind_mat[2,time] <- which.max(c(viterbi_mat[1,time-1] + log(tran[1,2]),
                                           viterbi_mat[2,time-1] + log(tran[2,2])))
  }

  decoded_mc <- c(which.max(viterbi_mat[,time]))
  for(time in day_length:2){
    decoded_mc <- c(viterbi_ind_mat[decoded_mc[1],time],decoded_mc)
  }

  return(decoded_mc-1)
}

#forward algorithm but only uses decoded data, not activity/light
ForwardAlt <- function(post_decode_collapsed,init,tran_list,vcovar_mat,
                       mix_num = nrow(init)){
  alpha_list <- list()

  for (ind in 1:dim(post_decode_collapsed)[2]){
    alpha_array <- array(NA,dim = c(dim(post_decode_collapsed)[1],2,mix_num))

    for (clust_i in 1:mix_num){
      alpha_array[,,clust_i] <- ForwardIndAltC(post_decode_collapsed[,ind],init[clust_i,],tran_list,clust_i-1,vcovar_mat[,ind])
      alpha_list[[ind]] <- alpha_array
    }
  }
  return(alpha_list)
}

#organizes forward output, mostly done in C
Forward <- function(act, light,init,tran_list,
                    emit_act, emit_light,
                    lod_act, lod_light, corr_mat, beta_vec, surv_coef, surv_covar_risk_vec,
                    event_vec, bline_vec, cbline_vec, lintegral_mat,
                    surv_covar, vcovar_mat, lambda_act_mat, lambda_light_mat, tobit, incl_surv,
                    beta_bool, mix_num = dim(emit_act)[3]){

  alpha_list <- list()
  day_length <- dim(act)[1]

  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))

  if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){
    adj_incl_surv <- 0
  } else {
    adj_incl_surv <- 1
  }

  for (ind in 1:dim(act)[2]){
    alpha_array <- array(NA,dim = c(day_length,2,mix_num))
    covar_risk <- surv_covar_risk_vec[ind]

    for (clust_i in 1:mix_num){
      alpha_array[,,clust_i] <- ForwardIndC(act[,ind], light[,ind], init[clust_i,], tran_list, emit_act_week, emit_light_week,
                                            emit_act_weekend, emit_light_weekend,clust_i-1, lod_act, lod_light, corr_mat,
                                            beta_vec,covar_risk,event_vec[ind], bline_vec[ind], cbline_vec[ind], lintegral_mat,
                                            vcovar_mat[,ind],lambda_act_mat, lambda_light_mat, tobit, adj_incl_surv*beta_bool)
      alpha_list[[ind]] <- alpha_array
    }
  }
  return(alpha_list)
}

#organizes backward output, mostly done in C
Backward <- function(act, light, tran_list,
                     emit_act, emit_light,
                     lod_act, lod_light, corr_mat, lintegral_mat, vcovar_mat,
                     lambda_act_mat, lambda_light_mat, tobit,
                     mix_num = dim(emit_act)[3],
                     day_length = dim(act)[1]){

  beta_list <- list()

  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))

  for (ind in 1:dim(act)[2]){
    beta_array <- array(NA,dim = c(day_length,2,mix_num))

    for (clust_i in 1:mix_num){
      beta_array[,,clust_i] <- BackwardIndC(act[,ind], light[,ind], tran_list, emit_act_week, emit_light_week,
                                            emit_act_weekend, emit_light_weekend,clust_i-1, lod_act, lod_light, corr_mat,
                                            lintegral_mat,
                                            vcovar_mat[,ind],lambda_act_mat, lambda_light_mat, tobit)
      beta_list[[ind]] <- beta_array
    }
  }
  return(beta_list)
}
