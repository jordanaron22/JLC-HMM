Params2TranVectorT <- function(re_ind,len,params_tran,period_len){
  return(t(sapply(c(2:(len)),FUN = Params2Tran,params_tran = params_tran,
                  re_ind = re_ind, period_len = period_len)))
}

#similar to below but transposes and compiles output into single matrix
#needed for faster derivations
Params2TranVectorTresid <- function(re_ind,len,params_tran,period_len){
  return(t(sapply(c(1:(len)),FUN = Params2Tran,params_tran = params_tran,
                  re_ind = re_ind, period_len = period_len)))
}

#Calculates list of transition probabilities over all times
TranByTimeVec <- function(re_ind, params_tran,time_vec,period_len){
  return(lapply(time_vec, Params2Tran, params_tran = params_tran,
                re_ind = re_ind, period_len = period_len))
}

Param2TranHelper <- function(p12,p21){
  tran <- matrix(0,2,2)
  tran[1,2] <- expit(p12)
  tran[1,1] <- 1- tran[1,2]
  tran[2,1] <- expit(p21)
  tran[2,2] <- 1 - tran[2,1]
  return(tran)
}

#takes vector of transition parameters and outputs transition matrix given current time and mixture
Params2Tran <- function(params_tran,time,re_ind,period_len){

  param_matrix <- matrix(params_tran[re_ind,],ncol=3,nrow=2, byrow = T)
  tran <- Param2TranHelper(param_matrix[1,1]+param_matrix[1,2]*cos(2*pi*time/period_len)+param_matrix[1,3]*sin(2*pi*time/period_len),
                           param_matrix[2,1]+param_matrix[2,2]*cos(2*pi*time/period_len)+param_matrix[2,3]*sin(2*pi*time/period_len))

  return(tran)
}

#Generates transition matrices across time as one large matrix for faster computation
#Outer index mixture and then by week/weekend
GenTranColVecList <- function(params_tran_array,len,mix_num,vcovar_num,period_len){

  mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()

  for (mixture_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      params_tran <- params_tran_array[,,vcovar_ind]
      if (mix_num == 1) {params_tran <- matrix(params_tran,nrow = 1)}
      vcovar_tran_list[[vcovar_ind]] <- Params2TranVectorT(mixture_ind,len,params_tran,
                                                           period_len = period_len)

    }
    mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
  }

  return(mixture_vcovar_tran_list)
}

#Generates big list of all transition matrices
#Outer index by mixture then by week/weekend then by time
GenTranList <- function(params_tran_array,time_vec,mix_num,vcovar_num,period_len){
  mixture_vcovar_tran_list <- list()
  vcovar_tran_list <- list()

  for (mixture_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){

      params_tran <- params_tran_array[,,vcovar_ind]
      if (mix_num == 1) {params_tran <- matrix(params_tran,nrow = 1)}
      vcovar_tran_list[[vcovar_ind]] <- TranByTimeVec(re_ind = mixture_ind,
                                                      params_tran = params_tran,
                                                      time_vec = time_vec,
                                                      period_len = period_len)

    }
    mixture_vcovar_tran_list[[mixture_ind]] <- vcovar_tran_list
  }

  return(mixture_vcovar_tran_list)
}

#Helper function for transition calculation
#Relies on helper in C but coded in R to organize
CalcTranHelper <- function(act, light, tran_list_mat, emit_act, emit_light,
                           ind_like_vec, alpha, beta, lod_act, lod_light,
                           corr_mat, lintegral_mat, pi_l,vcovar_mat,
                           lambda_act_mat, lambda_light_mat, tobit){

  num_people = dim(act)[2]
  len = dim(act)[1]
  num_re = dim(emit_act)[3]
  mix_num <- num_re

  tran_vals_re_array <- array(NA,c(2,2,len - 1,num_people,num_re))

  #cpp doesnt have 4d arrays so need to organize our 4d array 3d arrays for week/weekend
  emit_act_week <- array(emit_act[,,,1],dim = c(2,2,mix_num))
  emit_light_week <- array(emit_light[,,,1],dim = c(2,2,mix_num))
  emit_act_weekend <- array(emit_act[,,,2],dim = c(2,2,mix_num))
  emit_light_weekend <- array(emit_light[,,,2],dim = c(2,2,mix_num))

  for(init_state in 1:2){
    for(new_state in 1:2){
      for (clust_i in 1:num_re){
        tran_vals_re_array[init_state,new_state,,,clust_i] <- CalcTranHelperC(init_state = init_state-1, new_state = new_state-1,act = act,
                                                          light = light,tran_list_mat = tran_list_mat,
                                                          emit_act_week = emit_act_week,emit_light_week = emit_light_week,
                                                          emit_act_weekend = emit_act_weekend,emit_light_weekend = emit_light_weekend,
                                                          ind_like_vec = ind_like_vec,
                                                          alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                                          corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                                          clust_i = clust_i-1, vcovar_mat = vcovar_mat,
                                                          lambda_act_mat = lambda_act_mat,lambda_light_mat = lambda_light_mat,tobit = tobit)
      }
    }
  }

  return(tran_vals_re_array)
}

#Makes matrix symmetric
Symmetricize <- function(mat){
  mat <- mat + t(mat)
  diag(mat) <- diag(mat)/2
  return(mat)
}

#Calculates gradient and hessian for transition probabilities within transition helper
CalcGradHess <- function(gradient,hessian_vec,cos_part_vec,sin_part_vec,cos_sin_part){
  mix_num <- dim(gradient)[3]
  vcovar_num <- dim(gradient)[4]
  grad_array <- array(0,dim = c(6,mix_num,vcovar_num))
  hess_array <- array(0,dim = c(6,6,mix_num,vcovar_num))

  for(re_ind in 1:mix_num){
    for(vcovar_ind in 1:vcovar_num){

      grad_array[,re_ind,vcovar_ind] <- as.vector(t(gradient[,,re_ind,vcovar_ind]))

      hessian_vec_re <- hessian_vec[,,re_ind,vcovar_ind]
      hess_upper <- matrix(0,3,3)
      hess_lower <- matrix(0,3,3)

      #### HESS 1

      diag(hess_upper) <- c(hessian_vec_re[1,])
      hess_upper[1,2] <- cos_part_vec[1,re_ind,vcovar_ind]
      hess_upper[1,3] <- sin_part_vec[1,re_ind,vcovar_ind]
      hess_upper[2,3] <- cos_sin_part[1,re_ind,vcovar_ind]
      hess_upper <- Symmetricize(hess_upper)

      #### HESS 2

      diag(hess_lower) <- c(hessian_vec_re[2,])
      hess_lower[1,2] <- cos_part_vec[2,re_ind,vcovar_ind]
      hess_lower[1,3] <- sin_part_vec[2,re_ind,vcovar_ind]
      hess_lower[2,3] <- cos_sin_part[2,re_ind,vcovar_ind]
      hess_lower <- Symmetricize(hess_lower)

      #######

      hess_array[1:3,1:3,re_ind,vcovar_ind] <- hess_upper
      hess_array[4:6,4:6,re_ind,vcovar_ind] <- hess_lower
      # params_tran_working[re_ind,] <- params_tran_working[re_ind,] - solve(-hessian_re,-gradient_re)
    }
  }

  return(list(grad_array,hess_array))
}

#heavy lifting of LM for transition
CalcTranCHelper <- function(alpha,beta,act,light,params_tran_array,emit_act,emit_light,
                      corr_mat,pi_l,lod_act,lod_light,lintegral_mat, vcovar_mat,
                      lambda_act_mat, lambda_light_mat, tobit, check_tran,likelihood,
                      period_len, sweights_vec,
                      vcovar_num = dim(params_tran_array)[3]){

  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]

  gradient <- array(0,c(2,3,mix_num,vcovar_num))
  hessian_vec <- array(0,c(2,3,mix_num,vcovar_num))
  cos_part_vec <- array(0,c(2,mix_num,vcovar_num))
  sin_part_vec <- array(0,c(2,mix_num,vcovar_num))
  cos_sin_part <- array(0,c(2,mix_num,vcovar_num))

  # tran_list_mat <- lapply(c(1:mix_num),Params2TranVectorT, len = len, params_tran = params_tran)
  tran_list_mat <- GenTranColVecList(params_tran_array,len,mix_num,vcovar_num,
                                     period_len = period_len)
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))

  tran_vals_re_array <- CalcTranHelper(act = act,
                                       light = light,tran_list_mat = tran_list_mat,
                                       emit_act = emit_act,emit_light = emit_light,
                                       ind_like_vec = ind_like_vec,
                                       alpha = alpha,beta = beta,lod_act = lod_act,lod_light = lod_light,
                                       corr_mat = corr_mat,lintegral_mat = lintegral_mat,pi_l = pi_l,
                                       vcovar_mat = vcovar_mat[-1,],lambda_act_mat, lambda_light_mat, tobit)

  cos_vec <- cos(2*pi*c(2:(len))/period_len)
  sin_vec <- sin(2*pi*c(2:(len))/period_len)

  for (init_state in 1:2){
    for (new_state in 1:2){

      tran_vals <- tran_vals_re_array[init_state,new_state,,,]

      if (mix_num == 1){
        dim(tran_vals) <- c(dim(tran_vals), 1)
      }

      #tran_vals dimensions: time x participant x latent class
      #accounts for sample weights
      weighted_tran_vals <- sweep(
        tran_vals,
        MARGIN = 2,
        STATS = sweights_vec,
        FUN = "*"
      )

      for (ind in 1:length(alpha)){

        vcovar_vec <- vcovar_mat[-1,ind]
        #need to add 1 bc C is base 0
        vcovar_vecR <- vcovar_vec + 1

        for(re_ind in 1:mix_num){
          #calculates transition over week/weekend
          if (vcovar_num == 1){
            tran_mat <- tran_list_mat[[re_ind]][[1]]
          } else if (vcovar_num == 2){
            tran_mat_week <- tran_list_mat[[re_ind]][[1]]
            tran_mat_weekend <- tran_list_mat[[re_ind]][[2]]
            tran_mat <- tran_mat_week * (1-vcovar_vec) + tran_mat_weekend * vcovar_vec
          }

          if(init_state == 1 & new_state == 1){
            #Left these in for debugging
            #moved to putting all transition values into large matrix to vectorize
            # tran_prime <- -tran[1,2]
            # tran_prime_prime <- -tran[1,1] * tran[1,2]
            tran_prime <- -tran_mat[,3]
            tran_prime_prime <- -tran_mat[,3]*tran_mat[,1]

          } else if(init_state == 1 & new_state == 2){
            # tran_prime <- tran[1,1]
            # tran_prime_prime <- -tran[1,1] * tran[1,2]
            tran_prime <- tran_mat[,1]
            tran_prime_prime <- -tran_mat[,3]*tran_mat[,1]

          } else if(init_state == 2 & new_state == 2){
            # tran_prime <- -tran[2,1]
            # tran_prime_prime <- -tran[2,1] * tran[2,2]
            tran_prime <- -tran_mat[,2]
            tran_prime_prime <- -tran_mat[,2] * tran_mat[,4]

          } else if(init_state == 2 & new_state == 1){
            # tran_prime <- tran[2,2]
            # tran_prime_prime <- -tran[2,1] * tran[2,2]
            tran_prime <- tran_mat[,4]
            tran_prime_prime <- -tran_mat[,2] * tran_mat[,4]
          }

          for (vcovar_ind in 1:vcovar_num){

            #grad and hessian calculations
            gradient[init_state,1,re_ind,vcovar_ind] <- gradient[init_state,1,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime*(vcovar_vecR==vcovar_ind))
            gradient[init_state,2,re_ind,vcovar_ind] <- gradient[init_state,2,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            gradient[init_state,3,re_ind,vcovar_ind] <- gradient[init_state,3,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime*sin_vec*(vcovar_vecR==vcovar_ind))

            hessian_vec[init_state,1,re_ind,vcovar_ind] <- hessian_vec[init_state,1,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime_prime*(vcovar_vecR==vcovar_ind))
            hessian_vec[init_state,2,re_ind,vcovar_ind] <- hessian_vec[init_state,2,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec^2*(vcovar_vecR==vcovar_ind))
            hessian_vec[init_state,3,re_ind,vcovar_ind] <- hessian_vec[init_state,3,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec^2*(vcovar_vecR==vcovar_ind))

            cos_part_vec[init_state,re_ind,vcovar_ind] <- cos_part_vec[init_state,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*(vcovar_vecR==vcovar_ind))
            sin_part_vec[init_state,re_ind,vcovar_ind] <- sin_part_vec[init_state,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime_prime*sin_vec*(vcovar_vecR==vcovar_ind))

            cos_sin_part[init_state,re_ind,vcovar_ind] <- cos_sin_part[init_state,re_ind,vcovar_ind] + sum(weighted_tran_vals[,ind,re_ind]*tran_prime_prime*cos_vec*sin_vec*(vcovar_vecR==vcovar_ind))
          }
        }
      }
    }
  }

  grad_hess_list <- CalcGradHess(gradient,hessian_vec,cos_part_vec,sin_part_vec,cos_sin_part)
  grad_array <- grad_hess_list[[1]]
  hess_array <- grad_hess_list[[2]]

  return(list(grad_array,hess_array))
}

#levenberg marquardt
#interpolates btwn newton and gradient descent
LM <- function(grad_array,hess_array,params_tran_array,check_tran,likelihood,pi_l,
               mix_num = dim(params_tran_array)[1],
               vcovar_num = dim(params_tran_array)[3],
               tran_check_context = NULL,
               step_size = .01){
  params_tran_array_new <- params_tran_array
  new_likelihood <- -Inf

  for (re_ind in 1:mix_num){
    for (vcovar_ind in 1:vcovar_num){
      #returns all -1 if non-invertible
      inf_fact <- SolveCatch(hess_array[,,re_ind,vcovar_ind],grad_array[,re_ind,vcovar_ind])

      step_size <- 1
      #runs until step isnt too big
      #increases step size effectively increases hyperparam for gradient descent -> smaller grad step
      while(max(abs(inf_fact)) > 2 | all(inf_fact == -1)){
        step_fact <- matrix(0,6,6)
        diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) * step_size
        if(all(inf_fact == -1)){diag(step_fact) <- diag(hess_array[,,re_ind,vcovar_ind]) + step_size}

        inf_fact <- SolveCatch(hess_array[,,re_ind,vcovar_ind]+step_fact,grad_array[,re_ind,vcovar_ind])
        step_size <- step_size * 10
      }
      params_tran_array_new[re_ind,,vcovar_ind] <- params_tran_array_new[re_ind,,vcovar_ind] - inf_fact
    }
  }

  #Like decrease may happen here
  if (check_tran){
    if (is.null(tran_check_context)){
      stop("tran_check_context is required when check_tran is TRUE")
    }
    ctx <- tran_check_context
    tran_list <- GenTranList(params_tran_array_new,seq_len(ctx$day_length),
                             ctx$mix_num,ctx$vcovar_num,
                             period_len = ctx$period_len)
    alpha <- Forward(act = ctx$act,light = ctx$light,
                     init = ctx$init,tran_list = tran_list,
                     emit_act = ctx$emit_act,emit_light = ctx$emit_light,
                     lod_act = ctx$lod_act, lod_light = ctx$lod_light,
                     corr_mat = ctx$corr_mat, beta_vec = ctx$beta_vec,
                     surv_coef = ctx$surv_coef,
                     surv_covar_risk_vec = ctx$surv_covar_risk_vec,
                     event_vec = ctx$surv_event, bline_vec = ctx$bline_vec,
                     cbline_vec = ctx$cbline_vec,
                     lintegral_mat = ctx$lintegral_mat,
                     surv_covar = ctx$surv_covar, vcovar_mat = ctx$vcovar_mat,
                     lambda_act_mat = ctx$lambda_act_mat,
                     lambda_light_mat = ctx$lambda_light_mat,
                     tobit = ctx$tobit,incl_surv = ctx$incl_surv,
                     beta_bool = ctx$beta_bool, mix_num = ctx$mix_num)
    new_like <- CalcLikelihood(alpha,pi_l,ctx$sweights_vec)

    if (new_like < likelihood){
      # return(list(params_tran_array,alpha))
      print("Tran Like Dec")
    }
  }

  return(params_tran_array_new)
}

#inverts matrix, throws error if singular
SolveCatch <- function(hess,grad) {
  tryCatch(
    {
      solve(hess,grad,tol=1e-50)
    },
    error = function(cond) {
      # message("Non-Invertible Matrix")
      numeric(length(grad))-1
    },
    warning = function(cond) {
      NULL
    },
    finally = {}
  )
}

#turns transition matrices into dataframe for analyzing later
ParamsArray2DF <- function(params_tran_array,period_len,vcovar_num = dim(params_tran_array)[3]){

  tran_df <- data.frame(prob = c(),
                        type = c(),
                        time = c(),
                        age = c(),
                        weekend = c(),
                        mixture = c())

  for (re_ind in 1:dim(params_tran_array)[1]){
    for (vcovar_ind in 1:vcovar_num){

      params_tran <- params_tran_array[,,vcovar_ind]
      if (dim(params_tran_array)[1] == 1) {params_tran <- matrix(params_tran,nrow = 1)}

      tran_mat <- Params2TranVectorTresid(re_ind,period_len,params_tran,
                                          period_len = period_len)
      tosleep <- tran_mat[,3]
      towake <- tran_mat[,2]

      tran_df_working <- data.frame(prob = c(tosleep,towake),
                                    type = rep(c("Falling Asleep", "Waking"),each= period_len),
                                    time = rep(c(1:period_len)/4,2),
                                    weekend = vcovar_ind,
                                    mixture = re_ind)
      tran_df <- rbind(tran_df,tran_df_working)
    }
  }

  return(tran_df)
}

#Calculates prob of being in each mixture, just based on stationary and moderate activity
CalcPiHelper <- function(nu_mat,nu_covar_vec){
  pi_ind <- exp(colSums(nu_mat * nu_covar_vec))/sum(exp(colSums(nu_mat * nu_covar_vec)))
  if (any(is.na(pi_ind))){
    pi_ind[is.na(pi_ind)] <- 1
    pi_ind <- pi_ind/sum(pi_ind)
  }
  return(pi_ind)
}

CalcPi <- function(nu_mat,nu_covar_mat){

  pi_l_new <- matrix(NA,dim(nu_covar_mat)[1],dim(nu_mat)[2])
  for (ind in 1:dim(nu_covar_mat)[1]){
    pi_l_new[ind,] <- CalcPiHelper(nu_mat,nu_covar_mat[ind,])
  }

  return(pi_l_new)
}

#Calculates ordinal logistic regression coef for pi
#Uses LM
CalcNu <- function(nu_mat,re_prob,nu_covar_mat,alpha,sweights_vec,
                   mix_num = ncol(re_prob),
                   num_of_people = nrow(re_prob)){

  if (ncol(re_prob) == 1){
    return(nu_mat)
  }

  if (length(sweights_vec) != num_of_people){
    stop("sweights_vec must have one value per participant")
  }

  old_mlike <- CalcLikelihood(
    alpha,
    CalcPi(nu_mat,nu_covar_mat),
    sweights_vec
  )

  mlike_diff <- -1
  nnu_covar <- ncol(nu_covar_mat)
  gradient_nu <- numeric(mix_num * nnu_covar)
  hess_nu <- matrix(0,mix_num*nnu_covar,mix_num*nnu_covar)

  for (ind in seq_len(num_of_people)){
    age_ind_vec <- nu_covar_mat[ind,]
    age_ind_mat <- age_ind_vec %*% t(age_ind_vec)

    w_i <- sweights_vec[ind]

    #subtracts by max to avoid overflow
    #doesnt change probabilities
    eta <- colSums(nu_mat * nu_covar_mat[ind,])
    eta <- eta - max(eta)

    num <- exp(eta)
    p_vec <- num / sum(num)

    pvec_grad <- re_prob[ind,] - p_vec
    pvec_grad[1] <- 0

    gradient_nu <- gradient_nu +
      w_i * (pvec_grad %x% age_ind_vec)

    p_vec[1] <- 0
    p_mat <- p_vec %x% t(p_vec)
    diag(p_mat) <- -p_vec * (1-p_vec)

    hess_nu <- hess_nu +
      w_i * (p_mat %x% age_ind_mat)
  }

  reference_indices <- seq_len(nnu_covar)

  hess_nu <- hess_nu[
    -reference_indices,
    -reference_indices,
    drop = FALSE
  ]

  gradient_nu <- gradient_nu[-reference_indices]

  inf_mat <- numeric(length(gradient_nu)) - 1

  step_size <- .01

  while (mlike_diff < 0 || all(inf_mat == -1)){
    step_fact <- matrix(0,nrow(hess_nu),ncol(hess_nu))
    diag(step_fact) <- step_size - .01

    #this was accidentally + previously
    inf_mat <- SolveCatch(
      hess_nu - step_fact,
      gradient_nu
    )

    step_size <- step_size * 10

    nu_mat_lm <- cbind(
      rep(0,nnu_covar),
      matrix(inf_mat,nrow = nnu_covar,byrow = FALSE)
    )

    nu_mat_new <- nu_mat - nu_mat_lm
    pi_l_new <- CalcPi(nu_mat_new,nu_covar_mat)

    new_mlike <- CalcLikelihood(
      alpha,
      pi_l_new,
      sweights_vec
    )

    mlike_diff <- new_mlike - old_mlike
  }

  return(nu_mat_new)
}

CalcTranCHelperFast <- function(
  alpha,
  beta,
  act,
  light,
  params_tran_array,
  emit_act,
  emit_light,
  corr_mat,
  pi_l,
  lod_act,
  lod_light,
  lintegral_mat,
  vcovar_mat,
  lambda_act_mat,
  lambda_light_mat,
  tobit,
  check_tran,
  likelihood,
  period_len,
  sweights_vec,
  vcovar_num = dim(params_tran_array)[3]
) {
  len <- nrow(act)
  mix_num <- dim(emit_act)[3]

  if (nrow(vcovar_mat) != len ||
      ncol(vcovar_mat) != ncol(act)) {
    stop(
      paste(
        "vcovar_mat must have the same",
        "dimensions as act"
      )
    )
  }

  if (length(sweights_vec) != ncol(act)) {
    stop(
      "sweights_vec must contain one value per participant"
    )
  }

  if (vcovar_num < 1L || vcovar_num > 2L) {
    stop(
      paste(
        "CalcTranCHelperFast currently supports",
        "one or two transition day types"
      )
    )
  }

  tran_list_mat <- GenTranColVecList(
    params_tran_array = params_tran_array,
    len = len,
    mix_num = mix_num,
    vcovar_num = vcovar_num,
    period_len = period_len
  )

  ind_like_vec <- vapply(
    seq_along(alpha),
    FUN = IndLike,
    FUN.VALUE = numeric(1),
    alpha = alpha,
    pi_l = pi_l,
    len = len
  )

  emit_act_week <- array(
    emit_act[, , , 1L],
    dim = c(2L, 2L, mix_num)
  )

  emit_light_week <- array(
    emit_light[, , , 1L],
    dim = c(2L, 2L, mix_num)
  )

  if (vcovar_num == 2L) {
    emit_act_weekend <- array(
      emit_act[, , , 2L],
      dim = c(2L, 2L, mix_num)
    )

    emit_light_weekend <- array(
      emit_light[, , , 2L],
      dim = c(2L, 2L, mix_num)
    )
  } else {
    # The second cubes are not accessed when vcovar_num == 1,
    # but valid objects are still required by the C++ signature.
    emit_act_weekend <- emit_act_week
    emit_light_weekend <- emit_light_week
  }

  CalcTranGradHessFastC(
    act = act,
    light = light,
    tran_list_mat = tran_list_mat,
    emit_act_week = emit_act_week,
    emit_light_week = emit_light_week,
    emit_act_weekend = emit_act_weekend,
    emit_light_weekend = emit_light_weekend,
    ind_like_vec = ind_like_vec,
    alpha = alpha,
    beta = beta,
    lod_act = lod_act,
    lod_light = lod_light,
    corr_mat = corr_mat,
    lintegral_mat = lintegral_mat,
    pi_l = pi_l,
    vcovar_mat = vcovar_mat,
    sweights_vec = sweights_vec,
    lambda_act_mat = lambda_act_mat,
    lambda_light_mat = lambda_light_mat,
    tobit = tobit,
    vcovar_num = vcovar_num,
    period_len = period_len
  )
}