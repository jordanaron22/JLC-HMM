# Barebones Oakes H1 helpers.
#
# The Hessians returned here are Hessians of Q, so they are the objects that go
# into H1.

BuildThetaMapRows <- function(block,param_names,start_index,
                              vcovar_ind = NA_integer_,
                              re_ind = NA_integer_,
                              state_ind = NA_integer_,
                              covar_ind = NA_integer_,
                              param_ind = seq_along(param_names)){
  n <- length(param_names)
  if (n == 0){
    return(data.frame(index = integer(0),
                      block = character(0),
                      vcovar_ind = integer(0),
                      re_ind = integer(0),
                      state_ind = integer(0),
                      covar_ind = integer(0),
                      param_ind = integer(0),
                      param_name = character(0)))
  }

  expand_field <- function(x,field_name){
    if (length(x) == 1){
      return(rep(x,n))
    }
    if (length(x) == n){
      return(x)
    }
    stop(paste(field_name,"must have length 1 or",n))
  }

  data.frame(index = start_index + seq_len(n) - 1,
             block = rep(block,n),
             vcovar_ind = expand_field(vcovar_ind,"vcovar_ind"),
             re_ind = expand_field(re_ind,"re_ind"),
             state_ind = expand_field(state_ind,"state_ind"),
             covar_ind = expand_field(covar_ind,"covar_ind"),
             param_ind = expand_field(param_ind,"param_ind"),
             param_name = param_names)
}

PackOakesTheta <- function(init,params_tran_array,emit_act,emit_light,
                           corr_mat,nu_mat,beta_vec,surv_coef,vcovar_mat,
                           sleep_state = SLEEP_STATE,
                           fit_mix_num = length(beta_vec),
                           tran_vcovar_num = dim(params_tran_array)[3],
                           emit_vcovar_num = dim(emit_act)[4]){
  theta <- numeric(0)
  map_list <- list()

  add_block <- function(values,map_rows){
    force(map_rows)
    values <- as.numeric(values)
    theta <<- c(theta,values)
    map_list[[length(map_list) + 1]] <<- map_rows
  }

  next_index <- function(){
    length(theta) + 1
  }

  mix_num <- nrow(init)
  state_num <- dim(emit_act)[1]
  tran_param_num <- dim(params_tran_array)[2]
  nu_covar_num <- nrow(nu_mat)

  active_tran_vcovar_inds <- ActiveTransitionDayTypes(vcovar_mat,
                                                      tran_vcovar_num)
  active_emit_vcovar_inds <- ActiveEmissionDayTypes(vcovar_mat,
                                                    emit_vcovar_num)

  init_values <- init[,sleep_state]
  add_block(init_values,
            BuildThetaMapRows("initial",
                              rep("p_sleep",length(init_values)),
                              next_index(),
                              re_ind = seq_along(init_values),
                              param_ind = rep(sleep_state,
                                              length(init_values))))

  for (vcovar_ind in active_tran_vcovar_inds){
    for (re_ind in seq_len(mix_num)){
      tran_values <- params_tran_array[re_ind,,vcovar_ind]
      add_block(tran_values,
                BuildThetaMapRows("transition",
                                  paste0("tran_",seq_len(tran_param_num)),
                                  next_index(),
                                  vcovar_ind = vcovar_ind,
                                  re_ind = re_ind,
                                  param_ind = seq_len(tran_param_num)))
    }
  }

  emit_param_names <- c("mu_act","sig_act","mu_light","sig_light","corr")
  for (vcovar_ind in active_emit_vcovar_inds){
    for (re_ind in seq_len(mix_num)){
      for (state_ind in seq_len(state_num)){
        emit_values <- EmissionPsi(emit_act,emit_light,corr_mat,state_ind,
                                   re_ind,vcovar_ind)
        add_block(emit_values,
                  BuildThetaMapRows("emission",
                                    emit_param_names,
                                    next_index(),
                                    vcovar_ind = vcovar_ind,
                                    re_ind = re_ind,
                                    state_ind = state_ind,
                                    param_ind = seq_along(emit_param_names)))
      }
    }
  }

  if (ncol(nu_mat) > 1){
    for (re_ind in 2:ncol(nu_mat)){
      mix_values <- nu_mat[,re_ind]
      add_block(mix_values,
                BuildThetaMapRows("mixing",
                                  paste0("nu_",seq_len(nu_covar_num)),
                                  next_index(),
                                  re_ind = re_ind,
                                  covar_ind = seq_len(nu_covar_num),
                                  param_ind = seq_len(nu_covar_num)))
    }
  }

  surv_values <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
  class_names <- if (fit_mix_num > 1){
    paste0("class_",2:fit_mix_num)
  } else {
    character(0)
  }
  surv_covar_param_num <- length(surv_values) - fit_mix_num
  surv_covar_names <- if (surv_covar_param_num > 0){
    paste0("surv_covar_",seq_len(surv_covar_param_num))
  } else {
    character(0)
  }
  surv_names <- c("age",class_names,surv_covar_names)
  add_block(surv_values,
            BuildThetaMapRows("survival",
                              surv_names,
                              next_index(),
                              param_ind = seq_along(surv_values)))

  parameter_map <- do.call(rbind,map_list)
  rownames(parameter_map) <- NULL
  parameter_map$index <- seq_len(nrow(parameter_map))
  names(theta) <- paste0(parameter_map$block,".",parameter_map$param_name,
                         ".",parameter_map$index)

  list(theta = theta,
       parameter_map = parameter_map,
       active_tran_vcovar_inds = active_tran_vcovar_inds,
       active_emit_vcovar_inds = active_emit_vcovar_inds,
       sleep_state = sleep_state,
       fit_mix_num = fit_mix_num,
       surv_coef_len = unlist(lapply(surv_coef,length)),
       templates = list(init = init,
                        params_tran_array = params_tran_array,
                        emit_act = emit_act,
                        emit_light = emit_light,
                        corr_mat = corr_mat,
                        nu_mat = nu_mat,
                        beta_vec = beta_vec,
                        surv_coef = surv_coef))
}

UnpackOakesTheta <- function(theta,theta_pack){
  if (length(theta) != nrow(theta_pack$parameter_map)){
    stop("theta length must match theta_pack$parameter_map")
  }

  init <- theta_pack$templates$init
  params_tran_array <- theta_pack$templates$params_tran_array
  emit_act <- theta_pack$templates$emit_act
  emit_light <- theta_pack$templates$emit_light
  corr_mat <- theta_pack$templates$corr_mat
  nu_mat <- theta_pack$templates$nu_mat

  parameter_map <- theta_pack$parameter_map
  for (i in seq_along(theta)){
    map_row <- parameter_map[i,]
    value <- theta[i]

    if (map_row$block == "initial"){
      init[map_row$re_ind,theta_pack$sleep_state] <- value
      init[map_row$re_ind,-theta_pack$sleep_state] <- 1 - value
    } else if (map_row$block == "transition"){
      params_tran_array[map_row$re_ind,map_row$param_ind,
                        map_row$vcovar_ind] <- value
    } else if (map_row$block == "emission"){
      if (map_row$param_name == "mu_act"){
        emit_act[map_row$state_ind,1,map_row$re_ind,map_row$vcovar_ind] <-
          value
      } else if (map_row$param_name == "sig_act"){
        emit_act[map_row$state_ind,2,map_row$re_ind,map_row$vcovar_ind] <-
          value
      } else if (map_row$param_name == "mu_light"){
        emit_light[map_row$state_ind,1,map_row$re_ind,map_row$vcovar_ind] <-
          value
      } else if (map_row$param_name == "sig_light"){
        emit_light[map_row$state_ind,2,map_row$re_ind,map_row$vcovar_ind] <-
          value
      } else if (map_row$param_name == "corr"){
        corr_mat[map_row$re_ind,map_row$state_ind,map_row$vcovar_ind] <-
          value
      }
    } else if (map_row$block == "mixing"){
      nu_mat[map_row$covar_ind,map_row$re_ind] <- value
    }
  }

  survival_theta <- theta[parameter_map$block == "survival"]
  survival_params <- OutofBetaSurvCoef(survival_theta,
                                       theta_pack$surv_coef_len,
                                       theta_pack$fit_mix_num)

  list(init = init,
       params_tran_array = params_tran_array,
       emit_act = emit_act,
       emit_light = emit_light,
       corr_mat = corr_mat,
       nu_mat = nu_mat,
       beta_vec = survival_params[[1]],
       surv_coef = survival_params[[2]])
}

MakeInactiveDayTypesSafe <- function(params,vcovar_mat){
  safe_params <- params

  active_emit <- ActiveEmissionDayTypes(vcovar_mat,dim(safe_params$emit_act)[4])
  inactive_emit <- setdiff(seq_len(dim(safe_params$emit_act)[4]),active_emit)
  if (length(inactive_emit) > 0){
    source_emit <- active_emit[[1]]
    for (vcovar_ind in inactive_emit){
      safe_params$emit_act[,,,vcovar_ind] <-
        safe_params$emit_act[,,,source_emit]
      safe_params$emit_light[,,,vcovar_ind] <-
        safe_params$emit_light[,,,source_emit]
      safe_params$corr_mat[,,vcovar_ind] <-
        safe_params$corr_mat[,,source_emit]
    }
  }

  active_tran <- ActiveTransitionDayTypes(vcovar_mat,
                                          dim(safe_params$params_tran_array)[3])
  inactive_tran <- setdiff(seq_len(dim(safe_params$params_tran_array)[3]),
                           active_tran)
  if (length(inactive_tran) > 0){
    source_tran <- active_tran[[1]]
    for (vcovar_ind in inactive_tran){
      safe_params$params_tran_array[,,vcovar_ind] <-
        safe_params$params_tran_array[,,source_tran]
    }
  }

  safe_params
}

GetOakesContextValue <- function(posterior_context,name,default = NULL,
                                 required = TRUE){
  if (!is.null(posterior_context[[name]])){
    return(posterior_context[[name]])
  }
  if (!is.null(default) || !required){
    return(default)
  }
  stop(paste("posterior_context is missing required field:",name))
}

CalcOakesTransitionPosteriorWeights <- function(alpha,beta,pi_l,act,light,
                                                params_tran_array,
                                                emit_act,emit_light,corr_mat,
                                                lod_act,lod_light,
                                                lintegral_mat,vcovar_mat,
                                                lambda_act_mat,
                                                lambda_light_mat,tobit,
                                                period_len,
                                                vcovar_num =
                                                  dim(params_tran_array)[3]){
  len <- dim(act)[1]
  mix_num <- dim(emit_act)[3]
  tran_list_mat <- GenTranColVecList(params_tran_array,len,mix_num,
                                     vcovar_num,period_len = period_len)
  ind_like_vec <- unlist(lapply(seq_along(alpha),IndLike,alpha = alpha,
                                pi_l = pi_l,len = len))

  CalcTranHelper(act = act,
                 light = light,
                 tran_list_mat = tran_list_mat,
                 emit_act = emit_act,
                 emit_light = emit_light,
                 ind_like_vec = ind_like_vec,
                 alpha = alpha,
                 beta = beta,
                 lod_act = lod_act,
                 lod_light = lod_light,
                 corr_mat = corr_mat,
                 lintegral_mat = lintegral_mat,
                 pi_l = pi_l,
                 vcovar_mat = vcovar_mat[-1,,drop = FALSE],
                 lambda_act_mat = lambda_act_mat,
                 lambda_light_mat = lambda_light_mat,
                 tobit = tobit)
}

BuildOakesPosteriorContext <- function(alpha,beta,pi_l,log_sweights_vec,
                                       act,light,vcovar_mat,
                                       params_tran_array,
                                       emit_act,emit_light,corr_mat,
                                       lod_act,lod_light,
                                       lambda_act_mat,lambda_light_mat,tobit,
                                       period_len,lintegral_mat = NULL,
                                       nu_covar_mat = NULL,
                                       survival_context = NULL,
                                       re_prob = CalcProbRE(alpha,pi_l),
                                       bline_vec = NULL,
                                       cbline_vec = NULL,
                                       survival_baseline_mode = "fixed",
                                       surv_covar_risk_vec = NULL,
                                       profile_iterations = NA_integer_,
                                       profile_error = NA_real_){
  if (is.null(lintegral_mat)){
    lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,
                                      lod_act,lod_light)
  }

  if (!is.null(survival_context)){
    if (exists("update_survival_context_re_prob",mode = "function")){
      survival_context <- update_survival_context_re_prob(survival_context,
                                                          re_prob,ncol(pi_l))
    } else {
      survival_context$re_prob <- re_prob
    }
    if (!is.null(bline_vec)){
      survival_context$bline_vec <- bline_vec
    }
    if (!is.null(cbline_vec)){
      survival_context$cbline_vec <- cbline_vec
    }
  }

  list(alpha = alpha,
       beta = beta,
       pi_l = pi_l,
       log_sweights_vec = log_sweights_vec,
       init_counts = CalcOakesInitialCounts(alpha,beta,pi_l,
                                            log_sweights_vec),
       transition_weights =
         CalcOakesTransitionPosteriorWeights(alpha = alpha,
                                             beta = beta,
                                             pi_l = pi_l,
                                             act = act,
                                             light = light,
                                             params_tran_array =
                                               params_tran_array,
                                             emit_act = emit_act,
                                             emit_light = emit_light,
                                             corr_mat = corr_mat,
                                             lod_act = lod_act,
                                             lod_light = lod_light,
                                             lintegral_mat = lintegral_mat,
                                             vcovar_mat = vcovar_mat,
                                             lambda_act_mat =
                                               lambda_act_mat,
                                             lambda_light_mat =
                                               lambda_light_mat,
                                             tobit = tobit,
                                             period_len = period_len),
       re_prob = re_prob,
       nu_covar_mat = nu_covar_mat,
       survival_context = survival_context,
       bline_vec = bline_vec,
       cbline_vec = cbline_vec,
       survival_baseline_mode = survival_baseline_mode,
       surv_covar_risk_vec = surv_covar_risk_vec,
       profile_iterations = profile_iterations,
       profile_error = profile_error,
       act = act,
       light = light,
       vcovar_mat = vcovar_mat,
       lod_act = lod_act,
       lod_light = lod_light,
       lintegral_mat = lintegral_mat,
       lambda_act_mat = lambda_act_mat,
       lambda_light_mat = lambda_light_mat,
       tobit = tobit,
       period_len = period_len,
       params_tran_array = params_tran_array,
       emit_act = emit_act,
       emit_light = emit_light,
       corr_mat = corr_mat)
}

NormalizeLogRows <- function(log_mat){
  prob_mat <- matrix(NA_real_,nrow = nrow(log_mat),ncol = ncol(log_mat))

  for (i in seq_len(nrow(log_mat))){
    row_denom <- logSumExp(log_mat[i,])
    prob_mat[i,] <- exp(log_mat[i,] - row_denom)
  }

  prob_mat
}

CalcLongitudinalClassLogLik <- function(alpha_long,pi_l){
  len <- dim(alpha_long[[1]])[1]
  mix_num <- ncol(pi_l)
  log_like_mat <- matrix(NA_real_,nrow = length(alpha_long),ncol = mix_num)

  for (ind in seq_along(alpha_long)){
    for (re_ind in seq_len(mix_num)){
      log_like_mat[ind,re_ind] <-
        logSumExp(alpha_long[[ind]][len,,re_ind]) + log(pi_l[ind,re_ind])
    }
  }

  log_like_mat
}

CalcSurvivalClassLogLik <- function(beta_vec,surv_covar_risk_vec,
                                    surv_event,bline_vec,cbline_vec){
  n <- length(surv_event)
  mix_num <- length(beta_vec)
  eta_mat <- matrix(surv_covar_risk_vec,nrow = n,ncol = mix_num) +
    matrix(beta_vec,nrow = n,ncol = mix_num,byrow = TRUE)

  log_bline_event <- numeric(n)
  event_rows <- surv_event == 1
  log_bline_event[event_rows] <- log(bline_vec[event_rows])

  matrix(log_bline_event,nrow = n,ncol = mix_num) +
    matrix(surv_event,nrow = n,ncol = mix_num) * eta_mat -
    matrix(cbline_vec,nrow = n,ncol = mix_num) * exp(eta_mat)
}

AddClassLogOffsetsToAlpha <- function(alpha,log_offset_mat){
  alpha_offset <- alpha

  for (ind in seq_along(alpha_offset)){
    for (re_ind in seq_len(dim(alpha_offset[[ind]])[3])){
      alpha_offset[[ind]][,,re_ind] <-
        alpha_offset[[ind]][,,re_ind] + log_offset_mat[ind,re_ind]
    }
  }

  alpha_offset
}

RebuildProfiledPosteriorContext <- function(theta,theta_pack,data_context,
                                            make_inactive_day_types_safe =
                                              TRUE,
                                            maxit = 100L,
                                            tol = 1e-9,
                                            damping = 1.0){
  maxit <- as.integer(maxit)
  if (length(maxit) != 1 || is.na(maxit) || maxit < 1L){
    stop("maxit must be a positive integer")
  }
  if (length(tol) != 1 || is.na(tol) || tol < 0){
    stop("tol must be a non-negative scalar")
  }
  if (length(damping) != 1 || is.na(damping) ||
      damping < 0 || damping > 1){
    stop("damping must be in [0, 1]")
  }

  params <- UnpackOakesTheta(theta,theta_pack)
  fb_params <- params
  if (make_inactive_day_types_safe){
    fb_params <- MakeInactiveDayTypesSafe(fb_params,data_context$vcovar_mat)
  }

  survival_context <- data_context$survival_context
  pi_l <- CalcPi(params$nu_mat,data_context$nu_covar_mat)
  surv_covar_risk_vec <- SurvCovarRiskVec(survival_context$surv_covar,
                                          params$surv_coef)

  lintegral_mat <- CalcLintegralMat(fb_params$emit_act,fb_params$emit_light,
                                    fb_params$corr_mat,
                                    data_context$lod_act,
                                    data_context$lod_light)
  tran_list <- GenTranList(fb_params$params_tran_array,
                           seq_len(nrow(data_context$act)),
                           ncol(pi_l),
                           dim(fb_params$params_tran_array)[3],
                           period_len = data_context$period_len)

  alpha_long <- Forward(act = data_context$act,
                        light = data_context$light,
                        init = params$init,
                        tran_list = tran_list,
                        emit_act = fb_params$emit_act,
                        emit_light = fb_params$emit_light,
                        lod_act = data_context$lod_act,
                        lod_light = data_context$lod_light,
                        corr_mat = fb_params$corr_mat,
                        beta_vec = params$beta_vec,
                        surv_coef = params$surv_coef,
                        surv_covar_risk_vec = surv_covar_risk_vec,
                        event_vec = survival_context$surv_event,
                        bline_vec = rep(1,length(survival_context$surv_event)),
                        cbline_vec = rep(0,length(survival_context$surv_event)),
                        lintegral_mat = lintegral_mat,
                        log_sweights_vec = data_context$log_sweights_vec,
                        surv_covar = survival_context$surv_covar,
                        vcovar_mat = data_context$vcovar_mat,
                        lambda_act_mat = data_context$lambda_act_mat,
                        lambda_light_mat = data_context$lambda_light_mat,
                        tobit = data_context$tobit,
                        incl_surv = MODEL_TYPE_CODES[["two_stage"]],
                        beta_bool = 0,
                        mix_num = ncol(pi_l))

  log_long_mat <- CalcLongitudinalClassLogLik(alpha_long,pi_l)

  beta <- Backward(act = data_context$act,
                   light = data_context$light,
                   tran_list = tran_list,
                   emit_act = fb_params$emit_act,
                   emit_light = fb_params$emit_light,
                   lod_act = data_context$lod_act,
                   lod_light = data_context$lod_light,
                   corr_mat = fb_params$corr_mat,
                   lintegral_mat = lintegral_mat,
                   vcovar_mat = data_context$vcovar_mat,
                   lambda_act_mat = data_context$lambda_act_mat,
                   lambda_light_mat = data_context$lambda_light_mat,
                   tobit = data_context$tobit,
                   mix_num = ncol(pi_l),
                   day_length = nrow(data_context$act))

  re_prob_cur <- data_context$re_prob
  if (is.null(re_prob_cur) || is.null(dim(re_prob_cur)) ||
      any(dim(re_prob_cur) != dim(pi_l)) ||
      any(!is.finite(re_prob_cur))){
    re_prob_cur <- pi_l
  }

  profile_error <- NA_real_
  alpha <- NULL
  bline_vec <- NULL
  cbline_vec <- NULL
  log_surv_mat <- NULL

  for (iter in seq_len(maxit)){
    bhaz_vec <- CalcBLHaz(params$surv_coef,params$beta_vec,
                          re_prob_cur,surv_covar_risk_vec,
                          survival_context$surv_event,
                          survival_context$surv_time,
                          survival_context$surv_covar)
    bline_vec <- bhaz_vec[[1]]
    cbline_vec <- bhaz_vec[[2]]

    log_surv_mat <- CalcSurvivalClassLogLik(
      beta_vec = params$beta_vec,
      surv_covar_risk_vec = surv_covar_risk_vec,
      surv_event = survival_context$surv_event,
      bline_vec = bline_vec,
      cbline_vec = cbline_vec)
    re_prob_new <- NormalizeLogRows(log_long_mat + log_surv_mat)
    profile_error <- max(abs(re_prob_new - re_prob_cur),na.rm = TRUE)

    if (profile_error < tol || iter == maxit){
      alpha <- AddClassLogOffsetsToAlpha(alpha_long,log_surv_mat)

      return(BuildOakesPosteriorContext(
        alpha = alpha,
        beta = beta,
        pi_l = pi_l,
        log_sweights_vec = data_context$log_sweights_vec,
        act = data_context$act,
        light = data_context$light,
        vcovar_mat = data_context$vcovar_mat,
        params_tran_array = fb_params$params_tran_array,
        emit_act = fb_params$emit_act,
        emit_light = fb_params$emit_light,
        corr_mat = fb_params$corr_mat,
        lod_act = data_context$lod_act,
        lod_light = data_context$lod_light,
        lambda_act_mat = data_context$lambda_act_mat,
        lambda_light_mat = data_context$lambda_light_mat,
        tobit = data_context$tobit,
        period_len = data_context$period_len,
        lintegral_mat = lintegral_mat,
        nu_covar_mat = data_context$nu_covar_mat,
        survival_context = survival_context,
        re_prob = re_prob_new,
        bline_vec = bline_vec,
        cbline_vec = cbline_vec,
        survival_baseline_mode = "profiled",
        surv_covar_risk_vec = surv_covar_risk_vec,
        profile_iterations = iter,
        profile_error = profile_error))
    }

    re_prob_cur <- damping * re_prob_new + (1 - damping) * re_prob_cur
  }
}

RebuildOakesPosteriorContext <- function(theta,theta_pack,data_context,
                                         make_inactive_day_types_safe = TRUE){
  survival_baseline_mode <- data_context$survival_baseline_mode
  if (is.null(survival_baseline_mode)){
    survival_baseline_mode <- "fixed"
  }

  if (identical(survival_baseline_mode,"profiled")){
    profile_maxit <- data_context$profile_maxit
    if (is.null(profile_maxit)){
      profile_maxit <- 100L
    }
    profile_tol <- data_context$profile_tol
    if (is.null(profile_tol)){
      profile_tol <- 1e-9
    }
    profile_damping <- data_context$profile_damping
    if (is.null(profile_damping)){
      profile_damping <- 1.0
    }
    return(RebuildProfiledPosteriorContext(
      theta = theta,
      theta_pack = theta_pack,
      data_context = data_context,
      make_inactive_day_types_safe = make_inactive_day_types_safe,
      maxit = profile_maxit,
      tol = profile_tol,
      damping = profile_damping))
  }

  if (!identical(survival_baseline_mode,"fixed")){
    stop(paste("Unknown survival_baseline_mode:",
               survival_baseline_mode))
  }

  params <- UnpackOakesTheta(theta,theta_pack)
  fb_params <- params
  if (make_inactive_day_types_safe){
    fb_params <- MakeInactiveDayTypesSafe(fb_params,data_context$vcovar_mat)
  }

  pi_l <- CalcPi(params$nu_mat,data_context$nu_covar_mat)
  surv_covar_risk_vec <- SurvCovarRiskVec(data_context$survival_context$surv_covar,
                                          params$surv_coef)

  lintegral_mat <- CalcLintegralMat(fb_params$emit_act,fb_params$emit_light,
                                    fb_params$corr_mat,
                                    data_context$lod_act,
                                    data_context$lod_light)
  tran_list <- GenTranList(fb_params$params_tran_array,
                           seq_len(nrow(data_context$act)),
                           ncol(pi_l),
                           dim(fb_params$params_tran_array)[3],
                           period_len = data_context$period_len)

  survival_context <- data_context$survival_context
  bline_vec <- data_context$bline_vec
  cbline_vec <- data_context$cbline_vec
  if (is.null(bline_vec) || is.null(cbline_vec)){
    stop("survival_baseline_mode = 'fixed' requires bline_vec and cbline_vec")
  }

  alpha <- Forward(act = data_context$act,
                   light = data_context$light,
                   init = params$init,
                   tran_list = tran_list,
                   emit_act = fb_params$emit_act,
                   emit_light = fb_params$emit_light,
                   lod_act = data_context$lod_act,
                   lod_light = data_context$lod_light,
                   corr_mat = fb_params$corr_mat,
                   beta_vec = params$beta_vec,
                   surv_coef = params$surv_coef,
                   surv_covar_risk_vec = surv_covar_risk_vec,
                   event_vec = survival_context$surv_event,
                   bline_vec = bline_vec,
                   cbline_vec = cbline_vec,
                   lintegral_mat = lintegral_mat,
                   log_sweights_vec = data_context$log_sweights_vec,
                   surv_covar = survival_context$surv_covar,
                   vcovar_mat = data_context$vcovar_mat,
                   lambda_act_mat = data_context$lambda_act_mat,
                   lambda_light_mat = data_context$lambda_light_mat,
                   tobit = data_context$tobit,
                   incl_surv = if (is.null(data_context$incl_surv)){
                     MODEL_TYPE_CODES[["joint"]]
                   } else {
                     data_context$incl_surv
                   },
                   beta_bool = if (is.null(data_context$beta_bool)){
                     1
                   } else {
                     data_context$beta_bool
                   },
                   mix_num = ncol(pi_l))

  beta <- Backward(act = data_context$act,
                   light = data_context$light,
                   tran_list = tran_list,
                   emit_act = fb_params$emit_act,
                   emit_light = fb_params$emit_light,
                   lod_act = data_context$lod_act,
                   lod_light = data_context$lod_light,
                   corr_mat = fb_params$corr_mat,
                   lintegral_mat = lintegral_mat,
                   vcovar_mat = data_context$vcovar_mat,
                   lambda_act_mat = data_context$lambda_act_mat,
                   lambda_light_mat = data_context$lambda_light_mat,
                   tobit = data_context$tobit,
                   mix_num = ncol(pi_l),
                   day_length = nrow(data_context$act))

  re_prob <- CalcProbRE(alpha,pi_l)

  BuildOakesPosteriorContext(alpha = alpha,
                             beta = beta,
                             pi_l = pi_l,
                             log_sweights_vec =
                               data_context$log_sweights_vec,
                             act = data_context$act,
                             light = data_context$light,
                             vcovar_mat = data_context$vcovar_mat,
                             params_tran_array =
                               fb_params$params_tran_array,
                             emit_act = fb_params$emit_act,
                             emit_light = fb_params$emit_light,
                             corr_mat = fb_params$corr_mat,
                             lod_act = data_context$lod_act,
                             lod_light = data_context$lod_light,
                             lambda_act_mat = data_context$lambda_act_mat,
                             lambda_light_mat =
                               data_context$lambda_light_mat,
                             tobit = data_context$tobit,
                             period_len = data_context$period_len,
                             lintegral_mat = lintegral_mat,
                             nu_covar_mat = data_context$nu_covar_mat,
                             survival_context = survival_context,
                             re_prob = re_prob,
                             bline_vec = bline_vec,
                             cbline_vec = cbline_vec,
                             survival_baseline_mode =
                               survival_baseline_mode,
                             surv_covar_risk_vec = surv_covar_risk_vec,
                             profile_iterations = NA_integer_,
                             profile_error = NA_real_)
}

DiagnosePosteriorContext <- function(posterior_context){
  alpha_vec <- unlist(posterior_context$alpha)
  beta_vec <- unlist(posterior_context$beta)
  out <- data.frame(quantity = c("alpha","beta"),
                    length = c(length(alpha_vec),length(beta_vec)),
                    nonfinite = c(sum(!is.finite(alpha_vec)),
                                  sum(!is.finite(beta_vec))),
                    missing = c(sum(is.na(alpha_vec)),sum(is.na(beta_vec))))

  if (!is.null(posterior_context$transition_weights)){
    transition_weights <- posterior_context$transition_weights
    out <- rbind(out,
                 data.frame(quantity = "transition_weights",
                            length = length(transition_weights),
                            nonfinite = sum(!is.finite(transition_weights)),
                            missing = sum(is.na(transition_weights))))
  }

  if (!is.null(posterior_context$init_counts)){
    init_counts <- posterior_context$init_counts
    out <- rbind(out,
                 data.frame(quantity = "init_counts",
                            length = length(init_counts),
                            nonfinite = sum(!is.finite(init_counts)),
                            missing = sum(is.na(init_counts))))
  }

  if (!is.null(posterior_context$bline_vec)){
    bline_vec <- posterior_context$bline_vec
    out <- rbind(out,
                 data.frame(quantity = "bline_vec",
                            length = length(bline_vec),
                            nonfinite = sum(!is.finite(bline_vec)),
                            missing = sum(is.na(bline_vec))))
  }

  if (!is.null(posterior_context$cbline_vec)){
    cbline_vec <- posterior_context$cbline_vec
    out <- rbind(out,
                 data.frame(quantity = "cbline_vec",
                            length = length(cbline_vec),
                            nonfinite = sum(!is.finite(cbline_vec)),
                            missing = sum(is.na(cbline_vec))))
  }

  out
}

ValidateOakesThetaNaturalScale <- function(theta,theta_pack,
                                           label = "theta"){
  parameter_map <- theta_pack$parameter_map
  problems <- data.frame(index = integer(0),
                         block = character(0),
                         param_name = character(0),
                         value = numeric(0),
                         problem = character(0))

  add_problem <- function(idx,problem){
    if (length(idx) == 0){
      return()
    }
    problems <<- rbind(
      problems,
      data.frame(index = idx,
                 block = parameter_map$block[idx],
                 param_name = parameter_map$param_name[idx],
                 value = theta[idx],
                 problem = problem)
    )
  }

  init_idx <- which(parameter_map$block == "initial")
  add_problem(init_idx[theta[init_idx] <= 0 | theta[init_idx] >= 1],
              "initial probability must be in (0, 1)")

  sig_idx <- which(parameter_map$block == "emission" &
                    parameter_map$param_name %in% c("sig_act","sig_light"))
  add_problem(sig_idx[theta[sig_idx] <= 0],
              "emission standard deviation must be positive")

  corr_idx <- which(parameter_map$block == "emission" &
                      parameter_map$param_name == "corr")
  add_problem(corr_idx[abs(theta[corr_idx]) >= 1],
              "emission correlation must be in (-1, 1)")

  if (nrow(problems) > 0){
    print(problems)
    stop(paste(label,"has invalid natural-scale parameter values"))
  }

  invisible(TRUE)
}

CalcOakesH2 <- function(theta_pack,data_context,
                        theta = theta_pack$theta,
                        eps = 1e-4,
                        parameter_indices = seq_along(theta),
                        progress = TRUE,
                        make_inactive_day_types_safe = TRUE){
  ValidateOakesThetaNaturalScale(theta,theta_pack,"theta")

  param_num <- length(theta)
  h2 <- matrix(NA_real_,nrow = param_num,ncol = param_num)
  rownames(h2) <- names(theta)
  colnames(h2) <- names(theta)

  diagnostics <- theta_pack$parameter_map[parameter_indices,,drop = FALSE]
  diagnostics$eps <- eps
  diagnostics$plus_posterior_nonfinite <- NA_integer_
  diagnostics$minus_posterior_nonfinite <- NA_integer_
  diagnostics$plus_score_nonfinite <- NA_integer_
  diagnostics$minus_score_nonfinite <- NA_integer_
  diagnostics$column_abs_sum <- NA_real_

  base_context <- RebuildOakesPosteriorContext(
    theta = theta,
    theta_pack = theta_pack,
    data_context = data_context,
    make_inactive_day_types_safe = make_inactive_day_types_safe)
  base_score <- CalcOakesScore(theta = theta,
                               theta_pack = theta_pack,
                               posterior_context = base_context)

  for (k in seq_along(parameter_indices)){
    j <- parameter_indices[[k]]
    map_row <- theta_pack$parameter_map[j,,drop = FALSE]
    if (progress){
      message("H2 parameter ",k,"/",length(parameter_indices),
              " [",j,"] ",map_row$block,":",map_row$param_name)
    }

    theta_plus <- theta
    theta_minus <- theta
    theta_plus[j] <- theta_plus[j] + eps
    theta_minus[j] <- theta_minus[j] - eps

    ValidateOakesThetaNaturalScale(theta_plus,theta_pack,
                                   paste0("theta_plus[",j,"]"))
    ValidateOakesThetaNaturalScale(theta_minus,theta_pack,
                                   paste0("theta_minus[",j,"]"))

    posterior_plus <- RebuildOakesPosteriorContext(
      theta = theta_plus,
      theta_pack = theta_pack,
      data_context = data_context,
      make_inactive_day_types_safe = make_inactive_day_types_safe)
    score_plus <- CalcOakesScore(theta = theta,
                                 theta_pack = theta_pack,
                                 posterior_context = posterior_plus)

    posterior_minus <- RebuildOakesPosteriorContext(
      theta = theta_minus,
      theta_pack = theta_pack,
      data_context = data_context,
      make_inactive_day_types_safe = make_inactive_day_types_safe)
    score_minus <- CalcOakesScore(theta = theta,
                                  theta_pack = theta_pack,
                                  posterior_context = posterior_minus)

    h2[,j] <- (score_plus - score_minus) / (2 * eps)

    plus_post_diag <- DiagnosePosteriorContext(posterior_plus)
    minus_post_diag <- DiagnosePosteriorContext(posterior_minus)
    diagnostics$plus_posterior_nonfinite[k] <-
      sum(plus_post_diag$nonfinite)
    diagnostics$minus_posterior_nonfinite[k] <-
      sum(minus_post_diag$nonfinite)
    diagnostics$plus_score_nonfinite[k] <- sum(!is.finite(score_plus))
    diagnostics$minus_score_nonfinite[k] <- sum(!is.finite(score_minus))
    diagnostics$column_abs_sum[k] <- sum(abs(h2[,j]),na.rm = TRUE)
  }

  list(h2 = h2,
       hessian = h2,
       base_score = base_score,
       base_context = base_context,
       diagnostics = diagnostics,
       eps = eps,
       parameter_indices = parameter_indices,
       parameter_map = theta_pack$parameter_map)
}

CheckOakesH1AgainstScore <- function(theta_pack,posterior_context,h1,
                                     theta = theta_pack$theta,
                                     eps = 1e-5,
                                     parameter_indices = seq_along(theta),
                                     progress = TRUE){
  ValidateOakesThetaNaturalScale(theta,theta_pack,"theta")
  h1 <- as.matrix(h1)
  if (!all(dim(h1) == c(length(theta),length(theta)))){
    stop("h1 dimensions must match length(theta)")
  }

  numerical_h1 <- matrix(NA_real_,nrow = length(theta),ncol = length(theta))
  rownames(numerical_h1) <- names(theta)
  colnames(numerical_h1) <- names(theta)

  diagnostics <- theta_pack$parameter_map[parameter_indices,,drop = FALSE]
  diagnostics$eps <- eps
  diagnostics$max_abs_diff <- NA_real_
  diagnostics$max_abs_h1 <- NA_real_
  diagnostics$max_abs_numeric <- NA_real_
  diagnostics$cor_with_h1 <- NA_real_
  diagnostics$cor_with_negative_h1 <- NA_real_

  for (k in seq_along(parameter_indices)){
    j <- parameter_indices[[k]]
    map_row <- theta_pack$parameter_map[j,,drop = FALSE]
    if (progress){
      message("H1 check parameter ",k,"/",length(parameter_indices),
              " [",j,"] ",map_row$block,":",map_row$param_name)
    }

    theta_plus <- theta
    theta_minus <- theta
    theta_plus[j] <- theta_plus[j] + eps
    theta_minus[j] <- theta_minus[j] - eps

    ValidateOakesThetaNaturalScale(theta_plus,theta_pack,
                                   paste0("theta_plus[",j,"]"))
    ValidateOakesThetaNaturalScale(theta_minus,theta_pack,
                                   paste0("theta_minus[",j,"]"))

    score_plus <- CalcOakesScore(theta = theta_plus,
                                 theta_pack = theta_pack,
                                 posterior_context = posterior_context)
    score_minus <- CalcOakesScore(theta = theta_minus,
                                  theta_pack = theta_pack,
                                  posterior_context = posterior_context)

    numerical_h1[,j] <- (score_plus - score_minus) / (2 * eps)
    diff_col <- numerical_h1[,j] - h1[,j]

    diagnostics$max_abs_diff[k] <- max(abs(diff_col),na.rm = TRUE)
    diagnostics$max_abs_h1[k] <- max(abs(h1[,j]),na.rm = TRUE)
    diagnostics$max_abs_numeric[k] <- max(abs(numerical_h1[,j]),na.rm = TRUE)
    diagnostics$cor_with_h1[k] <- suppressWarnings(
      cor(numerical_h1[,j],h1[,j],use = "complete.obs"))
    diagnostics$cor_with_negative_h1[k] <- suppressWarnings(
      cor(numerical_h1[,j],-h1[,j],use = "complete.obs"))
  }

  list(numerical_h1 = numerical_h1,
       analytic_h1 = h1,
       diagnostics = diagnostics,
       eps = eps,
       parameter_indices = parameter_indices,
       parameter_map = theta_pack$parameter_map)
}

DiagnoseOakesInformation <- function(H1,H2){
  H1 <- as.matrix(H1)
  H2 <- as.matrix(H2)
  H_plus <- H1 + H2
  I_plus <- -(H_plus)
  I_plus_sym <- 0.5 * (I_plus + t(I_plus))

  H_minus <- H1 - H2
  I_minus <- -(H_minus)
  I_minus_sym <- 0.5 * (I_minus + t(I_minus))

  plus_eigen <- eigen(I_plus_sym,symmetric = TRUE,
                      only.values = TRUE)$values
  minus_eigen <- eigen(I_minus_sym,symmetric = TRUE,
                       only.values = TRUE)$values

  data.frame(formula = c("-(H1 + H2)","-(H1 - H2)"),
             max_abs_asymmetry = c(max(abs(I_plus - t(I_plus)),na.rm = TRUE),
                                   max(abs(I_minus - t(I_minus)),
                                       na.rm = TRUE)),
             min_eigen = c(min(plus_eigen),min(minus_eigen)),
             max_eigen = c(max(plus_eigen),max(minus_eigen)),
             nonpositive_eigen_count = c(sum(plus_eigen <= 0),
                                         sum(minus_eigen <= 0)))
}

CalcOakesInitialCounts <- function(alpha,beta,pi_l,log_sweights_vec){
  mix_num <- ncol(pi_l)
  state_num <- dim(alpha[[1]])[2]
  counts <- matrix(0,nrow = mix_num,ncol = state_num)

  if (length(log_sweights_vec) != length(alpha)){
    stop("log_sweights_vec must have one entry per subject")
  }

  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)

  for (ind in seq_along(alpha)){
    ind_like <- ind_like_vec[ind]

    for (g in seq_len(mix_num)){
      log_count <- alpha[[ind]][1,,g] +
        beta[[ind]][1,,g] +
        log(pi_l[ind,g]) -
        ind_like +
        log_sweights_vec[ind]

      counts[g,] <- counts[g,] + exp(log_count)
    }
  }

  counts
}

CalcOakesInitialH1 <- function(init,init_counts,
                               sleep_state = SLEEP_STATE){
  if (ncol(init) != 2){
    stop("CalcOakesInitialH1 currently assumes two Markov states")
  }
  if (!all(dim(init) == dim(init_counts))){
    stop("init and init_counts must have the same dimensions")
  }

  wake_state <- setdiff(seq_len(ncol(init)),sleep_state)
  p_sleep <- init[,sleep_state]

  if (any(p_sleep <= 0 | p_sleep >= 1)){
    stop("Initial sleep probabilities must be strictly between 0 and 1")
  }

  score <- init_counts[,sleep_state] / p_sleep -
    init_counts[,wake_state] / (1 - p_sleep)

  hessian_diag <- -init_counts[,sleep_state] / p_sleep^2 -
    init_counts[,wake_state] / (1 - p_sleep)^2

  list(score = score,
       hessian = diag(hessian_diag,nrow = length(hessian_diag)),
       counts = init_counts)
}

CalcOakesMixingH1 <- function(nu_mat,re_prob,nu_covar_mat){
  mix_num <- ncol(nu_mat)
  nu_covar_num <- nrow(nu_mat)
  num_people <- nrow(nu_covar_mat)

  if (ncol(re_prob) != mix_num){
    stop("re_prob and nu_mat must have the same number of classes")
  }
  if (nrow(re_prob) != num_people){
    stop("re_prob and nu_covar_mat must have the same number of subjects")
  }
  if (ncol(nu_covar_mat) != nu_covar_num){
    stop("nu_covar_mat columns must match nu_mat rows")
  }

  if (mix_num == 1){
    empty_hessian <- matrix(0,nrow = 0,ncol = 0)
    return(list(score = numeric(0),
                hessian = empty_hessian,
                pi_l = matrix(1,nrow = num_people,ncol = 1)))
  }

  full_score <- numeric(mix_num * nu_covar_num)
  full_hessian <- matrix(0,nrow = mix_num * nu_covar_num,
                         ncol = mix_num * nu_covar_num)
  pi_l <- matrix(NA,nrow = num_people,ncol = mix_num)

  for (ind in seq_len(num_people)){
    x_vec <- nu_covar_mat[ind,]
    x_outer <- x_vec %*% t(x_vec)

    eta <- colSums(nu_mat * x_vec)
    #soft max for numerical stability incase of overflow
    eta <- eta - max(eta)
    pi_vec <- exp(eta) / sum(exp(eta))
    pi_l[ind,] <- pi_vec

    class_score <- re_prob[ind,] - pi_vec
    class_score[1] <- 0
    full_score <- full_score + class_score %x% x_vec

    pi_work <- pi_vec
    pi_work[1] <- 0
    class_hessian <- pi_work %x% t(pi_work)
    diag(class_hessian) <- -pi_work * (1 - pi_work)
    full_hessian <- full_hessian + class_hessian %x% x_outer
  }

  reference_idx <- seq_len(nu_covar_num)

  list(score = full_score[-reference_idx],
       hessian = full_hessian[-reference_idx,-reference_idx,drop = FALSE],
       pi_l = pi_l)
}

BuildOakesSurvivalDesign <- function(surv_covar,surv_coef){
  age <- as.numeric(surv_covar[[1]])
  design <- matrix(age,ncol = 1)
  colnames(design) <- "age"

  if (length(surv_covar) > 1){
    for (covar_ind in 2:length(surv_covar)){
      covar_mat <- as.matrix(surv_covar[[covar_ind]])
      keep_cols <- seq_len(length(surv_coef[[covar_ind]]))
      keep_cols <- keep_cols[keep_cols <= ncol(covar_mat)]
      non_ref_cols <- keep_cols[keep_cols > 1]

      if (length(non_ref_cols) > 0){
        design <- cbind(design,
                        covar_mat[,non_ref_cols,drop = FALSE])
      }
    }
  }

  design
}


CalcOakesProfiledSurvivalH1 <- function(beta_vec,
                                        surv_coef,
                                        survival_context,
                                        fit_mix_num = length(beta_vec),
                                        surv_covar_risk_vec = NULL) {
  re_prob <- survival_context$re_prob
  if (is.null(re_prob)) {
    stop("CalcOakesProfiledSurvivalH1 requires survival_context$re_prob")
  }
  if (is.null(dim(re_prob))) {
    re_prob <- matrix(re_prob, ncol = fit_mix_num)
    survival_context$re_prob <- re_prob
  }
  if (ncol(re_prob) != fit_mix_num) {
    stop("survival_context$re_prob must have fit_mix_num columns")
  }

  beta_surv_coef <- IntoBetaSurvCoef(beta_vec = beta_vec,
                                    surv_coef = surv_coef,
                                    fit_mix_num = fit_mix_num)

  surv_coef_len <- unlist(lapply(surv_coef, length))

  out <- CalcSurvivalScoreInfo(beta_surv_coef = beta_surv_coef,
                              survival_context = survival_context,
                              surv_coef_len = surv_coef_len,
                              fit_mix_num = fit_mix_num,
                              surv_covar_risk_vec = surv_covar_risk_vec)

  list(score = out$score,
      info = out$info,
      hessian = -out$info,
      beta_vec = out$beta_vec,
      surv_coef = out$surv_coef,
      surv_covar_risk_vec = out$surv_covar_risk_vec,
      parameter_order = "age, class 2:G survival effects, non-reference survival covariate effects")
}




CalcOakesFixedBaselineSurvivalScoreInfo <- function(beta_vec,surv_coef,
                                                    survival_context,
                                                    cbline_vec,
                                                    fit_mix_num =
                                                      length(beta_vec)){
  re_prob <- survival_context$re_prob
  surv_event <- survival_context$surv_event
  surv_covar <- survival_context$surv_covar

  if (is.null(dim(re_prob))){
    re_prob <- matrix(re_prob,ncol = fit_mix_num)
  }

  if (ncol(re_prob) != fit_mix_num){
    stop("survival_context$re_prob must have fit_mix_num columns")
  }

  if (length(cbline_vec) != length(surv_event)){
    stop("cbline_vec must have one entry per survival observation")
  }

  surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
  design <- BuildOakesSurvivalDesign(surv_covar,surv_coef)

  common_param_num <- ncol(design)
  param_num <- fit_mix_num + length(unlist(surv_coef)) - length(surv_coef)
  score <- numeric(param_num)
  hessian <- matrix(0,nrow = param_num,ncol = param_num)

  common_index <- c(1)
  if (common_param_num > 1){
    common_index <- c(common_index,
                      fit_mix_num + seq_len(common_param_num - 1))
  }

  exp_beta <- exp(beta_vec)
  exp_surv_covar <- exp(surv_covar_risk_vec)
  weighted_class_risk <- sweep(re_prob,2,exp_beta,"*")
  class_risk <- rowSums(weighted_class_risk)
  baseline_risk <- cbline_vec * exp_surv_covar
  event_minus_risk <- surv_event - baseline_risk * class_risk

  score[common_index] <- colSums(design * event_minus_risk)
  hessian[common_index,common_index] <-
    -crossprod(design,design * (baseline_risk * class_risk))

  if (fit_mix_num > 1){
    for (re_ind in 2:fit_mix_num){
      class_risk_ind <- baseline_risk * weighted_class_risk[,re_ind]

      score[re_ind] <- sum(surv_event * re_prob[,re_ind] -
                             class_risk_ind)
      hessian[re_ind,re_ind] <- -sum(class_risk_ind)

      cross_val <- -colSums(design * class_risk_ind)
      hessian[common_index,re_ind] <- cross_val
      hessian[re_ind,common_index] <- cross_val
    }
  }

  hessian <- 0.5 * (hessian + t(hessian))

  list(score = score,
       info = -hessian,
       hessian = hessian,
       surv_covar_risk_vec = surv_covar_risk_vec,
       parameter_order =
         "age, class 2:G survival effects, non-reference survival covariate effects")
}

CalcOakesSurvivalH1 <- function(beta_vec,surv_coef,survival_context,
                                fit_mix_num = length(beta_vec),
                                cbline_vec = survival_context$cbline_vec){
  if (is.null(cbline_vec)){
    stop("CalcOakesSurvivalH1 requires fixed cbline_vec")
  }

  CalcOakesFixedBaselineSurvivalScoreInfo(beta_vec = beta_vec,
                                          surv_coef = surv_coef,
                                          survival_context =
                                            survival_context,
                                          cbline_vec = cbline_vec,
                                          fit_mix_num = fit_mix_num)
}

ActiveEmissionDayTypes <- function(vcovar_mat,vcovar_num){
  emission_day_types <- as.vector(vcovar_mat) + 1
  emission_day_types <- emission_day_types[!is.na(emission_day_types)]
  active_day_types <- sort(unique(emission_day_types))
  active_day_types <- active_day_types[active_day_types >= 1 &
                                         active_day_types <= vcovar_num]

  if (length(active_day_types) == 0){
    stop("No emission day types were observed in vcovar_mat")
  }

  active_day_types
}

EmissionPsi <- function(emit_act,emit_light,corr_mat,state_ind,
                        re_ind,vcovar_ind){
  c(emit_act[state_ind,1,re_ind,vcovar_ind],
    emit_act[state_ind,2,re_ind,vcovar_ind],
    emit_light[state_ind,1,re_ind,vcovar_ind],
    emit_light[state_ind,2,re_ind,vcovar_ind],
    corr_mat[re_ind,state_ind,vcovar_ind])
}

CleanEmissionInputs <- function(emit_inputs){
  finite_weight <- is.finite(emit_inputs$weights_vec)
  dropped_weight_count <- sum(!finite_weight)

  emit_inputs$weights_vec <- emit_inputs$weights_vec[finite_weight]
  if (!is.null(emit_inputs$emit_subset)){
    emit_inputs$emit_subset$index <- emit_inputs$emit_subset$index[finite_weight]
    emit_inputs$emit_subset$act <- emit_inputs$emit_subset$act[finite_weight]
    emit_inputs$emit_subset$light <- emit_inputs$emit_subset$light[finite_weight]
  }

  emit_inputs$dropped_weight_count <- dropped_weight_count
  emit_inputs
}

EmissionBlockObjective <- function(psi,act,light,lod_act,lod_light,
                                   vcovar_mat_emit,vcovar_ind,
                                   emit_inputs){
  if (any(!is.finite(psi)) || psi[2] <= 0 || psi[4] <= 0 ||
      abs(psi[5]) >= 0.999){
    return(1e100)
  }

  val <- tryCatch(
    withCallingHandlers(
      EmitLogLike(act = act,
                  light = light,
                  mu_act = psi[1],
                  sig_act = psi[2],
                  mu_light = psi[3],
                  sig_light = psi[4],
                  bivar_corr = psi[5],
                  lod_act = lod_act,
                  lod_light = lod_light,
                  vcovar_mat = vcovar_mat_emit,
                  vcovar_ind = vcovar_ind,
                  weights_mat = emit_inputs$weights_mat,
                  emit_subset = emit_inputs$emit_subset,
                  weights_vec = emit_inputs$weights_vec),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) 1e100
  )

  if (!is.finite(val)){
    return(1e100)
  }

  val
}

CalcOakesEmissionScore <- function(alpha,beta,pi_l,act,light,vcovar_mat,
                                   emit_act,emit_light,corr_mat,
                                   lod_act,lod_light,
                                   vcovar_num = dim(emit_act)[4]){
  if (!requireNamespace("numDeriv",quietly = TRUE)){
    stop("CalcOakesEmissionScore requires the numDeriv package")
  }

  state_num <- dim(emit_act)[1]
  mix_num <- dim(emit_act)[3]
  block_size <- 5

  active_vcovar_inds <- ActiveEmissionDayTypes(vcovar_mat,vcovar_num)
  vcovar_mat_emit <- vcovar_mat + 1
  emit_data <- PrepareEmitLogLikeData(act,light,vcovar_mat_emit)

  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array <- list(exp(weights_array_list[[1]]),
                        exp(weights_array_list[[2]]))

  block_count <- length(active_vcovar_inds) * mix_num * state_num
  score <- numeric(block_count * block_size)
  block_info <- data.frame(block = seq_len(block_count),
                           vcovar_ind = integer(block_count),
                           re_ind = integer(block_count),
                           state_ind = integer(block_count),
                           weight_sum = numeric(block_count),
                           dropped_weight_count = integer(block_count),
                           objective_at_value = numeric(block_count))

  block_ind <- 1
  pos <- 1
  for (vcovar_ind in active_vcovar_inds){
    for (re_ind in seq_len(mix_num)){
      for (state_ind in seq_len(state_num)){
        psi <- EmissionPsi(emit_act,emit_light,corr_mat,state_ind,
                           re_ind,vcovar_ind)
        emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                                     weights_array[[state_ind]],
                                                     re_ind)
        emit_inputs <- CleanEmissionInputs(emit_inputs)

        objective <- function(psi_work){
          EmissionBlockObjective(psi_work,act,light,lod_act,lod_light,
                                 vcovar_mat_emit,vcovar_ind,emit_inputs)
        }

        idx <- pos:(pos + block_size - 1)
        score[idx] <- -numDeriv::grad(objective,psi)
        block_info[block_ind,] <- list(block_ind,vcovar_ind,
                                       re_ind,state_ind,
                                       sum(emit_inputs$weights_vec),
                                       emit_inputs$dropped_weight_count,
                                       objective(psi))

        block_ind <- block_ind + 1
        pos <- pos + block_size
      }
    }
  }

  list(score = score,
       block_info = block_info,
       active_vcovar_inds = active_vcovar_inds,
       parameter_order = "mu_act, sig_act, mu_light, sig_light, corr")
}

CalcOakesEmissionH1 <- function(alpha,beta,pi_l,act,light,vcovar_mat,
                                emit_act,emit_light,corr_mat,
                                lod_act,lod_light,
                                vcovar_num = dim(emit_act)[4]){
  if (!requireNamespace("numDeriv",quietly = TRUE)){
    stop("CalcOakesEmissionH1 requires the numDeriv package")
  }

  state_num <- dim(emit_act)[1]
  mix_num <- dim(emit_act)[3]
  block_size <- 5

  active_vcovar_inds <- ActiveEmissionDayTypes(vcovar_mat,vcovar_num)
  vcovar_mat_emit <- vcovar_mat + 1
  emit_data <- PrepareEmitLogLikeData(act,light,vcovar_mat_emit)

  weights_array_list <- CondMarginalize(alpha,beta,pi_l)
  weights_array <- list(exp(weights_array_list[[1]]),
                        exp(weights_array_list[[2]]))

  block_count <- length(active_vcovar_inds) * mix_num * state_num
  hessian <- matrix(0,nrow = block_count * block_size,
                    ncol = block_count * block_size)
  score <- numeric(block_count * block_size)
  hessian_blocks <- vector(mode = "list",length = block_count)
  block_info <- data.frame(block = seq_len(block_count),
                           vcovar_ind = integer(block_count),
                           re_ind = integer(block_count),
                           state_ind = integer(block_count),
                           weight_sum = numeric(block_count),
                           dropped_weight_count = integer(block_count),
                           objective_at_hat = numeric(block_count),
                           hessian_abs_sum = numeric(block_count))

  block_ind <- 1
  pos <- 1
  for (vcovar_ind in active_vcovar_inds){
    for (re_ind in seq_len(mix_num)){
      for (state_ind in seq_len(state_num)){
        psi_hat <- EmissionPsi(emit_act,emit_light,corr_mat,state_ind,
                               re_ind,vcovar_ind)
        emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                                     weights_array[[state_ind]],
                                                     re_ind)
        emit_inputs <- CleanEmissionInputs(emit_inputs)

        objective <- function(psi){
          EmissionBlockObjective(psi,act,light,lod_act,lod_light,
                                 vcovar_mat_emit,vcovar_ind,emit_inputs)
        }

        objective_at_hat <- objective(psi_hat)
        block_score <- -numDeriv::grad(objective,psi_hat)
        block_hessian <- -numDeriv::hessian(objective,psi_hat)
        block_hessian <- 0.5 * (block_hessian + t(block_hessian))

        idx <- pos:(pos + block_size - 1)
        score[idx] <- block_score
        hessian[idx,idx] <- block_hessian
        hessian_blocks[[block_ind]] <- block_hessian
        block_info[block_ind,] <- list(block_ind,vcovar_ind,
                                       re_ind,state_ind,
                                       sum(emit_inputs$weights_vec),
                                       emit_inputs$dropped_weight_count,
                                       objective_at_hat,
                                       sum(abs(block_hessian)))

        block_ind <- block_ind + 1
        pos <- pos + block_size
      }
    }
  }

  list(score = score,
       hessian = hessian,
       hessian_blocks = hessian_blocks,
       block_info = block_info,
       active_vcovar_inds = active_vcovar_inds,
       parameter_order = "mu_act, sig_act, mu_light, sig_light, corr")
}

ActiveTransitionDayTypes <- function(vcovar_mat,vcovar_num){
  transition_day_types <- as.vector(vcovar_mat[-1,,drop = FALSE]) + 1
  transition_day_types <- transition_day_types[!is.na(transition_day_types)]
  active_day_types <- sort(unique(transition_day_types))
  active_day_types <- active_day_types[active_day_types >= 1 &
                                         active_day_types <= vcovar_num]

  if (length(active_day_types) == 0){
    stop("No transition day types were observed in vcovar_mat[-1, ]")
  }

  active_day_types
}

TransitionScoreArrayToVector <- function(score_array,active_vcovar_inds){
  block_size <- dim(score_array)[1]
  mix_num <- dim(score_array)[2]
  score <- numeric(block_size * mix_num * length(active_vcovar_inds))

  pos <- 1
  for (vcovar_ind in active_vcovar_inds){
    for (g in seq_len(mix_num)){
      idx <- pos:(pos + block_size - 1)
      score[idx] <- score_array[,g,vcovar_ind]
      pos <- pos + block_size
    }
  }

  score
}

CalcOakesTransitionScoreFromPosterior <- function(params_tran_array,
                                                  transition_weights,
                                                  vcovar_mat,
                                                  active_vcovar_inds,
                                                  period_len,
                                                  drop_nonfinite_weights =
                                                    TRUE){
  len <- dim(vcovar_mat)[1]
  mix_num <- dim(params_tran_array)[1]
  vcovar_num <- dim(params_tran_array)[3]

  transition_score <- array(0,dim = c(6,mix_num,vcovar_num))
  tran_list_mat <- GenTranColVecList(params_tran_array,len,mix_num,
                                     vcovar_num,period_len = period_len)
  cos_vec <- cos(2 * pi * seq(2,len) / period_len)
  sin_vec <- sin(2 * pi * seq(2,len) / period_len)
  dropped_weight_count <- 0

  for (init_state in 1:2){
    for (new_state in 1:2){
      tran_vals <- transition_weights[init_state,new_state,,,]
      if (mix_num == 1){
        dim(tran_vals) <- c(len - 1,dim(vcovar_mat)[2],1)
      }

      for (ind in seq_len(dim(vcovar_mat)[2])){
        vcovar_vec <- vcovar_mat[-1,ind]
        vcovar_vecR <- vcovar_vec + 1

        for (re_ind in seq_len(mix_num)){
          tran_mat <- matrix(NA,nrow = len - 1,ncol = 4)
          for (vcovar_ind in seq_len(vcovar_num)){
            vcovar_rows <- !is.na(vcovar_vecR) & vcovar_vecR == vcovar_ind
            if (any(vcovar_rows)){
              tran_mat[vcovar_rows,] <-
                tran_list_mat[[re_ind]][[vcovar_ind]][vcovar_rows,,
                                                       drop = FALSE]
            }
          }

          if(init_state == 1 & new_state == 1){
            tran_prime <- -tran_mat[,3]
          } else if(init_state == 1 & new_state == 2){
            tran_prime <- tran_mat[,1]
          } else if(init_state == 2 & new_state == 2){
            tran_prime <- -tran_mat[,2]
          } else if(init_state == 2 & new_state == 1){
            tran_prime <- tran_mat[,4]
          }

          transition_weight <- tran_vals[,ind,re_ind]
          finite_weight <- is.finite(transition_weight)
          dropped_weight_count <- dropped_weight_count + sum(!finite_weight)
          if (!drop_nonfinite_weights && any(!finite_weight)){
            stop("transition_weights contains non-finite values")
          }

          for (vcovar_ind in active_vcovar_inds){
            vcovar_rows <- !is.na(vcovar_vecR) & vcovar_vecR == vcovar_ind
            keep <- vcovar_rows & finite_weight
            if (!any(keep)){
              next
            }

            param_offset <- if (init_state == 1) 0 else 3
            transition_score[param_offset + 1,re_ind,vcovar_ind] <-
              transition_score[param_offset + 1,re_ind,vcovar_ind] +
              sum(transition_weight[keep] * tran_prime[keep])
            transition_score[param_offset + 2,re_ind,vcovar_ind] <-
              transition_score[param_offset + 2,re_ind,vcovar_ind] +
              sum(transition_weight[keep] * tran_prime[keep] *
                    cos_vec[keep])
            transition_score[param_offset + 3,re_ind,vcovar_ind] <-
              transition_score[param_offset + 3,re_ind,vcovar_ind] +
              sum(transition_weight[keep] * tran_prime[keep] *
                    sin_vec[keep])
          }
        }
      }
    }
  }

  list(score_array = transition_score,
       score = TransitionScoreArrayToVector(transition_score,
                                            active_vcovar_inds),
       dropped_weight_count = dropped_weight_count)
}

TransitionHessianArrayToMatrix <- function(hessian_array,active_vcovar_inds){
  block_size <- dim(hessian_array)[1]
  mix_num <- dim(hessian_array)[3]
  total_size <- block_size * mix_num * length(active_vcovar_inds)

  hessian <- matrix(0,nrow = total_size,ncol = total_size)

  # Order: active vcovar 1, class 1 params 1:6; active vcovar 1,
  # class 2 params 1:6; etc.
  pos <- 1
  for (vcovar_ind in active_vcovar_inds){
    for (g in seq_len(mix_num)){
      idx <- pos:(pos + block_size - 1)
      hessian[idx,idx] <- hessian_array[,,g,vcovar_ind]
      pos <- pos + block_size
    }
  }

  hessian
}

CalcOakesTransitionH1 <- function(alpha,beta,act,light,params_tran_array,
                                  emit_act,emit_light,corr_mat,pi_l,
                                  lod_act,lod_light,lintegral_mat,vcovar_mat,
                                  lambda_act_mat,lambda_light_mat,tobit,
                                  period_len,
                                  vcovar_num = dim(params_tran_array)[3]){
  active_vcovar_inds <- ActiveTransitionDayTypes(vcovar_mat,vcovar_num)

  tran_gradhess <- CalcTranCHelper(alpha = alpha,
                                   beta = beta,
                                   act = act,
                                   light = light,
                                   params_tran_array = params_tran_array,
                                   emit_act = emit_act,
                                   emit_light = emit_light,
                                   corr_mat = corr_mat,
                                   pi_l = pi_l,
                                   lod_act = lod_act,
                                   lod_light = lod_light,
                                   lintegral_mat = lintegral_mat,
                                   vcovar_mat = vcovar_mat,
                                   lambda_act_mat = lambda_act_mat,
                                   lambda_light_mat = lambda_light_mat,
                                   tobit = tobit,
                                   check_tran = FALSE,
                                   likelihood = NA,
                                   period_len = period_len,
                                   vcovar_num = vcovar_num)

  score_array <- tran_gradhess[[1]]
  hessian_array <- tran_gradhess[[2]]

  list(score_array = score_array,
       score = TransitionScoreArrayToVector(score_array,
                                            active_vcovar_inds),
       hessian_array = hessian_array,
       hessian = TransitionHessianArrayToMatrix(hessian_array,
                                                active_vcovar_inds),
       active_vcovar_inds = active_vcovar_inds)
}

CalcOakesScore <- function(theta,theta_pack,posterior_context,
                           return_components = FALSE,
                           drop_nonfinite_transition_weights = TRUE){
  params <- UnpackOakesTheta(theta,theta_pack)
  parameter_map <- theta_pack$parameter_map
  score <- numeric(length(theta))
  names(score) <- names(theta)
  components <- list()

  fill_score <- function(block,block_score){
    idx <- which(parameter_map$block == block)
    if (length(idx) != length(block_score)){
      stop(paste("Score length mismatch for block",block,
                 ": expected",length(idx),"got",length(block_score)))
    }
    score[idx] <<- block_score
  }

  init_idx <- which(parameter_map$block == "initial")
  if (length(init_idx) > 0){
    init_counts <- GetOakesContextValue(posterior_context,"init_counts",
                                        required = FALSE)
    if (is.null(init_counts)){
      init_counts <- CalcOakesInitialCounts(
        alpha = GetOakesContextValue(posterior_context,"alpha"),
        beta = GetOakesContextValue(posterior_context,"beta"),
        pi_l = GetOakesContextValue(posterior_context,"pi_l"),
        log_sweights_vec =
          GetOakesContextValue(posterior_context,"log_sweights_vec"))
    }

    init_score <- CalcOakesInitialH1(params$init,init_counts,
                                    theta_pack$sleep_state)
    fill_score("initial",init_score$score)
    components$initial <- init_score
  }

  tran_idx <- which(parameter_map$block == "transition")
  if (length(tran_idx) > 0){
    transition_weights <- GetOakesContextValue(posterior_context,
                                               "transition_weights",
                                               required = FALSE)
    if (is.null(transition_weights)){
      context_params_tran_array <- GetOakesContextValue(
        posterior_context,"params_tran_array",
        default = theta_pack$templates$params_tran_array,
        required = FALSE)
      context_emit_act <- GetOakesContextValue(
        posterior_context,"emit_act",
        default = theta_pack$templates$emit_act,
        required = FALSE)
      context_emit_light <- GetOakesContextValue(
        posterior_context,"emit_light",
        default = theta_pack$templates$emit_light,
        required = FALSE)
      context_corr_mat <- GetOakesContextValue(
        posterior_context,"corr_mat",
        default = theta_pack$templates$corr_mat,
        required = FALSE)
      context_lintegral_mat <- GetOakesContextValue(
        posterior_context,"lintegral_mat",required = FALSE)
      if (is.null(context_lintegral_mat)){
        context_lintegral_mat <- CalcLintegralMat(
          context_emit_act,context_emit_light,context_corr_mat,
          GetOakesContextValue(posterior_context,"lod_act"),
          GetOakesContextValue(posterior_context,"lod_light"))
      }

      transition_weights <- CalcOakesTransitionPosteriorWeights(
        alpha = GetOakesContextValue(posterior_context,"alpha"),
        beta = GetOakesContextValue(posterior_context,"beta"),
        pi_l = GetOakesContextValue(posterior_context,"pi_l"),
        act = GetOakesContextValue(posterior_context,"act"),
        light = GetOakesContextValue(posterior_context,"light"),
        params_tran_array = context_params_tran_array,
        emit_act = context_emit_act,
        emit_light = context_emit_light,
        corr_mat = context_corr_mat,
        lod_act = GetOakesContextValue(posterior_context,"lod_act"),
        lod_light = GetOakesContextValue(posterior_context,"lod_light"),
        lintegral_mat = context_lintegral_mat,
        vcovar_mat = GetOakesContextValue(posterior_context,"vcovar_mat"),
        lambda_act_mat =
          GetOakesContextValue(posterior_context,"lambda_act_mat"),
        lambda_light_mat =
          GetOakesContextValue(posterior_context,"lambda_light_mat"),
        tobit = GetOakesContextValue(posterior_context,"tobit"),
        period_len = GetOakesContextValue(posterior_context,"period_len"))
    }

    tran_score <- CalcOakesTransitionScoreFromPosterior(
      params_tran_array = params$params_tran_array,
      transition_weights = transition_weights,
      vcovar_mat = GetOakesContextValue(posterior_context,"vcovar_mat"),
      active_vcovar_inds = theta_pack$active_tran_vcovar_inds,
      period_len = GetOakesContextValue(posterior_context,"period_len"),
      drop_nonfinite_weights = drop_nonfinite_transition_weights)
    fill_score("transition",tran_score$score)
    components$transition <- tran_score
  }

  emit_idx <- which(parameter_map$block == "emission")
  if (length(emit_idx) > 0){
    emit_score <- CalcOakesEmissionScore(
      alpha = GetOakesContextValue(posterior_context,"alpha"),
      beta = GetOakesContextValue(posterior_context,"beta"),
      pi_l = GetOakesContextValue(posterior_context,"pi_l"),
      act = GetOakesContextValue(posterior_context,"act"),
      light = GetOakesContextValue(posterior_context,"light"),
      vcovar_mat = GetOakesContextValue(posterior_context,"vcovar_mat"),
      emit_act = params$emit_act,
      emit_light = params$emit_light,
      corr_mat = params$corr_mat,
      lod_act = GetOakesContextValue(posterior_context,"lod_act"),
      lod_light = GetOakesContextValue(posterior_context,"lod_light"))
    fill_score("emission",emit_score$score)
    components$emission <- emit_score
  }

  mix_idx <- which(parameter_map$block == "mixing")
  if (length(mix_idx) > 0){
    re_prob <- GetOakesContextValue(posterior_context,"re_prob",
                                    required = FALSE)
    if (is.null(re_prob)){
      re_prob <- CalcProbRE(GetOakesContextValue(posterior_context,"alpha"),
                            GetOakesContextValue(posterior_context,"pi_l"))
    }

    mix_score <- CalcOakesMixingH1(
      nu_mat = params$nu_mat,
      re_prob = re_prob,
      nu_covar_mat = GetOakesContextValue(posterior_context,
                                          "nu_covar_mat"))
    fill_score("mixing",mix_score$score)
    components$mixing <- mix_score
  }

  surv_idx <- which(parameter_map$block == "survival")
  if (length(surv_idx) > 0){
    survival_baseline_mode <- GetOakesContextValue(
      posterior_context,"survival_baseline_mode",
      default = "fixed",
      required = FALSE)

    survival_context <- GetOakesContextValue(posterior_context,
                                             "survival_context")
    re_prob <- GetOakesContextValue(posterior_context,"re_prob",
                                    required = FALSE)
    if (!is.null(re_prob)){
      if (exists("update_survival_context_re_prob",mode = "function")){
        survival_context <- update_survival_context_re_prob(
          survival_context,re_prob,theta_pack$fit_mix_num)
      } else {
        survival_context$re_prob <- re_prob
      }
    }

    if (identical(survival_baseline_mode, "fixed")){
      surv_score <- CalcOakesSurvivalH1(
        beta_vec = params$beta_vec,
        surv_coef = params$surv_coef,
        survival_context = survival_context,
        fit_mix_num = theta_pack$fit_mix_num,
        cbline_vec = GetOakesContextValue(posterior_context,"cbline_vec",
                                          required = FALSE)
      )
    } else if (identical(survival_baseline_mode, "profiled")){
      surv_covar_risk_vec_for_score <- SurvCovarRiskVec(
        survival_context$surv_covar,
        params$surv_coef
      )

      surv_score <- CalcOakesProfiledSurvivalH1(
        beta_vec = params$beta_vec,
        surv_coef = params$surv_coef,
        survival_context = survival_context,
        fit_mix_num = theta_pack$fit_mix_num,
        surv_covar_risk_vec = surv_covar_risk_vec_for_score
      )

    } else {
      stop(paste("Unknown survival_baseline_mode:", survival_baseline_mode))
    }

    fill_score("survival",surv_score$score)
    components$survival <- surv_score
  }

  if (!all(is.finite(score))){
    warning("CalcOakesScore produced non-finite score entries")
  }

  if (return_components){
    return(list(score = score,
                components = components,
                params = params,
                parameter_map = parameter_map))
  }

  score
}

DiagnoseOakesScore <- function(oakes_score){
  score <- if (is.list(oakes_score)) oakes_score$score else oakes_score
  parameter_map <- if (is.list(oakes_score)) oakes_score$parameter_map else NULL

  bad <- which(!is.finite(score))
  if (length(bad) == 0){
    return(data.frame())
  }

  if (!is.null(parameter_map)){
    out <- parameter_map[bad,,drop = FALSE]
  } else {
    out <- data.frame(index = bad)
  }
  out$score <- score[bad]
  rownames(out) <- NULL
  out
}
