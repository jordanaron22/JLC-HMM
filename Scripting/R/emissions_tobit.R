#Current implementation is a case-by-case basis that is pre-sorted
#Not nearly as easy to understand but much faster (one optimization call goes from 22 to 4, we do 10 calls)

#Calculates mu1|2
CalcCondMean <- function(mu1,sig1,mu2,sig2,bivar_corr,obs2){
  return(mu1 + bivar_corr*(sig1/sig2)*(obs2-mu2))
}

#Calculates sig1|2
CalcCondSig <- function(sig1,bivar_corr){
  return(sig1*sqrt(1-bivar_corr^2))
}

#Calcuates case where both activity and light are below LoD
Case4 <- function(act_obs,mu_act,sig_act,mu_light,sig_light,bivar_corr,light_LOD){

  mu_light_cond <- CalcCondMean(mu_light,sig_light,mu_act,sig_act,bivar_corr,act_obs)
  sig_light_cond <- CalcCondSig(sig_light,bivar_corr)

  lognorm_dens <- dnorm(act_obs,mu_act,sig_act) *
    pnorm(light_LOD,mu_light_cond,sig_light_cond)
  return(lognorm_dens)
}

#AI optimized Case4
#checked that is is equivalent to Case4 with integration
CalcCase4LogProbFast <- function(
    lower_act,
    upper_act,
    mu_act,
    sig_act,
    mu_light,
    sig_light,
    bivar_corr,
    upper_light
) {
  # Standardized upper limits.
  z_act_upper <- (upper_act - mu_act) / sig_act
  z_light_upper <- (upper_light - mu_light) / sig_light

  # P(A <= upper_act, L <= upper_light)
  upper_probability <- pbivnorm::pbivnorm(
    x = z_act_upper,
    y = z_light_upper,
    rho = bivar_corr
  )

  # P(A <= lower_act, L <= upper_light)
  lower_probability <- if (
    is.infinite(lower_act) && lower_act < 0
  ) {
    0
  } else {
    z_act_lower <- (lower_act - mu_act) / sig_act

    pbivnorm::pbivnorm(
      x = z_act_lower,
      y = z_light_upper,
      rho = bivar_corr
    )
  }

  rectangle_probability <-
    as.numeric(upper_probability - lower_probability)

  # Preserve the behavior of the current code:
  #
  # log(0) becomes -Inf and is then replaced with -9999.
  #
  # A very small negative value can also arise from floating-point
  # subtraction, so handle that the same way.
  if (
    !is.finite(rectangle_probability) ||
      rectangle_probability <= 0
  ) {
    return(-9999)
  }

  log(rectangle_probability)
}


##Calculates case 4 for given normal parameters and store in matrix for access later
# CalcLintegralMatOLD <- function(emit_act,emit_light,corr_mat,lod_act,lod_light){
#   mix_num <- dim(emit_act)[3]
#   if (is.na(mix_num)){mix_num <- 1}

#   lintegral_mat <- array(NA,dim = c(mix_num,2,2))
#   #j is week/weekend
#   for (j in 1:2){
#     for (i in 1:mix_num){

#       lintegral_mat[i,1,j] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
#                                           emit_act[1,1,i,j],emit_act[1,2,i,j],
#                                           emit_light[1,1,i,j],emit_light[1,2,i,j],
#                                           corr_mat[i,1,j],lod_light)[[1]])

#       lintegral_mat[i,2,j] <- log(integrate(Case4,lower = -Inf,upper = lod_act,
#                                           emit_act[2,1,i,j],emit_act[2,2,i,j],
#                                           emit_light[2,1,i,j],emit_light[2,2,i,j],
#                                           corr_mat[i,2,j],lod_light)[[1]])
#     }
#   }

#   #work on log scale so -9999 is effectively -Inf
#   lintegral_mat[lintegral_mat == -Inf] <- -9999

#   return(lintegral_mat)
# }

#AI optimized
CalcLintegralMat <- function(
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light
) {
  mix_num <- dim(emit_act)[3L]

  if (is.na(mix_num)) {
    mix_num <- 1L
  }

  state_num <- dim(emit_act)[1L]
  vcovar_num <- dim(emit_act)[4L]

  lintegral_mat <- array(
    NA_real_,
    dim = c(mix_num, state_num, vcovar_num)
  )

  for (vcovar_ind in seq_len(vcovar_num)) {
    for (re_ind in seq_len(mix_num)) {
      for (state_ind in seq_len(state_num)) {
        lintegral_mat[
          re_ind,
          state_ind,
          vcovar_ind
        ] <- CalcCase4LogProbFast(
          lower_act = -Inf,
          upper_act = lod_act,
          mu_act =
            emit_act[state_ind, 1L, re_ind, vcovar_ind],
          sig_act =
            emit_act[state_ind, 2L, re_ind, vcovar_ind],
          mu_light =
            emit_light[state_ind, 1L, re_ind, vcovar_ind],
          sig_light =
            emit_light[state_ind, 2L, re_ind, vcovar_ind],
          bivar_corr =
            corr_mat[re_ind, state_ind, vcovar_ind],
          upper_light = lod_light
        )
      }
    }
  }

  lintegral_mat
}
PrepareCase2Grouping <- function(
    act_subset,
    light_subset,
    case2_index,
    lod_light
) {
  if (length(case2_index) == 0L) {
    return(list(
      group = integer(0),
      unique_act = numeric(0),
      light_value = lod_light
    ))
  }

  case2_act <-
    act_subset[case2_index]

  case2_light <-
    light_subset[case2_index]

  # In the current NHANES data all censored-light values
  # are stored at one common value.
  case2_light_value <-
    case2_light[[1L]]

  if (any(case2_light != case2_light_value)) {
    stop(
      paste(
        "Case 2 contains more than one censored-light value.",
        "Grouping only by activity would not be exact."
      )
    )
  }

  case2_unique_act <-
    unique(case2_act)

  # One-based group number for each Case 2 observation.
  # This is created once before the EM loop.
  case2_group <-
    match(
      case2_act,
      case2_unique_act
    )

  if (anyNA(case2_group)) {
    stop(
      "Failed to match Case 2 observations to activity groups"
    )
  }

  list(
    group = as.integer(case2_group),
    unique_act = as.numeric(case2_unique_act),
    light_value = case2_light_value
  )
}

#above parameters are for debugging
#calculates likelihood of emission dist
#used in direct optimization
##### 9 cases
#1: Activity (A) above LoD (>) & Light (L) >
#2 A> & L below LoD (<)
#3 A< & L>
#4 A< & L<
#5 A> & L missing (m)
#6 A< & Lm
#7 Am & L>
#8 Am & L<
#9 Am & Lm
#case 9 is dropped from emission optimization as it does not give any information
PrepareEmitLogLikeData <- function(
    act,
    light,
    vcovar_mat,
    lod_act = NULL,
    lod_light = NULL
) {
  if (!identical(dim(act), dim(light))) {
    stop("act and light must have the same dimensions")
  }

  if (!identical(dim(act), dim(vcovar_mat))) {
    stop("act and vcovar_mat must have the same dimensions")
  }

  n_time <- nrow(act)
  n_subject <- ncol(act)

  act_vec <- as.vector(act)
  light_vec <- as.vector(light)
  vcovar_vec <- as.vector(vcovar_mat)

  vcovar_levels <- sort(
    unique(vcovar_vec[!is.na(vcovar_vec)])
  )

  vcovar_indices <- vector("list", length = 0L)
  vcovar_subject_indices <- vector("list", length = 0L)
  emit_subsets <- vector("list", length = 0L)

  for (vcovar_level in vcovar_levels) {
    vcovar_name <- as.character(vcovar_level)

    obs_index <- which(vcovar_vec == vcovar_level)

    subject_index <-
      ((obs_index - 1L) %/% n_time) + 1L

    act_subset <- act_vec[obs_index]
    light_subset <- light_vec[obs_index]

    case_indices <- NULL

    if (!is.null(lod_act) && !is.null(lod_light)) {
      act_missing <- is.na(act_subset)
      light_missing <- is.na(light_subset)

      both_present <- !act_missing & !light_missing

      
      case_indices <- list(
        case1 = which(
          both_present &
            act_subset > lod_act &
            light_subset > lod_light
        ),

        case2 = which(
          both_present &
            act_subset > lod_act &
            light_subset <= lod_light
        ),

        case3 = which(
          both_present &
            act_subset <= lod_act &
            light_subset > lod_light
        ),

        case4 = which(
          both_present &
            act_subset <= lod_act &
            light_subset <= lod_light
        ),

        case5 = which(
          !act_missing &
            light_missing &
            act_subset > lod_act
        ),

        case6 = which(
          !act_missing &
            light_missing &
            act_subset <= lod_act
        ),

        case7 = which(
          act_missing &
            !light_missing &
            light_subset > lod_light
        ),

        case8 = which(
          act_missing &
            !light_missing &
            light_subset <= lod_light
        )
      )

      case2_grouping <- PrepareCase2Grouping(
        act_subset = act_subset,
        light_subset = light_subset,
        case2_index = case_indices$case2,
        lod_light = lod_light
      )

      case_indices$case2_group <-
        case2_grouping$group

      case_indices$case2_unique_act <-
        case2_grouping$unique_act

      case_indices$case2_light_value <-
        case2_grouping$light_value


    }

    vcovar_indices[[vcovar_name]] <- obs_index
    vcovar_subject_indices[[vcovar_name]] <- subject_index

    emit_subsets[[vcovar_name]] <- list(
      index = obs_index,
      subject_index = subject_index,
      act = act_subset,
      light = light_subset,
      case_indices = case_indices
    )
  }

  list(
    act_vec = act_vec,
    light_vec = light_vec,
    vcovar_indices = vcovar_indices,
    vcovar_subject_indices = vcovar_subject_indices,
    emit_subsets = emit_subsets,
    n_time = n_time,
    n_subject = n_subject
  )
}

GetEmitLogLikeSubset <- function(emit_data, vcovar_ind) {
  if (is.null(emit_data)) {
    return(NULL)
  }

  vcovar_name <- as.character(vcovar_ind)

  emit_subset <- emit_data$emit_subsets[[vcovar_name]]

  # Preserve the current behavior if a day type has no observations.
  if (is.null(emit_subset)) {
    return(list(
      index = integer(0),
      subject_index = integer(0),
      act = numeric(0),
      light = numeric(0)
    ))
  }

  emit_subset
}



PrepareEmitOptimizationInputs <- function(
    emit_data,
    vcovar_ind,
    weights_array,
    re_ind
) {
  emit_subset <- GetEmitLogLikeSubset(
    emit_data = emit_data,
    vcovar_ind = vcovar_ind
  )

  # New cached path used by the main EM loop.
  #
  # UpdateNorm passes the cache for one Markov state, so:
  #
  # weights_array[[re_ind]][[vcovar_ind]]
  #
  # is already the exact vector needed by the optimizer.
  if (is.list(weights_array) && is.null(dim(weights_array))) {
    if (re_ind < 1L || re_ind > length(weights_array)) {
      stop("re_ind is outside the emission-weight cache")
    }

    class_cache <- weights_array[[re_ind]]

    if (
      vcovar_ind < 1L ||
      vcovar_ind > length(class_cache)
    ) {
      stop("vcovar_ind is outside the emission-weight cache")
    }

    weights_vec <- class_cache[[vcovar_ind]]

    if (length(weights_vec) != length(emit_subset$index)) {
      stop(
        paste(
          "Cached emission weights do not match",
          "the emission-data subset"
        )
      )
    }

    return(list(
      emit_subset = emit_subset,
      weights_mat = NULL,
      weights_vec = weights_vec
    ))
  }

  # Backward-compatible array path.
  #
  # This preserves uses elsewhere in the code, including the
  # Oakes-information calculations, which currently pass a full
  # time x subject x class array.
  if (is.null(emit_subset)) {
    return(list(
      emit_subset = NULL,
      weights_mat = as.vector(
        weights_array[, , re_ind]
      ),
      weights_vec = NULL
    ))
  }

  weights_vec <- ExtractEmitWeightSubset(
    weights_array = weights_array,
    re_ind = re_ind,
    obs_index = emit_subset$index
  )

  list(
    emit_subset = emit_subset,
    weights_mat = NULL,
    weights_vec = weights_vec
  )
}

#defunct now that we use case-by-case version
EmitLogLike <- function(act,light,mu_act,sig_act,mu_light,sig_light,bivar_corr,lod_act,lod_light,vcovar_mat,vcovar_ind,weights_mat,
                        emit_subset = NULL,weights_vec = NULL){

  #lower should theoretically be -Inf, but had some divergence issues
  lb <- min(mu_act - 5 * sig_act, -10)

  #old way of doing it, dont need to run a numerical integration
  # lintegral <- log(integrate(Case4,lower = lb,upper = lod_act,
  #                            mu_act,
  #                            sig_act,
  #                            mu_light,
  #                            sig_light,
  #                            bivar_corr,lod_light)[[1]])
  
  #AI improved way, validated against old way
  lintegral <- CalcCase4LogProbFast(
    lower_act = lb,
    upper_act = lod_act,
    mu_act = mu_act,
    sig_act = sig_act,
    mu_light = mu_light,
    sig_light = sig_light,
    bivar_corr = bivar_corr,
    upper_light = lod_light
  )

  if (lintegral == -Inf){
    lintegral <- -9999
  }

  if (is.null(emit_subset)){
    vcovar_vec <- as.vector(vcovar_mat)
    vcovar_vec_indicator <- vcovar_vec == vcovar_ind
    act_vec <- as.vector(act)[vcovar_vec_indicator]
    light_vec <- as.vector(light)[vcovar_vec_indicator]

    if (is.null(weights_vec)){
      weights_vec <- weights_mat[vcovar_vec_indicator]
    }
  } else {
    act_vec <- emit_subset$act
    light_vec <- emit_subset$light

    if (is.null(weights_vec)){
      weights_vec <- as.vector(weights_mat)[emit_subset$index]
    }
  }

  return(
    weightedLogClassificationCTobitC(
      act = act_vec,
      light = light_vec,
      weights = weights_vec,
      mu_act = mu_act,
      sig_act = sig_act,
      mu_light = mu_light,
      sig_light = sig_light,
      lod_act = lod_act,
      lod_light = lod_light,
      bivar_corr = bivar_corr,
      lintegral = lintegral
    )
  )
}

EmitLogLikeByCase <- function(
    mu_act,
    sig_act,
    mu_light,
    sig_light,
    bivar_corr,
    lod_act,
    lod_light,
    emit_subset,
    case_cache
) {
  lb <- min(
    mu_act - 5 * sig_act,
    -10
  )

  # Use the fast bivariate-normal probability helper
  lintegral <- CalcCase4LogProbFast(
    lower_act = lb,
    upper_act = lod_act,
    mu_act = mu_act,
    sig_act = sig_act,
    mu_light = mu_light,
    sig_light = sig_light,
    bivar_corr = bivar_corr,
    upper_light = lod_light
  )

  weightedTobitObjectiveByCaseC(
    act = emit_subset$act,
    light = emit_subset$light,
    case_indices = emit_subset$case_indices,
    case_cache = case_cache,
    mu_act = mu_act,
    sig_act = sig_act,
    mu_light = mu_light,
    sig_light = sig_light,
    lod_act = lod_act,
    lod_light = lod_light,
    bivar_corr = bivar_corr,
    lintegral = lintegral
  )
}

#optimizes activity mean
#all emission dist param calculated this way
#comment out which parameter currently being optimized
CalcActMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat,emit_data = NULL){
  emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                               weights_array,re_ind)

  mu_act <- optimize(EmitLogLike, c(-10,10), act = act, light = light,
                     # mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = emit_inputs$weights_mat,
                     emit_subset = emit_inputs$emit_subset,
                     weights_vec = emit_inputs$weights_vec,
                    tol = 1e-3)$minimum
  return(mu_act)
}

CalcActSig <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat,emit_data = NULL){
  emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                               weights_array,re_ind)

  mu_act <- optimize(EmitLogLike, c(0.1,10), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     # sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = emit_inputs$weights_mat,
                     emit_subset = emit_inputs$emit_subset,
                     weights_vec = emit_inputs$weights_vec,
                    tol = 1e-3)$minimum
  return(mu_act)
}

CalcLightMean <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat,emit_data = NULL){
  emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                               weights_array,re_ind)

  mu_act <- optimize(EmitLogLike, c(-30,10), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     # mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = emit_inputs$weights_mat,
                     emit_subset = emit_inputs$emit_subset,
                     weights_vec = emit_inputs$weights_vec,
                    tol = 1e-3)$minimum

  return(mu_act)
}

CalcLightSig <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat,emit_data = NULL){
  emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                               weights_array,re_ind)

  mu_act <- optimize(EmitLogLike, c(0.01,20), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     #sig_light = emit_light[mc_state,2,re_ind],
                     bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = emit_inputs$weights_mat,
                     emit_subset = emit_inputs$emit_subset,
                     weights_vec = emit_inputs$weights_vec,
                    tol = 1e-3)$minimum
  return(mu_act)
}

CalcBivarCorr <- function(mc_state,vcovar_ind,act,light,emit_act,emit_light,corr_mat,lod_act,lod_light,weights_array,re_ind,vcovar_mat,emit_data = NULL){
  emit_inputs <- PrepareEmitOptimizationInputs(emit_data,vcovar_ind,
                                               weights_array,re_ind)

  mu_act <- optimize(EmitLogLike, c(-.999,.999), act = act, light = light,
                     mu_act = emit_act[mc_state,1,re_ind,vcovar_ind],
                     sig_act = emit_act[mc_state,2,re_ind,vcovar_ind],
                     mu_light = emit_light[mc_state,1,re_ind,vcovar_ind],
                     sig_light = emit_light[mc_state,2,re_ind,vcovar_ind],
                     # bivar_corr = corr_mat[re_ind,mc_state,vcovar_ind],
                     lod_act = lod_act, lod_light = lod_light,
                     vcovar_mat = vcovar_mat, vcovar_ind = vcovar_ind,
                     weights_mat = emit_inputs$weights_mat,
                     emit_subset = emit_inputs$emit_subset,
                     weights_vec = emit_inputs$weights_vec,
                    tol = 1e-3)$minimum
  return(mu_act)
}

#Highest level for optimizing emission dist
#takes function as input, easier to process this way
UpdateNorm <- function(
    FUN,
    mc_state,
    act,
    light,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    emit_weight_cache,
    vcovar_mat,
    emit_data
) {
  mix_num <- dim(emit_act)[3L]
  vcovar_num <- dim(emit_act)[4L]

  opt_param_mat <- matrix(
    0,
    nrow = mix_num,
    ncol = vcovar_num
  )

  if (
    mc_state < 1L ||
    mc_state > length(emit_weight_cache)
  ) {
    stop("mc_state is outside emit_weight_cache")
  }

  # Select the state-specific cache once.
  #
  # Its structure is now:
  # weights_array[[class]][[day type]]
  weights_array <- emit_weight_cache[[mc_state]]

  for (re_ind in seq_len(mix_num)) {
    for (vcovar_ind in seq_len(vcovar_num)) {
      opt_param_mat[re_ind, vcovar_ind] <- FUN(
        mc_state,
        vcovar_ind,
        act,
        light,
        emit_act,
        emit_light,
        corr_mat,
        lod_act,
        lod_light,
        weights_array,
        re_ind,
        vcovar_mat,
        emit_data = emit_data
      )
    }
  }

  opt_param_mat
}


ExtractEmitWeightSubset <- function(
    weights_array,
    re_ind,
    obs_index
) {
  array_dim <- dim(weights_array)

  if (length(array_dim) != 3L) {
    stop(
      "weights_array must have dimensions time x subject x class"
    )
  }

  if (re_ind < 1L || re_ind > array_dim[3L]) {
    stop("re_ind is outside the class dimension")
  }

  observations_per_class <-
    array_dim[1L] * array_dim[2L]

  class_offset <-
    (re_ind - 1L) * observations_per_class

  weights_array[obs_index + class_offset]
}

BuildEmitWeightCache <- function(
    log_weights_by_state,
    emit_data,
    sweights_vec,
    mix_num,
    vcovar_num
) {
  state_num <- length(log_weights_by_state)

  if (state_num != NUM_MARKOV_STATES) {
    stop(
      "log_weights_by_state must contain one array per Markov state"
    )
  }

  reference_dim <- dim(log_weights_by_state[[1L]])

  if (length(reference_dim) != 3L) {
    stop(
      paste(
        "Each posterior-weight array must have dimensions",
        "time x subject x class"
      )
    )
  }

  if (reference_dim[2L] != length(sweights_vec)) {
    stop(
      "sweights_vec must contain one value per subject"
    )
  }

  if (reference_dim[3L] != mix_num) {
    stop(
      "The class dimension of the posterior weights does not match mix_num"
    )
  }

  for (state_ind in seq_len(state_num)) {
    if (!identical(
      dim(log_weights_by_state[[state_ind]]),
      reference_dim
    )) {
      stop(
        "All posterior-weight arrays must have the same dimensions"
      )
    }
  }

  # Structure:
  #
  # emit_weight_cache[[state_ind]][[re_ind]][[vcovar_ind]]
  #
  # For five classes and two day types, this creates:
  # 2 states x 5 classes x 2 day types = 20 vectors.
  emit_weight_cache <- vector(
    mode = "list",
    length = state_num
  )

  for (state_ind in seq_len(state_num)) {
    state_cache <- vector(
      mode = "list",
      length = mix_num
    )

    for (re_ind in seq_len(mix_num)) {
      class_cache <- vector(
        mode = "list",
        length = vcovar_num
      )

      for (vcovar_ind in seq_len(vcovar_num)) {
        emit_subset <- GetEmitLogLikeSubset(
          emit_data = emit_data,
          vcovar_ind = vcovar_ind
        )

        # Extract this state/class/day-type combination directly
        # from the log-posterior array.
        log_weights_vec <- ExtractEmitWeightSubset(
          weights_array =
            log_weights_by_state[[state_ind]],
          re_ind = re_ind,
          obs_index = emit_subset$index
        )

        # This is equivalent to:
        #
        # sweep(
        #   exp(weights_array),
        #   2,
        #   sweights_vec,
        #   "*"
        # )
        #
        # followed by subsetting, but it operates only on the
        # observations needed for this day type.
        weights_vec <-
          exp(log_weights_vec) *
          sweights_vec[emit_subset$subject_index]

        class_cache[[vcovar_ind]] <- weights_vec
      }

      state_cache[[re_ind]] <- class_cache
    }

    emit_weight_cache[[state_ind]] <- state_cache
  }

  names(emit_weight_cache) <- c("wake", "sleep")

  emit_weight_cache
}




BuildEmitCaseCache <- function(
    log_weights_by_state,
    emit_data,
    sweights_vec,
    mix_num,
    vcovar_num
) {
  state_num <- length(log_weights_by_state)

  reference_dim <- dim(
    log_weights_by_state[[1L]]
  )

  if (length(reference_dim) != 3L) {
    stop(
      paste(
        "Posterior weights must have dimensions",
        "time x subject x class"
      )
    )
  }

  if (reference_dim[2L] != length(sweights_vec)) {
    stop(
      "sweights_vec must contain one value per subject"
    )
  }

  emit_case_cache <- vector(
    "list",
    length = state_num
  )

  for (state_ind in seq_len(state_num)) {
    emit_case_cache[[state_ind]] <- vector(
      "list",
      length = mix_num
    )

    for (re_ind in seq_len(mix_num)) {
      emit_case_cache[[state_ind]][[re_ind]] <-
        vector(
          "list",
          length = vcovar_num
        )

      for (vcovar_ind in seq_len(vcovar_num)) {
        emit_subset <- GetEmitLogLikeSubset(
          emit_data = emit_data,
          vcovar_ind = vcovar_ind
        )

        if (is.null(emit_subset$case_indices)) {
          stop(
            paste(
              "Emission case indices are missing.",
              "Recreate emit_data with lod_act and lod_light."
            )
          )
        }

        log_weights_vec <- ExtractEmitWeightSubset(
          weights_array =
            log_weights_by_state[[state_ind]],
          re_ind = re_ind,
          obs_index = emit_subset$index
        )

        weights_vec <-
          exp(log_weights_vec) *
          sweights_vec[emit_subset$subject_index]

        emit_case_cache[[state_ind]][[re_ind]][[vcovar_ind]] <- buildTobitCaseCacheC(
          act = emit_subset$act,
          light = emit_subset$light,
          weights = weights_vec,
          case_indices =
            emit_subset$case_indices
        )

        rm(log_weights_vec, weights_vec)
      }
    }
  }

  names(emit_case_cache) <- c(
    "wake",
    "sleep"
  )

  emit_case_cache
}

CalcActMeanByCase <- function(
    mc_state,
    vcovar_ind,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    case_cache,
    re_ind,
    emit_subset
) {
  optimize(
    EmitLogLikeByCase,
    interval = c(-10, 10),

    sig_act =
      emit_act[mc_state, 2L, re_ind, vcovar_ind],

    mu_light =
      emit_light[mc_state, 1L, re_ind, vcovar_ind],

    sig_light =
      emit_light[mc_state, 2L, re_ind, vcovar_ind],

    bivar_corr =
      corr_mat[re_ind, mc_state, vcovar_ind],

    lod_act = lod_act,
    lod_light = lod_light,
    emit_subset = emit_subset,
    case_cache = case_cache,
    tol = 1e-3
  )$minimum
}


CalcActSigByCase <- function(
    mc_state,
    vcovar_ind,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    case_cache,
    re_ind,
    emit_subset
) {
  optimize(
    EmitLogLikeByCase,
    interval = c(0.1, 10),

    mu_act =
      emit_act[mc_state, 1L, re_ind, vcovar_ind],

    mu_light =
      emit_light[mc_state, 1L, re_ind, vcovar_ind],

    sig_light =
      emit_light[mc_state, 2L, re_ind, vcovar_ind],

    bivar_corr =
      corr_mat[re_ind, mc_state, vcovar_ind],

    lod_act = lod_act,
    lod_light = lod_light,
    emit_subset = emit_subset,
    case_cache = case_cache,
    tol = 1e-3
  )$minimum
}


CalcLightMeanByCase <- function(
    mc_state,
    vcovar_ind,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    case_cache,
    re_ind,
    emit_subset
) {
  optimize(
    EmitLogLikeByCase,
    interval = c(-30, 10),

    mu_act =
      emit_act[mc_state, 1L, re_ind, vcovar_ind],

    sig_act =
      emit_act[mc_state, 2L, re_ind, vcovar_ind],

    sig_light =
      emit_light[mc_state, 2L, re_ind, vcovar_ind],

    bivar_corr =
      corr_mat[re_ind, mc_state, vcovar_ind],

    lod_act = lod_act,
    lod_light = lod_light,
    emit_subset = emit_subset,
    case_cache = case_cache,
    tol = 1e-3
  )$minimum
}


CalcLightSigByCase <- function(
    mc_state,
    vcovar_ind,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    case_cache,
    re_ind,
    emit_subset
) {
  optimize(
    EmitLogLikeByCase,
    interval = c(0.01, 20),

    mu_act =
      emit_act[mc_state, 1L, re_ind, vcovar_ind],

    sig_act =
      emit_act[mc_state, 2L, re_ind, vcovar_ind],

    mu_light =
      emit_light[mc_state, 1L, re_ind, vcovar_ind],

    bivar_corr =
      corr_mat[re_ind, mc_state, vcovar_ind],

    lod_act = lod_act,
    lod_light = lod_light,
    emit_subset = emit_subset,
    case_cache = case_cache,
    tol = 1e-3
  )$minimum
}


CalcBivarCorrByCase <- function(
    mc_state,
    vcovar_ind,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    case_cache,
    re_ind,
    emit_subset
) {
  optimize(
    EmitLogLikeByCase,
    interval = c(-0.999, 0.999),

    mu_act =
      emit_act[mc_state, 1L, re_ind, vcovar_ind],

    sig_act =
      emit_act[mc_state, 2L, re_ind, vcovar_ind],

    mu_light =
      emit_light[mc_state, 1L, re_ind, vcovar_ind],

    sig_light =
      emit_light[mc_state, 2L, re_ind, vcovar_ind],

    lod_act = lod_act,
    lod_light = lod_light,
    emit_subset = emit_subset,
    case_cache = case_cache,
    tol = 1e-3
  )$minimum
}

UpdateNormByCase <- function(
    FUN,
    mc_state,
    emit_act,
    emit_light,
    corr_mat,
    lod_act,
    lod_light,
    emit_case_cache,
    emit_data
) {
  mix_num <- dim(emit_act)[3L]
  vcovar_num <- dim(emit_act)[4L]

  opt_param_mat <- matrix(
    NA_real_,
    nrow = mix_num,
    ncol = vcovar_num
  )

  for (re_ind in seq_len(mix_num)) {
    for (vcovar_ind in seq_len(vcovar_num)) {
      emit_subset <- GetEmitLogLikeSubset(
        emit_data,
        vcovar_ind
      )

      case_cache <- emit_case_cache[[mc_state]][[re_ind]][[vcovar_ind]]

      opt_param_mat[re_ind, vcovar_ind] <- FUN(
        mc_state = mc_state,
        vcovar_ind = vcovar_ind,
        emit_act = emit_act,
        emit_light = emit_light,
        corr_mat = corr_mat,
        lod_act = lod_act,
        lod_light = lod_light,
        case_cache = case_cache,
        re_ind = re_ind,
        emit_subset = emit_subset
      )
    }
  }

  opt_param_mat
}