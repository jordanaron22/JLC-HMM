################## Intro ################## 

jmhmm_start_time <- Sys.time()

read_memory_peak_gb <- function(){
  status_file <- "/proc/self/status"
  if (!file.exists(status_file)){
    return(NA_real_)
  }

  status_lines <- readLines(status_file,warn = FALSE)
  memory_line <- grep("^VmHWM:",status_lines,value = TRUE)
  if (length(memory_line) == 0){
    return(NA_real_)
  }

  memory_kb <- suppressWarnings(as.numeric(
    sub("^VmHWM:[[:space:]]*([0-9]+)[[:space:]]*kB.*$","\\1",
        memory_line[[1]])
  ))
  memory_kb / 1024^2
}

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)
library(dplyr)
library(numDeriv)
library(Matrix)
library(Hmisc)
library(survex)
#library(tidyverse)
library(pbivnorm)

source_jmhmm_module <- function(module_file){
  candidate_paths <- c(file.path("Scripting","R",module_file),
                       file.path("R",module_file),
                       file.path("..","Rcode","R",module_file))
  for (candidate_path in candidate_paths){
    if (file.exists(candidate_path)){
      source(candidate_path, local = parent.frame())
      return(invisible(candidate_path))
    }
  }
  stop(paste("Could not find module:",module_file))
}

# Load modules in dependency order.
for (module_file in c("constants.R","saved_results.R","validation.R",
                      "settings.R","params.R","transitions.R",
                      "emissions_tobit.R","forward_backward.R",
                      "oakes_info.R","data_simulation.R","helpers.R","data_nhanes.R",
                      "survival.R","diagnostics.R")){
  source_jmhmm_module(module_file)
}

settings <- build_settings()
settings$model_name <- build_model_name(settings)

cli_args <- settings$command_args
sim_num <- settings$sim_num

# Compatibility aliases: most of the legacy script still reads these names.
list2env(settings, envir = environment())
output_dir <- if (dir.exists("Routputs")){"Routputs"} else {"."}

print("Command line arguments:")
print(cli_args)
print("Run settings:")
print(settings)

if (set_seed){set.seed(sim_num)}

print(paste("Sim Seed:",sim_num,"Fit HMM Num:",fit_mix_num,"True HMM Num:",true_mix_num))
print(model_name)

################## EM Setup ################## 

readCpp( "Scripting/cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )

###### True Settings ###### 

#Sets up simulation sizing
day_length <- period_len * simulation_days
num_of_people <- num_people



true_param_list <- CreateDefaultParams(true_mix_num, vcovar_num)
if (!real_data){
  true_param_list <- ApplyEmissionOverlap(true_param_list,emission_overlap)
}
validate_param_list(true_param_list,true_mix_num,vcovar_num,"true_param_list")
init_true <- true_param_list$init
params_tran_array_true <- true_param_list$params_tran_array
emit_act_true <- true_param_list$emit_act
emit_light_true <- true_param_list$emit_light
corr_mat_true <- true_param_list$corr_mat
nu_mat_true <- true_param_list$nu_mat
beta_vec_true <- true_param_list$beta_vec
beta_age_true <- true_param_list$beta_age
lambda_act_mat_true <- true_param_list$lambda_act_mat
lambda_light_mat_true <- true_param_list$lambda_light_mat

#should only differ from in scenario where we simulate data to see if we can recover correct number of LC
fit_param_list <- CreateDefaultParams(fit_mix_num, vcovar_num)
validate_param_list(fit_param_list,fit_mix_num,vcovar_num,"fit_param_list")
init_start <- fit_param_list$init
params_tran_array_start <- fit_param_list$params_tran_array
emit_act_start <- fit_param_list$emit_act
emit_light_start <- fit_param_list$emit_light
corr_mat_start <- fit_param_list$corr_mat
nu_mat_start <- fit_param_list$nu_mat
beta_vec_start <- fit_param_list$beta_vec
lambda_act_mat_start <- fit_param_list$lambda_act_mat
lambda_light_mat_start <- fit_param_list$lambda_light_mat
#loads data in for a hot start
if (load_data){
  model_name_loadin <- "JMHMM"
  load_mix_num <- if (real_data){fit_mix_num} else {true_mix_num}
  folder_name <- paste(load_mix_num)
  if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){
    model_name_loadin <- paste0(model_name_loadin,"NoSurv")
    folder_name <- paste0("NS",folder_name)
  }


  legacy_model_name_loadin <- paste0(model_name_loadin,"Mix",load_mix_num,"Seed",".rda")
  model_name_loadin <- paste0(model_name_loadin,"FitMix",load_mix_num,"Seed",".rda")
  print(paste("Loading",model_name_loadin))
  setwd("Data")
  if (!file.exists(model_name_loadin)){
    model_name_loadin <- legacy_model_name_loadin
    print(paste("Loading legacy",model_name_loadin))
  }
  load(model_name_loadin)
  setwd("..")
  validate_saved_results(to_save,required_sections = c("est_params"),
                         source_name = model_name_loadin)

  # # model_name_loadin <- paste0(model_name_loadin,"Mix",mix_num,"Seed",sim_num,".rda")
  # # print(paste("Loading",model_name_loadin))
  #
  # model_name_loadin <- paste0("Inter",model_name)
  # print(paste("Loading",model_name_loadin))
  # # setwd(folder_name)
  # load(model_name_loadin)
  # setwd("..")
 
  loaded_est_params <- get_saved_section(to_save,"est_params")
  loaded_init <- get_saved_param(loaded_est_params,"init")
  loaded_params_tran_array <- get_saved_param(loaded_est_params,"params_tran_array")
  loaded_emit_act <- get_saved_param(loaded_est_params,"emit_act")
  loaded_emit_light <- get_saved_param(loaded_est_params,"emit_light")
  loaded_corr_mat <- get_saved_param(loaded_est_params,"corr_mat")
  loaded_nu_mat <- get_saved_param(loaded_est_params,"nu_mat")
  loaded_beta_vec <- get_saved_param(loaded_est_params,"beta_vec")
  loaded_lambda_act_mat <- lambda_act_mat_true
  loaded_lambda_light_mat <- lambda_light_mat_true
  surv_coef_true <- get_saved_param(loaded_est_params,"surv_coef")

  re_prob_true <- get_saved_param(loaded_est_params,"re_prob")
  re_prob <- re_prob_true

  loaded_lambda_act_mat <- get_saved_param(loaded_est_params,"lambda_act_mat",
                                           default = loaded_lambda_act_mat,
                                           required = FALSE)
  loaded_lambda_light_mat <- get_saved_param(loaded_est_params,"lambda_light_mat",
                                             default = loaded_lambda_light_mat,
                                             required = FALSE)

  if (!real_data){
    init_true <- loaded_init
    params_tran_array_true <- loaded_params_tran_array
    emit_act_true <- loaded_emit_act
    emit_light_true <- loaded_emit_light
    corr_mat_true <- loaded_corr_mat
    nu_mat_true <- loaded_nu_mat
    beta_vec_true <- loaded_beta_vec
    lambda_act_mat_true <- loaded_lambda_act_mat
    lambda_light_mat_true <- loaded_lambda_light_mat
  }

  if (dim(loaded_init)[1] == fit_mix_num){
    init_start <- loaded_init
    params_tran_array_start <- loaded_params_tran_array
    emit_act_start <- loaded_emit_act
    emit_light_start <- loaded_emit_light
    corr_mat_start <- loaded_corr_mat
    nu_mat_start <- loaded_nu_mat
    beta_vec_start <- loaded_beta_vec
    lambda_act_mat_start <- loaded_lambda_act_mat
    lambda_light_mat_start <- loaded_lambda_light_mat
  }
  
}

###### Simulate Data ###### 
if (!real_data){
  lod_act_true <- -5.809153
  lod_light_true <- -1.560658
  
  lod_act <- lod_act_true
  lod_light <- lod_light_true
  
  beta_covar_sim <- c(0,.6,-.5)
  
  simulated_hmm <- SimulateHMM(day_length,num_of_people,
                               init=init_true,params_tran_array = params_tran_array_true,
                               emit_act = emit_act_true,emit_light = emit_light_true,
                               corr_mat = corr_mat_true,
                               lod_act = lod_act_true,lod_light = lod_light_true,
                               nu_mat = nu_mat_true,
                               beta_age_true = beta_age_true,beta_covar_sim = beta_covar_sim,
                               missing_perc = missing_perc, beta_vec_true = beta_vec_true,
                               lambda_act_mat = lambda_act_mat_true,lambda_light_mat = lambda_light_mat_true,
                               true_mix_num = true_mix_num)
  mc <- simulated_hmm$mc
  act <- simulated_hmm$act
  light <- simulated_hmm$light
  mixture_mat <- simulated_hmm$mixture_mat
  age_vec <- simulated_hmm$age_vec
  nu_covar_mat <- simulated_hmm$nu_covar_mat
  vcovar_mat <-  simulated_hmm$vcovar_mat
  surv_list <- simulated_hmm$survival
  surv_covar_sim <- simulated_hmm$surv_covar_sim
  
  id_sim <- cbind(age_vec,surv_covar_sim-1)
  surv_covar <- list(age_vec,Vec2Mat(surv_covar_sim))
  surv_coef <- list(beta_age_true,beta_covar_sim)
  surv_coef_true <- surv_coef
  combined_covar_mat <- matrix(surv_covar_sim-1,nrow = num_of_people)
  combined_covar_mat <- as.factor(combined_covar_mat)
  
  surv_time <- surv_list$time
  surv_event <- surv_list$event
  
  #in simulated data sample weights are set to 1
  sweights_vec <- rep(1, dim(act)[2])

  
  # Independent test data for out-of-sample survival diagnostics.
  simulated_hmm_test <- SimulateHMM(
    day_length,num_of_people,
    init = init_true,params_tran_array = params_tran_array_true,
    emit_act = emit_act_true,emit_light = emit_light_true,
    corr_mat = corr_mat_true,
    lod_act = lod_act_true,lod_light = lod_light_true,
    nu_mat = nu_mat_true,
    beta_age_true = beta_age_true,beta_covar_sim = beta_covar_sim,
    missing_perc = missing_perc,beta_vec_true = beta_vec_true,
    lambda_act_mat = lambda_act_mat_true,
    lambda_light_mat = lambda_light_mat_true,
    true_mix_num = true_mix_num
  )
  test_data <- list(
    act = simulated_hmm_test$act,
    light = simulated_hmm_test$light,
    nu_covar_mat = simulated_hmm_test$nu_covar_mat,
    vcovar_mat = simulated_hmm_test$vcovar_mat,
    surv_time = simulated_hmm_test$survival$time,
    surv_event = simulated_hmm_test$survival$event,
    surv_covar = list(
      simulated_hmm_test$age_vec,
      Vec2Mat(simulated_hmm_test$surv_covar_sim)
    )
  )
  rm(simulated_hmm_test)
  
} 
###### Read in Data ###### 

if (real_data) {
  nhanes_data <- prepare_nhanes_data(
    period_len = period_len,
    bootstrap = bootstrap,
    leave_out = leave_out,
    sim_num = sim_num,
    single_day = single_day,
    weekend_only = weekend_only,
    load_data = load_data,
    surv_coef_true = get0("surv_coef_true",ifnotfound = NULL),
    data_dir = "Data"
  )
  invisible(list2env(nhanes_data,envir = environment()))
}

###### Initial Settings ###### 
##########
#if doing cv, load in full data values for hot start
if (leave_out){

  model_name_loadin <- "JMHMM"
  if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){model_name_loadin <- paste0(model_name_loadin,"NoSurv")}
  legacy_model_name_loadin <- paste0(model_name_loadin,"Mix",fit_mix_num,"Seed",".rda")
  model_name_loadin <- paste0(model_name_loadin,"FitMix",fit_mix_num,"Seed",".rda")

  full_data_file <- file.path("Data",model_name_loadin)
  if (!file.exists(full_data_file)){
    full_data_file <- file.path("Data",legacy_model_name_loadin)
    print(paste("Loading legacy",legacy_model_name_loadin))
  }
  if (!file.exists(full_data_file)){
    stop(paste("Could not find full-data reference file for CV:",
               full_data_file))
  }
  print(paste("Loading",full_data_file))
  load(full_data_file)
  full_data_to_save <- to_save
  validate_saved_results(to_save,required_sections = c("est_params"),
                         source_name = full_data_file)
  full_data_est_params <- get_saved_section(to_save,"est_params")
  full_data_re_prob <- get_saved_param(full_data_est_params,"re_prob")
  mix_assignment_true <- apply(full_data_re_prob,1,which.max)
  full_data_post_decode <- get_saved_param(
    full_data_est_params,
    "post_decode",
    required = FALSE
  )
  state_reference_available <-
    !is.null(full_data_post_decode) &&
    length(dim(full_data_post_decode)) == 2 &&
    ncol(full_data_post_decode) >= max(leave_out_inds)
  post_decode_collapsed_true <- NULL
  if (state_reference_available){
    post_decode_collapsed_true <-
      full_data_post_decode[,leave_out_inds,drop = FALSE]
  } else {
    warning(
      paste(
        "Full-data post_decode is unavailable or reduced.",
        "Skipping held-out state sensitivity/specificity diagnostics."
      )
    )
  }

  if (load_data){
    full_data_hot_start <- make_hot_start_variables(
      full_data_to_save,
      fit_mix_num = fit_mix_num,
      vcovar_num = vcovar_num,
      source_name = full_data_file,
      lambda_act_default = lambda_act_mat_start,
      lambda_light_default = lambda_light_mat_start
    )
    list2env(full_data_hot_start,envir = environment())
    pi_l_true <- CalcPi(nu_mat_start,nu_covar_mat)
  }
  
}

################## EM ##################

##### randomize starting parameters
# init <- matrix(rep(.5,mix_num*2),ncol = 2)
init <- init_start

params_tran_array <- params_tran_array_start + runif(unlist(length(params_tran_array_start)),-randomize_init*2,randomize_init*2)


emit_act <- emit_act_start + runif(length(unlist(emit_act_start)),-randomize_init,randomize_init)
emit_act[,2,,] <- abs(emit_act[,2,,])

emit_light <- emit_light_start + runif(length(unlist(emit_light_start)),-randomize_init*2,randomize_init*2)
emit_light[,2,,] <- abs(emit_light[,2,,])

#makes sure correlation makes sense
corr_mat <- corr_mat_start + runif(length(unlist(corr_mat_start)),-randomize_init/5,randomize_init/5)
corr_mat[corr_mat>.99] <- .99
corr_mat[corr_mat < -0.99] <- -0.99

#makes sure first val is always reference
beta_vec <- beta_vec_start + runif(mix_num,-randomize_init,randomize_init)
beta_vec[1] <- 0

##### Randomize survival starting values #####
for (i in 1:length(surv_coef)){
  surv_coef[[i]] <-surv_coef_true[[i]]  +  runif(length(surv_coef_true[[i]]),-randomize_init/10,randomize_init/10)
  if (length(surv_coef[[i]]) != 1){surv_coef[[i]][1] <- 0} 
}
surv_coef[[1]] <- surv_coef_true[[1]] + runif(1,-randomize_init/100,randomize_init/100)

#dont randomize as these are very sensitive
nu_mat <- nu_mat_start
lambda_act_mat <- lambda_act_mat_start
lambda_light_mat <- lambda_light_mat_start

start_params <- make_start_param_list(init = init,
                                      params_tran_array = params_tran_array,
                                      emit_act = emit_act,
                                      emit_light = emit_light,
                                      corr_mat = corr_mat,
                                      nu_mat = nu_mat,
                                      beta_vec = beta_vec,
                                      surv_coef = surv_coef,
                                      lambda_act_mat = lambda_act_mat,
                                      lambda_light_mat = lambda_light_mat)
validate_param_list(start_params,fit_mix_num,vcovar_num,"start_params")

time_vec <- c()
pi_l <- CalcPi(nu_mat,nu_covar_mat)
#sets some controls so matrix sizing lines up
if (!leave_out & !load_data){re_prob <- pi_l}
if (load_data & !real_data){re_prob <- pi_l}
if (load_data & period_len != 96){re_prob <- pi_l}
if (!is.null(dim(re_prob)) && ncol(re_prob) != fit_mix_num){re_prob <- pi_l}
if (!is.null(dim(re_prob)) && nrow(re_prob) != nrow(pi_l)){re_prob <- pi_l}
if (is.null(dim(re_prob))){re_prob <- matrix(re_prob,ncol = 1)}

validate_hmm_data(act,light,vcovar_mat)
validate_survival_inputs(surv_time,surv_event,surv_covar,num_of_people)
survival_context <- make_survival_context(surv_time,surv_event,surv_covar,
                                          re_prob,fit_mix_num,sweights_vec)

surv_coef_len <- unlist(lapply(surv_coef,length))
surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)

if(!incl_light){
  light <- matrix(NA,dim(light)[1],dim(light)[2])
  light_old <- matrix(NA,dim(light_old)[1],dim(light_old)[2])
}
if(!incl_act){
  act <- matrix(NA,dim(act)[1],dim(act)[2])
  act_old <- matrix(NA,dim(act_old)[1],dim(act_old)[2])
}

###needed for tobit emission optimization
vcovar_mat_emit <- vcovar_mat + 1
emit_data <- PrepareEmitLogLikeData(
  act = act,
  light = light,
  vcovar_mat = vcovar_mat_emit,
  lod_act = lod_act,
  lod_light = lod_light
)




#calculates baseline hazards
bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                      surv_covar_risk_vec,survival_context$surv_event,
                      survival_context$surv_time,survival_context$surv_covar, survival_context$sweights_vec)
bline_vec <- bhaz_vec[[1]]
cbline_vec <- bhaz_vec[[2]]

#caluclates case4 probabilities ahead of time and transition list
lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                         period_len = period_len)

print("Pre Forward-Backward")

fb_result <- ForwardBackward(
  act = act,
  light = light,
  init = init,
  tran_list = tran_list,
  emit_act = emit_act,
  emit_light = emit_light,
  lod_act = lod_act,
  lod_light = lod_light,
  corr_mat = corr_mat,
  beta_vec = beta_vec,
  surv_coef = surv_coef,
  surv_covar_risk_vec =
    surv_covar_risk_vec,
  event_vec = surv_event,
  bline_vec = bline_vec,
  cbline_vec = cbline_vec,
  lintegral_mat = lintegral_mat,
  surv_covar = surv_covar,
  vcovar_mat = vcovar_mat,
  lambda_act_mat =
    lambda_act_mat,
  lambda_light_mat =
    lambda_light_mat,
  tobit = tobit,
  incl_surv = incl_surv,
  beta_bool = beta_bool,
  mix_num = mix_num,
  vcovar_num = vcovar_num,
  period_len = period_len
)

alpha <- fb_result$alpha
beta <- fb_result$beta

rm(fb_result)

print("Post Forward-Backward")



new_likelihood <- CalcLikelihood(alpha,pi_l,sweights_vec)
if (incl_surv == MODEL_TYPE_CODES[["joint"]] & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,survival_context)}
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood 
#check to make sure all values are the same, simple sanity check
# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)
iter_count <- 1
stop_crit <- BASE_STOP_CRIT
if (!real_data){stop_crit <- stop_crit * SIM_STOP_CRIT_MULTIPLIER}
# if(mix_num > 8){stop_crit <- stop_crit * 10}
# if(mix_num > 12){stop_crit <- stop_crit * 10}
# if(mix_num > 15){stop_crit <- stop_crit * 5}

##### Named EM containers #####
# Fixed inputs are grouped separately from parameters that change during EM.
em_inputs <- list(
  longitudinal = list(
    act = act,
    light = light,
    vcovar_mat = vcovar_mat,
    lod_act = lod_act,
    lod_light = lod_light,
    sweights_vec = sweights_vec
  ),
  survival = list(
    time = surv_time,
    event = surv_event,
    covariates = surv_covar
  ),
  dimensions = list(
    day_length = day_length,
    num_of_people = num_of_people,
    mix_num = mix_num,
    vcovar_num = vcovar_num,
    period_len = period_len
  )
)

em_control <- list(
  convergence_tolerance = stop_crit,
  minimum_iterations = MIN_EM_ITERATIONS,
  interim_save_every = INTERIM_SAVE_EVERY,
  reorder_tolerance_multiplier = REORDER_STOP_CRIT_MULTIPLIER,
  survival_start_iteration = 3,
  run_only_survival = run_only_surv,
  check_transition_update = check_tran
)

# start_params already separates the randomized starting parameters from estimates.
em_initial_state <- list(
  parameters = start_params,
  mixture_probabilities = pi_l,
  random_effect_probabilities = re_prob,
  likelihood = new_likelihood
)

likelihood <- new_likelihood
while((abs(like_diff/likelihood) > em_control$convergence_tolerance |
       iter_count < em_control$minimum_iterations) &
      !em_control$run_only_survival){
  ##### EM iteration: bookkeeping #####
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### E-step: latent-class probabilities #####
  re_prob <- CalcProbRE(alpha,pi_l)
  survival_context <- update_survival_context_re_prob(survival_context,
                                                      re_prob,fit_mix_num)
  
  ##### M-step: mixing and survival parameters #####
  #need model to fit a bit first otherwise may run into some instability
  if(beta_bool){

    nu_mat  <- CalcNu(nu_mat,re_prob,nu_covar_mat,alpha = alpha, sweights_vec = sweights_vec,
                      mix_num = mix_num,num_of_people = num_of_people)
    pi_l <- CalcPi(nu_mat,nu_covar_mat)
    re_prob <- CalcProbRE(alpha,pi_l)
    survival_context <- update_survival_context_re_prob(survival_context,
                                                        re_prob,fit_mix_num)
    
    if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
      #calculates survival coef for JM 
      beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
      beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,
                                    surv_covar_risk_vec,incl_surv,
                                    survival_context,surv_coef_len,fit_mix_num)
      beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],
                                                    surv_coef_len,fit_mix_num)
      beta_vec <- beta_surv_coef_temp_list[[1]]
      surv_coef <- beta_surv_coef_temp_list[[2]]
      beta_se <- beta_surv_coef_se[[2]]

      surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)

      bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                            surv_covar_risk_vec,survival_context$surv_event,
                            survival_context$surv_time,survival_context$surv_covar, survival_context$sweights_vec)
      bline_vec <- bhaz_vec[[1]]
      cbline_vec <- bhaz_vec[[2]]
    }
      
  }
  
  ##### E-step: posterior wake/sleep weights #####
  #calculates wake/sleep probabilities, needed for emission dist estimation
  ##### E-step: posterior wake/sleep weights #####
  # Calculates log posterior wake/sleep probabilities.
  weights_array_list <- CondMarginalize(
    alpha,
    beta,
    pi_l
  )

  # Build each state/class/day-type weighted subset once.
  emit_case_cache <- BuildEmitCaseCache(
    log_weights_by_state = weights_array_list,
    emit_data = emit_data,
    sweights_vec = sweights_vec,
    mix_num = dim(emit_act)[3L],
    vcovar_num = dim(emit_act)[4L]
  )
  # The full log-posterior arrays are no longer needed by the
  # emission updates after the cache has been constructed.
  rm(weights_array_list)  

  
  
  ##### M-step: Tobit emission parameters #####
  ##### M-step: Tobit emission parameters #####

  # Update light-distribution parameters.
  if (incl_light) {

    # Wake-state light mean
    emit_light[1, 1, , ] <- UpdateNormByCase(
      FUN = CalcLightMeanByCase,
      mc_state = 1,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Sleep-state light mean
    emit_light[2, 1, , ] <- UpdateNormByCase(
      FUN = CalcLightMeanByCase,
      mc_state = 2,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Wake-state light standard deviation
    emit_light[1, 2, , ] <- UpdateNormByCase(
      FUN = CalcLightSigByCase,
      mc_state = 1,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Sleep-state light standard deviation
    emit_light[2, 2, , ] <- UpdateNormByCase(
      FUN = CalcLightSigByCase,
      mc_state = 2,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )
  }


  # Update activity-distribution parameters.
  if (incl_act) {

    # Wake-state activity standard deviation
    emit_act[1, 2, , ] <- UpdateNormByCase(
      FUN = CalcActSigByCase,
      mc_state = 1,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Sleep-state activity standard deviation
    emit_act[2, 2, , ] <- UpdateNormByCase(
      FUN = CalcActSigByCase,
      mc_state = 2,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Wake-state activity mean
    emit_act[1, 1, , ] <- UpdateNormByCase(
      FUN = CalcActMeanByCase,
      mc_state = 1,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Sleep-state activity mean
    emit_act[2, 1, , ] <- UpdateNormByCase(
      FUN = CalcActMeanByCase,
      mc_state = 2,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )
  }


  # Correlation is identifiable only when both emissions are included.
  if (incl_act && incl_light) {

    # Wake-state activity-light correlation
    corr_mat[, 1, ] <- UpdateNormByCase(
      FUN = CalcBivarCorrByCase,
      mc_state = 1,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )

    # Sleep-state activity-light correlation
    corr_mat[, 2, ] <- UpdateNormByCase(
      FUN = CalcBivarCorrByCase,
      mc_state = 2,
      emit_act = emit_act,
      emit_light = emit_light,
      corr_mat = corr_mat,
      lod_act = lod_act,
      lod_light = lod_light,
      emit_case_cache = emit_case_cache,
      emit_data = emit_data
    )
  }

  
  ##### M-step: initial-state probabilities #####
  #this only relies on normal parameters so calculate it now
  lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
  init <- CalcInit(alpha,beta,pi_l,sweights_vec)
  
  ##### M-step: transition parameters #####
  #saves old transition values in case likelihood decrease
  params_tran_array_old <- params_tran_array
  #gradient and hessian for tran parameters
  #Old approach
  # tran_gradhess_list <- CalcTranCHelper(alpha,beta,act,light,params_tran_array,
  #                                       emit_act,emit_light,corr_mat,
  #                                       pi_l,lod_act,lod_light,lintegral_mat,vcovar_mat,
  #                                       lambda_act_mat, lambda_light_mat, tobit,
  #                                       em_control$check_transition_update,likelihood,
  #                                       period_len = period_len, sweights_vec = sweights_vec,
  #                                       vcovar_num = vcovar_num)



  #AI improved approach
  #validated against old approach
  #260ish to <20 seconds
  tran_gradhess_list <- CalcTranCHelperFast(alpha,beta,act,light,params_tran_array,
                                        emit_act,emit_light,corr_mat,
                                        pi_l,lod_act,lod_light,lintegral_mat,vcovar_mat,
                                        lambda_act_mat, lambda_light_mat, tobit,
                                        em_control$check_transition_update,likelihood,
                                        period_len = period_len, sweights_vec = sweights_vec,
                                        vcovar_num = vcovar_num)
  

  tran_check_context <- list(day_length = em_inputs$dimensions$day_length,
                             period_len = em_inputs$dimensions$period_len,
                             mix_num = em_inputs$dimensions$mix_num,
                             vcovar_num = em_inputs$dimensions$vcovar_num,
                             act = em_inputs$longitudinal$act,
                             light = em_inputs$longitudinal$light,
                             init = init,
                             emit_act = emit_act,
                             emit_light = emit_light,
                             lod_act = em_inputs$longitudinal$lod_act,
                             lod_light = em_inputs$longitudinal$lod_light,
                             corr_mat = corr_mat,
                             beta_vec = beta_vec,
                             surv_coef = surv_coef,
                             surv_covar_risk_vec = surv_covar_risk_vec,
                             surv_event = em_inputs$survival$event,
                             bline_vec = bline_vec,
                             cbline_vec = cbline_vec,
                             lintegral_mat = lintegral_mat,
                             sweights_vec = em_inputs$longitudinal$sweights_vec,
                             surv_covar = em_inputs$survival$covariates,
                             vcovar_mat = em_inputs$longitudinal$vcovar_mat,
                             lambda_act_mat = lambda_act_mat,
                             lambda_light_mat = lambda_light_mat,
                             tobit = tobit,
                             incl_surv = incl_surv,
                             beta_bool = beta_bool)

  params_tran_array <- LM(tran_gradhess_list[[1]],tran_gradhess_list[[2]],
                          params_tran_array,
                          em_control$check_transition_update,likelihood,pi_l,
                          mix_num = mix_num,vcovar_num = vcovar_num,
                          tran_check_context = tran_check_context)

  tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                           period_len = period_len)

  ##### Likelihood refresh and transition rollback check #####
  fb_result <- ForwardBackward(
    act = act,
    light = light,
    init = init,
    tran_list = tran_list,
    emit_act = emit_act,
    emit_light = emit_light,
    lod_act = lod_act,
    lod_light = lod_light,
    corr_mat = corr_mat,
    beta_vec = beta_vec,
    surv_coef = surv_coef,
    surv_covar_risk_vec =
      surv_covar_risk_vec,
    event_vec = surv_event,
    bline_vec = bline_vec,
    cbline_vec = cbline_vec,
    lintegral_mat = lintegral_mat,
    surv_covar = surv_covar,
    vcovar_mat = vcovar_mat,
    lambda_act_mat =
      lambda_act_mat,
    lambda_light_mat =
      lambda_light_mat,
    tobit = tobit,
    incl_surv = incl_surv,
    beta_bool = beta_bool,
    mix_num = mix_num,
    vcovar_num = vcovar_num,
    period_len = period_len
  )

  alpha <- fb_result$alpha
  beta <- fb_result$beta

  rm(fb_result)
  
  new_likelihood <- CalcLikelihood(alpha,pi_l,sweights_vec)
  
  #if JM but during cold start, dont wan to actually include survial in likelihood yet
  if (incl_surv == MODEL_TYPE_CODES[["joint"]] & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,survival_context)}
  
  like_diff <- new_likelihood - likelihood
  
  if (like_diff < 0){
    #all other parameters are either
    #1) closed form
    #2) we can quickly calculate likelihood difference
    #tran likelihood requires forward algorithm and thus much slower
    #thus any like decrease is from transition parameters
    #effectively just doesnt optimize tran in this step of EM
    print("Transition Likelihood Decrease")

    params_tran_array <- params_tran_array_old
    tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                             period_len = period_len)

    fb_result <- ForwardBackward(
      act = act,
      light = light,
      init = init,
      tran_list = tran_list,
      emit_act = emit_act,
      emit_light = emit_light,
      lod_act = lod_act,
      lod_light = lod_light,
      corr_mat = corr_mat,
      beta_vec = beta_vec,
      surv_coef = surv_coef,
      surv_covar_risk_vec =
        surv_covar_risk_vec,
      event_vec = surv_event,
      bline_vec = bline_vec,
      cbline_vec = cbline_vec,
      lintegral_mat = lintegral_mat,
      surv_covar = surv_covar,
      vcovar_mat = vcovar_mat,
      lambda_act_mat =
        lambda_act_mat,
      lambda_light_mat =
        lambda_light_mat,
      tobit = tobit,
      incl_surv = incl_surv,
      beta_bool = beta_bool,
      mix_num = mix_num,
      vcovar_num = vcovar_num,
      period_len = period_len
    )

    alpha <- fb_result$alpha
    beta <- fb_result$beta

    rm(fb_result)

    new_likelihood <- CalcLikelihood(alpha,pi_l,sweights_vec)
    if (incl_surv == MODEL_TYPE_CODES[["joint"]] & beta_bool == 0){new_likelihood <- new_likelihood - SurvLike(beta_vec,surv_covar_risk_vec,surv_coef,survival_context)}
    like_diff <- new_likelihood - likelihood
  }
  
  print(paste("RE num:",mix_num,"Like:",round(abs(like_diff/likelihood),10)))
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
  # after a few iterations don't have any issues with survival estimation stability
  if (iter_count == em_control$survival_start_iteration){
    beta_bool <- T
    print("Starting to Est Survival and Age Mixing Effect")
  }
  
  iter_count <- iter_count + 1
  
  ##### Interim checkpoint #####
  #saves info every few iterations just in case
  if (iter_count %% em_control$interim_save_every == 0){
    
    tran_df <- ParamsArray2DF(params_tran_array,period_len = period_len,
                              vcovar_num = vcovar_num)
    if (!real_data & true_mix_num == fit_mix_num){
      tran_df_true <- ParamsArray2DF(params_tran_array_true,
                                     period_len = period_len,
                                     vcovar_num = vcovar_num)
      tran_df_truth <- tran_df_true[,1]
      
      tran_df <- tran_df %>% mutate(truth = tran_df_truth)
      tran_df <- tran_df %>% mutate(resid = prob - truth)
    }
    
    
    true_params <- make_true_param_list(init = init_true,
                                        params_tran_array = params_tran_array_true,
                                        emit_act = emit_act_true,
                                        emit_light = emit_light_true,
                                        corr_mat = corr_mat_true,
                                        nu_mat = nu_mat_true,
                                        beta_vec = beta_vec_true,
                                        beta_age = beta_age_true,
                                        lambda_act_mat = lambda_act_mat_true,
                                        lambda_light_mat = lambda_light_mat_true)

    est_params <- make_est_param_list(init = init,
                                      params_tran_array = params_tran_array,
                                      emit_act = emit_act,
                                      emit_light = emit_light,
                                      corr_mat = corr_mat,
                                      nu_mat = nu_mat,
                                      beta_vec = beta_vec,
                                      surv_coef = surv_coef,
                                      tran_df = tran_df,
                                      re_prob = re_prob,
                                      new_likelihood = new_likelihood,
                                      lambda_act_mat = lambda_act_mat,
                                      lambda_light_mat = lambda_light_mat)
    validate_param_list(est_params,fit_mix_num,vcovar_num,"est_params")

    to_save <- make_saved_results(true_params = true_params,
                                  est_params = est_params,
                                  settings = settings,
                                  start_params = start_params)
    
    if(!leave_out){
      save(to_save,file = file.path(output_dir,paste0("Inter",model_name)))
    }
    
      
  }
  
  
  ##### Latent-label alignment #####
  #reorders clusters from best to worst survival
  if ((abs(like_diff/likelihood) <
       em_control$convergence_tolerance *
       em_control$reorder_tolerance_multiplier) &
      !relabel_reset & !bootstrap & !leave_out & real_data){
    relabel_reset <- TRUE
    relabel_bool <- 0
    print("Relabelling")
    print("Potential Soft Reset")
    #### Reorder #####
    #Reorder to avoid label switching
    #Cluster means go from small to large by activity
    
    if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
      reord_inds <- order(beta_vec)
    } else if (incl_surv == MODEL_TYPE_CODES[["two_stage"]]){
      beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
      beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,
                                    surv_covar_risk_vec,incl_surv,
                                    survival_context,surv_coef_len,fit_mix_num)
      beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],
                                                    surv_coef_len,fit_mix_num)
      beta_vec_temp <- beta_surv_coef_temp_list[[1]]
      reord_inds <- order(beta_vec_temp)
    }
    
    # reord_inds <- c(0,rev(order(beta_vec[-1])))+1
    if (!all(reord_inds == c(1:mix_num)) & !leave_out){
      print("Swapping Labels")
      relabel_bool <- 1
      emit_act <- emit_act[,,reord_inds,]
      emit_light <- emit_light[,,reord_inds,]
      nu_mat <- nu_mat[,reord_inds]
      nu_mat <- nu_mat - nu_mat[,1]
      params_tran_array <- params_tran_array[reord_inds,,]
      corr_mat <- corr_mat[reord_inds,,]
      init <- init[reord_inds,]
      beta_vec <- beta_vec[reord_inds]
      beta_vec <- beta_vec-min(beta_vec)

      pi_l <- CalcPi(nu_mat,nu_covar_mat)
      re_prob <- re_prob[,reord_inds]
      survival_context <- update_survival_context_re_prob(survival_context,
                                                          re_prob,fit_mix_num)

      bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                            surv_covar_risk_vec,survival_context$surv_event,
                            survival_context$surv_time,survival_context$surv_covar, survival_context$sweights_vec)
      bline_vec <- bhaz_vec[[1]]
      cbline_vec <- bhaz_vec[[2]]

    }
    
    ##### Align wake/sleep labels within each class and day type #####
    for (re_ind in 1:mix_num){
      for (week_ind in 1:2){
        if (emit_act[2,1,re_ind,week_ind] > emit_act[1,1,re_ind,week_ind]){
          relabel_bool <- 1
          print(paste("Swapping wake/sleep for week_ind",week_ind,"Mixture",re_ind))

          if (week_ind == 1){
            temp <- init[re_ind,1]
            #NO WEEKEND INIT
            #SWAPPING MAY SLIGHTLEY DECREASE LIKE?
            init[re_ind,1] <- init[re_ind,2]
            init[re_ind,2] <- temp
          }

          temp <- emit_act[1,,re_ind,week_ind]
          emit_act[1,,re_ind,week_ind] <- emit_act[2,,re_ind,week_ind]
          emit_act[2,,re_ind,week_ind] <- temp

          temp <- emit_light[1,,re_ind,week_ind]
          emit_light[1,,re_ind,week_ind] <- emit_light[2,,re_ind,week_ind]
          emit_light[2,,re_ind,week_ind] <- temp

          #ISSUE HERE
          temp <- params_tran_array[re_ind,1:3,week_ind]
          params_tran_array[re_ind,1:3,week_ind] <- params_tran_array[re_ind,4:6,week_ind]
          params_tran_array[re_ind,4:6,week_ind] <- temp
          tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                                   period_len = period_len)

          temp <- corr_mat[re_ind,1,week_ind]
          corr_mat[re_ind,1,week_ind] <- corr_mat[re_ind,2,week_ind]
          corr_mat[re_ind,2,week_ind] <- temp
        }
      }
    }
    
    if (relabel_bool){
      tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                               period_len = period_len)
      lintegral_mat <- CalcLintegralMat(emit_act,emit_light,corr_mat,lod_act,lod_light)
      
      fb_result <- ForwardBackward(
        act = act,
        light = light,
        init = init,
        tran_list = tran_list,
        emit_act = emit_act,
        emit_light = emit_light,
        lod_act = lod_act,
        lod_light = lod_light,
        corr_mat = corr_mat,
        beta_vec = beta_vec,
        surv_coef = surv_coef,
        surv_covar_risk_vec =
          surv_covar_risk_vec,
        event_vec = surv_event,
        bline_vec = bline_vec,
        cbline_vec = cbline_vec,
        lintegral_mat = lintegral_mat,
        surv_covar = surv_covar,
        vcovar_mat = vcovar_mat,
        lambda_act_mat =
          lambda_act_mat,
        lambda_light_mat =
          lambda_light_mat,
        tobit = tobit,
        incl_surv = incl_surv,
        beta_bool = beta_bool,
        mix_num = mix_num,
        vcovar_num = vcovar_num,
        period_len = period_len
      )

      alpha <- fb_result$alpha
      beta <- fb_result$beta

      rm(fb_result)
      
      new_likelihood <- CalcLikelihood(alpha,pi_l,sweights_vec)
      like_diff <- em_control$convergence_tolerance * 2
    }
      
  }
} # end EM loop

######
#These probabilities are calculated in the beginning of the previous while loop
#need to recalculate them otherwise we'll get slightly stale values
re_prob <- CalcProbRE(alpha, pi_l)

survival_context <- update_survival_context_re_prob(
  survival_context,
  re_prob,
  fit_mix_num
)

weights_array_list <- CondMarginalize(alpha, beta, pi_l)

weights_array_wake <- exp(weights_array_list[[1]])
weights_array_sleep <- exp(weights_array_list[[2]])
######

likelihood_changes <- diff(likelihood_vec)
em_convergence <- list(
  em_was_run = !em_control$run_only_survival,
  converged = if (em_control$run_only_survival){
    NA
  } else {
    is.finite(like_diff) &&
      abs(like_diff/likelihood) <= em_control$convergence_tolerance &&
      iter_count >= em_control$minimum_iterations
  },
  completed_iterations = length(time_vec),
  iteration_counter = iter_count,
  final_likelihood = new_likelihood,
  final_likelihood_change = like_diff,
  convergence_tolerance = em_control$convergence_tolerance,
  minimum_iterations = em_control$minimum_iterations,
  likelihood_history = likelihood_vec,
  iteration_seconds = time_vec,
  total_iteration_seconds = sum(time_vec),
  likelihood_decrease_count = sum(likelihood_changes < 0,na.rm = TRUE)
)

#if 2-stage model, calculate survival here
if (incl_surv != MODEL_TYPE_CODES[["joint"]]){
  survival_context <- update_survival_context_re_prob(survival_context,
                                                      re_prob,fit_mix_num)
  beta_surv_coef <- IntoBetaSurvCoef(beta_vec,surv_coef,fit_mix_num)
  beta_surv_coef_se <- CalcBeta(beta_surv_coef,combined_covar_mat,
                                surv_covar_risk_vec,incl_surv,
                                survival_context,surv_coef_len,fit_mix_num)
  beta_surv_coef_temp_list <- OutofBetaSurvCoef(beta_surv_coef_se[[1]],
                                                surv_coef_len,fit_mix_num)
  beta_vec <- beta_surv_coef_temp_list[[1]]
  surv_coef <- beta_surv_coef_temp_list[[2]]
  beta_se <- beta_surv_coef_se[[2]]
  
  surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
  
  bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,survival_context$re_prob,
                        surv_covar_risk_vec,survival_context$surv_event,
                        survival_context$surv_time,survival_context$surv_covar, survival_context$sweights_vec)
  bline_vec <- bhaz_vec[[1]]
  cbline_vec <- bhaz_vec[[2]]
}

  
    



# Viterbi decoding is not needed for leave-out or class-selection diagnostics.
if(!leave_out && !class_selection_run){
  decoded_mat <- Viterbi(act,light,vcovar_mat)
  
  
  weights_array_wake_collapsed <- apply(weights_array_wake,c(1,2),sum)
  weights_array_sleep_collapsed <- apply(weights_array_sleep,c(1,2),sum)
  post_decode <- weights_array_wake_collapsed < .5
  
  # post_decode <- weights_array_wake
  # for (i in 1:mix_num){
  #   post_decode[,,i] <- weights_array_wake[,,i] < weights_array_sleep[,,i]
  # }
  
} else {
  decoded_mat <- matrix(NA,2,2)
  post_decode <- matrix(NA,2,2)
  
}

#transition parameters for later analysis
tran_df <- ParamsArray2DF(params_tran_array,period_len = period_len,
                          vcovar_num = vcovar_num)
if (!real_data & true_mix_num == fit_mix_num){
  tran_df_true <- ParamsArray2DF(params_tran_array_true,
                                 period_len = period_len,
                                 vcovar_num = vcovar_num)
  tran_df_truth <- tran_df_true[,1]
  
  tran_df <- tran_df %>% mutate(truth = tran_df_truth)
  tran_df <- tran_df %>% mutate(resid = prob - truth)
}

#concatenates true, starting, and estimated parameters into named lists to save
true_params <- make_true_param_list(init = init_true,
                                    params_tran_array = params_tran_array_true,
                                    emit_act = emit_act_true,
                                    emit_light = emit_light_true,
                                    corr_mat = corr_mat_true,
                                    nu_mat = nu_mat_true,
                                    beta_vec = beta_vec_true,
                                    beta_age = beta_age_true,
                                    lambda_act_mat = lambda_act_mat_true,
                                    lambda_light_mat = lambda_light_mat_true)

est_params <- make_est_param_list(init = init,
                                  params_tran_array = params_tran_array,
                                  emit_act = emit_act,
                                  emit_light = emit_light,
                                  corr_mat = corr_mat,
                                  nu_mat = nu_mat,
                                  beta_vec = beta_vec,
                                  surv_coef = surv_coef,
                                  tran_df = tran_df,
                                  re_prob = re_prob,
                                  new_likelihood = new_likelihood,
                                  decoded_mat = decoded_mat,
                                  lambda_act_mat = lambda_act_mat,
                                  lambda_light_mat = lambda_light_mat,
                                  bline_vec = bline_vec,
                                  cbline_vec = cbline_vec,
                                  beta_se = beta_se,
                                  post_decode = post_decode)
validate_param_list(est_params,fit_mix_num,vcovar_num,"est_params")

#if doing leave 100 out cross validation
#predict cluster assignment using varying levels of information
if (leave_out){
  new_act <- act_old[,leave_out_inds]
  new_light <- light_old[,leave_out_inds]
  new_vcovar_mat <- vcovar_mat_old[,leave_out_inds]
  len <- dim(new_act)[1]
  num_of_people <- dim(new_act)[2]
  new_surv_covar <- SubsetSurvCovar(surv_covar_old,leave_out_inds)
  new_pi_l <- CalcPi(nu_mat,nu_covar_mat_old[leave_out_inds,])

  surv_covar_risk_vec_new <- SurvCovarRiskVec(new_surv_covar,surv_coef)

  surv_event_new <- surv_event_old[leave_out_inds]
  surv_time_new <- surv_time_old[leave_out_inds]

  sweights_vec_new <- sweights_vec_old[leave_out_inds]


  ##### Initialize leave-out diagnostic containers #####
  leave_out_scenario_count <- length(CV_LEAVE_OUT_SCENARIOS)
  empty_list <- vector(mode = "list", length = leave_out_scenario_count)
  for (i in seq_len(leave_out_scenario_count)){
    empty_list[[i]] <- list()
  }
  
  empty_mat_list <- vector(mode = "list", length = leave_out_scenario_count)
  for (i in seq_len(leave_out_scenario_count)){
    empty_mat_list[[i]] <- matrix(0,mix_num,mix_num)
  }
  
  empty_vec_list <- vector(mode = "list", length = leave_out_scenario_count)
  
  empty_mat_sublist <- vector(mode = "list", length = 3)
  for (i in 1:3){
    empty_mat_sublist[[i]] <- matrix(0,2,2)
  }
  
  
  conf_mat_list <- empty_mat_list
  cindex_new_list <- empty_vec_list
  ibs_new_list <- empty_vec_list
  ibs2_new_list <- empty_vec_list
  senspec_list <- empty_list
  senspec_mix_list <- empty_list

  cv_longitudinal_loglik <- NULL
  cv_interval_survival_loglik <- NULL
  cv_cindex_data <- NULL
  cv_survival_predictions <- NULL
  
  #One is only cycle
  #Two is no light
  #Three is no activity
  #Four is no tran
  #Five is only act (no light\tran)
  #Six is standard
  leave_out_types <- if (class_selection_run){
    CV_LEAVE_OUT_SCENARIOS[["standard"]]
  } else {
    unname(CV_LEAVE_OUT_SCENARIOS)
  }
  
  ##### Evaluate each information-removal scenario #####
  for (leave_out_type in leave_out_types){
    ##### Scenario setup #####
    new_act_working <- new_act
    new_light_working <- new_light
    
    if (leave_out_type == 2 | leave_out_type == 5){
      new_light_working <- matrix(NA,nrow = dim(new_act)[1],ncol = dim(new_act)[2])
    } 
    if (leave_out_type == 3){
      new_act_working <- matrix(NA,nrow = dim(new_act)[1],ncol = dim(new_act)[2])
    } 
    
    if (leave_out_type == 4 | leave_out_type == 5){
      tran_list <- GenTranList(array(0,dim = dim(params_tran_array)),c(1:day_length),mix_num,vcovar_num,
                               period_len = period_len)
    } else {
      tran_list <- GenTranList(params_tran_array,c(1:day_length),mix_num,vcovar_num,
                               period_len = period_len)
    }
    
    ##### Forward-backward prediction #####
    # Match held-out observations to their own weekday/weekend sequence.

    fb_result <- ForwardBackward(
      act = new_act_working,
      light = new_light_working,
      init = init,
      tran_list = tran_list,
      emit_act = emit_act,
      emit_light = emit_light,
      lod_act = lod_act,
      lod_light = lod_light,
      corr_mat = corr_mat,
      beta_vec = beta_vec,
      surv_coef = surv_coef,
      surv_covar_risk_vec = surv_covar_risk_vec_new,
      event_vec = numeric(ncol(new_act_working)),
      bline_vec = numeric(ncol(new_act_working)),
      cbline_vec = numeric(ncol(new_act_working)),
      lintegral_mat = lintegral_mat,
      surv_covar = new_surv_covar,
      vcovar_mat = new_vcovar_mat,
      lambda_act_mat = lambda_act_mat,
      lambda_light_mat = lambda_light_mat,
      tobit = tobit,
      incl_surv = MODEL_TYPE_CODES[["two_stage"]],
      beta_bool = FALSE,
      mix_num = mix_num,
      vcovar_num = vcovar_num,
      period_len = period_len
    )

    alpha <- fb_result$alpha
    beta <- fb_result$beta
    rm(fb_result)



    ##### Posterior state and class decoding #####
    weights_array_list <- CondMarginalize(alpha,beta,new_pi_l)
    weights_array_wake <- exp(weights_array_list[[1]])
    weights_array_wake_collapsed <- apply(weights_array_wake,c(1,2),sum)
    post_decode_collapsed <- weights_array_wake_collapsed < .5
    
    if (leave_out_type == 1){
      alpha <- ForwardAlt(post_decode_collapsed,init,tran_list,new_vcovar_mat,
                          mix_num = mix_num)
    } 
    
    re_prob_new <- CalcProbRE(alpha,new_pi_l)
    mix_assignment_pred <- apply(re_prob_new,1,which.max)


    ###########################################################################
    # Primary model-selection diagnostics use standard LHOCV only
    ###########################################################################

    if (leave_out_type == 6) {

      cv_cindex_data <- cbind(
        data.frame(
          cv_fold_id = cv_fold_id,
          leave_out_ind = leave_out_inds,
          participant_id = id_old$SEQN[leave_out_inds],
          surv_time = surv_time_new,
          surv_event = surv_event_new,
          risk_score =
            log(drop(re_prob_new %*% exp(beta_vec))) +
            surv_covar_risk_vec_new,
          sweights_vec = sweights_vec_new,
          predicted_class = mix_assignment_pred,
          reference_class =
            mix_assignment_true[leave_out_inds],
          max_posterior =
            apply(re_prob_new, 1, max)
        ),
        setNames(
          as.data.frame(re_prob_new),
          paste0("class_probability_", seq_len(mix_num))
        )
      )

      #########################################################################
      # Held-out longitudinal likelihood
      #
      # alpha was calculated without survival information, so this is
      # log P(Y_i | fitted training parameters).
      #########################################################################

      cv_longitudinal_loglik <- CalcCVLongitudinalLogLik(
        alpha = alpha,
        pi_l = new_pi_l,
        sweights_vec = sweights_vec_new
      )

      #########################################################################
      # Held-out survival likelihood conditional on held-out actigraphy
      #
      # re_prob_new = P(class | held-out longitudinal data).
      # cbline_vec and surv_time come from the training fold.
      #########################################################################

      cv_interval_survival_loglik <-
        CalcCVIntervalSurvivalLogLik(
          surv_time = surv_time_new,
          surv_event = surv_event_new,
          cbline_vec = cbline_vec,
          beta_vec = beta_vec,
          re_prob = re_prob_new,
          surv_covar_risk_vec = surv_covar_risk_vec_new,
          sweights_vec = sweights_vec_new,

          # surv_time is the training-fold survival-time vector
          baseline_surv_time = surv_time,

          interval_breaks = CV_SURVIVAL_INTERVAL_BREAKS
        )

      cv_ibs_eval_times <- CV_IBS_EVAL_TIMES

      cv_ibs_cbline <- BaselineHazardAtTimes(
        baseline_surv_time = surv_time,
        cbline_vec = cbline_vec,
        prediction_times = cv_ibs_eval_times
      )

      cv_ibs_surv_prob <- CalcS(
        event_time = cv_ibs_eval_times,
        cbline_vec_new = cv_ibs_cbline,
        beta_vec = beta_vec,
        re_prob = re_prob_new,
        surv_covar_risk_vec = surv_covar_risk_vec_new
      )

      stopifnot(
        nrow(cv_ibs_surv_prob) == length(cv_ibs_eval_times),
        ncol(cv_ibs_surv_prob) == length(leave_out_inds)
      )

      cv_survival_predictions <- list(
        cv_fold_id = cv_fold_id,
        leave_out_ind = leave_out_inds,
        surv_time = surv_time_new,
        surv_event = surv_event_new,
        sweights_vec = sweights_vec_new,
        eval_times = cv_ibs_eval_times,
        surv_prob = cv_ibs_surv_prob
      )
    }



    senspec_list[[leave_out_type]] <- list()
    senspec_mix_list[[leave_out_type]] <- list()

    if (state_reference_available && leave_out_type != 3){

      post_decode_collapsed_true_vec <- as.vector(post_decode_collapsed_true)
      post_decode_collapsed_vec <- as.vector(post_decode_collapsed)
      senspec_list[[leave_out_type]] <- empty_mat_sublist
      senspec_mix_list[[leave_out_type]] <- vector(mode = "list", length = 3)
      ind_med <- apply(new_act_working,2,median,na.rm = T)

      ##### Activity-stratified state diagnostics #####
      #1 - below
      #2 - above
      #3 - total
      for (lohi_med in 1:3){
        if (lohi_med == 1){
          valid_inds <- t(t(new_act_working) < ind_med)
          pdc_colnames <-  c("Below-Med Pred Wake","Pred Sleep")
        } else if (lohi_med == 2) {
          valid_inds <- t(t(new_act_working) > ind_med)
          pdc_colnames <-  c("Above-Med Pred Wake","Pred Sleep")
        } else {
          valid_inds <- new_act_working > -Inf
          pdc_colnames <-  c("Pred Wake","Pred Sleep")
        }
        valid_inds_vec <- as.vector(valid_inds)
        
        valid_inds_vec[is.na(valid_inds_vec)] <- F
        
        coda <- c(T,F)
        
        pdc_tab <- table(c(post_decode_collapsed_vec[valid_inds_vec],coda),c(post_decode_collapsed_true_vec[valid_inds_vec],coda))
        diag(pdc_tab) <- diag(pdc_tab) - 1
        rownames(pdc_tab) <- pdc_colnames
        colnames(pdc_tab) <- c("True Wake","True Sleep")
        senspec_list[[leave_out_type]][[lohi_med]] <- pdc_tab
        
        senspec_mix_list[[leave_out_type]][[lohi_med]] <- vector(mode = "list", length = mix_num)
        ##### Class-specific state diagnostics #####
        for (curr_class in 1:mix_num){
          # senspec_mix_list[[leave_out_type]][[curr_class]] <- empty_mat_sublist
          valid_inds_vec_mix <- as.vector(valid_inds[,mix_assignment_pred == curr_class])
          valid_inds_vec_mix[is.na(valid_inds_vec_mix)] <- F
          post_decode_collapsed_true_vec_mix <- as.vector(post_decode_collapsed_true[,mix_assignment_pred == curr_class])
          post_decode_collapsed_vec_mix <- as.vector(post_decode_collapsed[,mix_assignment_pred == curr_class])
          pdc_tab_mix <- table(c(post_decode_collapsed_vec_mix[valid_inds_vec_mix],coda),c(post_decode_collapsed_true_vec_mix[valid_inds_vec_mix],coda))
          diag(pdc_tab_mix) <- diag(pdc_tab_mix) - 1
          senspec_mix_list[[leave_out_type]][[lohi_med]][[curr_class]] <- pdc_tab_mix
        }
         
        

      }
      
      
    }
    
    # senspec_list[[leave_out_type]] <- sens_spec_df
   
    
    ##### Class assignment and survival diagnostics #####
    mix_assignment_pred <- c(mix_assignment_pred,c(1:mix_num))
    mix_assignment_true_ind <- c(mix_assignment_true[leave_out_inds],c(1:mix_num))
    conf_mat_ind <- table(mix_assignment_pred,mix_assignment_true_ind)
    diag(conf_mat_ind) <- diag(conf_mat_ind) - 1
    
    cindex <- CalcCindex(
      surv_time = surv_time_new,
      surv_event = surv_event_new,
      beta_vec = beta_vec,
      re_prob = re_prob_new,
      surv_covar_risk_vec = surv_covar_risk_vec_new,
      sweights_vec = sweights_vec_new
    )


    ibs <- CalcIBS(
      surv_time = surv_time_new,
      surv_event = surv_event_new,
      cbline_vec = cbline_vec,
      beta_vec = beta_vec,
      surv_coef = surv_coef,
      surv_covar = new_surv_covar,
      re_prob = re_prob_new,
      incl_surv = incl_surv,
      mix_assignment = mix_assignment_pred,
      surv_covar_risk_vec = surv_covar_risk_vec_new,
      baseline_surv_time = surv_time
    )

    ibs2 <- CalcIBS2(
      surv_time = surv_time_new,
      surv_event = surv_event_new,
      cbline_vec = cbline_vec,
      beta_vec = beta_vec,
      re_prob = re_prob_new,
      surv_covar_risk_vec = surv_covar_risk_vec_new,
      baseline_surv_time = surv_time
    )
    
    conf_mat_list[[leave_out_type]] <- conf_mat_ind
    cindex_new_list[[leave_out_type]] <-cindex
    ibs_new_list[[leave_out_type]] <- ibs
    ibs2_new_list[[leave_out_type]] <- ibs2

  } # end information-removal scenarios
    
  
  
  leave_out_to_save <- make_leave_out_results(leave_out_inds = leave_out_inds,
                                              conf_mat_list = conf_mat_list,
                                              cindex_new_list = cindex_new_list,
                                              ibs_new_list = ibs_new_list,
                                              senspec_list = senspec_list,
                                              ibs2_new_list = ibs2_new_list,
                                              senspec_mix_list = senspec_mix_list)

  leave_out_to_save$cv_longitudinal_loglik <-
    cv_longitudinal_loglik

  leave_out_to_save$cv_interval_survival_loglik <-
    cv_interval_survival_loglik

  leave_out_to_save$cv_cindex_data <- cv_cindex_data

  leave_out_to_save$cv_survival_predictions <-
    cv_survival_predictions

  leave_out_to_save$cv_fold_id <- cv_fold_id

  leave_out_to_save$cv_fold_count <- if (exists("cv_fold_count")){
    cv_fold_count
  } else {
    NA_integer_
  }

  leave_out_to_save$cv_training_event_counts <-
    if (exists("cv_training_event_counts")){
      cv_training_event_counts
    } else {
      NULL
    }

  leave_out_to_save$leave_out_scenarios <- leave_out_types
  leave_out_to_save$cv_interval_breaks <- CV_SURVIVAL_INTERVAL_BREAKS
  leave_out_to_save$state_reference_available <-
    if (exists("state_reference_available")){
      state_reference_available
    } else {
      FALSE
    }

} else {
  leave_out_to_save <- list()
}

#if not using simulated data or want to save space
if (real_data | save_space){
  simulated_hmm <- list()
}


#mixture predections
mix_assignment <- apply(re_prob,1,which.max)
if (!real_data){
  true_class <- factor(as.vector(mixture_mat),levels = seq_len(true_mix_num))
  fitted_class <- factor(mix_assignment,levels = seq_len(fit_mix_num))
  tab <- table(true_class,fitted_class)
} else {
  fitted_class <- factor(c(mix_assignment,seq_len(fit_mix_num)), levels = seq_len(fit_mix_num))
  tab <- table(fitted_class,fitted_class)
  diag(tab) <- diag(tab) - 1
}
  
#removes some data from saving
if (save_space){
  est_params$decoded_mat <- 0
  est_params$bline_vec <- 0
  est_params$cbline_vec <- 0
  est_params$post_decode <- 0
}

#diagnostics
ibs2 <- CalcIBS2(surv_time,surv_event,cbline_vec,beta_vec,re_prob,surv_covar_risk_vec)
ibs <- CalcIBS(surv_time,surv_event,cbline_vec,beta_vec,surv_coef,surv_covar,
               re_prob,incl_surv,mix_assignment,surv_covar_risk_vec)
cindex <- CalcCindex(
  surv_time,surv_event,beta_vec,re_prob,surv_covar_risk_vec,sweights_vec
)
diagnostics <- make_diagnostics_list(cindex = cindex,
                                     ibs = ibs,
                                     confusion_table = tab,
                                     ibs2 = ibs2)
diagnostics$training <- list(cindex = cindex,ibs = ibs)
diagnostics$convergence <- em_convergence
diagnostics$viterbi_skipped <- leave_out || class_selection_run

if (!real_data){
  final_re_prob <- CalcProbRE(alpha,pi_l)
  diagnostics$class_entropy <- CalcClassEntropy(final_re_prob)

  test_act <- test_data$act
  test_light <- test_data$light
  if (!incl_act){test_act[,] <- NA}
  if (!incl_light){test_light[,] <- NA}

  test_pi_l <- CalcPi(nu_mat,test_data$nu_covar_mat)
  test_tran_list <- GenTranList(
    params_tran_array,seq_len(nrow(test_act)),mix_num,vcovar_num,
    period_len = period_len
  )
  test_surv_covar_risk_vec <- SurvCovarRiskVec(
    test_data$surv_covar,surv_coef
  )
  test_alpha <- Forward(
    act = test_act,light = test_light,
    init = init,tran_list = test_tran_list,
    emit_act = emit_act,emit_light = emit_light,
    lod_act = lod_act,lod_light = lod_light,
    corr_mat = corr_mat,beta_vec = beta_vec,surv_coef = surv_coef,
    surv_covar_risk_vec = test_surv_covar_risk_vec,
    event_vec = numeric(ncol(test_act)),
    bline_vec = numeric(ncol(test_act)),
    cbline_vec = numeric(ncol(test_act)),
    lintegral_mat = lintegral_mat,
    surv_covar = test_data$surv_covar,
    vcovar_mat = test_data$vcovar_mat,
    lambda_act_mat = lambda_act_mat,
    lambda_light_mat = lambda_light_mat,
    tobit = tobit,incl_surv = MODEL_TYPE_CODES[["two_stage"]],
    beta_bool = FALSE,mix_num = mix_num
  )
  test_re_prob <- CalcProbRE(test_alpha,test_pi_l)
  test_mix_assignment <- apply(test_re_prob,1,which.max)
  test_cindex <- CalcCindex(
    test_data$surv_time,test_data$surv_event,beta_vec,test_re_prob,
    test_surv_covar_risk_vec, sweights_vec
  )
  test_ibs <- CalcIBS(
    test_data$surv_time,test_data$surv_event,cbline_vec,beta_vec,
    surv_coef,test_data$surv_covar,test_re_prob,incl_surv,
    test_mix_assignment,test_surv_covar_risk_vec,
    baseline_surv_time = surv_time
  )
  diagnostics$test <- list(cindex = test_cindex,ibs = test_ibs)
}

#save everything
bic_parameter_count <- CountModelParameters(
  mix_num = mix_num,
  vcovar_num = vcovar_num,
  nu_covar_num = ncol(nu_covar_mat),
  incl_act = incl_act,
  incl_light = incl_light,
  joint_model = incl_surv == MODEL_TYPE_CODES[["joint"]],
  surv_coef = surv_coef
)
diagnostics$bic_parameter_count <- bic_parameter_count[["total"]]
diagnostics$bic_parameter_breakdown <- bic_parameter_count
bic_sample_size <- ncol(act)
diagnostics$bic_sample_size <- bic_sample_size
diagnostics$bic_likelihood <- if (incl_surv == MODEL_TYPE_CODES[["joint"]]){
  "joint"
} else {
  "longitudinal_only"
}
bic <- CalcBIC(new_likelihood,bic_sample_size,
               bic_parameter_count[["total"]])
aic <- CalcAIC(new_likelihood,bic_parameter_count[["total"]])
diagnostics$aic_parameter_count <- bic_parameter_count[["total"]]
diagnostics$aic_likelihood <- diagnostics$bic_likelihood
diagnostics$runtime_seconds <- as.numeric(difftime(Sys.time(),
                                                   jmhmm_start_time,
                                                   units = "secs"))
diagnostics$memory_peak_gb <- read_memory_peak_gb()
to_save <- make_saved_results(true_params = true_params,
                              est_params = est_params,
                              bic = bic,
                              leave_out = leave_out_to_save,
                              simulated_hmm = simulated_hmm,
                              diagnostics = diagnostics,
                              settings = settings,
                              start_params = start_params,
                              aic = aic)
# model_name <- paste0("ReRun",model_name)
final_file <- file.path(output_dir,model_name)
inter_file <- file.path(output_dir,paste0("Inter",model_name))
save(to_save,file = final_file)
#deletes inter file if final file exists
if (file.exists(final_file) && file.exists(inter_file)){
  unlink(inter_file)
}
