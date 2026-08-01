jlc_hmmse_two_stage_start_time <- Sys.time()

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

add_local_lib <- function(path){
  if (dir.exists(path)){
    .libPaths(unique(c(normalizePath(path,winslash = "/",mustWork = TRUE),
                       .libPaths())))
  }
}

for (lib_path in c("Rlib",file.path("..","Rlib"),
                   file.path("..","Rcode","Rlib"))){
  add_local_lib(lib_path)
}

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(survival)
library(numDeriv)
library(Matrix)
library(Hmisc)

parameter_wiggle <- 1e-5
parameter_wiggle <- .005

source_jmhmm_module <- function(module_file){
  candidate_paths <- c(file.path("Scripting","R",module_file),
                       file.path("R",module_file),
                       file.path("..","Rcode","R",module_file))
  for (candidate_path in candidate_paths){
    if (file.exists(candidate_path)){
      source(candidate_path,local = parent.frame())
      return(invisible(candidate_path))
    }
  }
  stop(paste("Could not find module:",module_file))
}

for (module_file in c("constants.R","saved_results.R","validation.R",
                      "settings.R","params.R","transitions.R",
                      "emissions_tobit.R","forward_backward.R",
                      "oakes_info.R","data_simulation.R","helpers.R",
                      "data_nhanes.R","analysis_data.R",
                      "survival.R","diagnostics.R")){
  source_jmhmm_module(module_file)
}

compile_cpp_helpers <- function(){
  candidate_paths <- c(file.path("Scripting","cFunctions.cpp"),
                       "cFunctions.cpp",
                       file.path("..","Rcode","cFunctions.cpp"))
  found_path <- NULL
  for (candidate_path in candidate_paths){
    if (file.exists(candidate_path)){
      found_path <- candidate_path
      break
    }
  }
  if (is.null(found_path)){
    stop("Could not find cFunctions.cpp")
  }
  readCpp(found_path)
}

find_saved_model_file <- function(model_name){
  candidate_paths <- c(model_name,
                       file.path("Routputs",model_name),
                       file.path("Output","Routputs","RoutputsFull",
                                 model_name),
                       file.path("..","Routputs",model_name))
  for (candidate_path in candidate_paths){
    if (file.exists(candidate_path)){
      return(candidate_path)
    }
  }
  stop(paste("Could not find saved JMHMM output:",model_name))
}

parse_eps_arg <- function(cli_args){
  eps <- if (length(cli_args) >= 16) {
    suppressWarnings(as.numeric(cli_args[[16]]))
  } else {
    1e-5
  }
  if (length(eps) != 1 || !is.finite(eps) || eps <= 0){
    stop("Murphy-Topel finite-difference epsilon must be positive")
  }
  eps
}

parse_retain_scores_arg <- function(cli_args){
  if (length(cli_args) < 17){
    return(FALSE)
  }
  value <- tolower(trimws(as.character(cli_args[[17]])))
  if (value %in% c("1","true","t","yes","y")){
    return(TRUE)
  }
  if (value %in% c("0","false","f","no","n","")){
    return(FALSE)
  }
  stop("retain score matrix argument must be TRUE/FALSE or 1/0")
}

format_eps_label <- function(eps){
  label <- format(eps,scientific = TRUE,trim = TRUE)
  gsub("[^0-9A-Za-z-]+","",label)
}

make_murphy_topel_output_file <- function(input_file,eps){
  suffix <- paste0("_two_stage_murphy_topel_eps",
                   format_eps_label(eps),"_se.rda")
  sub("\\.rda$",suffix,input_file,ignore.case = TRUE)
}

load_to_save <- function(input_file){
  load_env <- new.env(parent = emptyenv())
  load(input_file,envir = load_env)
  if (!exists("to_save",envir = load_env,inherits = FALSE)){
    stop(paste("Loaded file does not contain to_save:",input_file))
  }
  get("to_save",envir = load_env,inherits = FALSE)
}

get_saved_or_stop <- function(saved_params,param_name){
  get_saved_param(saved_params,param_name,required = TRUE)
}

regenerate_simulated_hmm <- function(to_save,settings){
  if (!settings$use_seed){
    stop("Reduced-output regeneration requires settings$use_seed = TRUE")
  }
  if (is.na(settings$sim_num)){
    stop("Reduced-output regeneration requires settings$sim_num")
  }

  true_params <- get_saved_section(to_save,"true_params",required = TRUE)
  init_true <- get_saved_or_stop(true_params,"init")
  params_tran_array_true <- get_saved_or_stop(true_params,
                                              "params_tran_array")
  emit_act_true <- get_saved_or_stop(true_params,"emit_act")
  emit_light_true <- get_saved_or_stop(true_params,"emit_light")
  corr_mat_true <- get_saved_or_stop(true_params,"corr_mat")
  nu_mat_true <- get_saved_or_stop(true_params,"nu_mat")
  beta_vec_true <- get_saved_or_stop(true_params,"beta_vec")
  beta_age_true <- get_saved_or_stop(true_params,"beta_age")
  lambda_act_mat_true <- get_saved_or_stop(true_params,"lambda_act_mat")
  lambda_light_mat_true <- get_saved_or_stop(true_params,"lambda_light_mat")

  set.seed(settings$sim_num)
  SimulateHMM(
    day_length = settings$period_len * settings$simulation_days,
    num_of_people = settings$num_people,
    init = init_true,
    params_tran_array = params_tran_array_true,
    emit_act = emit_act_true,
    emit_light = emit_light_true,
    corr_mat = corr_mat_true,
    lod_act = -5.809153,
    lod_light = -1.560658,
    nu_mat = nu_mat_true,
    beta_vec_true = beta_vec_true,
    beta_age_true = beta_age_true,
    beta_covar_sim = c(0,.6,-.5),
    missing_perc = settings$missing_perc,
    lambda_act_mat = lambda_act_mat_true,
    lambda_light_mat = lambda_light_mat_true,
    true_mix_num = settings$true_mix_num
  )
}

get_simulated_hmm <- function(to_save,settings){
  simulated_hmm <- get_saved_section(to_save,"simulated_hmm",
                                     required = FALSE)
  if (!is.null(simulated_hmm) && length(simulated_hmm) > 0 &&
      !is.null(simulated_hmm$act)){
    return(list(simulated_hmm = simulated_hmm,
                regenerated = FALSE))
  }

  list(simulated_hmm = regenerate_simulated_hmm(to_save,settings),
       regenerated = TRUE)
}

get_se_analysis_data <- function(to_save,est_params,settings){
  if (settings$data_source == DATA_SOURCE[["simulation"]]){
    simulated <- get_simulated_hmm(to_save,settings)
    return(list(
      analysis_data = MakeSimulationSEAnalysisData(
        simulated$simulated_hmm
      ),
      data_regenerated = simulated$regenerated,
      data_reconstructed = simulated$regenerated
    ))
  }

  if (settings$data_source == DATA_SOURCE[["nhanes"]]){
    return(list(
      analysis_data = PrepareNHANESSEAnalysisData(settings,est_params),
      data_regenerated = FALSE,
      data_reconstructed = TRUE
    ))
  }

  stop(paste("Unsupported SE data source:",settings$data_source))
}

safe_solve <- function(A,B = NULL){
  tryCatch(
    if (is.null(B)) solve(A) else solve(A,B),
    error = function(e) NULL
  )
}

repair_information_for_vcov <- function(I_obs_sym,
                                        parameter_map = NULL,
                                        rel_tol = 1e-5,
                                        floor_rel = 1e-7) {
  I_obs_sym <- 0.5 * (I_obs_sym + t(I_obs_sym))

  eig_raw <- eigen(I_obs_sym,symmetric = TRUE)
  max_eval <- max(abs(eig_raw$values))
  min_eval <- min(eig_raw$values)
  rel_min <- min_eval / max_eval
  num_negative <- sum(eig_raw$values <= 0)

  raw_vcov <- safe_solve(I_obs_sym)
  raw_vcov_diag <- if (is.null(raw_vcov)) {
    rep(NA_real_,nrow(I_obs_sym))
  } else {
    diag(raw_vcov)
  }

  raw_negative_diag_idx <- which(raw_vcov_diag < 0)

  raw_negative_diag_names <- if (!is.null(parameter_map) &&
                                 length(raw_negative_diag_idx) > 0) {
    paste(parameter_map$block[raw_negative_diag_idx],
          parameter_map$param_name[raw_negative_diag_idx],
          raw_negative_diag_idx,
          sep = ".")
  } else {
    names(raw_vcov_diag)[raw_negative_diag_idx]
  }

  if (rel_min < -rel_tol) {
    return(list(
      I_obs_for_vcov = I_obs_sym,
      vcov = NULL,
      se = rep(NA_real_,nrow(I_obs_sym)),
      repaired = FALSE,
      failed = TRUE,
      repair_method = "failed_rel_min_too_negative",
      rel_tol = rel_tol,
      floor_rel = floor_rel,
      min_eigen_raw = min_eval,
      max_eigen_raw = max_eval,
      rel_min_raw = rel_min,
      num_negative_raw = num_negative,
      eigen_floor = NA_real_,
      min_eigen_repaired = NA_real_,
      max_eigen_repaired = NA_real_,
      rel_min_repaired = NA_real_,
      raw_vcov_diag = raw_vcov_diag,
      raw_negative_diag_idx = raw_negative_diag_idx,
      raw_negative_diag_names = raw_negative_diag_names
    ))
  }

  eigen_floor <- floor_rel * max_eval

  if (min_eval <= eigen_floor) {
    values_repaired <- pmax(eig_raw$values,eigen_floor)

    I_obs_for_vcov <- eig_raw$vectors %*%
      diag(values_repaired,nrow = length(values_repaired)) %*%
      t(eig_raw$vectors)

    I_obs_for_vcov <- 0.5 * (I_obs_for_vcov + t(I_obs_for_vcov))
    repaired <- TRUE
    repair_method <- "eigen_floor"
  } else {
    I_obs_for_vcov <- I_obs_sym
    repaired <- FALSE
    repair_method <- "none"
  }

  eig_repaired <- eigen(I_obs_for_vcov,symmetric = TRUE,
                        only.values = TRUE)$values

  vcov <- solve(I_obs_for_vcov)
  se <- sqrt(diag(vcov))

  list(
    I_obs_for_vcov = I_obs_for_vcov,
    vcov = vcov,
    se = se,
    repaired = repaired,
    failed = FALSE,
    repair_method = repair_method,
    rel_tol = rel_tol,
    floor_rel = floor_rel,
    min_eigen_raw = min_eval,
    max_eigen_raw = max_eval,
    rel_min_raw = rel_min,
    num_negative_raw = num_negative,
    eigen_floor = eigen_floor,
    min_eigen_repaired = min(eig_repaired),
    max_eigen_repaired = max(eig_repaired),
    rel_min_repaired = min(eig_repaired) / max(abs(eig_repaired)),
    raw_vcov_diag = raw_vcov_diag,
    raw_negative_diag_idx = raw_negative_diag_idx,
    raw_negative_diag_names = raw_negative_diag_names
  )
}

repair_covariance_if_tiny_negative <- function(V,rel_tol = 1e-6,
                                               floor_rel = 1e-10){
  V <- 0.5 * (V + t(V))
  eig_raw <- eigen(V,symmetric = TRUE)
  max_eval <- max(abs(eig_raw$values))
  min_eval <- min(eig_raw$values)
  rel_min <- if (max_eval == 0) NA_real_ else min_eval / max_eval

  if (min_eval >= 0){
    return(list(vcov = V,
                raw_vcov = V,
                repaired = FALSE,
                failed = FALSE,
                repair_method = "none",
                min_eigen_raw = min_eval,
                max_eigen_raw = max_eval,
                rel_min_raw = rel_min,
                min_eigen_repaired = min_eval,
                max_eigen_repaired = max_eval,
                rel_min_repaired = rel_min))
  }

  if (!is.na(rel_min) && rel_min < -rel_tol){
    return(list(vcov = V,
                raw_vcov = V,
                repaired = FALSE,
                failed = TRUE,
                repair_method = "failed_rel_min_too_negative",
                min_eigen_raw = min_eval,
                max_eigen_raw = max_eval,
                rel_min_raw = rel_min,
                min_eigen_repaired = NA_real_,
                max_eigen_repaired = NA_real_,
                rel_min_repaired = NA_real_))
  }

  eigen_floor <- floor_rel * max_eval
  values_repaired <- pmax(eig_raw$values,eigen_floor)
  V_repaired <- eig_raw$vectors %*%
    diag(values_repaired,nrow = length(values_repaired)) %*%
    t(eig_raw$vectors)
  V_repaired <- 0.5 * (V_repaired + t(V_repaired))
  eig_repaired <- eigen(V_repaired,symmetric = TRUE,
                        only.values = TRUE)$values

  list(vcov = V_repaired,
       raw_vcov = V,
       repaired = TRUE,
       failed = FALSE,
       repair_method = "eigen_floor",
       min_eigen_raw = min_eval,
       max_eigen_raw = max_eval,
       rel_min_raw = rel_min,
       min_eigen_repaired = min(eig_repaired),
       max_eigen_repaired = max(eig_repaired),
       rel_min_repaired = min(eig_repaired) / max(abs(eig_repaired)))
}

matrix_max_asymmetry <- function(mat){
  mat <- as.matrix(mat)
  vals <- abs(mat - t(mat))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0){
    return(NA_real_)
  }
  max(vals)
}

matrix_eigenvalues <- function(mat){
  mat <- as.matrix(mat)
  if (any(!is.finite(mat))){
    return(rep(NA_real_,nrow(mat)))
  }
  eigen(0.5 * (mat + t(mat)),symmetric = TRUE,
        only.values = TRUE)$values
}

score_sum_diagnostics <- function(score_mat,name){
  score_mat <- as.matrix(score_mat)
  sums <- colSums(score_mat)
  norms <- sqrt(colSums(score_mat^2))
  data.frame(
    score_matrix = name,
    max_abs_col_sum = max(abs(sums),na.rm = TRUE),
    max_relative_col_sum =
      max(abs(sums) / pmax(norms,.Machine$double.eps),na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

build_two_stage_mt_data <- function(analysis_data,est_params,settings){
  init <- get_saved_or_stop(est_params,"init")
  params_tran_array <- get_saved_or_stop(est_params,"params_tran_array")
  emit_act <- get_saved_or_stop(est_params,"emit_act")
  emit_light <- get_saved_or_stop(est_params,"emit_light")
  corr_mat <- get_saved_or_stop(est_params,"corr_mat")
  nu_mat <- get_saved_or_stop(est_params,"nu_mat")
  beta_vec <- get_saved_or_stop(est_params,"beta_vec")
  surv_coef <- get_saved_or_stop(est_params,"surv_coef")
  re_prob <- get_saved_or_stop(est_params,"re_prob")
  lambda_act_mat <- get_saved_or_stop(est_params,"lambda_act_mat")
  lambda_light_mat <- get_saved_or_stop(est_params,"lambda_light_mat")

  act <- analysis_data$act
  light <- analysis_data$light
  vcovar_mat <- analysis_data$vcovar_mat
  nu_covar_mat <- analysis_data$nu_covar_mat
  surv_time <- analysis_data$surv_time
  surv_event <- analysis_data$surv_event
  surv_covar <- analysis_data$surv_covar
  combined_covar_mat <- analysis_data$combined_covar_mat
  sweights_vec <- analysis_data$sweights_vec
  lod_act <- analysis_data$lod_act
  lod_light <- analysis_data$lod_light

  if (!settings$include_light){
    light <- matrix(NA,nrow = nrow(light),ncol = ncol(light))
  }
  if (!settings$include_activity){
    act <- matrix(NA,nrow = nrow(act),ncol = ncol(act))
  }

  validate_hmm_data(act,light,vcovar_mat)
  validate_survival_inputs(surv_time,surv_event,surv_covar,ncol(act))
  validate_re_prob(re_prob,ncol(act),settings$fit_mix_num)
  ValidateSEAnalysisData(analysis_data,settings$fit_mix_num)

  survival_context <- make_survival_context(surv_time,surv_event,surv_covar,
                                            re_prob,settings$fit_mix_num,
                                            sweights_vec = sweights_vec)

  list(
    params = list(init = init,
                  params_tran_array = params_tran_array,
                  emit_act = emit_act,
                  emit_light = emit_light,
                  corr_mat = corr_mat,
                  nu_mat = nu_mat,
                  beta_vec = beta_vec,
                  surv_coef = surv_coef,
                  lambda_act_mat = lambda_act_mat,
                  lambda_light_mat = lambda_light_mat),
    data_context = list(act = act,
                        light = light,
                        vcovar_mat = vcovar_mat,
                         lod_act = lod_act,
                         lod_light = lod_light,
                        sweights_vec = sweights_vec,
                        lambda_act_mat = lambda_act_mat,
                        lambda_light_mat = lambda_light_mat,
                        tobit = settings$tobit,
                        period_len = settings$period_len,
                        nu_covar_mat = nu_covar_mat,
                        survival_context = survival_context,
                        re_prob = re_prob,
                        incl_surv = settings$incl_surv,
                        beta_bool = settings$beta_bool),
    combined_covar_mat = combined_covar_mat
  )
}

build_simulation_gamma_truth <- function(coef_names,true_params,fit_mix_num){
  true_beta <- as.numeric(get_saved_or_stop(true_params,"beta_vec"))
  beta_age <- as.numeric(get_saved_or_stop(true_params,"beta_age"))
  truth <- rep(NA_real_,length(coef_names))
  names(truth) <- coef_names

  if (length(truth) >= 1){
    truth[[1]] <- beta_age[[1]]
  }
  if (fit_mix_num > 1){
    class_pos <- 2:fit_mix_num
    truth[class_pos] <- true_beta[class_pos]
  }
  if (length(truth) > fit_mix_num){
    extra_n <- length(truth) - fit_mix_num
    beta_covar_sim <- c(.6,-.5)
    if (extra_n > length(beta_covar_sim)){
      stop("More simulated Cox covariate coefficients than known truth values")
    }
    truth[(fit_mix_num + 1):length(truth)] <- beta_covar_sim[seq_len(extra_n)]
  }

  truth
}

check_dimensions <- function(I1,V1,D,S1,S2,C12,V_MT){
  p <- nrow(I1)
  q <- nrow(D)
  n <- nrow(S1)
  data.frame(
    quantity = c("I1","V1","D","S1","S2","C12","V_MT"),
    observed = c(paste(dim(I1),collapse = "x"),
                 paste(dim(V1),collapse = "x"),
                 paste(dim(D),collapse = "x"),
                 paste(dim(S1),collapse = "x"),
                 paste(dim(S2),collapse = "x"),
                 paste(dim(C12),collapse = "x"),
                 paste(dim(V_MT),collapse = "x")),
    expected = c(paste(c(p,p),collapse = "x"),
                 paste(c(p,p),collapse = "x"),
                 paste(c(q,p),collapse = "x"),
                 paste(c(n,p),collapse = "x"),
                 paste(c(n,q),collapse = "x"),
                 paste(c(p,q),collapse = "x"),
                 paste(c(q,q),collapse = "x")),
    passed = c(all(dim(I1) == c(p,p)),
               all(dim(V1) == c(p,p)),
               all(dim(D) == c(q,p)),
               all(dim(S1) == c(n,p)),
               all(dim(S2) == c(n,q)),
               all(dim(C12) == c(p,q)),
               all(dim(V_MT) == c(q,q))),
    stringsAsFactors = FALSE
  )
}

settings <- build_settings()
settings$model_name <- build_model_name(settings)

cli_args <- settings$command_args
list2env(settings,envir = environment())

print("Command line arguments:")
print(cli_args)
print("Run settings:")
print(settings)
print(paste("JLC-HMMse-two-stage input model:",settings$model_name))

compile_cpp_helpers()

eps <- parse_eps_arg(cli_args)
retain_score_matrices <- parse_retain_scores_arg(cli_args)
print(paste("Two-stage Murphy-Topel epsilon:",eps))
print(paste("Retain S1/S2 matrices:",retain_score_matrices))

input_file <- find_saved_model_file(settings$model_name)
output_file <- make_murphy_topel_output_file(input_file,eps)
to_save <- load_to_save(input_file)

if (settings$model_type != "two_stage"){
  stop("JLC-HMMse-two-stage supports only two-stage models")
}

validate_saved_results(to_save,required_sections = c("est_params","settings"),
                        source_name = input_file)

saved_settings <- get_saved_section(to_save,"settings",required = TRUE)
settings_fields <- c("sim_num","model_name","fit_mix_num","model_type",
                     "data_source","period_len","run_bootstrap",
                     "run_leave_one_out_cv","target_weekday","weekend_only",
                     "include_activity","include_light")
if (settings$data_source == DATA_SOURCE[["simulation"]]){
  settings_fields <- c(settings_fields,"true_mix_num","simulation_days",
                       "num_people","emission_overlap")
}
for (field in settings_fields){
  if (!identical(settings[[field]],saved_settings[[field]])){
    stop(paste("Command-line settings do not match saved setting:",field))
  }
}

est_params <- get_saved_section(to_save,"est_params",required = TRUE)
true_params <- if (settings$data_source == DATA_SOURCE[["simulation"]]){
  get_saved_section(to_save,"true_params",required = TRUE)
} else {
  NULL
}
analysis <- get_se_analysis_data(to_save,est_params,settings)
analysis_data <- analysis$analysis_data
mt_data <- build_two_stage_mt_data(analysis_data,est_params,
                                   settings)
params <- mt_data$params
data_context <- mt_data$data_context
combined_covar_mat <- mt_data$combined_covar_mat

theta_pack <- PackOakesTheta(
  init = params$init,
  params_tran_array = params$params_tran_array,
  emit_act = params$emit_act,
  emit_light = params$emit_light,
  corr_mat = params$corr_mat,
  nu_mat = params$nu_mat,
  beta_vec = params$beta_vec,
  surv_coef = params$surv_coef,
  vcovar_mat = data_context$vcovar_mat,
  fit_mix_num = settings$fit_mix_num)

long_idx <- GetLongitudinalIndex(theta_pack)
long_map <- theta_pack$parameter_map[long_idx,,drop = FALSE]

posterior_context <- RebuildLongitudinalPosteriorContext(
  theta = theta_pack$theta,
  theta_pack = theta_pack,
  data_context = data_context)

posterior_max_abs_diff <- max(abs(posterior_context$re_prob -
                                    data_context$re_prob),na.rm = TRUE)
if (!is.finite(posterior_max_abs_diff) ||
    posterior_max_abs_diff > parameter_wiggle){
  stop(paste("Rebuilt longitudinal re_prob differs from saved re_prob;",
             "max abs diff =",posterior_max_abs_diff))
}

cox_fit_base <- FitTwoStageCox(
  re_prob = posterior_context$re_prob,
  survival_context = data_context$survival_context,
  combined_covar_mat = combined_covar_mat,
  ties = "breslow")

gamma_hat <- stats::coef(cox_fit_base)
saved_gamma <- IntoBetaSurvCoef(params$beta_vec,params$surv_coef,
                                settings$fit_mix_num)
if (length(saved_gamma) != length(gamma_hat)){
  stop("Saved survival coefficient vector length does not match Cox fit")
}
names(saved_gamma) <- names(gamma_hat)
cox_coef_max_abs_diff <- max(abs(gamma_hat - saved_gamma),na.rm = TRUE)
if (!is.finite(cox_coef_max_abs_diff) || cox_coef_max_abs_diff > parameter_wiggle){
  stop(paste("Base Cox coefficients differ from saved two-stage estimates;",
             "max abs diff =",cox_coef_max_abs_diff))
}

oakes_params <- MakeInactiveDayTypesSafe(
  UnpackOakesTheta(theta_pack$theta,theta_pack),
  data_context$vcovar_mat)

h1 <- BuildLongitudinalH1(theta_pack = theta_pack,
                          posterior_context = posterior_context,
                          oakes_params = oakes_params,
                          data_context = data_context,
                          long_idx = long_idx)
H1_long <- h1$H1

mt_perturb <- CalcTwoStageOakesMurphyTopel(
  theta_pack = theta_pack,
  data_context = data_context,
  long_idx = long_idx,
  cox_fit_base = cox_fit_base,
  combined_covar_mat = combined_covar_mat,
  eps = eps,
  progress = TRUE)
H2_long <- mt_perturb$H2_long
S1 <- mt_perturb$S1
D <- mt_perturb$D

H_obs_long <- H1_long + H2_long
I1_raw <- -0.5 * (H_obs_long + t(H_obs_long))
rownames(I1_raw) <- names(theta_pack$theta)[long_idx]
colnames(I1_raw) <- names(theta_pack$theta)[long_idx]

I1_repair <- repair_information_for_vcov(
  I_obs_sym = I1_raw,
  parameter_map = long_map,
  rel_tol = 1e-5,
  floor_rel = 1e-7
)
if (isTRUE(I1_repair$failed) || is.null(I1_repair$vcov)){
  stop("Longitudinal observed information repair failed")
}
I1_for_vcov <- I1_repair$I_obs_for_vcov
V1 <- I1_repair$vcov

V_naive <- as.matrix(cox_fit_base$var)
I2 <- safe_solve(V_naive)
if (is.null(I2)){
  stop("Could not invert naive Cox covariance")
}
rownames(V_naive) <- names(gamma_hat)
colnames(V_naive) <- names(gamma_hat)
rownames(I2) <- names(gamma_hat)
colnames(I2) <- names(gamma_hat)

S2 <- residuals(cox_fit_base,type = "score",weighted = TRUE)
S2 <- as.matrix(S2)
if (ncol(S2) != length(gamma_hat) && nrow(S2) == length(gamma_hat)){
  S2 <- t(S2)
}
if (!all(dim(S2) == c(ncol(data_context$act),length(gamma_hat)))){
  stop("Cox score residual dimensions do not match participants and coefficients")
}
if (any(!is.finite(S2))){
  stop("Cox score residuals contain nonfinite entries")
}
colnames(S2) <- names(gamma_hat)

C12 <- crossprod(S1,S2)
rownames(C12) <- colnames(S1)
colnames(C12) <- colnames(S2)

C_if <- V1 %*% C12 %*% V_naive
rownames(C_if) <- colnames(S1)
colnames(C_if) <- names(gamma_hat)

delta_increment <- D %*% V1 %*% t(D)
cross_increment <- D %*% C_if + t(C_if) %*% t(D)
rownames(delta_increment) <- names(gamma_hat)
colnames(delta_increment) <- names(gamma_hat)
rownames(cross_increment) <- names(gamma_hat)
colnames(cross_increment) <- names(gamma_hat)

V_delta <- V_naive + delta_increment
V_MT_raw <- V_delta + cross_increment
V_delta <- 0.5 * (V_delta + t(V_delta))
V_MT_raw <- 0.5 * (V_MT_raw + t(V_MT_raw))
rownames(V_delta) <- names(gamma_hat)
colnames(V_delta) <- names(gamma_hat)
rownames(V_MT_raw) <- names(gamma_hat)
colnames(V_MT_raw) <- names(gamma_hat)

V_MT_repair <- repair_covariance_if_tiny_negative(V_MT_raw)
if (isTRUE(V_MT_repair$failed)){
  stop("Murphy-Topel covariance repair failed")
}
V_MT <- V_MT_repair$vcov
rownames(V_MT) <- names(gamma_hat)
colnames(V_MT) <- names(gamma_hat)

se_naive <- sqrt(diag(V_naive))
se_delta <- sqrt(diag(V_delta))
se_MT <- sqrt(diag(V_MT))
if (any(!is.finite(se_MT))){
  stop("Murphy-Topel standard errors contain nonfinite values")
}

gamma_truth <- if (settings$data_source == DATA_SOURCE[["simulation"]]){
  build_simulation_gamma_truth(names(gamma_hat),true_params,
                               settings$fit_mix_num)
} else {
  stats::setNames(rep(NA_real_,length(gamma_hat)),names(gamma_hat))
}

V_MT_zero_C12 <- V_delta
V_delta_zero_D <- V_naive
V_MT_zero_D <- V_naive

dimension_checks <- check_dimensions(I1_for_vcov,V1,D,S1,S2,C12,V_MT)
symmetry_checks <- data.frame(
  matrix = c("I1","V_naive","V_delta","V_MT"),
  max_abs_asymmetry = c(matrix_max_asymmetry(I1_for_vcov),
                        matrix_max_asymmetry(V_naive),
                        matrix_max_asymmetry(V_delta),
                        matrix_max_asymmetry(V_MT)),
  stringsAsFactors = FALSE
)
score_checks <- rbind(score_sum_diagnostics(S1,"S1"),
                      score_sum_diagnostics(S2,"S2"))
formula_checks <- data.frame(
  check = c("zero_C12_V_MT_equals_V_delta",
            "zero_D_V_delta_equals_V_naive",
            "zero_D_V_MT_equals_V_naive"),
  max_abs_diff = c(max(abs(V_MT_zero_C12 - V_delta),na.rm = TRUE),
                   max(abs(V_delta_zero_D - V_naive),na.rm = TRUE),
                   max(abs(V_MT_zero_D - V_naive),na.rm = TRUE)),
  stringsAsFactors = FALSE
)

to_save$murphy_topel <- list(
  variance_method = "weighted_murphy_topel",
  settings = list(h2_eps = eps,
                  cox_ties = "breslow",
                  cox_score_residuals_weighted = TRUE,
                  variance_method = "weighted_murphy_topel",
                  input_file = input_file,
                  output_file = output_file,
                  data_source = settings$data_source,
                  data_regenerated = analysis$data_regenerated,
                  data_reconstructed = analysis$data_reconstructed,
                  weight_summary =
                    SummarizeSEWeights(data_context$sweights_vec),
                  retain_score_matrices = retain_score_matrices),
  coefficient_names = names(gamma_hat),
  gamma_hat = gamma_hat,
  coefficient_truth = gamma_truth,
  saved_gamma = saved_gamma,
  saved_gamma_max_abs_diff = cox_coef_max_abs_diff,
  posterior_re_prob_max_abs_diff = posterior_max_abs_diff,
  longitudinal_parameter_map = long_map,
  h2_eps = eps,
  theta_pack = theta_pack,
  H1_long = H1_long,
  H2_long = H2_long,
  H_obs_long = H_obs_long,
  I1_raw = I1_raw,
  I1_for_vcov = I1_for_vcov,
  V1 = V1,
  I2 = I2,
  V_naive = V_naive,
  D = D,
  S1 = if (retain_score_matrices) S1 else NULL,
  S2 = if (retain_score_matrices) S2 else NULL,
  C12 = C12,
  C_if = C_if,
  delta_increment = delta_increment,
  cross_increment = cross_increment,
  V_delta = V_delta,
  V_MT_raw = V_MT_raw,
  V_MT = V_MT,
  se_naive = se_naive,
  se_delta = se_delta,
  se_MT = se_MT,
  components = list(H1 = h1$components),
  perturbation_summary = list(
    diagnostics = mt_perturb$diagnostics,
    score_base = mt_perturb$base_score
  ),
  diagnostics = list(
    posterior = DiagnosePosteriorContext(posterior_context),
    information = DiagnoseOakesInformation(H1_long,H2_long),
    I1_repair = I1_repair,
    V_MT_repair = V_MT_repair,
    dimension_checks = dimension_checks,
    symmetry_checks = symmetry_checks,
    score_checks = score_checks,
    formula_checks = formula_checks,
    perturbation = mt_perturb$diagnostics,
    I1_eigenvalues = matrix_eigenvalues(I1_raw),
    V_delta_eigenvalues = matrix_eigenvalues(V_delta),
    V_MT_raw_eigenvalues = matrix_eigenvalues(V_MT_raw),
    V_MT_eigenvalues = matrix_eigenvalues(V_MT),
    runtime_seconds = as.numeric(difftime(Sys.time(),
                                          jlc_hmmse_two_stage_start_time,
                                          units = "secs")),
    memory_peak_gb = read_memory_peak_gb()
  )
)

save(to_save,file = output_file)
print(paste("Saved two-stage Murphy-Topel SE output:",output_file))
