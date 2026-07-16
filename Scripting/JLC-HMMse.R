jlc_hmmse_start_time <- Sys.time()

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
                       file.path("..","Routputs",model_name))
  for (candidate_path in candidate_paths){
    if (file.exists(candidate_path)){
      return(candidate_path)
    }
  }
  stop(paste("Could not find saved JMHMM output:",model_name))
}

get_survival_baseline_mode <- function(cli_args){
  mode <- if (length(cli_args) >= 16) cli_args[[16]] else "profiled"
  mode <- tolower(trimws(mode))
  if (!mode %in% c("fixed","profiled")){
    stop(paste("survival_baseline_mode must be fixed or profiled;",
               "got:",mode))
  }
  mode
}

make_oakes_output_file <- function(input_file,survival_baseline_mode){
  suffix <- paste0("_oakes_",survival_baseline_mode,"_se.rda")
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

is_missing_saved_matrix <- function(x){
  is.null(x) || (length(x) == 1 && is.numeric(x) && x == 0)
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
  params_tran_array_true <- get_saved_or_stop(true_params,"params_tran_array")
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

repair_information_for_vcov <- function(I_obs_sym,
                                        parameter_map = NULL,
                                        rel_tol = 1e-5,
                                        floor_rel = 1e-7) {
  I_obs_sym <- 0.5 * (I_obs_sym + t(I_obs_sym))

  eig_raw <- eigen(I_obs_sym, symmetric = TRUE)
  max_eval <- max(abs(eig_raw$values))
  min_eval <- min(eig_raw$values)
  rel_min <- min_eval / max_eval
  num_negative <- sum(eig_raw$values <= 0)

  raw_vcov <- safe_solve(I_obs_sym)
  raw_vcov_diag <- if (is.null(raw_vcov)) {
    rep(NA_real_, nrow(I_obs_sym))
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
      se = rep(NA_real_, nrow(I_obs_sym)),
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
    values_repaired <- pmax(eig_raw$values, eigen_floor)

    I_obs_for_vcov <- eig_raw$vectors %*%
      diag(values_repaired, nrow = length(values_repaired)) %*%
      t(eig_raw$vectors)

    I_obs_for_vcov <- 0.5 * (I_obs_for_vcov + t(I_obs_for_vcov))
    repaired <- TRUE
    repair_method <- "eigen_floor"
  } else {
    I_obs_for_vcov <- I_obs_sym
    repaired <- FALSE
    repair_method <- "none"
  }

  eig_repaired <- eigen(I_obs_for_vcov, symmetric = TRUE,
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



build_oakes_data <- function(simulated_hmm,est_params,settings,
                             survival_baseline_mode = "profiled"){
  if (!survival_baseline_mode %in% c("fixed","profiled")){
    stop("survival_baseline_mode must be fixed or profiled")
  }

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
  bline_vec <- get_saved_param(est_params,"bline_vec",required = FALSE)
  cbline_vec <- get_saved_param(est_params,"cbline_vec",required = FALSE)

  act <- simulated_hmm$act
  light <- simulated_hmm$light
  vcovar_mat <- simulated_hmm$vcovar_mat
  nu_covar_mat <- simulated_hmm$nu_covar_mat
  surv_time <- simulated_hmm$survival$time
  surv_event <- simulated_hmm$survival$event
  surv_covar <- list(simulated_hmm$age_vec,
                     Vec2Mat(simulated_hmm$surv_covar_sim))
  log_sweights_vec <- numeric(ncol(act))

  if (!settings$include_light){
    light <- matrix(NA,nrow = nrow(light),ncol = ncol(light))
  }
  if (!settings$include_activity){
    act <- matrix(NA,nrow = nrow(act),ncol = ncol(act))
  }

  validate_hmm_data(act,light,vcovar_mat)
  validate_survival_inputs(surv_time,surv_event,surv_covar,ncol(act))
  validate_re_prob(re_prob,ncol(act),settings$fit_mix_num)

  survival_context <- make_survival_context(surv_time,surv_event,surv_covar,
                                            re_prob,settings$fit_mix_num)

  bline_vec <- if (is_missing_saved_matrix(bline_vec)) NULL else bline_vec
  cbline_vec <- if (is_missing_saved_matrix(cbline_vec)) NULL else cbline_vec
  if (identical(survival_baseline_mode,"fixed") &&
      (is.null(bline_vec) || is.null(cbline_vec))){
    surv_covar_risk_vec <- SurvCovarRiskVec(surv_covar,surv_coef)
    bhaz_vec <- CalcBLHaz(surv_coef,beta_vec,re_prob,
                          surv_covar_risk_vec,surv_event,
                          surv_time,surv_covar)
    bline_vec <- bhaz_vec[[1]]
    cbline_vec <- bhaz_vec[[2]]
  }

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
                        lod_act = -5.809153,
                        lod_light = -1.560658,
                        log_sweights_vec = log_sweights_vec,
                        lambda_act_mat = lambda_act_mat,
                        lambda_light_mat = lambda_light_mat,
                        tobit = settings$tobit,
                        period_len = settings$period_len,
                        nu_covar_mat = nu_covar_mat,
                        survival_context = survival_context,
                        re_prob = re_prob,
                        bline_vec = bline_vec,
                        cbline_vec = cbline_vec,
                        incl_surv = settings$incl_surv,
                        beta_bool = settings$beta_bool,
                        survival_baseline_mode = survival_baseline_mode,
                        profile_maxit = 100L,
                        profile_tol = 1e-9,
                        profile_damping = 1.0)
  )
}

build_h1 <- function(theta_pack,posterior_context,oakes_params,
                     data_context){
  init_h1 <- CalcOakesInitialH1(init = oakes_params$init,
                                init_counts =
                                  posterior_context$init_counts,
                                sleep_state = SLEEP_STATE)

  tran_h1 <- CalcOakesTransitionH1(
    alpha = posterior_context$alpha,
    beta = posterior_context$beta,
    act = data_context$act,
    light = data_context$light,
    params_tran_array = oakes_params$params_tran_array,
    emit_act = oakes_params$emit_act,
    emit_light = oakes_params$emit_light,
    corr_mat = oakes_params$corr_mat,
    pi_l = posterior_context$pi_l,
    lod_act = data_context$lod_act,
    lod_light = data_context$lod_light,
    lintegral_mat = posterior_context$lintegral_mat,
    vcovar_mat = data_context$vcovar_mat,
    lambda_act_mat = data_context$lambda_act_mat,
    lambda_light_mat = data_context$lambda_light_mat,
    tobit = data_context$tobit,
    period_len = data_context$period_len)

  emit_h1 <- CalcOakesEmissionH1(
    alpha = posterior_context$alpha,
    beta = posterior_context$beta,
    pi_l = posterior_context$pi_l,
    act = data_context$act,
    light = data_context$light,
    vcovar_mat = data_context$vcovar_mat,
    emit_act = oakes_params$emit_act,
    emit_light = oakes_params$emit_light,
    corr_mat = oakes_params$corr_mat,
    lod_act = data_context$lod_act,
    lod_light = data_context$lod_light)

  mix_h1 <- CalcOakesMixingH1(
    nu_mat = oakes_params$nu_mat,
    re_prob = posterior_context$re_prob,
    nu_covar_mat = data_context$nu_covar_mat)

  if (identical(data_context$survival_baseline_mode,"fixed")){
    surv_h1 <- CalcOakesSurvivalH1(
      beta_vec = oakes_params$beta_vec,
      surv_coef = oakes_params$surv_coef,
      survival_context = posterior_context$survival_context,
      fit_mix_num = theta_pack$fit_mix_num,
      cbline_vec = posterior_context$cbline_vec)
  } else {
    surv_h1 <- CalcOakesProfiledSurvivalH1(
      beta_vec = oakes_params$beta_vec,
      surv_coef = oakes_params$surv_coef,
      survival_context = posterior_context$survival_context,
      fit_mix_num = theta_pack$fit_mix_num,
      surv_covar_risk_vec = posterior_context$surv_covar_risk_vec)
  }

  H1 <- bdiag(init_h1$hessian,tran_h1$hessian,emit_h1$hessian,
              mix_h1$hessian,surv_h1$hessian)
  H1 <- as.matrix(H1)
  rownames(H1) <- names(theta_pack$theta)
  colnames(H1) <- names(theta_pack$theta)

  list(H1 = H1,
       components = list(initial = init_h1,
                         transition = tran_h1,
                         emission = emit_h1,
                         mixing = mix_h1,
                         survival = surv_h1))
}

safe_solve <- function(A,B = NULL){
  tryCatch(
    if (is.null(B)) solve(A) else solve(A,B),
    error = function(e) NULL
  )
}

matrix_diagnostics <- function(mat,name){
  mat <- as.matrix(mat)
  sym <- 0.5 * (mat + t(mat))
  eigen_values <- tryCatch(
    eigen(sym,symmetric = TRUE,only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  data.frame(
    matrix = name,
    nrow = nrow(mat),
    ncol = ncol(mat),
    finite = all(is.finite(mat)),
    max_abs = max(abs(mat),na.rm = TRUE),
    max_abs_asymmetry = max(abs(mat - t(mat)),na.rm = TRUE),
    min_symmetric_eigen = min(eigen_values,na.rm = TRUE),
    max_symmetric_eigen = max(eigen_values,na.rm = TRUE),
    nonpositive_symmetric_eigen_count = sum(eigen_values <= 0,na.rm = TRUE)
  )
}

schur_survival <- function(I_obs,score,parameter_map){
  surv_idx <- which(parameter_map$block == "survival")
  nuis_idx <- setdiff(seq_len(nrow(I_obs)),surv_idx)
  I_ss <- I_obs[surv_idx,surv_idx,drop = FALSE]
  I_nn <- I_obs[nuis_idx,nuis_idx,drop = FALSE]
  I_ns <- I_obs[nuis_idx,surv_idx,drop = FALSE]

  solve_I_nn_I_ns <- safe_solve(I_nn,I_ns)
  used_near_pd <- FALSE
  near_pd_norm_f <- NA_real_
  near_pd_eigenvalues <- NA

  if (is.null(solve_I_nn_I_ns)){
    I_nn_near <- as.matrix(Matrix::nearPD(I_nn,corr = FALSE)$mat)
    solve_I_nn_I_ns <- solve(I_nn_near,I_ns)
  }

  S_surv <- I_ss - t(I_ns) %*% solve_I_nn_I_ns
  S_surv <- 0.5 * (S_surv + t(S_surv))
  S_surv_for_solve <- S_surv
  S_surv_eigen <- eigen(S_surv,symmetric = TRUE,only.values = TRUE)$values

  vcov_surv <- safe_solve(S_surv_for_solve)
  if (is.null(vcov_surv) || any(!is.finite(vcov_surv)) ||
      any(S_surv_eigen <= 0)){
    near_pd <- Matrix::nearPD(S_surv,corr = FALSE)
    S_surv_for_solve <- as.matrix(near_pd$mat)
    vcov_surv <- solve(S_surv_for_solve)
    used_near_pd <- TRUE
    near_pd_norm_f <- near_pd$normF
    near_pd_eigenvalues <- eigen(S_surv_for_solve,symmetric = TRUE,
                                 only.values = TRUE)$values
  }

  se_surv <- sqrt(diag(vcov_surv))
  s_surv <- score[surv_idx]
  delta_surv <- safe_solve(S_surv_for_solve,s_surv)

  score_check <- data.frame(
    parameter = parameter_map$param_name[surv_idx],
    score = s_surv,
    newton_step = as.numeric(delta_surv),
    se = se_surv,
    step_over_se = as.numeric(delta_surv) / se_surv
  )

  list(surv_idx = surv_idx,
       nuisance_idx = nuis_idx,
       S_surv = S_surv,
       S_surv_for_solve = S_surv_for_solve,
       vcov_surv = vcov_surv,
       se_surv = se_surv,
       score_check = score_check,
       diagnostics = list(used_near_pd = used_near_pd,
                          near_pd_norm_f = near_pd_norm_f,
                          S_surv_eigenvalues = S_surv_eigen,
                          S_surv_for_solve_eigenvalues =
                            if (used_near_pd) near_pd_eigenvalues else
                              S_surv_eigen))
}



settings <- build_settings()
settings$model_name <- build_model_name(settings)

cli_args <- settings$command_args
list2env(settings,envir = environment())

print("Command line arguments:")
print(cli_args)
print("Run settings:")
print(settings)
print(paste("JLC-HMMse input model:",settings$model_name))

compile_cpp_helpers()

survival_baseline_mode <- get_survival_baseline_mode(cli_args)
print(paste("JLC-HMMse survival baseline mode:",
            survival_baseline_mode))

input_file <- find_saved_model_file(settings$model_name)
output_file <- make_oakes_output_file(input_file,survival_baseline_mode)
to_save <- load_to_save(input_file)










if (settings$data_source != DATA_SOURCE[["simulation"]]){
  stop("JLC-HMMse currently supports only simulated data")
}
if (settings$model_type != "joint"){
  stop("JLC-HMMse currently supports only joint simulated models")
}

validate_saved_results(to_save,required_sections = c("true_params",
                                                      "est_params",
                                                      "settings"),
                        source_name = input_file)

saved_settings <- get_saved_section(to_save,"settings",required = TRUE)
for (field in c("sim_num","model_name","fit_mix_num","true_mix_num",
                "simulation_days","num_people","emission_overlap")){
  if (!identical(settings[[field]],saved_settings[[field]])){
    stop(paste("Command-line settings do not match saved setting:",field))
  }
}

est_params <- get_saved_section(to_save,"est_params",required = TRUE)
simulated <- get_simulated_hmm(to_save,settings)
oakes_data <- build_oakes_data(simulated$simulated_hmm,est_params,
                                settings,survival_baseline_mode)
params <- oakes_data$params
data_context <- oakes_data$data_context

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

posterior_context <- RebuildOakesPosteriorContext(
  theta = theta_pack$theta,
  theta_pack = theta_pack,
  data_context = data_context)

oakes_params <- MakeInactiveDayTypesSafe(
  UnpackOakesTheta(theta_pack$theta,theta_pack),
  data_context$vcovar_mat)

h1 <- build_h1(theta_pack,posterior_context,oakes_params,data_context)
H1 <- h1$H1

h2_eps <- 1e-5
h2 <- CalcOakesH2(theta_pack = theta_pack,
                  data_context = data_context,
                  eps = h2_eps)
H2 <- h2$hessian

H_obs <- H1 + H2
I_obs <- -H_obs
I_obs_sym <- 0.5 * (I_obs + t(I_obs))
rownames(I_obs_sym) <- names(theta_pack$theta)
colnames(I_obs_sym) <- names(theta_pack$theta)

vcov_repair <- repair_information_for_vcov(
  I_obs_sym = I_obs_sym,
  parameter_map = theta_pack$parameter_map,
  rel_tol = 1e-5,
  floor_rel = 1e-7
)

I_obs_for_vcov <- vcov_repair$I_obs_for_vcov
vcov <- vcov_repair$vcov
se <- vcov_repair$se

oakes_score <- CalcOakesScore(theta = theta_pack$theta,
                              theta_pack = theta_pack,
                              posterior_context = posterior_context,
                              return_components = TRUE)
posterior_diagnostics <- DiagnosePosteriorContext(posterior_context)
score_diagnostics <- DiagnoseOakesScore(oakes_score)
information_diagnostics <- DiagnoseOakesInformation(H1,H2)
matrix_checks <- rbind(matrix_diagnostics(H1,"H1"),
                        matrix_diagnostics(H2,"H2"),
                        matrix_diagnostics(H_obs,"H_obs"),
                        matrix_diagnostics(I_obs_sym,"I_obs"))
schur <- schur_survival(I_obs_sym,oakes_score$score,
                        theta_pack$parameter_map)

to_save$oakes <- list(
  settings = list(h2_eps = h2_eps,
                  survival_baseline_mode = survival_baseline_mode,
                  input_file = input_file,
                  output_file = output_file,
                  data_regenerated = simulated$regenerated,
                  profile_maxit = data_context$profile_maxit,
                  profile_tol = data_context$profile_tol,
                  profile_damping = data_context$profile_damping),
  theta_pack = theta_pack,
  H1 = H1,
  H2 = H2,
  H_obs = H_obs,
  I_obs = I_obs_sym,
  I_obs_for_vcov = I_obs_for_vcov,
  vcov = vcov,
  se = se,
  parameter_map = theta_pack$parameter_map,
  components = list(H1 = h1$components),
  H2_summary = list(base_score = h2$base_score,
                    eps = h2$eps,
                    parameter_indices = h2$parameter_indices,
                    parameter_map = h2$parameter_map),
  posterior_context_summary = list(
    profile_iterations = posterior_context$profile_iterations,
    profile_error = posterior_context$profile_error),
  survival_schur = schur,
  diagnostics = list(
    posterior = posterior_diagnostics,
    score = score_diagnostics,
    matrix_checks = matrix_checks,
    h2 = h2$diagnostics,
    information = information_diagnostics,
    observed_information_eigenvalues =
      eigen(I_obs_sym,symmetric = TRUE,only.values = TRUE)$values,
    vcov_eigen_repair_used = vcov_repair$repaired,
    vcov_eigen_repair_failed = vcov_repair$failed,
    vcov_repair_method = vcov_repair$repair_method,
    schur = schur$diagnostics,
    vcov_repair = list(
      repaired = vcov_repair$repaired,
      failed = vcov_repair$failed,
      repair_method = vcov_repair$repair_method,
      rel_tol = vcov_repair$rel_tol,
      floor_rel = vcov_repair$floor_rel,
      min_eigen_raw = vcov_repair$min_eigen_raw,
      max_eigen_raw = vcov_repair$max_eigen_raw,
      rel_min_raw = vcov_repair$rel_min_raw,
      num_negative_raw = vcov_repair$num_negative_raw,
      eigen_floor = vcov_repair$eigen_floor,
      min_eigen_repaired = vcov_repair$min_eigen_repaired,
      max_eigen_repaired = vcov_repair$max_eigen_repaired,
      rel_min_repaired = vcov_repair$rel_min_repaired,
      raw_negative_diag_idx = vcov_repair$raw_negative_diag_idx,
      raw_negative_diag_names = vcov_repair$raw_negative_diag_names
    ),
    raw_vcov_diag = vcov_repair$raw_vcov_diag,
    runtime_seconds = as.numeric(difftime(Sys.time(),
                                          jlc_hmmse_start_time,
                                          units = "secs")),
    memory_peak_gb = read_memory_peak_gb()))




save(to_save,file = output_file)
print(paste("Saved Oakes SE output:",output_file))
