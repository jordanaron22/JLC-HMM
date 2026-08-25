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
                       file.path("Output","Routputs","Routputs",model_name),
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

# Repair the longitudinal information before it is used in either Murphy-Topel
# target solve.  The diagonal scaling makes the near-PD projection insensitive
# to the units of the longitudinal parameters.  The unscaled raw and repaired
# matrices are retained for diagnostics.
solve_information_for_targets <- function(
    I1,D,C12){
  I1 <- 0.5 * (as.matrix(I1) + t(as.matrix(I1)))
  D <- as.matrix(D)
  C12 <- as.matrix(C12)

  p <- nrow(I1)
  if (ncol(I1) != p || ncol(D) != p || nrow(C12) != p ||
      ncol(C12) != nrow(D)){
    stop("Target information solve received incompatible dimensions")
  }
  if (any(!is.finite(I1)) || any(!is.finite(D)) ||
      any(!is.finite(C12))){
    stop("Target information solve received nonfinite values")
  }

  eig_raw <- eigen(I1,symmetric = TRUE)
  max_eval <- max(abs(eig_raw$values))
  if (!is.finite(max_eval) || max_eval <= 0){
    return(list(
      failed = TRUE,
      failure_reason = "longitudinal_information_has_no_finite_scale",
      I1_for_solve = NULL,
      solve_D = NULL,
      solve_C12 = NULL,
      diagnostics = list(method = "failed",
                         used_target_restriction = FALSE)
    ))
  }

  near_pd_formals <- formals(Matrix::nearPD)
  near_pd_eig_tol <- as.numeric(near_pd_formals$eig.tol)
  near_pd_conv_tol <- as.numeric(near_pd_formals$conv.tol)
  near_pd_posd_tol <- as.numeric(near_pd_formals$posd.tol)

  max_positive_eval <- max(eig_raw$values)
  eigen_cutoff <- if (max_positive_eval > 0){
    near_pd_eig_tol * max_positive_eval
  } else {
    0
  }
  retained_raw <- eig_raw$values > eigen_cutoff

  raw_condition_number <- if (min(abs(eig_raw$values)) == 0){
      Inf
    } else {
      max_eval / min(abs(eig_raw$values))
    }

  diag_raw <- diag(I1)
  use_diagonal_scaling <- all(is.finite(diag_raw)) && all(diag_raw > 0)
  information_scale <- if (use_diagonal_scaling){
    sqrt(diag_raw)
  } else {
    rep(1,nrow(I1))
  }
  I1_scaled <- sweep(sweep(I1,1,information_scale,"/"),
                     2,information_scale,"/")
  I1_scaled <- 0.5 * (I1_scaled + t(I1_scaled))
  eig_scaled_raw <- eigen(I1_scaled,symmetric = TRUE)
  max_scaled_positive_eval <- max(eig_scaled_raw$values)
  scaled_eigen_cutoff <- if (max_scaled_positive_eval > 0){
    near_pd_eig_tol * max_scaled_positive_eval
  } else {
    0
  }
  retained_scaled <- eig_scaled_raw$values > scaled_eigen_cutoff
  discarded_scaled <- !retained_scaled

  # R is the information in standardized parameter coordinates.  D must use
  # the corresponding coordinates for this diagnostic projection.
  D_scaled <- sweep(D,2,information_scale,"/")
  target_norm <- sqrt(sum(D_scaled^2))
  discarded_target_norm <- if (any(discarded_scaled)){
    sqrt(sum((D_scaled %*%
                eig_scaled_raw$vectors[,discarded_scaled,drop = FALSE])^2))
  } else {
    0
  }
  target_projection_ratio <- if (target_norm == 0){
    0
  } else {
    discarded_target_norm / target_norm
  }

  # Use the same relative classification as nearPD(). This also regularizes
  # positive but numerically null directions that nearPD() treats as zero.
  repair_needed <- min(eig_scaled_raw$values) <= scaled_eigen_cutoff

  near_pd_converged <- NA
  near_pd_iterations <- NA_integer_
  near_pd_norm_f_scaled <- 0
  I1_scaled_for_solve <- I1_scaled
  repair_method <- "none"

  if (repair_needed){
    near_pd <- tryCatch(
      Matrix::nearPD(
        I1_scaled,
        corr = use_diagonal_scaling,
        keepDiag = use_diagonal_scaling
      ),
      error = function(e) NULL
    )
    if (is.null(near_pd) || !isTRUE(near_pd$converged)){
      return(list(
        failed = TRUE,
        failure_reason = "longitudinal_information_nearPD_failed",
        I1_for_solve = NULL,
        solve_D = NULL,
        solve_C12 = NULL,
        diagnostics = list(
          method = "nearPD_failed",
          repaired = FALSE,
          full_information_positive = all(eig_raw$values > 0),
          min_eigen_raw = min(eig_raw$values),
          max_eigen_raw = max(eig_raw$values),
          rel_min_raw = min(eig_raw$values) / max_eval,
          target_projection_ratio = target_projection_ratio,
          near_pd_eig_tol = near_pd_eig_tol,
          near_pd_conv_tol = near_pd_conv_tol,
          near_pd_posd_tol = near_pd_posd_tol
        )
      ))
    }
    I1_scaled_for_solve <- as.matrix(near_pd$mat)
    I1_scaled_for_solve <- 0.5 *
      (I1_scaled_for_solve + t(I1_scaled_for_solve))
    near_pd_converged <- near_pd$converged
    near_pd_iterations <- near_pd$iterations
    near_pd_norm_f_scaled <- near_pd$normF
    repair_method <- if (use_diagonal_scaling){
      "scaled_nearPD_information_solve"
    } else {
      "nearPD_information_solve"
    }
  }

  I1_for_solve <- sweep(sweep(I1_scaled_for_solve,1,information_scale,"*"),
                        2,information_scale,"*")
  I1_for_solve <- 0.5 * (I1_for_solve + t(I1_for_solve))
  dimnames(I1_for_solve) <- dimnames(I1)
  eig_repaired <- eigen(I1_for_solve,symmetric = TRUE,
                        only.values = TRUE)$values
  eig_scaled_repaired <- eigen(I1_scaled_for_solve,symmetric = TRUE,
                               only.values = TRUE)$values
  if (any(eig_scaled_repaired <= 0)){
    return(list(
      failed = TRUE,
      failure_reason = "longitudinal_information_not_positive_after_repair",
      I1_for_solve = I1_for_solve,
      I1_scaled_for_solve = I1_scaled_for_solve,
      information_scale = information_scale,
      solve_D = NULL,
      solve_C12 = NULL,
      diagnostics = list(
        method = if (repair_needed) repair_method else "direct_target_solve",
        repaired = repair_needed,
        min_eigen_raw = min(eig_raw$values),
        min_eigen_scaled_raw = min(eig_scaled_raw$values),
        min_eigen_repaired = min(eig_repaired),
        min_eigen_scaled_repaired = min(eig_scaled_repaired),
        near_pd_eig_tol = near_pd_eig_tol,
        near_pd_conv_tol = near_pd_conv_tol,
        near_pd_posd_tol = near_pd_posd_tol
      )
    ))
  }
  repaired_max_eval <- max(abs(eig_repaired))
  repaired_condition_number <- if (min(eig_repaired) <= 0){
    Inf
  } else {
    repaired_max_eval / min(eig_repaired)
  }
  repaired_condition_number_scaled <-
    max(eig_scaled_repaired) / min(eig_scaled_repaired)
  repair_norm_f <- sqrt(sum((I1_for_solve - I1)^2))
  raw_norm_f <- sqrt(sum(I1^2))

  diagnostics <- list(
    method = if (repair_needed) repair_method else "direct_target_solve",
    repaired = repair_needed,
    repair_method = repair_method,
    used_target_restriction = FALSE,
    full_information_positive = all(eig_raw$values > 0),
    condition_number = raw_condition_number,
    condition_number_repaired = repaired_condition_number,
    condition_number_scaled_repaired = repaired_condition_number_scaled,
    near_pd_eig_tol = near_pd_eig_tol,
    near_pd_conv_tol = near_pd_conv_tol,
    near_pd_posd_tol = near_pd_posd_tol,
    eigen_cutoff = eigen_cutoff,
    raw_retained_rank = sum(retained_raw),
    scaled_retained_rank = sum(retained_scaled),
    retained_rank = sum(retained_scaled),
    discarded_rank = sum(discarded_scaled),
    target_norm = target_norm,
    discarded_target_norm = discarded_target_norm,
    target_projection_ratio = target_projection_ratio,
    min_eigen_raw = min(eig_raw$values),
    max_eigen_raw = max(eig_raw$values),
    rel_min_raw = min(eig_raw$values) / max_eval,
    use_diagonal_scaling = use_diagonal_scaling,
    scaled_eigen_cutoff = scaled_eigen_cutoff,
    min_eigen_scaled_raw = min(eig_scaled_raw$values),
    max_eigen_scaled_raw = max(eig_scaled_raw$values),
    min_eigen_scaled_repaired = min(eig_scaled_repaired),
    max_eigen_scaled_repaired = max(eig_scaled_repaired),
    min_eigen_repaired = min(eig_repaired),
    max_eigen_repaired = max(eig_repaired),
    rel_min_repaired = min(eig_repaired) / repaired_max_eval,
    near_pd_converged = near_pd_converged,
    near_pd_iterations = near_pd_iterations,
    near_pd_norm_f_scaled = near_pd_norm_f_scaled,
    repair_norm_f = repair_norm_f,
    relative_repair_norm_f = if (raw_norm_f == 0){
      NA_real_
    } else {
      repair_norm_f / raw_norm_f
    }
  )

  target_rhs <- cbind(t(D),C12)
  target_rhs_scaled <- sweep(target_rhs,1,information_scale,"/")
  scaled_solution <- safe_solve(I1_scaled_for_solve,target_rhs_scaled)
  direct_solution <- if (is.null(scaled_solution)){
    NULL
  } else {
    sweep(scaled_solution,1,information_scale,"/")
  }
  if (!is.null(direct_solution) && all(is.finite(direct_solution))){
    q <- nrow(D)
    diagnostics$relative_solve_residual_scaled <-
      max(abs(I1_scaled_for_solve %*% scaled_solution -
                target_rhs_scaled)) /
      max(1,max(abs(target_rhs_scaled)))
    diagnostics$relative_solve_residual_original_scale <-
      max(abs(I1_for_solve %*% direct_solution - target_rhs)) /
      max(1,max(abs(target_rhs)))
    diagnostics$relative_solve_residual <-
      diagnostics$relative_solve_residual_scaled
    diagnostics$solve_residual_tolerance <- 1e-6
    if (!is.finite(diagnostics$relative_solve_residual) ||
        diagnostics$relative_solve_residual >
          diagnostics$solve_residual_tolerance){
      return(list(
        failed = TRUE,
        failure_reason = "longitudinal_information_solve_residual_too_large",
        I1_for_solve = I1_for_solve,
        solve_D = NULL,
        solve_C12 = NULL,
        diagnostics = diagnostics
      ))
    }
    return(list(
      failed = FALSE,
      failure_reason = NA_character_,
      I1_for_solve = I1_for_solve,
      I1_scaled_for_solve = I1_scaled_for_solve,
      information_scale = information_scale,
      solve_D = direct_solution[,seq_len(q),drop = FALSE],
      solve_C12 = direct_solution[,q + seq_len(q),drop = FALSE],
      diagnostics = diagnostics
    ))
  }

  list(
    failed = TRUE,
    failure_reason = "repaired_longitudinal_information_solve_failed",
    I1_for_solve = I1_for_solve,
    I1_scaled_for_solve = I1_scaled_for_solve,
    information_scale = information_scale,
    solve_D = NULL,
    solve_C12 = NULL,
    diagnostics = diagnostics
  )
}

repair_target_covariance_near_pd <- function(V){
  V_raw <- 0.5 * (as.matrix(V) + t(as.matrix(V)))

  if (any(!is.finite(V_raw))){
    stop("Target covariance contains nonfinite values")
  }

  eig_raw <- eigen(V_raw,symmetric = TRUE,only.values = TRUE)$values
  max_eval <- max(abs(eig_raw))
  min_eval <- min(eig_raw)
  rel_min <- if (max_eval == 0) NA_real_ else min_eval / max_eval
  raw_norm_f <- sqrt(sum(V_raw^2))
  near_pd_formals <- formals(Matrix::nearPD)
  near_pd_eig_tol <- as.numeric(near_pd_formals$eig.tol)
  near_pd_conv_tol <- as.numeric(near_pd_formals$conv.tol)
  near_pd_posd_tol <- as.numeric(near_pd_formals$posd.tol)

  if (all(eig_raw > 0)){
    return(list(
      vcov = V_raw,
      raw_vcov = V_raw,
      repaired = FALSE,
      failed = FALSE,
      repair_method = "none",
      near_pd_norm_f = 0,
      relative_repair_norm_f = 0,
      near_pd_converged = NA,
      near_pd_iterations = NA_integer_,
      near_pd_eig_tol = near_pd_eig_tol,
      near_pd_conv_tol = near_pd_conv_tol,
      near_pd_posd_tol = near_pd_posd_tol,
      min_eigen_raw = min_eval,
      max_eigen_raw = max_eval,
      rel_min_raw = rel_min,
      min_eigen_repaired = min_eval,
      max_eigen_repaired = max_eval,
      rel_min_repaired = rel_min
    ))
  }

  near_pd <- tryCatch(
    Matrix::nearPD(V_raw,corr = FALSE),
    error = function(e) NULL
  )

  if (is.null(near_pd) || !isTRUE(near_pd$converged)){
    return(list(
      vcov = NULL,
      raw_vcov = V_raw,
      repaired = FALSE,
      failed = TRUE,
      repair_method = "nearPD_failed",
      near_pd_norm_f = NA_real_,
      relative_repair_norm_f = NA_real_,
      near_pd_converged = if (is.null(near_pd)) NA else near_pd$converged,
      near_pd_iterations = if (is.null(near_pd)) {
        NA_integer_
      } else {
        near_pd$iterations
      },
      near_pd_eig_tol = near_pd_eig_tol,
      near_pd_conv_tol = near_pd_conv_tol,
      near_pd_posd_tol = near_pd_posd_tol,
      min_eigen_raw = min_eval,
      max_eigen_raw = max_eval,
      rel_min_raw = rel_min,
      min_eigen_repaired = NA_real_,
      max_eigen_repaired = NA_real_,
      rel_min_repaired = NA_real_
    ))
  }

  V_repaired <- as.matrix(near_pd$mat)
  V_repaired <- 0.5 * (V_repaired + t(V_repaired))
  eig_repaired <- eigen(V_repaired,symmetric = TRUE,
                        only.values = TRUE)$values
  if (any(eig_repaired <= 0)){
    return(list(
      vcov = NULL,
      raw_vcov = V_raw,
      repaired = FALSE,
      failed = TRUE,
      repair_method = "nearPD_not_positive_definite",
      near_pd_norm_f = near_pd$normF,
      relative_repair_norm_f = if (raw_norm_f == 0) {
        NA_real_
      } else {
        near_pd$normF / raw_norm_f
      },
      near_pd_converged = near_pd$converged,
      near_pd_iterations = near_pd$iterations,
      near_pd_eig_tol = near_pd_eig_tol,
      near_pd_conv_tol = near_pd_conv_tol,
      near_pd_posd_tol = near_pd_posd_tol,
      min_eigen_raw = min_eval,
      max_eigen_raw = max_eval,
      rel_min_raw = rel_min,
      min_eigen_repaired = min(eig_repaired),
      max_eigen_repaired = max(eig_repaired),
      rel_min_repaired = min(eig_repaired) / max(abs(eig_repaired))
    ))
  }

  list(
    vcov = V_repaired,
    raw_vcov = V_raw,
    repaired = TRUE,
    failed = FALSE,
    repair_method = "nearPD",
    near_pd_norm_f = near_pd$normF,
    relative_repair_norm_f = if (raw_norm_f == 0) {
      NA_real_
    } else {
      near_pd$normF / raw_norm_f
    },
    near_pd_converged = near_pd$converged,
    near_pd_iterations = near_pd$iterations,
    near_pd_eig_tol = near_pd_eig_tol,
    near_pd_conv_tol = near_pd_conv_tol,
    near_pd_posd_tol = near_pd_posd_tol,
    min_eigen_raw = min_eval,
    max_eigen_raw = max_eval,
    rel_min_raw = rel_min,
    min_eigen_repaired = min(eig_repaired),
    max_eigen_repaired = max(eig_repaired),
    rel_min_repaired = min(eig_repaired) / max(abs(eig_repaired))
  )
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

check_dimensions <- function(I1,D,I1_solve_D,I1_solve_C12,
                             S1,S2,C12,V_MT){
  p <- nrow(I1)
  q <- nrow(D)
  n <- nrow(S1)
  data.frame(
    quantity = c("I1","D","I1_solve_D","I1_solve_C12",
                 "S1","S2","C12","V_MT"),
    observed = c(paste(dim(I1),collapse = "x"),
                 paste(dim(D),collapse = "x"),
                 paste(dim(I1_solve_D),collapse = "x"),
                 paste(dim(I1_solve_C12),collapse = "x"),
                 paste(dim(S1),collapse = "x"),
                 paste(dim(S2),collapse = "x"),
                 paste(dim(C12),collapse = "x"),
                 paste(dim(V_MT),collapse = "x")),
    expected = c(paste(c(p,p),collapse = "x"),
                 paste(c(q,p),collapse = "x"),
                 paste(c(p,q),collapse = "x"),
                 paste(c(p,q),collapse = "x"),
                 paste(c(n,p),collapse = "x"),
                 paste(c(n,q),collapse = "x"),
                 paste(c(p,q),collapse = "x"),
                 paste(c(q,q),collapse = "x")),
    passed = c(all(dim(I1) == c(p,p)),
               all(dim(D) == c(q,p)),
               all(dim(I1_solve_D) == c(p,q)),
               all(dim(I1_solve_C12) == c(p,q)),
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

V_naive_raw <- as.matrix(cox_fit_base$var)
V_naive_raw <- 0.5 * (V_naive_raw + t(V_naive_raw))
V_naive_repair <- repair_target_covariance_near_pd(V_naive_raw)
if (isTRUE(V_naive_repair$failed)){
  stop("Could not repair naive Cox covariance")
}
V_naive <- V_naive_repair$vcov
I2 <- safe_solve(V_naive)
if (is.null(I2) || any(!is.finite(I2))){
  stop("Could not invert naive Cox covariance")
}
rownames(V_naive) <- names(gamma_hat)
colnames(V_naive) <- names(gamma_hat)
rownames(V_naive_raw) <- names(gamma_hat)
colnames(V_naive_raw) <- names(gamma_hat)
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

I1_target_solve <- solve_information_for_targets(
  I1 = I1_raw,
  D = D,
  C12 = C12
)
if (isTRUE(I1_target_solve$failed)){
  stop(paste("Longitudinal information repair/solve failed:",
             I1_target_solve$failure_reason))
}
print(paste("Longitudinal information method:",
            I1_target_solve$diagnostics$method,
            "scaled retained rank:",I1_target_solve$diagnostics$retained_rank,
            "of",nrow(I1_raw)))
if (isTRUE(I1_target_solve$diagnostics$repaired)){
  print(paste("Longitudinal information nearPD relative Frobenius change:",
              signif(I1_target_solve$diagnostics$relative_repair_norm_f,6)))
}
I1_for_solve <- I1_target_solve$I1_for_solve
I1_scaled_for_solve <- I1_target_solve$I1_scaled_for_solve
I1_information_scale <- I1_target_solve$information_scale
I1_solve_D <- I1_target_solve$solve_D
I1_solve_C12 <- I1_target_solve$solve_C12

C_if <- I1_solve_C12 %*% V_naive
rownames(C_if) <- colnames(S1)
colnames(C_if) <- names(gamma_hat)

delta_increment_raw <- D %*% I1_solve_D
delta_increment_raw <- 0.5 * (delta_increment_raw + t(delta_increment_raw))
rownames(delta_increment_raw) <- names(gamma_hat)
colnames(delta_increment_raw) <- names(gamma_hat)

delta_increment_eigenvalues <- matrix_eigenvalues(delta_increment_raw)
delta_increment_min_eigen <- min(delta_increment_eigenvalues)
delta_increment_max_abs_eigen <- max(abs(delta_increment_eigenvalues))
delta_increment_rel_min <- if (delta_increment_max_abs_eigen == 0){
  NA_real_
} else {
  delta_increment_min_eigen / delta_increment_max_abs_eigen
}
delta_increment_negative_tol <-
  1e-10 * max(1,delta_increment_max_abs_eigen)

# The stage-one propagated term is theoretically positive semidefinite.  If
# finite-precision arithmetic or an I1 repair makes it indefinite, retain the
# raw result for diagnostics and use its nearest positive-definite projection
# in the downstream delta-method and Murphy-Topel covariances.
delta_increment_repair <-
  repair_target_covariance_near_pd(delta_increment_raw)
if (isTRUE(delta_increment_repair$failed)){
  stop("Could not repair stage-one propagated uncertainty using nearPD")
}

delta_increment <- delta_increment_repair$vcov
rownames(delta_increment) <- names(gamma_hat)
colnames(delta_increment) <- names(gamma_hat)

if (isTRUE(delta_increment_repair$repaired)){
  warning(paste(
    "Stage-one propagated uncertainty repaired using nearPD;",
    "raw relative minimum eigenvalue:",
    signif(delta_increment_rel_min,6),
    "relative Frobenius change:",
    signif(delta_increment_repair$relative_repair_norm_f,6)
  ))
}

cross_increment <- D %*% C_if + t(C_if) %*% t(D)
rownames(cross_increment) <- names(gamma_hat)
colnames(cross_increment) <- names(gamma_hat)

V_delta_raw <- V_naive + delta_increment
V_delta_raw <- 0.5 * (V_delta_raw + t(V_delta_raw))
rownames(V_delta_raw) <- names(gamma_hat)
colnames(V_delta_raw) <- names(gamma_hat)

V_MT_raw <- V_delta_raw + cross_increment
V_MT_raw <- 0.5 * (V_MT_raw + t(V_MT_raw))
rownames(V_MT_raw) <- names(gamma_hat)
colnames(V_MT_raw) <- names(gamma_hat)

V_delta_repair <- repair_target_covariance_near_pd(V_delta_raw)
if (isTRUE(V_delta_repair$failed)){
  stop("Could not repair delta-method target covariance")
}
V_delta <- V_delta_repair$vcov
rownames(V_delta) <- names(gamma_hat)
colnames(V_delta) <- names(gamma_hat)

V_MT_repair <- repair_target_covariance_near_pd(V_MT_raw)
if (isTRUE(V_MT_repair$failed)){
  stop("Could not repair Murphy-Topel target covariance")
}
V_MT <- V_MT_repair$vcov
rownames(V_MT) <- names(gamma_hat)
colnames(V_MT) <- names(gamma_hat)

se_naive <- sqrt(diag(V_naive))
se_delta <- sqrt(diag(V_delta))
se_MT <- sqrt(diag(V_MT))
if (any(!is.finite(c(se_naive,se_delta,se_MT)))){
  stop("Two-stage standard errors contain nonfinite values")
}

gamma_truth <- if (settings$data_source == DATA_SOURCE[["simulation"]]){
  build_simulation_gamma_truth(names(gamma_hat),true_params,
                               settings$fit_mix_num)
} else {
  stats::setNames(rep(NA_real_,length(gamma_hat)),names(gamma_hat))
}

zero_C_if <- matrix(0,nrow = nrow(C_if),ncol = ncol(C_if))
dimnames(zero_C_if) <- dimnames(C_if)
V_MT_zero_C12_raw <- V_naive + delta_increment +
  D %*% zero_C_if + t(zero_C_if) %*% t(D)
V_MT_zero_C12_raw <- 0.5 *
  (V_MT_zero_C12_raw + t(V_MT_zero_C12_raw))

zero_D <- matrix(0,nrow = nrow(D),ncol = ncol(D))
dimnames(zero_D) <- dimnames(D)
V_delta_zero_D <- V_naive + zero_D %*% I1_solve_D
V_delta_zero_D <- 0.5 * (V_delta_zero_D + t(V_delta_zero_D))
V_MT_zero_D <- V_delta_zero_D +
  zero_D %*% C_if + t(C_if) %*% t(zero_D)
V_MT_zero_D <- 0.5 * (V_MT_zero_D + t(V_MT_zero_D))

dimension_checks <- check_dimensions(I1_raw,D,I1_solve_D,I1_solve_C12,
                                     S1,S2,C12,V_MT)
if (!all(dimension_checks$passed)){
  stop("Murphy-Topel dimension checks failed")
}
symmetry_checks <- data.frame(
  matrix = c("I1","V_naive","V_delta_raw","V_delta",
             "V_MT_raw","V_MT"),
  max_abs_asymmetry = c(matrix_max_asymmetry(I1_raw),
                        matrix_max_asymmetry(V_naive),
                        matrix_max_asymmetry(V_delta_raw),
                        matrix_max_asymmetry(V_delta),
                        matrix_max_asymmetry(V_MT_raw),
                        matrix_max_asymmetry(V_MT)),
  stringsAsFactors = FALSE
)
score_checks <- rbind(score_sum_diagnostics(S1,"S1"),
                      score_sum_diagnostics(S2,"S2"))
stage1_score_base <- as.numeric(mt_perturb$base_score)
names(stage1_score_base) <- colnames(S1)
stage1_score_sum_numeric <- colSums(S1)
stage1_score_column_norm <- sqrt(colSums(S1^2))
stage1_score_relative_imbalance <-
  abs(stage1_score_sum_numeric) /
  pmax(stage1_score_column_norm,.Machine$double.eps)
stage1_score_rhs_scaled <- stage1_score_base / I1_information_scale
stage1_newton_step_scaled <- safe_solve(I1_scaled_for_solve,
                                        stage1_score_rhs_scaled)
if (is.null(stage1_newton_step_scaled) ||
    any(!is.finite(stage1_newton_step_scaled))){
  stop("Could not calculate stage-one stationarity diagnostics")
}
stage1_newton_step <- stage1_newton_step_scaled / I1_information_scale
names(stage1_newton_step) <- names(stage1_score_base)
stage1_newton_decrement_squared <-
  sum(stage1_score_base * stage1_newton_step)
stage1_stationarity <- list(
  max_abs_analytic_score = max(abs(stage1_score_base)),
  max_abs_numeric_score_sum = max(abs(stage1_score_sum_numeric)),
  max_abs_analytic_numeric_score_diff =
    max(abs(stage1_score_base - stage1_score_sum_numeric)),
  max_relative_analytic_numeric_score_diff =
    max(abs(stage1_score_base - stage1_score_sum_numeric)) /
    max(1,max(abs(stage1_score_base)),
        max(abs(stage1_score_sum_numeric))),
  max_relative_score_imbalance =
    max(stage1_score_relative_imbalance),
  max_abs_repaired_information_newton_step =
    max(abs(stage1_newton_step)),
  repaired_information_newton_decrement =
    sqrt(max(0,stage1_newton_decrement_squared)),
  analytic_score = stage1_score_base,
  numeric_score_sum = stage1_score_sum_numeric,
  relative_score_imbalance = stage1_score_relative_imbalance,
  repaired_information_newton_step = stage1_newton_step
)
print(paste("Stage-one maximum relative score imbalance:",
            signif(stage1_stationarity$max_relative_score_imbalance,6)))
formula_checks <- data.frame(
  check = c("zero_C12_V_MT_equals_V_delta",
            "zero_D_V_delta_equals_V_naive",
            "zero_D_V_MT_equals_V_naive"),
  max_abs_diff = c(max(abs(V_MT_zero_C12_raw - V_delta_raw),na.rm = TRUE),
                   max(abs(V_delta_zero_D - V_naive),na.rm = TRUE),
                   max(abs(V_MT_zero_D - V_naive),na.rm = TRUE)),
  stringsAsFactors = FALSE
)

to_save$murphy_topel <- list(
  variance_method = "weighted_target_murphy_topel",
  settings = list(h2_eps = eps,
                  cox_ties = "breslow",
                  cox_score_residuals_weighted = TRUE,
                  variance_method = "weighted_target_murphy_topel",
                  longitudinal_information_method =
                    I1_target_solve$diagnostics$method,
                  longitudinal_information_repaired =
                    I1_target_solve$diagnostics$repaired,
                  near_pd_eig_tol =
                    I1_target_solve$diagnostics$near_pd_eig_tol,
                  near_pd_conv_tol =
                    I1_target_solve$diagnostics$near_pd_conv_tol,
                  near_pd_posd_tol =
                    I1_target_solve$diagnostics$near_pd_posd_tol,
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
  I1_for_solve = I1_for_solve,
  I1_for_vcov = I1_for_solve,
  I1_scaled_for_solve = I1_scaled_for_solve,
  I1_information_scale = I1_information_scale,
  I1_solve_D = I1_solve_D,
  I1_solve_C12 = I1_solve_C12,
  I2 = I2,
  V_naive_raw = V_naive_raw,
  V_naive = V_naive,
  D = D,
  S1 = if (retain_score_matrices) S1 else NULL,
  S2 = if (retain_score_matrices) S2 else NULL,
  C12 = C12,
  C_if = C_if,
  delta_increment = delta_increment,
  delta_increment_raw = delta_increment_raw,
  cross_increment = cross_increment,
  V_delta_raw = V_delta_raw,
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
    I1_target_solve = I1_target_solve$diagnostics,
    I1_repair = I1_target_solve$diagnostics,
    delta_increment = list(
      min_eigen_raw = delta_increment_min_eigen,
      max_abs_eigen_raw = delta_increment_max_abs_eigen,
      rel_min_raw = delta_increment_rel_min,
      negative_eigen_warning_tolerance = delta_increment_negative_tol,
      eigenvalues = delta_increment_eigenvalues
    ),
    delta_increment_repair = delta_increment_repair,
    V_naive_repair = V_naive_repair,
    V_delta_repair = V_delta_repair,
    V_MT_repair = V_MT_repair,
    dimension_checks = dimension_checks,
    symmetry_checks = symmetry_checks,
    score_checks = score_checks,
    stage1_stationarity = stage1_stationarity,
    formula_checks = formula_checks,
    perturbation = mt_perturb$diagnostics,
    I1_eigenvalues = matrix_eigenvalues(I1_raw),
    I1_for_solve_eigenvalues = matrix_eigenvalues(I1_for_solve),
    I1_scaled_for_solve_eigenvalues =
      matrix_eigenvalues(I1_scaled_for_solve),
    V_delta_raw_eigenvalues = matrix_eigenvalues(V_delta_raw),
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
