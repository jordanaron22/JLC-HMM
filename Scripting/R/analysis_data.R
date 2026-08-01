FindNHANESDataDir <- function(
    candidate_paths = c("Data",file.path("..","Data"),
                        file.path("..","Rcode","Data"))){
  required_file <- "NHANES_2011_2012_2013_2014.rda"
  for (candidate_path in candidate_paths){
    if (dir.exists(candidate_path) &&
        file.exists(file.path(candidate_path,required_file))){
      return(normalizePath(candidate_path,winslash = "/",mustWork = TRUE))
    }
  }

  stop(paste("Could not find the NHANES data directory; searched:",
             paste(candidate_paths,collapse = ", ")))
}

SummarizeSEWeights <- function(sweights_vec){
  list(
    n = length(sweights_vec),
    sum = sum(sweights_vec),
    mean = mean(sweights_vec),
    min = min(sweights_vec),
    max = max(sweights_vec)
  )
}

ValidateSEAnalysisData <- function(analysis_data,fit_mix_num = NULL){
  required_fields <- c("act","light","vcovar_mat","nu_covar_mat",
                       "surv_time","surv_event","surv_covar",
                       "combined_covar_mat","sweights_vec",
                       "lod_act","lod_light","data_source")
  missing_fields <- setdiff(required_fields,names(analysis_data))
  if (length(missing_fields) > 0){
    stop("SE analysis data are missing fields: ",
         paste(missing_fields,collapse = ", "))
  }

  act <- as.matrix(analysis_data$act)
  light <- as.matrix(analysis_data$light)
  n <- ncol(act)
  if (!all(dim(light) == dim(act))){
    stop("SE analysis activity and light matrices must have equal dimensions")
  }
  if (!all(dim(analysis_data$vcovar_mat) == dim(act))){
    stop("SE analysis vcovar_mat dimensions must match activity")
  }
  if (nrow(as.matrix(analysis_data$nu_covar_mat)) != n){
    stop("SE analysis nu_covar_mat must have one row per participant")
  }
  if (length(analysis_data$surv_time) != n ||
      length(analysis_data$surv_event) != n){
    stop("SE analysis survival vectors must have one value per participant")
  }

  sweights_vec <- as.numeric(analysis_data$sweights_vec)
  if (length(sweights_vec) != n){
    stop("SE analysis weights must contain one value per participant")
  }
  if (any(!is.finite(sweights_vec))){
    stop("SE analysis weights contain nonfinite values")
  }
  if (any(sweights_vec <= 0)){
    stop("SE analysis weights must be strictly positive")
  }

  if (!is.null(fit_mix_num) && !is.null(analysis_data$re_prob)){
    validate_re_prob(analysis_data$re_prob,n,fit_mix_num)
  }

  invisible(analysis_data)
}

MakeSimulationSEAnalysisData <- function(simulated_hmm){
  act <- simulated_hmm$act
  surv_covar_sim <- simulated_hmm$surv_covar_sim
  combined_covar_mat <- matrix(surv_covar_sim - 1,nrow = ncol(act))
  combined_covar_mat <- as.factor(combined_covar_mat)

  analysis_data <- list(
    act = act,
    light = simulated_hmm$light,
    vcovar_mat = simulated_hmm$vcovar_mat,
    nu_covar_mat = simulated_hmm$nu_covar_mat,
    surv_time = simulated_hmm$survival$time,
    surv_event = simulated_hmm$survival$event,
    surv_covar = list(simulated_hmm$age_vec,Vec2Mat(surv_covar_sim)),
    combined_covar_mat = combined_covar_mat,
    sweights_vec = rep(1,ncol(act)),
    lod_act = -5.809153,
    lod_light = -1.560658,
    subject_id = seq_len(ncol(act)),
    data_source = DATA_SOURCE[["simulation"]]
  )
  ValidateSEAnalysisData(analysis_data)
  analysis_data
}

PrepareNHANESSEAnalysisData <- function(settings,est_params,data_dir = NULL){
  if (is.null(data_dir)){
    data_dir <- FindNHANESDataDir()
  }

  bootstrap <- isTRUE(settings$run_bootstrap)
  leave_out <- isTRUE(settings$run_leave_one_out_cv)
  if (is.null(settings$run_bootstrap)){
    bootstrap <- isTRUE(as.logical(settings$bootstrap))
  }
  if (is.null(settings$run_leave_one_out_cv)){
    leave_out <- isTRUE(as.logical(settings$leave_out))
  }

  # The original NHANES fit seeds before data preparation. Reproduce that
  # ordering so bootstrap samples are reconstructed exactly.
  if (bootstrap && isTRUE(settings$use_seed)){
    set.seed(settings$sim_num)
  }

  surv_coef <- get_saved_param(est_params,"surv_coef",required = TRUE)
  nhanes_data <- prepare_nhanes_data(
    period_len = settings$period_len,
    bootstrap = bootstrap,
    leave_out = leave_out,
    sim_num = settings$sim_num,
    single_day = settings$single_day,
    weekend_only = settings$weekend_only,
    load_data = TRUE,
    surv_coef_true = surv_coef,
    data_dir = data_dir
  )

  analysis_data <- list(
    act = nhanes_data$act,
    light = nhanes_data$light,
    vcovar_mat = nhanes_data$vcovar_mat,
    nu_covar_mat = nhanes_data$nu_covar_mat,
    surv_time = nhanes_data$surv_time,
    surv_event = nhanes_data$surv_event,
    surv_covar = nhanes_data$surv_covar,
    combined_covar_mat = nhanes_data$combined_covar_mat,
    sweights_vec = nhanes_data$sweights_vec,
    lod_act = nhanes_data$lod_act,
    lod_light = nhanes_data$lod_light,
    subject_id = as.character(nhanes_data$id$SEQN),
    data_source = DATA_SOURCE[["nhanes"]],
    data_dir = data_dir
  )
  ValidateSEAnalysisData(analysis_data)
  if (abs(mean(analysis_data$sweights_vec) - 1) > 1e-10){
    stop("Reconstructed NHANES weights are not normalized to mean one")
  }
  analysis_data
}
