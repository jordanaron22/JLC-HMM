strict_logical_setting <- function(value, setting_name){
  if (is.logical(value)){
    return(value)
  }
  normalized <- tolower(as.character(value))
  if (normalized == "true"){
    return(TRUE)
  }
  if (normalized == "false"){
    return(FALSE)
  }
  stop(paste(setting_name,"must be exactly true or false; got:",value))
}

validate_count_setting <- function(value, setting_name, min_value = 1){
  if (is.na(value) || value != floor(value) || value < min_value){
    stop(paste(setting_name,"must be an integer >=",min_value,"; got:",value))
  }
}

normalize_target_weekday <- function(target_weekday, weekday_codes){
  if (is.character(target_weekday) && target_weekday %in% names(weekday_codes)){
    target_weekday <- unname(weekday_codes[[target_weekday]])
  }
  target_weekday <- suppressWarnings(as.numeric(target_weekday))
  if (is.na(target_weekday) || !target_weekday %in% unique(weekday_codes)){
    stop(paste("target_weekday must be 0-7 or one of:",
               paste(names(weekday_codes),collapse = ", ")))
  }
  target_weekday
}

validate_settings <- function(settings, model_type_codes, data_source_values,
                              weekday_codes, emission_overlap_factor){
  settings$fit_mix_num <- as.numeric(settings$fit_mix_num)
  settings$true_mix_num <- as.numeric(settings$true_mix_num)
  validate_count_setting(settings$fit_mix_num,"fit_mix_num")
  if (!is.na(settings$true_mix_num)){
    validate_count_setting(settings$true_mix_num,"true_mix_num")
  }

  settings$model_type <- as.character(settings$model_type)
  if (!settings$model_type %in% names(model_type_codes)){
    stop(paste("model_type must be one of:",paste(names(model_type_codes),collapse = ", ")))
  }
  settings$data_source <- as.character(settings$data_source)
  if (!settings$data_source %in% data_source_values){
    stop(paste("data_source must be one of:",paste(data_source_values,collapse = ", ")))
  }

  logical_setting_names <- c("run_bootstrap","run_leave_one_out_cv",
                             "use_hot_start","use_seed","include_activity",
                             "include_light","check_transition_update",
                             "weekend_only","save_reduced_output",
                             "survival_only","include_tlogmims",
                             "class_selection_run")
  for (setting_name in logical_setting_names){
    settings[[setting_name]] <- strict_logical_setting(settings[[setting_name]],setting_name)
  }

  settings$target_weekday <- normalize_target_weekday(settings$target_weekday, weekday_codes)
  settings$periods_per_day <- as.numeric(settings$periods_per_day)
  settings$num_day_type_levels <- as.numeric(settings$num_day_type_levels)
  settings$init_jitter_scale <- as.numeric(settings$init_jitter_scale)
  settings$simulation_days <- as.numeric(settings$simulation_days)
  settings$num_people <- as.numeric(settings$num_people)
  settings$missing_perc <- as.numeric(settings$missing_perc)

  validate_count_setting(settings$periods_per_day,"periods_per_day")
  validate_count_setting(settings$num_day_type_levels,"num_day_type_levels")
  validate_count_setting(settings$simulation_days,"simulation_days")
  validate_count_setting(settings$num_people,"num_people")
  if (is.na(settings$init_jitter_scale) || settings$init_jitter_scale < 0){
    stop(paste("init_jitter_scale must be nonnegative; got:",settings$init_jitter_scale))
  }
  if (is.na(settings$missing_perc) ||
      settings$missing_perc < 0 || settings$missing_perc > 1){
    stop(paste("missing_perc must be between 0 and 1; got:",
               settings$missing_perc))
  }
  settings$emission_overlap <- as.character(settings$emission_overlap)
  if (!settings$emission_overlap %in% names(emission_overlap_factor)){
    stop(paste("emission_overlap must be one of:",
               paste(names(emission_overlap_factor),collapse = ", ")))
  }

  settings
}

validate_saved_results <- function(to_save, required_sections = c("est_params"),
                                   source_name = "saved results"){
  if (is.null(to_save)){
    stop(paste("Could not load",source_name))
  }
  for (section_name in required_sections){
    get_saved_section(to_save,section_name,required = TRUE)
  }
  invisible(TRUE)
}

validate_dimensions <- function(value, expected_dim, value_name){
  actual_dim <- dim(value)
  if (is.null(actual_dim) || length(actual_dim) != length(expected_dim) ||
      any(actual_dim != expected_dim)){
    stop(paste(value_name,"must have dimensions",
               paste(expected_dim,collapse = " x "),
               "but has",
               paste(actual_dim,collapse = " x ")))
  }
}

validate_param_list <- function(params, expected_mix_num, expected_vcovar_num,
                                param_label = "params"){
  required_names <- c("init","params_tran_array","emit_act","emit_light",
                      "corr_mat","nu_mat","beta_vec","lambda_act_mat",
                      "lambda_light_mat")
  missing_names <- setdiff(required_names,names(params))
  if (length(missing_names) > 0){
    stop(paste(param_label,"is missing required fields:",
               paste(missing_names,collapse = ", ")))
  }

  validate_dimensions(params$init,c(expected_mix_num,NUM_MARKOV_STATES),
                      paste0(param_label,"$init"))
  validate_dimensions(params$params_tran_array,
                      c(expected_mix_num,6,expected_vcovar_num),
                      paste0(param_label,"$params_tran_array"))
  validate_dimensions(params$emit_act,
                      c(NUM_MARKOV_STATES,2,expected_mix_num,expected_vcovar_num),
                      paste0(param_label,"$emit_act"))
  validate_dimensions(params$emit_light,
                      c(NUM_MARKOV_STATES,2,expected_mix_num,expected_vcovar_num),
                      paste0(param_label,"$emit_light"))
  validate_dimensions(params$corr_mat,
                      c(expected_mix_num,NUM_MARKOV_STATES,expected_vcovar_num),
                      paste0(param_label,"$corr_mat"))
  if (ncol(params$nu_mat) != expected_mix_num){
    stop(paste(param_label,"$nu_mat must have",expected_mix_num,"columns"))
  }
  if (length(params$beta_vec) != expected_mix_num){
    stop(paste(param_label,"$beta_vec must have length",expected_mix_num))
  }
  validate_dimensions(params$lambda_act_mat,
                      c(expected_mix_num,NUM_MARKOV_STATES,expected_vcovar_num),
                      paste0(param_label,"$lambda_act_mat"))
  validate_dimensions(params$lambda_light_mat,
                      c(expected_mix_num,NUM_MARKOV_STATES,expected_vcovar_num),
                      paste0(param_label,"$lambda_light_mat"))
  invisible(TRUE)
}

validate_hmm_data <- function(act, light, vcovar_mat){
  if (is.null(dim(act)) || is.null(dim(light))){
    stop("act and light must both be matrices")
  }
  if (any(dim(act) != dim(light))){
    stop("act and light must have identical dimensions")
  }
  if (length(vcovar_mat) != length(act)){
    stop("vcovar_mat must have one entry for each act/light observation")
  }
  invisible(TRUE)
}

validate_survival_inputs <- function(surv_time, surv_event, surv_covar,
                                     expected_n = length(surv_time)){
  if (length(surv_time) != expected_n || length(surv_event) != expected_n){
    stop("surv_time and surv_event must match the number of individuals")
  }
  if (any(is.na(surv_time)) || any(surv_time < 0)){
    stop("surv_time must be non-missing and nonnegative")
  }
  if (any(is.na(surv_event)) || any(!surv_event %in% c(0,1))){
    stop("surv_event must contain only 0/1 indicators")
  }
  for (covar_ind in seq_along(surv_covar)){
    covar_n <- if (is.null(dim(surv_covar[[covar_ind]]))){
      length(surv_covar[[covar_ind]])
    } else {
      nrow(surv_covar[[covar_ind]])
    }
    if (covar_n != expected_n){
      stop(paste("surv_covar[[",covar_ind,"]] has",covar_n,
                 "rows but expected",expected_n))
    }
  }
  invisible(TRUE)
}

validate_re_prob <- function(re_prob, expected_n, expected_mix_num){
  if (is.null(dim(re_prob))){
    stop("re_prob must be a matrix")
  }
  validate_dimensions(re_prob,c(expected_n,expected_mix_num),"re_prob")
  invisible(TRUE)
}

make_survival_context <- function(surv_time, surv_event, surv_covar, re_prob,
                                  expected_mix_num){
  validate_survival_inputs(surv_time,surv_event,surv_covar)
  validate_re_prob(re_prob,length(surv_time),expected_mix_num)
  list(surv_time = surv_time,
       surv_event = surv_event,
       surv_covar = surv_covar,
       re_prob = re_prob)
}

update_survival_context_re_prob <- function(survival_context, re_prob,
                                            expected_mix_num){
  validate_re_prob(re_prob,length(survival_context$surv_time),expected_mix_num)
  survival_context$re_prob <- re_prob
  survival_context
}
