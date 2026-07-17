get_named_or_positional <- function(x, name, slot, default = NULL, required = FALSE){
  if (!is.null(x) && !is.null(names(x)) && name %in% names(x)){
    return(x[[name]])
  }
  if (!is.null(x) && length(x) >= slot){
    return(x[[slot]])
  }
  if (required){
    stop(paste("Missing saved field:",name))
  }
  default
}

get_saved_section <- function(to_save, section_name, required = TRUE){
  get_named_or_positional(to_save, section_name, SAVED_SECTION_SLOT[[section_name]],
                          required = required)
}

get_saved_param <- function(saved_params, param_name, default = NULL, required = TRUE){
  get_named_or_positional(saved_params, param_name, PARAM_SLOT[[param_name]],
                          default = default, required = required)
}

load_to_save_from_path <- function(path, object_name = "to_save"){
  if (!file.exists(path)){
    stop(paste("Hot-start file does not exist:",path))
  }

  if (grepl("\\.rds$",path,ignore.case = TRUE)){
    return(readRDS(path))
  }

  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(path,envir = load_env)

  if (object_name %in% loaded_names){
    return(get(object_name,envir = load_env))
  }
  if (length(loaded_names) == 1){
    return(get(loaded_names[[1]],envir = load_env))
  }

  stop(paste("Could not identify",object_name,"in",path,
             "loaded objects:",paste(loaded_names,collapse = ", ")))
}

make_hot_start_variables <- function(to_save, fit_mix_num = NULL,
                                     vcovar_num = NULL,
                                     source_name = "saved results",
                                     lambda_act_default = NULL,
                                     lambda_light_default = NULL,
                                     assign_true_params = FALSE){
  validate_saved_results(to_save,required_sections = c("est_params"),
                         source_name = source_name)

  loaded_est_params <- get_saved_section(to_save,"est_params")
  loaded_init <- get_saved_param(loaded_est_params,"init")
  loaded_params_tran_array <- get_saved_param(loaded_est_params,
                                              "params_tran_array")
  loaded_emit_act <- get_saved_param(loaded_est_params,"emit_act")
  loaded_emit_light <- get_saved_param(loaded_est_params,"emit_light")
  loaded_corr_mat <- get_saved_param(loaded_est_params,"corr_mat")
  loaded_nu_mat <- get_saved_param(loaded_est_params,"nu_mat")
  loaded_beta_vec <- get_saved_param(loaded_est_params,"beta_vec")
  loaded_surv_coef <- get_saved_param(loaded_est_params,"surv_coef")
  loaded_re_prob <- get_saved_param(loaded_est_params,"re_prob",
                                    required = FALSE)
  loaded_lambda_act_mat <- get_saved_param(loaded_est_params,
                                           "lambda_act_mat",
                                           default = lambda_act_default,
                                           required = FALSE)
  loaded_lambda_light_mat <- get_saved_param(loaded_est_params,
                                             "lambda_light_mat",
                                             default = lambda_light_default,
                                             required = FALSE)

  hot_start_params <- make_start_param_list(
    init = loaded_init,
    params_tran_array = loaded_params_tran_array,
    emit_act = loaded_emit_act,
    emit_light = loaded_emit_light,
    corr_mat = loaded_corr_mat,
    nu_mat = loaded_nu_mat,
    beta_vec = loaded_beta_vec,
    surv_coef = loaded_surv_coef,
    lambda_act_mat = loaded_lambda_act_mat,
    lambda_light_mat = loaded_lambda_light_mat
  )

  if (!is.null(fit_mix_num) && !is.null(vcovar_num)){
    validate_param_list(hot_start_params,fit_mix_num,vcovar_num,
                        "hot_start_params")
  }

  hot_start_vars <- list(
    loaded_est_params = loaded_est_params,
    loaded_init = loaded_init,
    loaded_params_tran_array = loaded_params_tran_array,
    loaded_emit_act = loaded_emit_act,
    loaded_emit_light = loaded_emit_light,
    loaded_corr_mat = loaded_corr_mat,
    loaded_nu_mat = loaded_nu_mat,
    loaded_beta_vec = loaded_beta_vec,
    loaded_lambda_act_mat = loaded_lambda_act_mat,
    loaded_lambda_light_mat = loaded_lambda_light_mat,
    init_start = loaded_init,
    params_tran_array_start = loaded_params_tran_array,
    emit_act_start = loaded_emit_act,
    emit_light_start = loaded_emit_light,
    corr_mat_start = loaded_corr_mat,
    nu_mat_start = loaded_nu_mat,
    beta_vec_start = loaded_beta_vec,
    lambda_act_mat_start = loaded_lambda_act_mat,
    lambda_light_mat_start = loaded_lambda_light_mat,
    surv_coef_true = loaded_surv_coef,
    re_prob_true = loaded_re_prob,
    re_prob = loaded_re_prob,
    hot_start_params = hot_start_params
  )

  if (assign_true_params){
    hot_start_vars <- c(
      hot_start_vars,
      list(init_true = loaded_init,
           params_tran_array_true = loaded_params_tran_array,
           emit_act_true = loaded_emit_act,
           emit_light_true = loaded_emit_light,
           corr_mat_true = loaded_corr_mat,
           nu_mat_true = loaded_nu_mat,
           beta_vec_true = loaded_beta_vec,
           lambda_act_mat_true = loaded_lambda_act_mat,
           lambda_light_mat_true = loaded_lambda_light_mat)
    )
  }

  hot_start_vars
}

load_hot_start_from_path <- function(path, fit_mix_num = NULL,
                                     vcovar_num = NULL,
                                     assign_to = parent.frame(),
                                     assign_true_params = FALSE,
                                     lambda_act_default = NULL,
                                     lambda_light_default = NULL){
  normalized_path <- normalizePath(path,mustWork = TRUE)
  to_save <- load_to_save_from_path(normalized_path)
  hot_start_vars <- make_hot_start_variables(
    to_save,
    fit_mix_num = fit_mix_num,
    vcovar_num = vcovar_num,
    source_name = normalized_path,
    lambda_act_default = lambda_act_default,
    lambda_light_default = lambda_light_default,
    assign_true_params = assign_true_params
  )
  hot_start_vars$hot_start_path <- normalized_path
  hot_start_vars$hot_start_to_save <- to_save

  if (!is.null(assign_to)){
    list2env(hot_start_vars,envir = assign_to)
  }

  invisible(hot_start_vars)
}

make_true_param_list <- function(init, params_tran_array, emit_act, emit_light,
                                 corr_mat, nu_mat, beta_vec, beta_age,
                                 lambda_act_mat = NULL, lambda_light_mat = NULL){
  list(init = init,
       params_tran_array = params_tran_array,
       emit_act = emit_act,
       emit_light = emit_light,
       corr_mat = corr_mat,
       nu_mat = nu_mat,
       beta_vec = beta_vec,
       beta_age = beta_age,
       lambda_act_mat = lambda_act_mat,
       lambda_light_mat = lambda_light_mat)
}

make_start_param_list <- function(init, params_tran_array, emit_act, emit_light,
                                  corr_mat, nu_mat, beta_vec, surv_coef,
                                  lambda_act_mat = NULL, lambda_light_mat = NULL){
  list(init = init,
       params_tran_array = params_tran_array,
       emit_act = emit_act,
       emit_light = emit_light,
       corr_mat = corr_mat,
       nu_mat = nu_mat,
       beta_vec = beta_vec,
       surv_coef = surv_coef,
       lambda_act_mat = lambda_act_mat,
       lambda_light_mat = lambda_light_mat)
}

make_est_param_list <- function(init, params_tran_array, emit_act, emit_light,
                                corr_mat, nu_mat, beta_vec, surv_coef,
                                tran_df = NULL, re_prob = NULL,
                                new_likelihood = NULL, decoded_mat = NULL,
                                lambda_act_mat = NULL, lambda_light_mat = NULL,
                                bline_vec = NULL, cbline_vec = NULL,
                                beta_se = NULL, post_decode = NULL){
  list(init = init,
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
}

make_saved_results <- function(true_params, est_params, bic = NULL,
                               leave_out = list(), simulated_hmm = list(),
                               diagnostics = list(), settings = NULL,
                               start_params = NULL, aic = NULL){
  list(true_params = true_params,
       est_params = est_params,
       bic = bic,
       leave_out = leave_out,
       simulated_hmm = simulated_hmm,
       diagnostics = diagnostics,
       settings = settings,
       start_params = start_params,
       aic = aic)
}

make_simulated_hmm_list <- function(mc, act, light, mixture_mat, age_vec,
                                    nu_covar_mat, vcovar_mat, survival,
                                    surv_covar_sim){
  list(mc = mc,
       act = act,
       light = light,
       mixture_mat = mixture_mat,
       age_vec = age_vec,
       nu_covar_mat = nu_covar_mat,
       vcovar_mat = vcovar_mat,
       survival = survival,
       surv_covar_sim = surv_covar_sim)
}

make_leave_out_results <- function(leave_out_inds, conf_mat_list, cindex_new_list,
                                   ibs_new_list, senspec_list, ibs2_new_list,
                                   senspec_mix_list){
  list(leave_out_inds = leave_out_inds,
       conf_mat_list = conf_mat_list,
       cindex_new_list = cindex_new_list,
       ibs_new_list = ibs_new_list,
       senspec_list = senspec_list,
       ibs2_new_list = ibs2_new_list,
       senspec_mix_list = senspec_mix_list)
}

make_diagnostics_list <- function(cindex, ibs, confusion_table, ibs2){
  list(cindex = cindex,
       ibs = ibs,
       confusion_table = confusion_table,
       ibs2 = ibs2)
}
