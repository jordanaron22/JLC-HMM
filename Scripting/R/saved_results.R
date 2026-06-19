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
