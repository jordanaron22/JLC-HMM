#ToDo better load in, right now I just run JMHMM first

theta_pack <- PackOakesTheta(init = init,
                             params_tran_array = params_tran_array,
                             emit_act = emit_act,
                             emit_light = emit_light,
                             corr_mat = corr_mat,
                             nu_mat = nu_mat,
                             beta_vec = beta_vec,
                             surv_coef = surv_coef,
                             vcovar_mat = vcovar_mat,
                             fit_mix_num = mix_num)

oakes_data_context <- list(act = act,
                           light = light,
                           vcovar_mat = vcovar_mat,
                           lod_act = lod_act,
                           lod_light = lod_light,
                           log_sweights_vec = log_sweights_vec,
                           lambda_act_mat = lambda_act_mat,
                           lambda_light_mat = lambda_light_mat,
                           tobit = tobit,
                           period_len = period_len,
                           nu_covar_mat = nu_covar_mat,
                           survival_context = survival_context,
                           re_prob = re_prob,
                           bline_vec = bline_vec,
                           cbline_vec = cbline_vec,
                           incl_surv = incl_surv,
                           beta_bool = beta_bool)

posterior_context <- RebuildOakesPosteriorContext(
  theta = theta_pack$theta,
  theta_pack = theta_pack,
  data_context = oakes_data_context)

oakes_params <- MakeInactiveDayTypesSafe(
  UnpackOakesTheta(theta_pack$theta,theta_pack),
  vcovar_mat)

init_counts <- posterior_context$init_counts

init_h1 <- CalcOakesInitialH1(init = oakes_params$init,
                              init_counts = init_counts,
                              sleep_state = SLEEP_STATE)

tran_h1 <- CalcOakesTransitionH1(alpha = posterior_context$alpha,
                                 beta = posterior_context$beta,
                                 act = act,
                                 light = light,
                                 params_tran_array =
                                   oakes_params$params_tran_array,
                                 emit_act = oakes_params$emit_act,
                                 emit_light = oakes_params$emit_light,
                                 corr_mat = oakes_params$corr_mat,
                                 pi_l = posterior_context$pi_l,
                                 lod_act = lod_act,
                                 lod_light = lod_light,
                                 lintegral_mat =
                                   posterior_context$lintegral_mat,
                                 vcovar_mat = vcovar_mat,
                                 lambda_act_mat = lambda_act_mat,
                                 lambda_light_mat = lambda_light_mat,
                                 tobit = tobit,
                                 period_len = period_len)

mix_h1 <- CalcOakesMixingH1(nu_mat = oakes_params$nu_mat,
                            re_prob = posterior_context$re_prob,
                            nu_covar_mat = nu_covar_mat)

surv_h1 <- CalcOakesSurvivalH1(
  beta_vec = oakes_params$beta_vec,
  surv_coef = oakes_params$surv_coef,
  survival_context = posterior_context$survival_context,
  fit_mix_num = mix_num,
  cbline_vec = posterior_context$cbline_vec)

emit_h1 <- CalcOakesEmissionH1(alpha = posterior_context$alpha,
                               beta = posterior_context$beta,
                               pi_l = posterior_context$pi_l,
                               act = act,
                               light = light,
                               vcovar_mat = vcovar_mat,
                               emit_act = oakes_params$emit_act,
                               emit_light = oakes_params$emit_light,
                               corr_mat = oakes_params$corr_mat,
                               lod_act = lod_act,
                               lod_light = lod_light)

oakes_score <- CalcOakesScore(theta = theta_pack$theta,
                              theta_pack = theta_pack,
                              posterior_context = posterior_context,
                              return_components = TRUE)

posterior_diagnostics <- DiagnosePosteriorContext(posterior_context)
score_diagnostics <- DiagnoseOakesScore(oakes_score)

init_h1$hessian
tran_h1$hessian
mix_h1$hessian
surv_h1$hessian
emit_h1$hessian
oakes_score$score
posterior_diagnostics
score_diagnostics

H1 <- bdiag(init_h1$hessian,tran_h1$hessian,emit_h1$hessian,
           mix_h1$hessian,surv_h1$hessian)
# solve(-H1)

tic()
h2_1e5 <- CalcOakesH2(theta_pack = theta_pack,
                  data_context = oakes_data_context,
                  eps = 1e-5)
toc()

H2 <- h2_1e5$hessian

H_obs <- as.matrix(H1) + H2
I_obs <- -H_obs
I_obs_sym <- (I_obs + t(I_obs)) / 2

info_diagnostics <- list(
  h1_dim = dim(H1),
  h2_dim = dim(H2),
  finite_h1 = all(is.finite(H1)),
  finite_h2 = all(is.finite(H2)),
  max_abs_info_asymmetry = max(abs(I_obs - t(I_obs)),na.rm = TRUE),
  observed_info_eigen_range =
    range(eigen(I_obs_sym,symmetric = TRUE,only.values = TRUE)$values),
  nonpositive_observed_info_eigen_count =
    sum(eigen(I_obs_sym,symmetric = TRUE,only.values = TRUE)$values <= 0)
)

vcov = solve(I_obs_sym)

bad <- which(diag(vcov) < 0)

se = sqrt(diag(vcov))

# DiagnoseOakesInformation(H1, H2)


h2_1e4_new <- CalcOakesH2(
  theta_pack = theta_pack,
  data_context = oakes_data_context,
  eps = 1e-4
)

H2_1e4 <- h2_1e4_new$hessian
DiagnoseOakesInformation(H1, H2_1e4)

eig <- eigen(I_obs_sym, symmetric = TRUE)

tol <- 1e-8 * max(eig$values)
vals_fixed <- pmax(eig$values, tol)

I_obs_pd <- eig$vectors %*% diag(vals_fixed) %*% t(eig$vectors)
I_obs_pd <- 0.5 * (I_obs_pd + t(I_obs_pd))

vcov_pd <- solve(I_obs_pd)
se_pd <- sqrt(diag(vcov_pd))



pm <- theta_pack$parameter_map

surv_class_idx <- which(
  pm$block == "survival" &
    grepl("^class_", pm$param_name)
)

pm[surv_class_idx, c("index", "block", "param_name")]
se_pd[surv_class_idx]


I_h1_only <- -as.matrix(H1)
vcov_h1_only <- solve(I_h1_only)
se_h1_only <- sqrt(diag(vcov_h1_only))

compare_surv_class_se <- data.frame(
  index = surv_class_idx,
  parameter = pm$param_name[surv_class_idx],
  se_h1_only = se_h1_only[surv_class_idx],
  se_oakes_pd = se_pd[surv_class_idx],
  ratio_oakes_to_h1 = se_pd[surv_class_idx] / se_h1_only[surv_class_idx]
)

compare_surv_class_se
