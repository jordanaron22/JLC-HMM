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
                           beta_bool = beta_bool,
                           survival_baseline_mode = "profiled",
                           profile_maxit = 100L,
                           profile_tol = 1e-9,
                           profile_damping = 1.0)

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

if (oakes_data_context$survival_baseline_mode == "fixed") {
  surv_h1 <- CalcOakesSurvivalH1(
    beta_vec = oakes_params$beta_vec,
    surv_coef = oakes_params$surv_coef,
    survival_context = posterior_context$survival_context,
    fit_mix_num = mix_num,
    cbline_vec = posterior_context$cbline_vec)
} else {
  surv_h1 <- CalcOakesProfiledSurvivalH1(
    beta_vec = oakes_params$beta_vec,
    surv_coef = oakes_params$surv_coef,
    survival_context = posterior_context$survival_context,
    fit_mix_num = mix_num,
    surv_covar_risk_vec = posterior_context$surv_covar_risk_vec
  )
}


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



h2_1e5 <- CalcOakesH2(
  theta_pack = theta_pack,
  data_context = oakes_data_context,
  eps = 1e-5
)

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



I <- I_obs_sym
pm <- theta_pack$parameter_map

surv_idx <- which(pm$block == "survival")
nuis_idx <- setdiff(seq_len(nrow(I)), surv_idx)

I_ss <- I[surv_idx, surv_idx, drop = FALSE]
I_nn <- I[nuis_idx, nuis_idx, drop = FALSE]
I_ns <- I[nuis_idx, surv_idx, drop = FALSE]

S_surv <- I_ss - t(I_ns) %*% solve(I_nn, I_ns)
S_surv <- 0.5 * (S_surv + t(S_surv))

eigen(S_surv, symmetric = TRUE)$values

vcov_surv <- solve(S_surv)
se_surv <- sqrt(diag(vcov_surv))

data.frame(
  parameter = pm$param_name[surv_idx],
  se = se_surv
)



s_surv <- oakes_score$score[surv_idx]
delta_surv <- solve(S_surv, s_surv)
vcov_surv <- solve(S_surv)
se_surv <- sqrt(diag(vcov_surv))

score_check <- data.frame(
  parameter = pm$param_name[surv_idx],
  score = s_surv,
  newton_step = delta_surv,
  se = se_surv,
  step_over_se = delta_surv / se_surv
)

score_check
max(abs(score_check$step_over_se))
max(abs(delta_surv / se_surv)) < 0.05
