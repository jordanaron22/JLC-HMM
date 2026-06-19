CLI_ARG <- list(
  fit_mix_num = 1,
  model_type = 2,
  data_source = 3,
  run_bootstrap = 4,
  init_jitter_scale = 5,
  run_leave_one_out_cv = 6,
  use_hot_start = 7,
  sim_scenario = 8,
  true_mix_num = 9,
  save_reduced_output = 10,
  class_selection_run = 11,
  emission_overlap = 12
)

MODEL_TYPE_CODES <- c(two_stage = 0, joint = 2)
DATA_SOURCE <- c(simulation = "simulation", nhanes = "nhanes")
WEEKDAY_CODES <- c(all = 0, sunday = 1, monday = 2, tuesday = 3,
                   wednesday = 4, thursday = 5, friday = 6, saturday = 7)

TOBIT_EMISSION <- 1
NUM_MARKOV_STATES <- 2
WAKE_STATE <- 1
SLEEP_STATE <- 2
WEEKDAY_LEVELS <- 2
DEFAULT_PERIODS_PER_DAY <- 96
HOURLY_PERIODS_PER_DAY <- 24
MINUTE_PERIODS_PER_DAY <- 1440
LOD_OFFSET <- 1e-5
NHANES_NUM_WAVES <- 2
BASE_STOP_CRIT <- 1e-3
SIM_STOP_CRIT_MULTIPLIER <- 10
MIN_EM_ITERATIONS <- 5
INTERIM_SAVE_EVERY <- 10
REORDER_STOP_CRIT_MULTIPLIER <- 10
# TODO: After the outer EM loop converges, consider rerunning the joint survival
# update with a stricter tolerance to polish survival coefficients, then
# recompute baseline hazards, forward/backward probabilities, and likelihood.
JOINT_BETA_STOP_CRIT <- 100

EMISSION_OVERLAP_FACTOR <- c(low = 1, high = 0.25)

SIM_SCENARIOS <- list(
  `-1` = list(days = 1, num_people = 600, missing_perc = 0.2),
  `0` = list(days = 1, num_people = 6000, missing_perc = 0.2),
  `1` = list(days = 3, num_people = 6000, missing_perc = 0.2),
  `2` = list(days = 7, num_people = 6000, missing_perc = 0.2),
  `3` = list(days = 1, num_people = 1000, missing_perc = 0.2),
  `4` = list(days = 3, num_people = 1000, missing_perc = 0.2),
  `5` = list(days = 7, num_people = 1000, missing_perc = 0.2)
)

SAVED_SECTION_SLOT <- c(true_params = 1, est_params = 2, bic = 3,
                        leave_out = 4, simulated_hmm = 5, diagnostics = 6,
                        settings = 7, start_params = 8, aic = 9)

PARAM_SLOT <- c(init = 1, params_tran_array = 2, emit_act = 3,
                emit_light = 4, corr_mat = 5, nu_mat = 6,
                beta_vec = 7, surv_coef = 8, beta_age = 8,
                tran_df = 9, re_prob = 10, new_likelihood = 11,
                decoded_mat = 12, lambda_act_mat = 13,
                lambda_light_mat = 14, bline_vec = 15, cbline_vec = 16,
                beta_se = 17, post_decode = 18)
