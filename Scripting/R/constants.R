CLI_ARG <- list(
  fit_mix_num = 1,
  model_type = 2,
  data_source = 3,
  run_bootstrap = 4,
  init_jitter_scale = 5,
  run_leave_one_out_cv = 6,
  use_hot_start = 7,
  simulation_days = 8,
  num_people = 9,
  true_mix_num = 10,
  save_reduced_output = 11,
  class_selection_run = 12,
  emission_overlap = 13,
  time_limit_hours = 14,
  memory_limit_gb = 15
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
NHANES_CV_FOLD_COUNT <- 20
NHANES_CV_FOLD_SEED <- 999
CV_SURVIVAL_INTERVAL_BREAKS <- seq(0,102,by = 6)
CV_IBS_EVAL_TIMES <- seq(0,102,by = 1)
CV_LEAVE_OUT_SCENARIOS <- c(sequence_only = 1, no_light = 2,
                            no_activity = 3, no_transition = 4,
                            activity_only = 5, standard = 6)
BASE_STOP_CRIT <- 1e-8
SIM_STOP_CRIT_MULTIPLIER <- 10
MIN_EM_ITERATIONS <- 5
INTERIM_SAVE_EVERY <- 10
REORDER_STOP_CRIT_MULTIPLIER <- 10
# TODO: After the outer EM loop converges, consider rerunning the joint survival
# update with a stricter tolerance to polish survival coefficients, then
# recompute baseline hazards, forward/backward probabilities, and likelihood.
JOINT_BETA_STOP_CRIT <- 100

EMISSION_OVERLAP_FACTOR <- c(low = 1, mid = 0.75, high = 0.5)
DEFAULT_SIMULATION_DAYS <- 1
DEFAULT_SIMULATION_PEOPLE <- 600
DEFAULT_MISSING_PERC <- 0.2

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
