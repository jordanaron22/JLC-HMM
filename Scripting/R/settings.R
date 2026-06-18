default_settings <- list(
  use_seed = TRUE,
  include_activity = TRUE,
  include_light = TRUE,
  check_transition_update = FALSE,
  weekend_only = FALSE,
  save_reduced_output = FALSE,
  # 0 -> all days, 1-7 -> subset to that day of the week only.
  target_weekday = 0,
  survival_only = FALSE,
  include_tlogmims = FALSE,
  # 24 hourly, 96 fifteen min, 1440 single min.
  periods_per_day = DEFAULT_PERIODS_PER_DAY,
  num_day_type_levels = WEEKDAY_LEVELS,
  # Command-line defaults used when not running through bash/slurm.
  fit_mix_num = 5,
  true_mix_num = NA_real_, # set to  NA_real_ to match fit_mix_num by default
  model_type = "joint", # "joint", "two_stage"
  data_source = "simulation", #"nhanes" or "simulation"
  run_bootstrap = FALSE,
  init_jitter_scale = 0,
  run_leave_one_out_cv = FALSE,
  use_hot_start = FALSE,
  # Compact simulation scenario code used by the shell grid.
  sim_scenario = -1
)

get_numeric_arg <- function(cli_args, arg_index){
  if (length(cli_args) < arg_index){
    return(NA_real_)
  }
  as.numeric(cli_args[arg_index])
}

get_arg <- function(cli_args, arg_index){
  if (length(cli_args) < arg_index){
    return(NA_character_)
  }
  cli_args[arg_index]
}

is_missing_setting <- function(value){
  length(value) == 0 || is.na(value) || as.character(value) == ""
}

add_legacy_settings <- function(settings){
  # Legacy aliases: the rest of this script still reads these historical names.
  settings$set_seed <- as.numeric(settings$use_seed)
  # Internal constant for legacy helper signatures; Tobit is now the only emission model.
  settings$tobit <- TOBIT_EMISSION
  settings$incl_act <- as.numeric(settings$include_activity)
  settings$incl_light <- as.numeric(settings$include_light)
  settings$check_tran <- as.numeric(settings$check_transition_update)
  settings$save_space <- as.numeric(settings$save_reduced_output)
  settings$single_day <- settings$target_weekday
  settings$run_only_surv <- as.numeric(settings$survival_only)
  settings$tlogmims_bool <- as.numeric(settings$include_tlogmims)
  settings$period_len <- settings$periods_per_day
  settings$vcovar_num <- settings$num_day_type_levels
  settings$incl_surv <- MODEL_TYPE_CODES[[settings$model_type]]
  settings$real_data <- as.numeric(settings$data_source == DATA_SOURCE[["nhanes"]])
  settings$bootstrap <- as.numeric(settings$run_bootstrap)
  settings$randomize_init <- settings$init_jitter_scale
  settings$leave_out <- as.numeric(settings$run_leave_one_out_cv)
  settings$load_data <- as.numeric(settings$use_hot_start)
  settings$sim_size <- settings$sim_scenario
  # Legacy alias: most fitted-model calculations still read mix_num internally.
  settings$mix_num <- settings$fit_mix_num

  # Start by not estimating survival coef for stability, unless doing hot start.
  settings$beta_bool <- as.numeric(settings$use_hot_start | settings$run_leave_one_out_cv)
  if (settings$run_bootstrap){settings$beta_bool <- 0}

  # Switch to order clusters from best to worst survival.
  settings$relabel_reset <- FALSE

  settings
}

build_settings <- function(cli_args = commandArgs(TRUE),
                           sim_num = suppressWarnings(as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID')))){
  if (is.na(sim_num)){sim_num <- 1}

  settings <- default_settings
  settings$sim_num <- sim_num
  settings$command_args <- cli_args

  command_settings <- list(
    fit_mix_num = get_numeric_arg(cli_args,CLI_ARG$fit_mix_num),
    model_type = get_arg(cli_args,CLI_ARG$model_type),
    data_source = get_arg(cli_args,CLI_ARG$data_source),
    run_bootstrap = get_arg(cli_args,CLI_ARG$run_bootstrap),
    init_jitter_scale = get_numeric_arg(cli_args,CLI_ARG$init_jitter_scale),
    run_leave_one_out_cv = get_arg(cli_args,CLI_ARG$run_leave_one_out_cv),
    use_hot_start = get_arg(cli_args,CLI_ARG$use_hot_start),
    sim_scenario = get_numeric_arg(cli_args,CLI_ARG$sim_scenario),
    true_mix_num = get_numeric_arg(cli_args,CLI_ARG$true_mix_num)
  )

  for (setting_name in names(command_settings)){
    if (!is_missing_setting(command_settings[[setting_name]])){
      settings[[setting_name]] <- command_settings[[setting_name]]
    }
  }

  settings <- validate_settings(settings, MODEL_TYPE_CODES, DATA_SOURCE,
                                WEEKDAY_CODES, SIM_SCENARIOS)

  if (is.na(settings$true_mix_num)){
    settings$true_mix_num <- settings$fit_mix_num
  }

  add_legacy_settings(settings)
}

build_model_name <- function(settings){
  model_name <- "JMHMM"
  if (!settings$real_data){model_name <- paste0(model_name,"SimSize",settings$sim_size)}
  if (settings$bootstrap){model_name <- paste0(model_name,"Bootstrap")}
  if (settings$leave_out){model_name <- paste0(model_name,"LeaveOut")}
  if (settings$incl_surv == MODEL_TYPE_CODES[["two_stage"]]){model_name <- paste0(model_name,"NoSurv")}
  if (!settings$incl_light){model_name <- paste0(model_name,"NoLight")}
  if (!settings$incl_act){model_name <- paste0(model_name,"NoAct")}
  if (settings$weekend_only){model_name <- paste0(model_name,"WeekendOnly")}
  if (settings$single_day){model_name <- paste0(model_name,"SingleDay",settings$single_day)}
  if (settings$randomize_init){model_name <- paste0(model_name,"RandInit")}
  if (settings$load_data){model_name <- paste0(model_name,"LoadIn")}

  if (settings$real_data){
    paste0(model_name,"FitMix",settings$fit_mix_num,"Seed",
           settings$sim_num,"len",settings$period_len,".rda")
  } else {
    paste0(model_name,"TrueMix",settings$true_mix_num,"FitMix",
           settings$fit_mix_num,"Seed",settings$sim_num,
           "len",settings$period_len,".rda")
  }
}
