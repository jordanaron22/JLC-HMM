load_rda_into_env <- function(file_path){
  if (!file.exists(file_path)){
    stop(paste("Required data file not found:",file_path))
  }

  data_env <- new.env(parent = emptyenv())
  load(file_path,envir = data_env)
  data_env
}

make_nhanes_survival_covariates <- function(id){
  list(age = id$age,
       gender = Vec2Mat(id$gender+1),
       race = Vec2Mat(id$race+1),
       overall_health = Vec2Mat(id$overall_health+1),
       education = Vec2Mat(id$education+1),
       bmi = Vec2Mat(id$bmi_disc+1),
       diabetes = Vec2Mat(id$diabetes+1),
       chd = Vec2Mat(id$CHD+1),
       chf = Vec2Mat(id$CHF+1),
       heart_attack = Vec2Mat(id$heart_attack+1),
       stroke = Vec2Mat(id$stroke+1),
       alcohol = Vec2Mat(id$alcohol+1),
       smoking = Vec2Mat(id$smoking+1),
       physical_function = Vec2Mat(id$phyfunc+1))
}

assign_indices_to_folds <- function(indices,fold_count){
  if (length(indices) == 0){
    return(data.frame(subject_index = integer(),fold_id = integer()))
  }

  shuffled_indices <- indices[sample.int(length(indices))]
  data.frame(
    subject_index = shuffled_indices,
    fold_id = rep(seq_len(fold_count),length.out = length(shuffled_indices))
  )
}

ValidateNHANESCVFolds <- function(cv_folds,fold_count,interval_breaks){
  required_cols <- c("subject_index","fold_id","surv_event",
                     "event_interval_index")
  missing_cols <- setdiff(required_cols,names(cv_folds))
  if (length(missing_cols) > 0){
    stop("cv_folds is missing required columns: ",
         paste(missing_cols,collapse = ", "))
  }

  if (!all(seq_len(fold_count) %in% cv_folds$fold_id)){
    stop("Every NHANES CV fold must contain at least one participant")
  }

  event_interval_levels <- seq_len(length(interval_breaks) - 1L)
  training_event_counts <- sapply(seq_len(fold_count),function(current_fold){
    training_rows <- cv_folds$fold_id != current_fold
    tabulate(
      cv_folds$event_interval_index[
        training_rows & cv_folds$surv_event == 1
      ],
      nbins = length(event_interval_levels)
    )
  })

  rownames(training_event_counts) <- paste0(
    interval_breaks[-length(interval_breaks)],"-",
    interval_breaks[-1]
  )
  colnames(training_event_counts) <- paste0("fold_",seq_len(fold_count))

  if (any(training_event_counts <= 0)){
    failing_positions <- which(training_event_counts <= 0,arr.ind = TRUE)
    failing_labels <- paste0(
      rownames(training_event_counts)[failing_positions[,1]],
      " in ",
      colnames(training_event_counts)[failing_positions[,2]]
    )
    stop(
      paste0(
        "Cannot create NHANES CV folds: at least one training fold ",
        "has no death in interval(s): ",
        paste(failing_labels,collapse = ", ")
      )
    )
  }

  invisible(training_event_counts)
}

MakeNHANESCVFolds <- function(surv_time,
                              surv_event,
                              fold_count = NHANES_CV_FOLD_COUNT,
                              interval_breaks = CV_SURVIVAL_INTERVAL_BREAKS,
                              seed = NHANES_CV_FOLD_SEED){
  num_people <- length(surv_time)

  if (length(surv_event) != num_people){
    stop("surv_time and surv_event must have equal lengths")
  }
  if (length(fold_count) != 1 ||
      is.na(fold_count) ||
      fold_count != floor(fold_count) ||
      fold_count < 2){
    stop("fold_count must be an integer >= 2")
  }
  if (any(!surv_event %in% c(0,1))){
    stop("surv_event must contain only 0 and 1")
  }
  if (length(interval_breaks) < 2 ||
      any(!is.finite(interval_breaks)) ||
      any(diff(interval_breaks) <= 0)){
    stop("interval_breaks must be a strictly increasing finite vector")
  }

  event_indicator <- surv_event == 1
  event_interval_index <- rep(NA_integer_,num_people)
  event_interval_index[event_indicator] <- cut(
    surv_time[event_indicator],
    breaks = interval_breaks,
    include.lowest = TRUE,
    right = TRUE,
    labels = FALSE
  )

  if (anyNA(event_interval_index[event_indicator])){
    stop("Every observed death must fall within CV_SURVIVAL_INTERVAL_BREAKS")
  }

  event_interval_counts <- tabulate(
    event_interval_index[event_indicator],
    nbins = length(interval_breaks) - 1L
  )
  if (any(event_interval_counts < 2)){
    interval_labels <- paste0(
      interval_breaks[-length(interval_breaks)],"-",
      interval_breaks[-1],
      ": ",
      event_interval_counts
    )
    stop(
      paste0(
        "Cannot create NHANES CV folds with at least one training death ",
        "per interval. Each interval needs at least two deaths; counts are ",
        paste(interval_labels,collapse = ", "),
        "."
      )
    )
  }

  old_seed <- if (exists(".Random.seed",envir = .GlobalEnv,
                         inherits = FALSE)){
    get(".Random.seed",envir = .GlobalEnv,inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)){
      if (exists(".Random.seed",envir = .GlobalEnv,inherits = FALSE)){
        rm(".Random.seed",envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed",old_seed,envir = .GlobalEnv)
    }
  },add = TRUE)
  set.seed(seed)

  fold_id <- integer(num_people)
  for (current_interval in seq_len(length(interval_breaks) - 1L)){
    interval_indices <- which(
      event_indicator & event_interval_index == current_interval
    )
    interval_assignment <- assign_indices_to_folds(
      interval_indices,fold_count
    )
    fold_id[interval_assignment$subject_index] <-
      interval_assignment$fold_id
  }

  censored_indices <- which(!event_indicator)
  censored_assignment <- assign_indices_to_folds(censored_indices,fold_count)
  fold_id[censored_assignment$subject_index] <- censored_assignment$fold_id

  if (any(fold_id == 0)){
    stop("Internal error: at least one participant was not assigned a fold")
  }

  cv_folds <- data.frame(
    subject_index = seq_len(num_people),
    fold_id = fold_id,
    surv_time = surv_time,
    surv_event = surv_event,
    event_interval_index = event_interval_index
  )

  ValidateNHANESCVFolds(cv_folds,fold_count,interval_breaks)
  cv_folds
}

prepare_nhanes_data <- function(period_len,bootstrap,leave_out,sim_num,
                                single_day,weekend_only,load_data,
                                surv_coef_true = NULL,data_dir = "Data"){
  mortality_env <- load_rda_into_env(
    file.path(data_dir,"NHANES_2011_2012_2013_2014.rda")
  )
  if (!exists("NHANES_mort_list",envir = mortality_env,inherits = FALSE)){
    stop("NHANES mortality file does not contain NHANES_mort_list")
  }

  nhanes_mort_list <- get("NHANES_mort_list",envir = mortality_env)
  nhanes1 <- nhanes_mort_list[[1]] %>% filter(eligstat == 1)
  nhanes2 <- nhanes_mort_list[[2]] %>% filter(eligstat == 1)
  lmf_data <- rbind(nhanes1,nhanes2)

  wave_files <- if (period_len == HOURLY_PERIODS_PER_DAY){
    c("Wavedata24_G.rda","Wavedata24_H.rda")
  } else if (period_len == DEFAULT_PERIODS_PER_DAY){
    c("Wavedata_G.rda","Wavedata_H.rda")
  } else if (period_len == MINUTE_PERIODS_PER_DAY){
    c("Wavedata1440_G.rda","Wavedata1440_H.rda")
  } else {
    stop(paste("Unsupported NHANES period_len:",period_len))
  }

  wave_g_env <- load_rda_into_env(file.path(data_dir,wave_files[1]))
  wave_h_env <- load_rda_into_env(file.path(data_dir,wave_files[2]))
  if (!exists("wave_data_G",envir = wave_g_env,inherits = FALSE) ||
      !exists("wave_data_H",envir = wave_h_env,inherits = FALSE)){
    stop("NHANES wave files do not contain wave_data_G and wave_data_H")
  }

  wave_data_G <- get("wave_data_G",envir = wave_g_env)
  wave_data_H <- get("wave_data_H",envir = wave_h_env)

  act <- t(rbind(wave_data_G[[1]],wave_data_H[[1]])[,-1])
  act0 <- act == 0
  act <- log(act)
  lod_act <- min(act[act != -Inf],na.rm = TRUE) - LOD_OFFSET
  act[act0] <- lod_act

  light <- t(rbind(wave_data_G[[2]],wave_data_H[[2]])[,-1])
  light0 <- light == 0
  light <- log(light)
  lod_light <- min(light[light != -Inf],na.rm = TRUE) - LOD_OFFSET
  light[light0] <- lod_light

  id <- rbind(wave_data_G[[3]],wave_data_H[[3]])
  mims <- rbind(wave_data_G[[4]],wave_data_H[[4]])

  seqn_com_id <- id$SEQN %in% lmf_data$seqn
  seqn_com_lmf <- lmf_data$seqn %in% id$SEQN

  id <- id[seqn_com_id,]
  act <- act[,seqn_com_id,drop = FALSE]
  light <- light[,seqn_com_id,drop = FALSE]
  mims <- mims[seqn_com_id,]
  lmf_data <- lmf_data[seqn_com_lmf,]

  if (sum(id$SEQN - lmf_data$seqn) != 0){
    stop("LMF data are not linked to NHANES actigraphy records correctly")
  }

  sweights_vec_raw <- id$sweights/NHANES_NUM_WAVES
  sweights_vec <- sweights_vec_raw/mean(sweights_vec_raw)

  id <- id %>% mutate(age_disc = case_when(age <= 30 ~ 1,
                                           age <= 50 & age > 30 ~ 2,
                                           age <= 65 & age > 50 ~ 3,
                                           age > 65 ~ 4))
  id <- id %>% mutate(pov_disc = floor(poverty)+1)
  id$modact <- id$modact - 1

  surv_event <- lmf_data$mortstat
  surv_time <- lmf_data$permth_exm

  if (bootstrap){
    boot_inds <- sample(ncol(act),ncol(act),replace = TRUE)
    act <- act[,boot_inds,drop = FALSE]
    light <- light[,boot_inds,drop = FALSE]
    id <- id[boot_inds,]
    surv_event <- surv_event[boot_inds]
    surv_time <- surv_time[boot_inds]
    sweights_vec <- id$sweights/NHANES_NUM_WAVES
    sweights_vec <- sweights_vec/mean(sweights_vec)
  }

  act_old <- act
  light_old <- light

  leave_out_data <- list()
  if (leave_out){
    cv_fold_count <- NHANES_CV_FOLD_COUNT
    cv_fold_id <- as.integer(sim_num)
    if (is.na(cv_fold_id) ||
        cv_fold_id < 1 ||
        cv_fold_id > cv_fold_count){
      stop(
        paste0(
          "For NHANES CV, sim_num is the fold_id and must be between 1 and ",
          cv_fold_count,
          "; got ",
          sim_num
        )
      )
    }

    cv_folds <- MakeNHANESCVFolds(
      surv_time = surv_time,
      surv_event = surv_event,
      fold_count = cv_fold_count,
      interval_breaks = CV_SURVIVAL_INTERVAL_BREAKS,
      seed = NHANES_CV_FOLD_SEED
    )
    cv_training_event_counts <- ValidateNHANESCVFolds(
      cv_folds = cv_folds,
      fold_count = cv_fold_count,
      interval_breaks = CV_SURVIVAL_INTERVAL_BREAKS
    )
    leave_out_inds <- cv_folds$subject_index[
      cv_folds$fold_id == cv_fold_id
    ]

    id_old <- id
    surv_event_old <- surv_event
    surv_time_old <- surv_time
    sweights_vec_old <- sweights_vec

    first_day_vec_old <- as.numeric(id_old$PAXDAYWM)
    vcovar_mat_old <- sapply(first_day_vec_old,FirstDay2WeekInd)
    surv_covar_old <- make_nhanes_survival_covariates(id_old)

    age_vec_old <- id_old$age
    statact_vec_old <- id_old$statact
    nu_covar_mat_old <- cbind(age_vec_old/10,(age_vec_old/10)^2,
                              statact_vec_old,statact_vec_old^2)

    act <- act[,-leave_out_inds,drop = FALSE]
    light <- light[,-leave_out_inds,drop = FALSE]
    id <- id[-leave_out_inds,,drop = FALSE]
    surv_event <- surv_event[-leave_out_inds]
    surv_time <- surv_time[-leave_out_inds]
    sweights_vec <- id$sweights/NHANES_NUM_WAVES
    sweights_vec <- sweights_vec/mean(sweights_vec)

    leave_out_data <- list(leave_out_inds = leave_out_inds,
                           cv_fold_id = cv_fold_id,
                           cv_fold_count = cv_fold_count,
                           cv_folds = cv_folds,
                           cv_training_event_counts =
                             cv_training_event_counts,
                           id_old = id_old,
                           surv_event_old = surv_event_old,
                           surv_time_old = surv_time_old,
                           sweights_vec_old = sweights_vec_old,
                           vcovar_mat_old = vcovar_mat_old,
                           surv_covar_old = surv_covar_old,
                           nu_covar_mat_old = nu_covar_mat_old)
  }

  first_day_vec <- as.numeric(id$PAXDAYWM)
  vcovar_mat <- sapply(first_day_vec,FirstDay2WeekInd)

  if (single_day != WEEKDAY_CODES[["all"]]){
    single_day_mat <- sapply(first_day_vec,FirstDay2SingleDay,
                             target_day = single_day)
    new_act <- matrix(NA,period_len,ncol(act))
    new_light <- matrix(NA,period_len,ncol(light))
    vcovar_mat <- matrix(0,period_len,ncol(light))

    # Extract the requested weekday for each participant.
    for (i in seq_len(ncol(act))){
      new_act[,i] <- act[,i][single_day_mat[,i] == 1]
      new_light[,i] <- light[,i][single_day_mat[,i] == 1]
    }

    act <- new_act
    light <- new_light
  }

  day_length <- nrow(act)
  num_of_people <- ncol(act)

  age_vec <- id$age
  modact_vec <- id$modact
  statact_vec <- id$statact
  nu_covar_mat <- cbind(age_vec/10,(age_vec/10)^2,
                        statact_vec,statact_vec^2)
  surv_covar <- make_nhanes_survival_covariates(id)

  if (!load_data){
    surv_coef_true <- lapply(surv_covar[-1],SurvCovar2Coef)
    surv_coef_true <- append(list(.05),surv_coef_true)
  } else if (is.null(surv_coef_true)){
    stop("surv_coef_true is required when loading NHANES starting values")
  }

  surv_coef_len <- unlist(lapply(surv_coef_true,length))
  surv_coef <- surv_coef_true

  combined_covar_mat <- id %>%
    dplyr::select(gender,race,overall_health,education,bmi_disc,diabetes,
                  race,CHD,CHF,heart_attack,stroke,alcohol,smoking,phyfunc)
  combined_covar_mat <- lapply(combined_covar_mat,factor)

  if (weekend_only){
    act[vcovar_mat == 0] <- NA
    light[vcovar_mat == 0] <- NA
  }

  c(list(act = act,
         light = light,
         id = id,
         mims = mims,
         lod_act = lod_act,
         lod_light = lod_light,
         sweights_vec = sweights_vec,
         surv_event = surv_event,
         surv_time = surv_time,
         act_old = act_old,
         light_old = light_old,
         first_day_vec = first_day_vec,
         vcovar_mat = vcovar_mat,
         day_length = day_length,
         num_of_people = num_of_people,
         age_vec = age_vec,
         modact_vec = modact_vec,
         statact_vec = statact_vec,
         nu_covar_mat = nu_covar_mat,
         surv_covar = surv_covar,
         surv_coef_true = surv_coef_true,
         surv_coef_len = surv_coef_len,
         surv_coef = surv_coef,
         combined_covar_mat = combined_covar_mat),
    leave_out_data)
}
