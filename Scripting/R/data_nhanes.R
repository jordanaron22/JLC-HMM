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
    leave_out_env <- load_rda_into_env(file.path(data_dir,"LeaveOutMat.rda"))
    if (!exists("leave_out_mat",envir = leave_out_env,inherits = FALSE)){
      stop("LeaveOutMat.rda does not contain leave_out_mat")
    }

    leave_out_mat <- get("leave_out_mat",envir = leave_out_env)
    leave_out_inds <- leave_out_mat[sim_num,]
    leave_out_inds <- leave_out_inds[!is.na(leave_out_inds)]

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
    id <- id[-leave_out_inds,]
    surv_event <- surv_event[-leave_out_inds]
    surv_time <- surv_time[-leave_out_inds]
    sweights_vec <- id$sweights/NHANES_NUM_WAVES
    sweights_vec <- sweights_vec/mean(sweights_vec)

    leave_out_data <- list(leave_out_inds = leave_out_inds,
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
