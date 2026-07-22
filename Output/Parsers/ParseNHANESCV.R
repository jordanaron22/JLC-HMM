#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

source(file.path("Scripting","R","constants.R"))
source(file.path("Scripting","R","diagnostics.R"))

input_dir <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Output","Routputs","RoutputsNHANESCV")
}

output_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Output","parse_nhanes_cv.rds")
}

manifest_file <- if (length(args) >= 3){
  args[[3]]
} else {
  candidates <- c("expected_jobs_nhanes.tsv",
                  file.path("Output","expected_jobs_nhanes.tsv"))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) existing[[1]] else NA_character_
}

SAVED_SECTION_SLOT <- c(true_params = 1, est_params = 2, bic = 3,
                        leave_out = 4, simulated_hmm = 5,
                        diagnostics = 6, settings = 7)

LEAVE_OUT_SLOT <- c(leave_out_inds = 1, conf_mat_list = 2,
                    cindex_new_list = 3, ibs_new_list = 4,
                    senspec_list = 5, ibs2_new_list = 6,
                    senspec_mix_list = 7)

CV_STANDARD_SCENARIO <- 6L

regex_value <- function(pattern,text,default = NA_character_){
  match <- regexec(pattern,text,perl = TRUE)
  found <- regmatches(text,match)[[1]]
  if (length(found) < 2){
    return(default)
  }
  found[[2]]
}

scalar_or_na <- function(x){
  if (is.null(x) || length(x) == 0){
    return(NA)
  }
  x[[1]]
}

numeric_or_na <- function(x){
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0){
    return(NA_real_)
  }
  out[[1]]
}

finite_or_na <- function(x){
  out <- numeric_or_na(x)
  if (!is.finite(out)){
    return(NA_real_)
  }
  out
}

as_true <- function(x){
  if (is.null(x) || length(x) == 0){
    return(FALSE)
  }
  normalized <- tolower(as.character(x[[1]]))
  normalized %in% c("true","t","1","yes")
}

mean_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  mean(x)
}

median_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  median(x)
}

sd_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2){
    return(NA_real_)
  }
  sd(x)
}

min_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  min(x)
}

max_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  max(x)
}

q95_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  unname(quantile(x,0.95,names = FALSE,type = 7))
}

collapse_values <- function(x){
  if (is.null(x) || length(x) == 0){
    return(NA_character_)
  }
  if (all(is.na(x))){
    return(NA_character_)
  }
  paste(x,collapse = ",")
}

rbind_fill <- function(data_list){
  data_list <- Filter(function(x) !is.null(x) && nrow(x) > 0,data_list)
  if (length(data_list) == 0){
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(data_list,names),use.names = FALSE))
  aligned <- lapply(data_list,function(current){
    missing_names <- setdiff(all_names,names(current))
    for (missing_name in missing_names){
      current[[missing_name]] <- NA
    }
    current[,all_names,drop = FALSE]
  })

  out <- do.call(rbind,aligned)
  row.names(out) <- NULL
  out
}

get_section <- function(to_save,name,slot,required = FALSE){
  if (!is.null(to_save) && !is.null(names(to_save)) && name %in% names(to_save)){
    return(to_save[[name]])
  }
  if (!is.null(to_save) && length(to_save) >= slot){
    return(to_save[[slot]])
  }
  if (required){
    stop("Missing saved section: ",name)
  }
  NULL
}

get_setting <- function(settings,name,default = NA){
  if (!is.null(settings) && !is.null(settings[[name]]) &&
      length(settings[[name]]) > 0){
    return(settings[[name]][[1]])
  }
  default
}

get_saved_scalar <- function(x,name,default = NA){
  if (!is.null(x) && !is.null(x[[name]]) && length(x[[name]]) > 0){
    return(x[[name]][[1]])
  }
  default
}

get_leave_out_field <- function(leave_out,name,slot,default = NULL){
  if (!is.null(leave_out) &&
      !is.null(names(leave_out)) &&
      name %in% names(leave_out)){
    return(leave_out[[name]])
  }
  if (!is.null(leave_out) && length(leave_out) >= slot){
    return(leave_out[[slot]])
  }
  default
}

load_to_save <- function(file){
  if (grepl("[.]rds$",file,ignore.case = TRUE)){
    return(readRDS(file))
  }

  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file,envir = load_env)
  if ("to_save" %in% loaded_names && exists("to_save",envir = load_env)){
    return(get("to_save",envir = load_env))
  }
  if (length(loaded_names) == 1){
    return(get(loaded_names[[1]],envir = load_env))
  }
  stop("File does not contain object named to_save")
}

model_type_from_file <- function(file_name){
  if (grepl("NoSurv",file_name,fixed = TRUE)){
    return("two_stage")
  }
  "joint"
}

scenario_value <- function(x,scenario_index = CV_STANDARD_SCENARIO){
  if (is.null(x) || length(x) < scenario_index){
    return(NA_real_)
  }
  finite_or_na(x[[scenario_index]])
}

cv_loglik_fields <- function(leave_out,diagnostic_name){
  empty <- list(
    weighted_mean = NA_real_,
    weighted_sum = NA_real_,
    weight_sum = NA_real_,
    individual_count = NA_integer_
  )

  if (is.null(leave_out) ||
      is.null(leave_out[[diagnostic_name]]) ||
      !is.list(leave_out[[diagnostic_name]])){
    return(empty)
  }

  diagnostic <- leave_out[[diagnostic_name]]
  empty$weighted_mean <- finite_or_na(diagnostic[["weighted_mean"]])
  empty$weighted_sum <- finite_or_na(diagnostic[["weighted_sum"]])
  empty$weight_sum <- finite_or_na(diagnostic[["weight_sum"]])
  if (!is.null(diagnostic[["individual_log_score"]])){
    empty$individual_count <- length(diagnostic[["individual_log_score"]])
  }
  empty
}

interval_survival_counts <- function(leave_out){
  empty <- list(
    heldout_event_count = NA_integer_,
    heldout_censored_count = NA_integer_,
    heldout_event_interval_counts = NA_character_
  )

  diagnostic <- leave_out[["cv_interval_survival_loglik"]]
  if (is.null(diagnostic) || !is.list(diagnostic)){
    return(empty)
  }

  surv_event <- diagnostic[["surv_event"]]
  event_interval_index <- diagnostic[["event_interval_index"]]
  interval_breaks <- diagnostic[["interval_breaks"]]
  if (is.null(surv_event) || length(surv_event) == 0){
    return(empty)
  }

  empty$heldout_event_count <- sum(surv_event == 1,na.rm = TRUE)
  empty$heldout_censored_count <- sum(surv_event == 0,na.rm = TRUE)

  if (!is.null(event_interval_index) &&
      !is.null(interval_breaks) &&
      length(interval_breaks) >= 2){
    counts <- tabulate(
      event_interval_index[surv_event == 1],
      nbins = length(interval_breaks) - 1L
    )
    labels <- paste0(
      interval_breaks[-length(interval_breaks)],
      "-",
      interval_breaks[-1],
      "=",
      counts
    )
    empty$heldout_event_interval_counts <- paste(labels,collapse = ";")
  }

  empty
}

extract_cv_cindex_data <- function(leave_out,file,file_name,fit_mix_num,
                                   model_type,cv_fold_id){
  if (is.null(leave_out[["cv_cindex_data"]])){
    return(data.frame())
  }

  cindex_data <- as.data.frame(leave_out[["cv_cindex_data"]])
  if (nrow(cindex_data) == 0){
    return(data.frame())
  }

  required_cols <- c("surv_time","surv_event","risk_score","sweights_vec")
  missing_cols <- setdiff(required_cols,names(cindex_data))
  if (length(missing_cols) > 0){
    stop("cv_cindex_data is missing required columns: ",
         paste(missing_cols,collapse = ", "))
  }

  if (!"cv_fold_id" %in% names(cindex_data)){
    cindex_data$cv_fold_id <- cv_fold_id
  }

  cindex_data$file <- normalizePath(file,mustWork = TRUE)
  cindex_data$file_name <- file_name
  cindex_data$fit_mix_num <- fit_mix_num
  cindex_data$model_type <- model_type
  cindex_data$surv_time <- suppressWarnings(
    as.numeric(cindex_data$surv_time)
  )
  cindex_data$surv_event <- suppressWarnings(
    as.numeric(cindex_data$surv_event)
  )
  cindex_data$risk_score <- suppressWarnings(
    as.numeric(cindex_data$risk_score)
  )
  cindex_data$sweights_vec <- suppressWarnings(
    as.numeric(cindex_data$sweights_vec)
  )
  cindex_data$cv_fold_id <- suppressWarnings(
    as.numeric(cindex_data$cv_fold_id)
  )

  cindex_data
}

extract_cv_survival_predictions <- function(leave_out,file_name,fit_mix_num,
                                            model_type){
  prediction <- leave_out[["cv_survival_predictions"]]
  if (is.null(prediction) || !is.list(prediction)){
    stop("Saved object does not contain cv_survival_predictions")
  }

  required_fields <- c(
    "cv_fold_id",
    "leave_out_ind",
    "surv_time",
    "surv_event",
    "sweights_vec",
    "eval_times",
    "surv_prob"
  )
  missing_fields <- setdiff(required_fields,names(prediction))
  if (length(missing_fields) > 0){
    stop("cv_survival_predictions is missing fields: ",
         paste(missing_fields,collapse = ", "))
  }

  leave_out_ind <- suppressWarnings(
    as.integer(prediction[["leave_out_ind"]])
  )
  surv_time <- suppressWarnings(as.numeric(prediction[["surv_time"]]))
  surv_event <- suppressWarnings(as.numeric(prediction[["surv_event"]]))
  sweights_vec <- suppressWarnings(as.numeric(prediction[["sweights_vec"]]))
  eval_times <- suppressWarnings(as.numeric(prediction[["eval_times"]]))
  surv_prob <- as.matrix(prediction[["surv_prob"]])
  cv_fold_id <- numeric_or_na(prediction[["cv_fold_id"]])

  participant_count <- length(leave_out_ind)
  if (participant_count == 0){
    stop("cv_survival_predictions contains no held-out participants")
  }

  if (length(surv_time) != participant_count ||
      length(surv_event) != participant_count ||
      length(sweights_vec) != participant_count){
    stop(
      paste(
        "leave_out_ind, surv_time, surv_event, and sweights_vec",
        "must have equal lengths in cv_survival_predictions"
      )
    )
  }

  if (is.null(dim(surv_prob)) || length(dim(surv_prob)) != 2 ||
      nrow(surv_prob) != length(eval_times) ||
      ncol(surv_prob) != participant_count){
    stop(
      paste(
        "cv_survival_predictions$surv_prob must have rows equal",
        "to eval_times and columns equal to held-out participants"
      )
    )
  }

  if (!is.finite(cv_fold_id)){
    stop("cv_survival_predictions$cv_fold_id must be finite")
  }
  if (any(is.na(leave_out_ind)) || any(leave_out_ind < 1)){
    stop("cv_survival_predictions$leave_out_ind must be positive integers")
  }
  if (any(!is.finite(surv_time)) || any(surv_time < 0)){
    stop("cv_survival_predictions$surv_time must be finite and nonnegative")
  }
  if (any(is.na(surv_event)) || any(!surv_event %in% c(0,1))){
    stop("cv_survival_predictions$surv_event must contain only 0 and 1")
  }
  if (any(!is.finite(sweights_vec)) ||
      any(sweights_vec < 0) ||
      sum(sweights_vec) <= 0){
    stop(
      paste(
        "cv_survival_predictions$sweights_vec must be finite,",
        "nonnegative, and have positive total weight"
      )
    )
  }
  if (length(eval_times) < 2 ||
      any(!is.finite(eval_times)) ||
      any(diff(eval_times) <= 0)){
    stop("cv_survival_predictions$eval_times must be strictly increasing")
  }
  if (any(!is.finite(surv_prob)) ||
      any(surv_prob < -1e-12) ||
      any(surv_prob > 1 + 1e-12)){
    stop("cv_survival_predictions$surv_prob contains invalid probabilities")
  }

  list(
    cv_fold_id = cv_fold_id,
    leave_out_ind = leave_out_ind,
    surv_time = surv_time,
    surv_event = surv_event,
    sweights_vec = sweights_vec,
    eval_times = eval_times,
    surv_prob = pmin(pmax(surv_prob,0),1),
    fit_mix_num = fit_mix_num,
    model_type = model_type,
    file_name = file_name
  )
}

parse_one_file <- function(file){
  to_save <- load_to_save(file)
  file_name <- basename(file)

  settings <- get_section(to_save,"settings",
                          SAVED_SECTION_SLOT[["settings"]])
  diagnostics <- get_section(to_save,"diagnostics",
                             SAVED_SECTION_SLOT[["diagnostics"]])
  leave_out <- get_section(to_save,"leave_out",
                           SAVED_SECTION_SLOT[["leave_out"]])

  if (is.null(leave_out) || length(leave_out) == 0){
    stop("Saved object does not contain leave_out results")
  }

  convergence <- diagnostics[["convergence"]]
  leave_out_inds <- get_leave_out_field(
    leave_out,"leave_out_inds",LEAVE_OUT_SLOT[["leave_out_inds"]]
  )
  cindex_new_list <- get_leave_out_field(
    leave_out,"cindex_new_list",LEAVE_OUT_SLOT[["cindex_new_list"]]
  )
  ibs_new_list <- get_leave_out_field(
    leave_out,"ibs_new_list",LEAVE_OUT_SLOT[["ibs_new_list"]]
  )
  ibs2_new_list <- get_leave_out_field(
    leave_out,"ibs2_new_list",LEAVE_OUT_SLOT[["ibs2_new_list"]]
  )

  cv_longitudinal_loglik <- cv_loglik_fields(
    leave_out,"cv_longitudinal_loglik"
  )
  cv_interval_survival_loglik <- cv_loglik_fields(
    leave_out,"cv_interval_survival_loglik"
  )
  interval_counts <- interval_survival_counts(leave_out)

  cv_fold_id <- numeric_or_na(get_saved_scalar(
    leave_out,
    "cv_fold_id",
    get_setting(settings,"sim_num",regex_value("Seed([0-9]+)",file_name))
  ))

  fit_mix_num <- numeric_or_na(get_setting(
    settings,"fit_mix_num",regex_value("FitMix([0-9]+)",file_name)
  ))
  model_type <- as.character(get_setting(
    settings,"model_type",model_type_from_file(file_name)
  ))

  parsed_row <- data.frame(
    file = normalizePath(file,mustWork = TRUE),
    file_name = file_name,
    fit_mix_num = fit_mix_num,
    model_type = model_type,
    data_source = as.character(get_setting(settings,"data_source","nhanes")),
    cv_fold_id = cv_fold_id,
    sim_num = numeric_or_na(get_setting(
      settings,"sim_num",regex_value("Seed([0-9]+)",file_name)
    )),
    cv_fold_count = numeric_or_na(get_saved_scalar(
      leave_out,"cv_fold_count",NA
    )),
    heldout_n = length(leave_out_inds),
    heldout_event_count = interval_counts$heldout_event_count,
    heldout_censored_count = interval_counts$heldout_censored_count,
    heldout_event_interval_counts =
      interval_counts$heldout_event_interval_counts,
    run_leave_one_out_cv = as.character(get_setting(
      settings,"run_leave_one_out_cv",NA
    )),
    class_selection_run = as.character(get_setting(
      settings,"class_selection_run",NA
    )),
    use_hot_start = as.character(get_setting(settings,"use_hot_start",NA)),
    init_jitter_scale = numeric_or_na(get_setting(
      settings,"init_jitter_scale",NA
    )),
    save_reduced_output = as.character(get_setting(
      settings,"save_reduced_output",NA
    )),
    leave_out_scenarios = collapse_values(
      leave_out[["leave_out_scenarios"]]
    ),
    state_reference_available = scalar_or_na(get_saved_scalar(
      leave_out,"state_reference_available",NA
    )),
    runtime_seconds = finite_or_na(diagnostics[["runtime_seconds"]]),
    runtime_hours = finite_or_na(diagnostics[["runtime_seconds"]]) / 3600,
    memory_peak_gb = finite_or_na(diagnostics[["memory_peak_gb"]]),
    converged = scalar_or_na(convergence[["converged"]]),
    completed_iterations = finite_or_na(
      convergence[["completed_iterations"]]
    ),
    final_likelihood = finite_or_na(convergence[["final_likelihood"]]),
    cindex_scenario6 = scenario_value(cindex_new_list),
    ibs_scenario6 = scenario_value(ibs_new_list),
    ibs2_scenario6 = scenario_value(ibs2_new_list),
    cv_longitudinal_loglik_weighted_mean =
      cv_longitudinal_loglik$weighted_mean,
    cv_longitudinal_loglik_weighted_sum =
      cv_longitudinal_loglik$weighted_sum,
    cv_longitudinal_loglik_weight_sum =
      cv_longitudinal_loglik$weight_sum,
    cv_longitudinal_loglik_individual_count =
      cv_longitudinal_loglik$individual_count,
    cv_interval_survival_loglik_weighted_mean =
      cv_interval_survival_loglik$weighted_mean,
    cv_interval_survival_loglik_weighted_sum =
      cv_interval_survival_loglik$weighted_sum,
    cv_interval_survival_loglik_weight_sum =
      cv_interval_survival_loglik$weight_sum,
    cv_interval_survival_loglik_individual_count =
      cv_interval_survival_loglik$individual_count,
    stringsAsFactors = FALSE
  )

  cv_cindex_data <- extract_cv_cindex_data(
    leave_out = leave_out,
    file = file,
    file_name = file_name,
    fit_mix_num = fit_mix_num,
    model_type = model_type,
    cv_fold_id = cv_fold_id
  )

  cv_survival_predictions <- extract_cv_survival_predictions(
    leave_out = leave_out,
    file_name = file_name,
    fit_mix_num = fit_mix_num,
    model_type = model_type
  )

  list(
    fold_summary = parsed_row,
    cv_cindex_data = cv_cindex_data,
    cv_survival_predictions = cv_survival_predictions
  )
}

pool_cv_loglik <- function(data,group_cols,score_prefix){
  weighted_sum_col <- paste0(score_prefix,"_weighted_sum")
  weight_sum_col <- paste0(score_prefix,"_weight_sum")

  if (nrow(data) == 0 ||
      !all(c(weighted_sum_col,weight_sum_col) %in% names(data))){
    return(data.frame())
  }

  valid_rows <- is.finite(data[[weighted_sum_col]]) &
    is.finite(data[[weight_sum_col]]) &
    data[[weight_sum_col]] > 0
  if (!any(valid_rows)){
    return(data.frame())
  }

  valid_data <- data[valid_rows,,drop = FALSE]
  pooled <- aggregate(
    valid_data[,c(weighted_sum_col,weight_sum_col),drop = FALSE],
    by = valid_data[,group_cols,drop = FALSE],
    FUN = sum
  )

  names(pooled)[names(pooled) == weighted_sum_col] <-
    paste0(score_prefix,"_pooled_weighted_sum")
  names(pooled)[names(pooled) == weight_sum_col] <-
    paste0(score_prefix,"_pooled_weight_sum")

  pooled[[paste0(score_prefix,"_pooled_weighted_mean")]] <-
    pooled[[paste0(score_prefix,"_pooled_weighted_sum")]] /
    pooled[[paste0(score_prefix,"_pooled_weight_sum")]]

  folds_pooled <- aggregate(
    valid_data$cv_fold_id,
    by = valid_data[,group_cols,drop = FALSE],
    FUN = function(x) length(unique(x[is.finite(x)]))
  )
  names(folds_pooled)[ncol(folds_pooled)] <-
    paste0(score_prefix,"_folds_pooled")

  merge(pooled,folds_pooled,by = group_cols,all.x = TRUE)
}

calc_pooled_cindex_value <- function(data){
  required_cols <- c("surv_time","surv_event","risk_score","sweights_vec")
  missing_cols <- setdiff(required_cols,names(data))
  if (length(missing_cols) > 0){
    stop("Pooled C-index data is missing columns: ",
         paste(missing_cols,collapse = ", "))
  }

  valid_rows <- is.finite(data$surv_time) &
    data$surv_event %in% c(0,1) &
    is.finite(data$risk_score) &
    is.finite(data$sweights_vec) &
    data$sweights_vec >= 0
  data <- data[valid_rows,,drop = FALSE]

  if (nrow(data) == 0 ||
      sum(data$sweights_vec) <= 0 ||
      sum(data$surv_event == 1) == 0 ||
      length(unique(data$risk_score)) < 2){
    return(NA_real_)
  }

  result <- survival::concordance(
    survival::Surv(surv_time, surv_event) ~ risk_score,
    data = data,
    weights = sweights_vec,
    reverse = TRUE,
    timewt = "n"
  )

  unname(result$concordance)
}

pool_oof_cindex <- function(data,group_cols,
                            output_col = "cindex_scenario6_pooled_oof"){
  if (nrow(data) == 0){
    return(data.frame())
  }

  missing_group_cols <- setdiff(group_cols,names(data))
  if (length(missing_group_cols) > 0){
    stop("cv_cindex_data is missing grouping columns: ",
         paste(missing_group_cols,collapse = ", "))
  }

  split_key <- interaction(data[,group_cols,drop = FALSE],
                           drop = TRUE,lex.order = TRUE)
  rows <- lapply(split(data,split_key),function(group_data){
    group_values <- group_data[1,group_cols,drop = FALSE]
    group_values[[output_col]] <- calc_pooled_cindex_value(group_data)
    group_values$cindex_scenario6_pooled_oof_n <- nrow(group_data)
    group_values$cindex_scenario6_pooled_oof_events <-
      sum(group_data$surv_event == 1,na.rm = TRUE)
    group_values$cindex_scenario6_pooled_oof_folds <-
      length(unique(group_data$cv_fold_id[is.finite(group_data$cv_fold_id)]))
    group_values
  })

  out <- do.call(rbind,rows)
  row.names(out) <- NULL
  out
}

pool_oof_partial_loglik <- function(data,group_cols){
  if (nrow(data) == 0){
    return(data.frame())
  }

  required_cols <- c(
    group_cols,
    "cv_fold_id",
    "leave_out_ind",
    "surv_time",
    "surv_event",
    "risk_score",
    "sweights_vec"
  )
  missing_cols <- setdiff(required_cols,names(data))
  if (length(missing_cols) > 0){
    stop("Partial likelihood data is missing columns: ",
         paste(missing_cols,collapse = ", "))
  }

  split_key <- interaction(data[,group_cols,drop = FALSE],
                           drop = TRUE,lex.order = TRUE)
  rows <- lapply(split(data,split_key),function(group_data){
    group_values <- group_data[1,group_cols,drop = FALSE]
    group_label <- paste0(
      "fit_mix_num=",group_values$fit_mix_num,
      ", model_type=",group_values$model_type,
      ": "
    )

    fold_ids <- sort(unique(as.integer(
      group_data$cv_fold_id[is.finite(group_data$cv_fold_id)]
    )))
    if (!identical(fold_ids,seq_len(NHANES_CV_FOLD_COUNT))){
      stop(
        group_label,
        "expected folds 1:",
        NHANES_CV_FOLD_COUNT,
        " but found folds ",
        paste(fold_ids,collapse = ", ")
      )
    }

    leave_out_ind <- suppressWarnings(as.integer(group_data$leave_out_ind))
    if (anyNA(leave_out_ind) || any(leave_out_ind < 1)){
      stop(group_label,"leave_out_ind must contain positive integers")
    }

    duplicated_indices <- unique(leave_out_ind[duplicated(leave_out_ind)])
    if (length(duplicated_indices) > 0){
      stop(
        group_label,
        "participants appear in more than one fold: ",
        paste(head(duplicated_indices,10),collapse = ", ")
      )
    }

    expected_indices <- seq_len(max(leave_out_ind))
    if (!identical(sort(leave_out_ind),expected_indices)){
      missing_indices <- setdiff(expected_indices,leave_out_ind)
      stop(
        group_label,
        "at least one participant is missing from pooled folds: ",
        paste(head(missing_indices,10),collapse = ", ")
      )
    }

    ordered_data <- group_data[order(leave_out_ind),,drop = FALSE]
    pooled_result <- CalcSurveyWeightedPartialLogLik(
      surv_time = ordered_data$surv_time,
      surv_event = ordered_data$surv_event,
      risk_score = ordered_data$risk_score,
      sweights_vec = ordered_data$sweights_vec
    )

    group_values$cv_partial_loglik_pooled_weighted <-
      pooled_result$partial_loglik
    group_values$cv_partial_loglik_event_weighted_mean <-
      pooled_result$event_weighted_mean
    group_values$cv_partial_loglik_event_weight_sum <-
      pooled_result$event_weight_sum
    group_values$cv_partial_loglik_weight_sum <-
      pooled_result$survey_weight_sum
    group_values$cv_partial_loglik_n <-
      pooled_result$participant_count
    group_values$cv_partial_loglik_events <-
      pooled_result$event_count
    group_values$cv_partial_loglik_unique_event_times <-
      pooled_result$unique_event_time_count
    group_values$cv_partial_loglik_folds <-
      length(fold_ids)
    group_values$cv_partial_loglik_ties_method <-
      pooled_result$ties_method

    group_values
  })

  out <- do.call(rbind,rows)
  row.names(out) <- NULL
  out
}

pool_cv_weighted_ibs <- function(prediction_list){
  prediction_list <- Filter(Negate(is.null),prediction_list)
  if (length(prediction_list) == 0){
    stop(
      paste(
        "No cv_survival_predictions were parsed.",
        "Rerun the CV jobs with the current JMHMM.R before pooling IBS."
      )
    )
  }

  metadata <- data.frame(
    list_index = seq_along(prediction_list),
    fit_mix_num = vapply(
      prediction_list,
      function(x) as.numeric(x$fit_mix_num),
      numeric(1)
    ),
    model_type = vapply(
      prediction_list,
      function(x) as.character(x$model_type),
      character(1)
    ),
    stringsAsFactors = FALSE
  )

  split_key <- interaction(
    metadata[,c("fit_mix_num","model_type"),drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )

  summary_rows <- list()
  curve_rows <- list()

  for (group_indices in split(metadata$list_index,split_key)){
    group_predictions <- prediction_list[group_indices]
    fit_mix_num <- group_predictions[[1]]$fit_mix_num
    model_type <- group_predictions[[1]]$model_type
    group_label <- paste0(
      "fit_mix_num=",fit_mix_num,
      ", model_type=",model_type,
      ": "
    )

    eval_times <- group_predictions[[1]]$eval_times
    for (prediction in group_predictions){
      if (!identical(prediction$eval_times,eval_times)){
        stop(group_label,"evaluation grids differ between folds")
      }
    }

    fold_ids <- sort(unique(vapply(
      group_predictions,
      function(x) as.integer(x$cv_fold_id),
      integer(1)
    )))
    if (!identical(fold_ids,seq_len(NHANES_CV_FOLD_COUNT))){
      stop(
        group_label,
        "expected folds 1:",
        NHANES_CV_FOLD_COUNT,
        " but found folds ",
        paste(fold_ids,collapse = ", ")
      )
    }

    all_indices <- unlist(
      lapply(group_predictions,function(x) x$leave_out_ind),
      use.names = FALSE
    )
    duplicated_indices <- unique(all_indices[duplicated(all_indices)])
    if (length(duplicated_indices) > 0){
      stop(
        group_label,
        "participants appear in more than one fold: ",
        paste(head(duplicated_indices,10),collapse = ", ")
      )
    }

    expected_indices <- seq_len(max(all_indices))
    if (!identical(sort(all_indices),expected_indices)){
      missing_indices <- setdiff(expected_indices,all_indices)
      stop(
        group_label,
        "at least one participant is missing from pooled folds: ",
        paste(head(missing_indices,10),collapse = ", ")
      )
    }

    participant_count <- length(expected_indices)
    pooled_surv_time <- rep(NA_real_,participant_count)
    pooled_surv_event <- rep(NA_real_,participant_count)
    pooled_sweights <- rep(NA_real_,participant_count)
    pooled_surv_prob <- matrix(
      NA_real_,
      nrow = length(eval_times),
      ncol = participant_count
    )

    for (prediction in group_predictions){
      current_indices <- prediction$leave_out_ind
      pooled_surv_time[current_indices] <- prediction$surv_time
      pooled_surv_event[current_indices] <- prediction$surv_event
      pooled_sweights[current_indices] <- prediction$sweights_vec
      pooled_surv_prob[,current_indices] <- prediction$surv_prob
    }

    if (anyNA(pooled_surv_time) ||
        anyNA(pooled_surv_event) ||
        anyNA(pooled_sweights) ||
        anyNA(pooled_surv_prob)){
      stop(group_label,"pooled out-of-fold IBS inputs contain missing values")
    }

    pooled_result <- CalcSurveyWeightedIBS(
      surv_time = pooled_surv_time,
      surv_event = pooled_surv_event,
      surv_prob = pooled_surv_prob,
      eval_times = eval_times,
      sweights_vec = pooled_sweights
    )

    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      fit_mix_num = fit_mix_num,
      model_type = model_type,
      cv_ibs_pooled_weighted = pooled_result$ibs,
      cv_ibs_pooled_n = pooled_result$participant_count,
      cv_ibs_pooled_events = pooled_result$event_count,
      cv_ibs_pooled_weight_sum = pooled_result$survey_weight_sum,
      cv_ibs_pooled_folds = length(fold_ids),
      cv_ibs_eval_time_min = min(eval_times),
      cv_ibs_eval_time_max = max(eval_times),
      cv_ibs_min_censoring_survival =
        pooled_result$minimum_censoring_survival,
      stringsAsFactors = FALSE
    )

    curve_rows[[length(curve_rows) + 1L]] <- data.frame(
      fit_mix_num = fit_mix_num,
      model_type = model_type,
      eval_time = pooled_result$eval_times,
      brier_score = pooled_result$brier_score,
      censoring_survival =
        pooled_result$censoring_survival_at_eval_times,
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind,summary_rows)
  curve <- do.call(rbind,curve_rows)
  row.names(summary) <- NULL
  row.names(curve) <- NULL

  list(summary = summary, curve = curve)
}

summarize_metric <- function(data,group_cols,metric_name){
  if (nrow(data) == 0 || !metric_name %in% names(data)){
    return(data.frame())
  }

  metric_data <- data[,c(group_cols,metric_name),drop = FALSE]

  summaries <- list(
    mean = mean_finite,
    median = median_finite,
    sd = sd_finite,
    min = min_finite,
    max = max_finite
  )

  out <- NULL
  for (summary_name in names(summaries)){
    current <- aggregate(
      metric_data[[metric_name]],
      by = metric_data[,group_cols,drop = FALSE],
      FUN = summaries[[summary_name]]
    )
    names(current)[ncol(current)] <- paste0(metric_name,"_",summary_name)
    out <- if (is.null(out)){
      current
    } else {
      merge(out,current,by = group_cols,all = TRUE)
    }
  }

  out
}

complete_cv_rows <- function(data){
  required_cols <- c("completed_folds","expected_folds","missing_folds")
  if (!all(required_cols %in% names(data))){
    return(rep(TRUE,nrow(data)))
  }

  completed_folds <- suppressWarnings(as.numeric(data$completed_folds))
  expected_folds <- suppressWarnings(as.numeric(data$expected_folds))
  missing_folds <- suppressWarnings(as.numeric(data$missing_folds))

  is.finite(completed_folds) &
    is.finite(expected_folds) &
    is.finite(missing_folds) &
    expected_folds > 0 &
    completed_folds == expected_folds &
    missing_folds == 0
}

add_metric_rank <- function(data,value_col,direction,prefix,
                            rank_group_cols = "model_type",
                            eligible = NULL){
  data[[paste0(prefix,"_rank")]] <- NA_integer_
  data[[paste0(prefix,"_selected")]] <- FALSE

  if (!value_col %in% names(data)){
    return(data)
  }

  if (is.null(eligible)){
    eligible <- rep(TRUE,nrow(data))
  } else {
    eligible <- as.logical(eligible)
    if (length(eligible) != nrow(data)){
      stop("eligible must have one value per row")
    }
    eligible <- !is.na(eligible) & eligible
  }

  rank_group_cols <- rank_group_cols[rank_group_cols %in% names(data)]
  rank_groups <- if (length(rank_group_cols) == 0){
    list(all = seq_len(nrow(data)))
  } else {
    split(
      seq_len(nrow(data)),
      interaction(data[,rank_group_cols,drop = FALSE],
                  drop = TRUE,lex.order = TRUE)
    )
  }

  for (group_indices in rank_groups){
    values <- suppressWarnings(as.numeric(data[[value_col]][group_indices]))
    valid <- is.finite(values) & eligible[group_indices]
    if (!any(valid)){
      next
    }

    valid_indices <- group_indices[valid]
    ordered_indices <- if (direction == "max"){
      valid_indices[order(values[valid],decreasing = TRUE)]
    } else {
      valid_indices[order(values[valid],decreasing = FALSE)]
    }

    data[[paste0(prefix,"_rank")]][ordered_indices] <-
      seq_along(ordered_indices)
    data[[paste0(prefix,"_selected")]][ordered_indices[[1]]] <- TRUE
  }
  data
}

rank_cv_models <- function(summary_data){
  if (nrow(summary_data) == 0){
    return(summary_data)
  }

  summary_data$selection_group <- if ("model_type" %in% names(summary_data)){
    summary_data$model_type
  } else {
    "all"
  }
  summary_data$cindex_scenario6_rank_source <- "fold_mean"
  summary_data$eligible_for_cv_selection <- complete_cv_rows(summary_data)

  summary_data <- add_metric_rank(
    summary_data,
    "cv_longitudinal_loglik_pooled_weighted_mean",
    "max",
    "cv_longitudinal_loglik",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "cv_interval_survival_loglik_pooled_weighted_mean",
    "max",
    "cv_interval_survival_loglik",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "cv_partial_loglik_pooled_weighted",
    "max",
    "cv_partial_loglik",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "cindex_scenario6_mean",
    "max",
    "cindex_scenario6",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "cindex_scenario6_pooled_oof",
    "max",
    "cindex_scenario6_pooled_oof",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "ibs_scenario6_mean",
    "min",
    "ibs_scenario6",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "cv_ibs_pooled_weighted",
    "min",
    "cv_ibs_pooled_weighted",
    eligible = summary_data$eligible_for_cv_selection
  )
  summary_data <- add_metric_rank(
    summary_data,
    "ibs2_scenario6_mean",
    "min",
    "ibs2_scenario6",
    eligible = summary_data$eligible_for_cv_selection
  )

  summary_data[order(summary_data$model_type,summary_data$fit_mix_num),,
               drop = FALSE]
}

read_manifest <- function(path){
  if (is.na(path) || !file.exists(path)){
    return(data.frame())
  }

  manifest <- read.delim(path,stringsAsFactors = FALSE,check.names = FALSE)
  required_cols <- c("sim_num","fit_mix_num","model_type","expected_file")
  missing_cols <- setdiff(required_cols,names(manifest))
  if (length(missing_cols) > 0){
    stop("Manifest is missing required columns: ",
         paste(missing_cols,collapse = ", "))
  }

  if ("run_leave_one_out_cv" %in% names(manifest)){
    manifest <- manifest[vapply(manifest$run_leave_one_out_cv,
                                as_true,logical(1)),,drop = FALSE]
  }

  if ("run_id" %in% names(manifest) && nrow(manifest) > 0){
    latest_run_id <- sort(unique(manifest$run_id))[[length(unique(manifest$run_id))]]
    manifest <- manifest[manifest$run_id == latest_run_id,,drop = FALSE]
  }

  manifest
}

missing_from_manifest <- function(manifest,per_fold){
  if (nrow(manifest) == 0){
    return(data.frame())
  }

  manifest_key <- paste(manifest$fit_mix_num,manifest$model_type,
                        manifest$sim_num,sep = "|")
  completed_key <- if (nrow(per_fold) > 0){
    paste(per_fold$fit_mix_num,per_fold$model_type,
          per_fold$cv_fold_id,sep = "|")
  } else {
    character()
  }

  manifest$completed <- manifest_key %in% completed_key
  manifest[!manifest$completed,,drop = FALSE]
}

missing_from_observed_folds <- function(per_fold){
  if (nrow(per_fold) == 0 ||
      !"cv_fold_count" %in% names(per_fold) ||
      !any(is.finite(per_fold$cv_fold_count))){
    return(data.frame())
  }

  group_cols <- c("fit_mix_num","model_type")
  split_key <- interaction(per_fold[,group_cols,drop = FALSE],
                           drop = TRUE,lex.order = TRUE)
  rows <- lapply(split(per_fold,split_key),function(group_data){
    fold_count <- max_finite(group_data$cv_fold_count)
    if (!is.finite(fold_count)){
      return(NULL)
    }
    expected <- seq_len(as.integer(fold_count))
    observed <- unique(group_data$cv_fold_id[is.finite(group_data$cv_fold_id)])
    missing <- setdiff(expected,observed)
    if (length(missing) == 0){
      return(NULL)
    }
    data.frame(
      fit_mix_num = group_data$fit_mix_num[[1]],
      model_type = group_data$model_type[[1]],
      cv_fold_id = missing,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null),rows)
  if (length(rows) == 0){
    return(data.frame())
  }
  out <- do.call(rbind,rows)
  row.names(out) <- NULL
  out
}

if (!dir.exists(input_dir)){
  stop("Input directory does not exist: ",input_dir)
}

files <- list.files(input_dir,pattern = "[.]rda$",full.names = TRUE,
                    recursive = TRUE,ignore.case = TRUE)
files <- files[!grepl("(^|[/\\\\])Inter",files)]
files <- files[grepl("LeaveOut",basename(files),fixed = TRUE)]
files <- files[!grepl("^LeaveOutMat[.]rda$",basename(files),
                      ignore.case = TRUE)]

rows <- list()
cindex_rows <- list()
survival_prediction_rows <- list()
errors <- list()
for (file in files){
  parsed <- tryCatch(
    parse_one_file(file),
    error = function(e){
      errors[[length(errors) + 1L]] <<- data.frame(
        file = file,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )
  if (!is.null(parsed)){
    rows[[length(rows) + 1L]] <- parsed$fold_summary
    cindex_data <- parsed$cv_cindex_data
    if (!is.null(cindex_data) && nrow(cindex_data) > 0){
      cindex_rows[[length(cindex_rows) + 1L]] <- cindex_data
    }
    if (!is.null(parsed$cv_survival_predictions)){
      survival_prediction_rows[[length(survival_prediction_rows) + 1L]] <-
        parsed$cv_survival_predictions
    }
  }
}

per_fold <- if (length(rows) > 0){
  do.call(rbind,rows)
} else {
  data.frame()
}
row.names(per_fold) <- NULL

cv_cindex_data <- if (length(cindex_rows) > 0){
  rbind_fill(cindex_rows)
} else {
  data.frame()
}
row.names(cv_cindex_data) <- NULL

pooled_ibs_result <- pool_cv_weighted_ibs(survival_prediction_rows)
cv_ibs_summary <- pooled_ibs_result$summary
cv_ibs_curve <- pooled_ibs_result$curve
cv_partial_loglik_summary <- pool_oof_partial_loglik(
  cv_cindex_data,
  c("fit_mix_num","model_type")
)

group_cols <- c("fit_mix_num","model_type")
if (nrow(per_fold) > 0){
  cv_summary_by_model <- aggregate(
    cv_fold_id ~ fit_mix_num + model_type,
    data = per_fold,
    FUN = function(x) length(unique(x[is.finite(x)]))
  )
  names(cv_summary_by_model)[names(cv_summary_by_model) == "cv_fold_id"] <-
    "completed_folds"

  heldout_summary <- aggregate(
    heldout_n ~ fit_mix_num + model_type,
    data = per_fold,
    FUN = sum
  )
  names(heldout_summary)[names(heldout_summary) == "heldout_n"] <-
    "total_heldout_n"
  cv_summary_by_model <- merge(cv_summary_by_model,heldout_summary,
                               by = group_cols,all.x = TRUE)

  expected_fold_summary <- aggregate(
    cv_fold_count ~ fit_mix_num + model_type,
    data = per_fold,
    FUN = max_finite
  )
  names(expected_fold_summary)[
    names(expected_fold_summary) == "cv_fold_count"
  ] <- "expected_folds"
  cv_summary_by_model <- merge(cv_summary_by_model,expected_fold_summary,
                               by = group_cols,all.x = TRUE)

  convergence_summary <- aggregate(
    converged ~ fit_mix_num + model_type,
    data = per_fold,
    FUN = function(x) sum(vapply(x,as_true,logical(1)))
  )
  names(convergence_summary)[names(convergence_summary) == "converged"] <-
    "converged_folds"
  cv_summary_by_model <- merge(cv_summary_by_model,convergence_summary,
                               by = group_cols,all.x = TRUE)

  numeric_summaries <- list(
    mean_runtime_seconds = aggregate(runtime_seconds ~ fit_mix_num + model_type,
                                     data = per_fold,FUN = mean_finite),
    q95_runtime_seconds = aggregate(runtime_seconds ~ fit_mix_num + model_type,
                                    data = per_fold,FUN = q95_finite),
    max_runtime_seconds = aggregate(runtime_seconds ~ fit_mix_num + model_type,
                                    data = per_fold,FUN = max_finite),
    mean_memory_peak_gb = aggregate(memory_peak_gb ~ fit_mix_num + model_type,
                                    data = per_fold,FUN = mean_finite),
    q95_memory_peak_gb = aggregate(memory_peak_gb ~ fit_mix_num + model_type,
                                   data = per_fold,FUN = q95_finite),
    max_memory_peak_gb = aggregate(memory_peak_gb ~ fit_mix_num + model_type,
                                   data = per_fold,FUN = max_finite),
    mean_completed_iterations = aggregate(
      completed_iterations ~ fit_mix_num + model_type,
      data = per_fold,FUN = mean_finite
    )
  )

  for (summary_name in names(numeric_summaries)){
    summary_data <- numeric_summaries[[summary_name]]
    names(summary_data)[ncol(summary_data)] <- summary_name
    cv_summary_by_model <- merge(cv_summary_by_model,summary_data,
                                 by = group_cols,all.x = TRUE)
  }

  for (metric_name in c("cindex_scenario6","ibs_scenario6",
                        "ibs2_scenario6")){
    metric_summary <- summarize_metric(per_fold,group_cols,metric_name)
    if (nrow(metric_summary) > 0){
      cv_summary_by_model <- merge(cv_summary_by_model,metric_summary,
                                   by = group_cols,all.x = TRUE)
    }
  }

  cv_cindex_summary <- pool_oof_cindex(cv_cindex_data,group_cols)
  if (nrow(cv_cindex_summary) > 0){
    cv_summary_by_model <- merge(cv_summary_by_model,cv_cindex_summary,
                                 by = group_cols,all.x = TRUE)
  }

  if (nrow(cv_ibs_summary) > 0){
    cv_summary_by_model <- merge(cv_summary_by_model,cv_ibs_summary,
                                 by = group_cols,all.x = TRUE)
  }

  if (nrow(cv_partial_loglik_summary) > 0){
    cv_summary_by_model <- merge(cv_summary_by_model,
                                 cv_partial_loglik_summary,
                                 by = group_cols,all.x = TRUE)
  }

  cv_longitudinal_loglik_summary <- pool_cv_loglik(
    per_fold,group_cols,"cv_longitudinal_loglik"
  )
  if (nrow(cv_longitudinal_loglik_summary) > 0){
    cv_summary_by_model <- merge(cv_summary_by_model,
                                 cv_longitudinal_loglik_summary,
                                 by = group_cols,all.x = TRUE)
  }

  cv_interval_survival_loglik_summary <- pool_cv_loglik(
    per_fold,group_cols,"cv_interval_survival_loglik"
  )
  if (nrow(cv_interval_survival_loglik_summary) > 0){
    cv_summary_by_model <- merge(cv_summary_by_model,
                                 cv_interval_survival_loglik_summary,
                                 by = group_cols,all.x = TRUE)
  }

  cv_summary_by_model$missing_folds <-
    cv_summary_by_model$expected_folds -
    cv_summary_by_model$completed_folds
  cv_summary_by_model$mean_runtime_hours <-
    cv_summary_by_model$mean_runtime_seconds / 3600
  cv_summary_by_model$q95_runtime_hours <-
    cv_summary_by_model$q95_runtime_seconds / 3600
  cv_summary_by_model$max_runtime_hours <-
    cv_summary_by_model$max_runtime_seconds / 3600
} else {
  cv_summary_by_model <- data.frame()
  cv_cindex_summary <- data.frame()
  cv_ibs_summary <- data.frame()
  cv_ibs_curve <- data.frame()
  cv_partial_loglik_summary <- data.frame()
  cv_longitudinal_loglik_summary <- data.frame()
  cv_interval_survival_loglik_summary <- data.frame()
}

manifest <- read_manifest(manifest_file)
missing_cv_jobs <- missing_from_manifest(manifest,per_fold)
if (nrow(missing_cv_jobs) == 0){
  missing_cv_jobs <- missing_from_observed_folds(per_fold)
}

cv_model_selection <- rank_cv_models(cv_summary_by_model)

out <- list(
  cv_per_fold = per_fold,
  cv_cindex_data = cv_cindex_data,
  cv_summary_by_model = cv_summary_by_model,
  cv_model_selection = cv_model_selection,
  cv_cindex_summary = cv_cindex_summary,
  cv_ibs_summary = cv_ibs_summary,
  cv_ibs_curve = cv_ibs_curve,
  cv_partial_loglik_summary = cv_partial_loglik_summary,
  cv_longitudinal_loglik_summary = cv_longitudinal_loglik_summary,
  cv_interval_survival_loglik_summary =
    cv_interval_survival_loglik_summary,
  manifest = manifest,
  missing_cv_jobs = missing_cv_jobs,
  parse_errors = if (length(errors) > 0) do.call(rbind,errors) else data.frame(),
  input_dir = normalizePath(input_dir,mustWork = TRUE),
  manifest_file = if (!is.na(manifest_file) && file.exists(manifest_file)){
    normalizePath(manifest_file,mustWork = TRUE)
  } else {
    NA_character_
  }
)

dir.create(dirname(output_file),recursive = TRUE,showWarnings = FALSE)
saveRDS(out,output_file)

message("Read ",length(files)," leave-out files from: ",input_dir)
message("Saved parsed NHANES CV data to: ",output_file)

if (length(errors) > 0){
  message("Parse errors: ",length(errors))
}
if (nrow(missing_cv_jobs) > 0){
  message("Missing CV jobs: ",nrow(missing_cv_jobs))
}
