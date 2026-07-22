#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Output","Routputs","RoutputsNHANES")
}

output_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Output","parse_nhanes_runtime.rds")
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
EST_PARAM_SLOT <- c(re_prob = 10)
MIN_CLASS_PROP_THRESHOLD <- 0.03

get_optional_class_size_weights <- function(){
  if (!exists("sweights_vec",envir = .GlobalEnv,inherits = FALSE)){
    return(NULL)
  }

  weights_vec <- suppressWarnings(
    as.numeric(get("sweights_vec",envir = .GlobalEnv,inherits = FALSE))
  )
  if (length(weights_vec) == 0 || any(!is.finite(weights_vec))){
    return(NULL)
  }
  weights_vec
}

CLASS_SIZE_WEIGHTS <- get_optional_class_size_weights()

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

mean_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  mean(x)
}

max_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  max(x)
}

min_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  min(x)
}

q95_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  unname(quantile(x,0.95,names = FALSE,type = 7))
}

select_best_rows <- function(data,group_cols,value_col,direction = "min"){
  if (nrow(data) == 0 || !value_col %in% names(data)){
    return(data.frame())
  }

  split_key <- interaction(data[,group_cols,drop = FALSE],
                           drop = TRUE,lex.order = TRUE)
  selected_rows <- lapply(split(data,split_key),function(group_data){
    values <- suppressWarnings(as.numeric(group_data[[value_col]]))
    valid <- is.finite(values)
    if (!any(valid)){
      return(NULL)
    }

    valid_indices <- which(valid)
    best_offset <- if (direction == "max"){
      which.max(values[valid])
    } else {
      which.min(values[valid])
    }
    group_data[valid_indices[[best_offset]],,drop = FALSE]
  })
  selected_rows <- Filter(Negate(is.null),selected_rows)

  if (length(selected_rows) == 0){
    return(data.frame())
  }

  out <- do.call(rbind,selected_rows)
  row.names(out) <- NULL
  out
}

add_bic_ranking <- function(model_selection){
  if (nrow(model_selection) == 0 || !"bic" %in% names(model_selection)){
    return(model_selection)
  }

  bic <- suppressWarnings(as.numeric(model_selection$bic))
  valid <- is.finite(bic)
  model_selection$bic_rank <- NA_integer_
  model_selection$delta_bic <- NA_real_
  model_selection$selected_by_bic <- FALSE

  if (!any(valid)){
    return(model_selection)
  }

  valid_order <- order(bic[valid])
  valid_indices <- which(valid)
  ranked_indices <- valid_indices[valid_order]
  model_selection$bic_rank[ranked_indices] <- seq_along(ranked_indices)
  model_selection$delta_bic[valid] <- bic[valid] - min(bic[valid])
  model_selection$selected_by_bic[ranked_indices[[1]]] <- TRUE
  model_selection[order(model_selection$bic_rank,na.last = TRUE),,
                  drop = FALSE]
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

get_param <- function(params,name,slot,default = NULL){
  if (!is.null(params) && !is.null(names(params)) && name %in% names(params)){
    return(params[[name]])
  }
  if (!is.null(params) && length(params) >= slot){
    return(params[[slot]])
  }
  default
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

class_size_summary <- function(re_prob,
                               weights_vec = CLASS_SIZE_WEIGHTS,
                               min_prop_threshold = MIN_CLASS_PROP_THRESHOLD){
  empty_summary <- list(
    min_class_size = NA_real_,
    min_class_prop = NA_real_,
    passes_min_class_prop = NA,
    class_size_weighted = FALSE,
    class_size_weight_source = "unweighted"
  )

  if (is.null(re_prob) || length(re_prob) == 0){
    return(empty_summary)
  }

  re_prob <- as.matrix(re_prob)
  storage.mode(re_prob) <- "numeric"
  if (nrow(re_prob) == 0 || ncol(re_prob) == 0){
    return(empty_summary)
  }

  weights_used <- FALSE
  weight_source <- "unweighted"
  if (!is.null(weights_vec)){
    weights_vec <- suppressWarnings(as.numeric(weights_vec))
    if (length(weights_vec) == nrow(re_prob) && all(is.finite(weights_vec))){
      weights_used <- TRUE
      weight_source <- "sweights_vec"
    } else {
      weight_source <- "unweighted_weight_length_mismatch"
    }
  }

  if (weights_used){
    class_size <- colSums(
      re_prob *
        matrix(weights_vec,nrow = nrow(re_prob),ncol = ncol(re_prob))
    )
  } else {
    class_size <- colSums(re_prob,na.rm = TRUE)
  }
  total_size <- sum(class_size)
  if (!is.finite(total_size) || total_size <= 0){
    empty_summary$class_size_weighted <- weights_used
    empty_summary$class_size_weight_source <- weight_source
    return(empty_summary)
  }

  class_prop <- class_size / total_size
  min_class_prop <- min(class_prop)
  list(
    min_class_size = min(class_size),
    min_class_prop = min_class_prop,
    passes_min_class_prop = is.finite(min_class_prop) &&
      min_class_prop >= min_prop_threshold,
    class_size_weighted = weights_used,
    class_size_weight_source = weight_source
  )
}

model_type_from_file <- function(file_name){
  if (grepl("NoSurv",file_name,fixed = TRUE)){
    return("two_stage")
  }
  "joint"
}

load_to_save <- function(file){
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

parse_one_file <- function(file){
  to_save <- load_to_save(file)
  file_name <- basename(file)
  settings <- get_section(to_save,"settings",
                          SAVED_SECTION_SLOT[["settings"]])
  est_params <- get_section(to_save,"est_params",
                            SAVED_SECTION_SLOT[["est_params"]])
  leave_out <- get_section(to_save,"leave_out",
                           SAVED_SECTION_SLOT[["leave_out"]])
  diagnostics <- get_section(to_save,"diagnostics",
                             SAVED_SECTION_SLOT[["diagnostics"]])
  convergence <- diagnostics[["convergence"]]

  runtime_seconds <- finite_or_na(diagnostics[["runtime_seconds"]])
  memory_peak_gb <- finite_or_na(diagnostics[["memory_peak_gb"]])
  bic <- finite_or_na(to_save[["bic"]])
  aic <- finite_or_na(to_save[["aic"]])
  class_summary <- class_size_summary(
    get_param(est_params,"re_prob",EST_PARAM_SLOT[["re_prob"]])
  )
  cv_longitudinal_loglik <- cv_loglik_fields(
    leave_out,"cv_longitudinal_loglik"
  )
  cv_interval_survival_loglik <- cv_loglik_fields(
    leave_out,"cv_interval_survival_loglik"
  )

  data.frame(
    file = normalizePath(file,mustWork = TRUE),
    file_name = file_name,
    sim_num = numeric_or_na(get_setting(
      settings,"sim_num",regex_value("Seed([0-9]+)",file_name)
    )),
    fit_mix_num = numeric_or_na(get_setting(
      settings,"fit_mix_num",regex_value("FitMix([0-9]+)",file_name)
    )),
    model_type = as.character(get_setting(
      settings,"model_type",model_type_from_file(file_name)
    )),
    data_source = as.character(get_setting(settings,"data_source","nhanes")),
    run_bootstrap = as.character(get_setting(settings,"run_bootstrap",NA)),
    init_jitter_scale = numeric_or_na(get_setting(settings,
                                                  "init_jitter_scale",NA)),
    run_leave_one_out_cv = as.character(get_setting(
      settings,"run_leave_one_out_cv",NA
    )),
    use_hot_start = as.character(get_setting(settings,"use_hot_start",NA)),
    save_reduced_output = as.character(get_setting(
      settings,"save_reduced_output",NA
    )),
    class_selection_run = as.character(get_setting(
      settings,"class_selection_run",NA
    )),
    period_len = numeric_or_na(get_setting(settings,"period_len",
                                           regex_value("len([0-9]+)",file_name))),
    runtime_seconds = runtime_seconds,
    runtime_hours = runtime_seconds / 3600,
    memory_peak_gb = memory_peak_gb,
    completed_iterations = finite_or_na(
      convergence[["completed_iterations"]]
    ),
    total_iteration_seconds = finite_or_na(
      convergence[["total_iteration_seconds"]]
    ),
    final_likelihood = finite_or_na(convergence[["final_likelihood"]]),
    bic = bic,
    aic = aic,
    min_class_size = class_summary$min_class_size,
    min_class_prop = class_summary$min_class_prop,
    min_class_prop_threshold = MIN_CLASS_PROP_THRESHOLD,
    passes_min_class_prop = class_summary$passes_min_class_prop,
    class_size_weighted = class_summary$class_size_weighted,
    class_size_weight_source = class_summary$class_size_weight_source,
    cv_fold_id = numeric_or_na(get_saved_scalar(leave_out,
                                                "cv_fold_id",NA)),
    cv_fold_count = numeric_or_na(get_saved_scalar(leave_out,
                                                   "cv_fold_count",NA)),
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
    converged = scalar_or_na(convergence[["converged"]]),
    time_limit_hours = numeric_or_na(get_setting(settings,
                                                 "time_limit_hours",NA)),
    memory_limit_gb = numeric_or_na(get_setting(settings,
                                                "memory_limit_gb",NA)),
    stringsAsFactors = FALSE
  )
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
  manifest
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

  fold_source <- if ("cv_fold_id" %in% names(valid_data) &&
                     any(is.finite(valid_data$cv_fold_id))){
    "cv_fold_id"
  } else {
    "sim_num"
  }

  folds_pooled <- aggregate(
    valid_data[[fold_source]],
    by = valid_data[,group_cols,drop = FALSE],
    FUN = function(x) length(unique(x[is.finite(x)]))
  )
  names(folds_pooled)[ncol(folds_pooled)] <-
    paste0(score_prefix,"_folds_pooled")

  merge(pooled,folds_pooled,by = group_cols,all.x = TRUE)
}

if (!dir.exists(input_dir)){
  stop("Input directory does not exist: ",input_dir)
}

files <- list.files(input_dir,pattern = "[.]rda$",full.names = TRUE,
                    recursive = TRUE,ignore.case = TRUE)
files <- files[!grepl("(^|[/\\\\])Inter",files)]

rows <- list()
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
    rows[[length(rows) + 1L]] <- parsed
  }
}

per_run <- if (length(rows) > 0){
  do.call(rbind,rows)
} else {
  data.frame()
}
row.names(per_run) <- NULL

scenario_cols <- c("fit_mix_num","model_type")
eligible_per_run <- if (nrow(per_run) > 0 &&
                        "passes_min_class_prop" %in% names(per_run)){
  per_run[per_run$passes_min_class_prop %in% TRUE,,drop = FALSE]
} else {
  data.frame()
}
cv_longitudinal_loglik_summary <- pool_cv_loglik(
  per_run,scenario_cols,"cv_longitudinal_loglik"
)
cv_interval_survival_loglik_summary <- pool_cv_loglik(
  per_run,scenario_cols,"cv_interval_survival_loglik"
)

if (nrow(per_run) > 0){
  scenario_summary <- aggregate(
    sim_num ~ fit_mix_num + model_type,
    data = per_run,
    FUN = function(x) length(unique(x))
  )
  names(scenario_summary)[names(scenario_summary) == "sim_num"] <-
    "completed_runs"

  numeric_summaries <- list(
    mean_runtime_seconds = aggregate(runtime_seconds ~ fit_mix_num + model_type,
                                     data = per_run,FUN = mean_finite),
    q95_runtime_seconds = aggregate(runtime_seconds ~ fit_mix_num + model_type,
                                    data = per_run,FUN = q95_finite),
    max_runtime_seconds = aggregate(runtime_seconds ~ fit_mix_num + model_type,
                                    data = per_run,FUN = max_finite),
    mean_memory_peak_gb = aggregate(memory_peak_gb ~ fit_mix_num + model_type,
                                    data = per_run,FUN = mean_finite),
    q95_memory_peak_gb = aggregate(memory_peak_gb ~ fit_mix_num + model_type,
                                   data = per_run,FUN = q95_finite),
    max_memory_peak_gb = aggregate(memory_peak_gb ~ fit_mix_num + model_type,
                                   data = per_run,FUN = max_finite)
  )

  for (summary_name in names(numeric_summaries)){
    summary_data <- numeric_summaries[[summary_name]]
    names(summary_data)[ncol(summary_data)] <- summary_name
    scenario_summary <- merge(scenario_summary,summary_data,
                              by = scenario_cols,all.x = TRUE)
  }

  scenario_summary$mean_runtime_hours <-
    scenario_summary$mean_runtime_seconds / 3600
  scenario_summary$q95_runtime_hours <-
    scenario_summary$q95_runtime_seconds / 3600
  scenario_summary$max_runtime_hours <-
    scenario_summary$max_runtime_seconds / 3600

  if ("final_likelihood" %in% names(per_run)){
    likelihood_summary <- aggregate(
      final_likelihood ~ fit_mix_num + model_type,
      data = per_run,
      FUN = max_finite
    )
    names(likelihood_summary)[names(likelihood_summary) == "final_likelihood"] <-
      "best_final_likelihood"
    scenario_summary <- merge(scenario_summary,likelihood_summary,
                              by = scenario_cols,all.x = TRUE)
  }

  if ("bic" %in% names(per_run)){
    bic_summary <- aggregate(
      bic ~ fit_mix_num + model_type,
      data = per_run,
      FUN = min_finite
    )
    names(bic_summary)[names(bic_summary) == "bic"] <- "best_bic"
    scenario_summary <- merge(scenario_summary,bic_summary,
                              by = scenario_cols,all.x = TRUE)

    if (nrow(eligible_per_run) > 0){
      screened_bic_summary <- aggregate(
        bic ~ fit_mix_num + model_type,
        data = eligible_per_run,
        FUN = min_finite
      )
      names(screened_bic_summary)[names(screened_bic_summary) == "bic"] <-
        "best_screened_bic"
      scenario_summary <- merge(scenario_summary,screened_bic_summary,
                                by = scenario_cols,all.x = TRUE)
    } else {
      scenario_summary$best_screened_bic <- NA_real_
    }
  }

  if ("aic" %in% names(per_run)){
    aic_summary <- aggregate(
      aic ~ fit_mix_num + model_type,
      data = per_run,
      FUN = min_finite
    )
    names(aic_summary)[names(aic_summary) == "aic"] <- "best_aic"
    scenario_summary <- merge(scenario_summary,aic_summary,
                              by = scenario_cols,all.x = TRUE)
  }

  if ("passes_min_class_prop" %in% names(per_run)){
    eligible_summary <- aggregate(
      passes_min_class_prop ~ fit_mix_num + model_type,
      data = per_run,
      FUN = function(x) sum(x %in% TRUE)
    )
    names(eligible_summary)[
      names(eligible_summary) == "passes_min_class_prop"
    ] <- "eligible_runs"
    scenario_summary <- merge(scenario_summary,eligible_summary,
                              by = scenario_cols,all.x = TRUE)
  }

  if (nrow(cv_longitudinal_loglik_summary) > 0){
    scenario_summary <- merge(
      scenario_summary,
      cv_longitudinal_loglik_summary,
      by = scenario_cols,
      all.x = TRUE
    )
  }

  if (nrow(cv_interval_survival_loglik_summary) > 0){
    scenario_summary <- merge(
      scenario_summary,
      cv_interval_survival_loglik_summary,
      by = scenario_cols,
      all.x = TRUE
    )
  }
} else {
  scenario_summary <- data.frame()
}

manifest <- read_manifest(manifest_file)
missing_jobs <- data.frame()

if (nrow(manifest) > 0){
  manifest_key <- paste(manifest$fit_mix_num,manifest$model_type,
                        manifest$sim_num,sep = "|")
  completed_key <- if (nrow(per_run) > 0){
    paste(per_run$fit_mix_num,per_run$model_type,per_run$sim_num,sep = "|")
  } else {
    character()
  }
  manifest$completed <- manifest_key %in% completed_key
  missing_jobs <- manifest[!manifest$completed,,drop = FALSE]

  submitted_counts <- aggregate(
    sim_num ~ fit_mix_num + model_type,
    data = manifest,
    FUN = function(x) length(unique(x))
  )
  names(submitted_counts)[names(submitted_counts) == "sim_num"] <-
    "submitted_runs"

  if (nrow(scenario_summary) == 0){
    scenario_summary <- submitted_counts
    scenario_summary$completed_runs <- 0L
  } else {
    scenario_summary <- merge(scenario_summary,submitted_counts,
                              by = scenario_cols,all = TRUE)
    scenario_summary$completed_runs[is.na(scenario_summary$completed_runs)] <- 0L
  }
  scenario_summary$missing_runs <-
    scenario_summary$submitted_runs - scenario_summary$completed_runs
}

scenario_summary <- scenario_summary[order(scenario_summary$fit_mix_num,
                                           scenario_summary$model_type),,
                                     drop = FALSE]

best_bic_by_scenario <- select_best_rows(per_run,scenario_cols,"bic","min")
best_bic_screened_by_scenario <- select_best_rows(eligible_per_run,
                                                  scenario_cols,"bic","min")
best_likelihood_by_scenario <- select_best_rows(per_run,scenario_cols,
                                                "final_likelihood","max")
model_selection <- add_bic_ranking(best_bic_by_scenario)
model_selection_screened <- add_bic_ranking(best_bic_screened_by_scenario)

out <- list(
  per_run = per_run,
  scenario_summary = scenario_summary,
  best_bic_by_scenario = best_bic_by_scenario,
  best_bic_screened_by_scenario = best_bic_screened_by_scenario,
  best_likelihood_by_scenario = best_likelihood_by_scenario,
  cv_longitudinal_loglik_summary = cv_longitudinal_loglik_summary,
  cv_interval_survival_loglik_summary =
    cv_interval_survival_loglik_summary,
  model_selection = model_selection,
  model_selection_screened = model_selection_screened,
  manifest = manifest,
  missing_jobs = missing_jobs,
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


if (length(errors) > 0){
  message("Parse errors: ",length(errors))
}
if (nrow(missing_jobs) > 0){
  message("Missing submitted jobs: ",nrow(missing_jobs))
}

print(scenario_summary)
print(model_selection)
print(model_selection_screened)
