#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

find_repo_root <- function(){
  candidates <- c(
    getwd(),
    file.path(getwd(),".."),
    file.path(getwd(),"..","..")
  )
  candidates <- normalizePath(candidates,mustWork = FALSE)
  for (candidate in candidates){
    if (file.exists(file.path(candidate,"Scripting","R","constants.R"))){
      return(candidate)
    }
  }
  stop("Could not find repository root containing Scripting/R/constants.R")
}

is_absolute_path <- function(path){
  grepl("^([A-Za-z]:)?[\\/]",path)
}

resolve_path <- function(path,repo_root){
  if (is_absolute_path(path)){
    return(normalizePath(path,mustWork = FALSE))
  }
  normalizePath(file.path(repo_root,path),mustWork = FALSE)
}

repo_root <- find_repo_root()
input_dir <- resolve_path(if (length(args) >= 1) args[[1]] else "Routputs",
                          repo_root)
output_dir <- resolve_path(if (length(args) >= 2) args[[2]] else
                             file.path("Output","parsed_results"),
                           repo_root)

source(file.path(repo_root,"Scripting","R","constants.R"))
source(file.path(repo_root,"Scripting","R","saved_results.R"))

dir.create(output_dir,recursive = TRUE,showWarnings = FALSE)

empty_run_level <- function(){
  data.frame(
    file = character(),
    file_name = character(),
    sim_num = integer(),
    true_mix_num = integer(),
    fit_mix_num = integer(),
    model_type = character(),
    data_source = character(),
    sim_scenario = integer(),
    days = integer(),
    num_people = integer(),
    missing_perc = numeric(),
    emission_overlap = character(),
    init_jitter_scale = numeric(),
    class_selection_run = logical(),
    save_reduced_output = logical(),
    bic = numeric(),
    aic = numeric(),
    new_likelihood = numeric(),
    parameter_count = numeric(),
    likelihood_type = character(),
    training_ibs = numeric(),
    training_cindex = numeric(),
    test_ibs = numeric(),
    test_cindex = numeric(),
    ibs2 = numeric(),
    class_entropy_mean = numeric(),
    class_entropy_median = numeric(),
    class_entropy_q25 = numeric(),
    class_entropy_q75 = numeric(),
    class_entropy_mean_max_posterior = numeric(),
    survival_beta_comparable = logical(),
    survival_beta_rmse = numeric(),
    survival_beta_rmse_non_reference = numeric(),
    survival_beta_sse = numeric(),
    file_size = numeric(),
    file_modified_time = as.POSIXct(character()),
    stringsAsFactors = FALSE
  )
}

empty_confusion <- function(){
  data.frame(file = character(),
             sim_num = integer(),
             true_mix_num = integer(),
             fit_mix_num = integer(),
             model_type = character(),
             sim_scenario = integer(),
             days = integer(),
             emission_overlap = character(),
             true_class = character(),
             fitted_class = character(),
             n = numeric(),
             stringsAsFactors = FALSE)
}

empty_parameter_counts <- function(){
  data.frame(file = character(),
             sim_num = integer(),
             true_mix_num = integer(),
             fit_mix_num = integer(),
             model_type = character(),
             sim_scenario = integer(),
             days = integer(),
             emission_overlap = character(),
             parameter = character(),
             count = numeric(),
             stringsAsFactors = FALSE)
}

empty_survival_coefficients <- function(){
  data.frame(file = character(),
             sim_num = integer(),
             true_mix_num = integer(),
             fit_mix_num = integer(),
             model_type = character(),
             sim_scenario = integer(),
             days = integer(),
             emission_overlap = character(),
             class = integer(),
             true_beta = numeric(),
             estimated_beta = numeric(),
             error = numeric(),
             squared_error = numeric(),
             stringsAsFactors = FALSE)
}

empty_class_selection <- function(){
  data.frame(true_mix_num = integer(),
             model_type = character(),
             sim_scenario = integer(),
             days = integer(),
             num_people = integer(),
             missing_perc = numeric(),
             emission_overlap = character(),
             sim_num = integer(),
             criterion = character(),
             direction = character(),
             selected_fit_mix_num = integer(),
             selected_value = numeric(),
             correct_true_mix_num = logical(),
             num_fits_available = integer(),
             selected_file = character(),
             stringsAsFactors = FALSE)
}

empty_file_issue <- function(){
  data.frame(file = character(),
             file_name = character(),
             reason = character(),
             stringsAsFactors = FALSE)
}

bind_rows <- function(rows,empty_table){
  rows <- Filter(Negate(is.null),rows)
  if (length(rows) == 0){
    return(empty_table())
  }
  do.call(rbind,rows)
}

regex_value <- function(pattern,text,default = NA_character_){
  match <- regexec(pattern,text,perl = TRUE)
  found <- regmatches(text,match)[[1]]
  if (length(found) < 2){
    return(default)
  }
  found[[2]]
}

as_number_or_na <- function(value){
  if (is.null(value) || length(value) == 0){
    return(NA_real_)
  }
  suppressWarnings(as.numeric(value[[1]]))
}

as_character_or_na <- function(value){
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])){
    return(NA_character_)
  }
  as.character(value[[1]])
}

as_logical_or_na <- function(value){
  if (is.null(value) || length(value) == 0 || is.na(value[[1]])){
    return(NA)
  }
  if (is.logical(value)){
    return(value[[1]])
  }
  if (is.numeric(value)){
    return(value[[1]] != 0)
  }
  normalized <- tolower(as.character(value[[1]]))
  if (normalized %in% c("true","t","1","yes")){
    return(TRUE)
  }
  if (normalized %in% c("false","f","0","no")){
    return(FALSE)
  }
  NA
}

coalesce_value <- function(...){
  values <- list(...)
  for (value in values){
    if (!is.null(value) && length(value) > 0 && !is.na(value[[1]])){
      return(value[[1]])
    }
  }
  NA
}

section_safe <- function(to_save,section_name){
  tryCatch(
    get_saved_section(to_save,section_name,required = FALSE),
    error = function(e) NULL
  )
}

param_safe <- function(params,param_name){
  tryCatch(
    get_saved_param(params,param_name,default = NULL,required = FALSE),
    error = function(e) NULL
  )
}

diag_value <- function(diagnostics,name,slot = NULL,default = NA_real_){
  if (is.null(diagnostics)){
    return(default)
  }
  if (!is.null(names(diagnostics)) && name %in% names(diagnostics)){
    value <- diagnostics[[name]]
    if (length(value) > 0){
      return(value[[1]])
    }
  }
  if (!is.null(slot) && length(diagnostics) >= slot){
    value <- diagnostics[[slot]]
    if (length(value) > 0 && is.atomic(value)){
      return(value[[1]])
    }
  }
  default
}

list_value <- function(x,name,default = NA_real_){
  if (is.null(x) || is.null(names(x)) || !name %in% names(x)){
    return(default)
  }
  value <- x[[name]]
  if (length(value) == 0){
    return(default)
  }
  value[[1]]
}

extract_filename_metadata <- function(file_name){
  list(
    sim_num = as_number_or_na(regex_value("Seed([0-9]+)",file_name)),
    true_mix_num = as_number_or_na(regex_value("TrueMix([0-9]+)",file_name)),
    fit_mix_num = as_number_or_na(regex_value("FitMix([0-9]+)",file_name)),
    sim_scenario = as_number_or_na(regex_value("SimSize(-?[0-9]+)",file_name)),
    period_len = as_number_or_na(regex_value("len([0-9]+)",file_name)),
    emission_overlap = tolower(regex_value("Overlap([A-Za-z]+)",file_name)),
    model_type = if (grepl("NoSurv",file_name)) "two_stage" else "joint",
    class_selection_run = grepl("ClassSelection",file_name),
    save_reduced_output = grepl("Reduced",file_name),
    init_jitter_scale = if (grepl("RandInit",file_name)) NA_real_ else 0
  )
}

setting_or_file <- function(settings,file_meta,name,default = NA){
  setting_value <- if (!is.null(settings) && !is.null(settings[[name]])){
    settings[[name]]
  } else {
    NULL
  }
  file_value <- if (!is.null(file_meta[[name]])){
    file_meta[[name]]
  } else {
    NULL
  }
  coalesce_value(setting_value,file_value,default)
}

scenario_lookup <- function(sim_scenario){
  scenario_name <- as.character(sim_scenario)
  if (is.na(sim_scenario) || !scenario_name %in% names(SIM_SCENARIOS)){
    return(list(days = NA_integer_,num_people = NA_integer_,
                missing_perc = NA_real_))
  }
  SIM_SCENARIOS[[scenario_name]]
}

metric_from_nested <- function(diagnostics,section_name,metric_name,default){
  if (is.null(diagnostics) ||
      is.null(diagnostics[[section_name]]) ||
      is.null(diagnostics[[section_name]][[metric_name]])){
    return(default)
  }
  diagnostics[[section_name]][[metric_name]][[1]]
}

extract_confusion_table <- function(file,run_meta,diagnostics){
  confusion_table <- NULL
  if (!is.null(diagnostics)){
    if (!is.null(names(diagnostics)) &&
        "confusion_table" %in% names(diagnostics)){
      confusion_table <- diagnostics$confusion_table
    } else if (length(diagnostics) >= 3){
      confusion_table <- diagnostics[[3]]
    }
  }
  if (is.null(confusion_table) || is.null(dim(confusion_table))){
    return(NULL)
  }

  if (!inherits(confusion_table,"table")){
    confusion_table <- as.table(confusion_table)
  }
  confusion_df <- as.data.frame(confusion_table,stringsAsFactors = FALSE)
  names(confusion_df) <- c("true_class","fitted_class","n")
  data.frame(file = file,
             sim_num = run_meta$sim_num,
             true_mix_num = run_meta$true_mix_num,
             fit_mix_num = run_meta$fit_mix_num,
             model_type = run_meta$model_type,
             sim_scenario = run_meta$sim_scenario,
             days = run_meta$days,
             emission_overlap = run_meta$emission_overlap,
             confusion_df,
             stringsAsFactors = FALSE)
}

extract_parameter_counts <- function(file,run_meta,diagnostics){
  if (is.null(diagnostics) ||
      is.null(diagnostics$bic_parameter_breakdown)){
    return(NULL)
  }
  counts <- diagnostics$bic_parameter_breakdown
  data.frame(file = file,
             sim_num = run_meta$sim_num,
             true_mix_num = run_meta$true_mix_num,
             fit_mix_num = run_meta$fit_mix_num,
             model_type = run_meta$model_type,
             sim_scenario = run_meta$sim_scenario,
             days = run_meta$days,
             emission_overlap = run_meta$emission_overlap,
             parameter = names(counts),
             count = as.numeric(counts),
             stringsAsFactors = FALSE)
}

extract_survival_coefficients <- function(file,run_meta,true_params,est_params){
  true_beta <- param_safe(true_params,"beta_vec")
  estimated_beta <- param_safe(est_params,"beta_vec")
  if (is.null(true_beta) || is.null(estimated_beta) ||
      length(true_beta) == 0 || length(estimated_beta) == 0 ||
      length(true_beta) != length(estimated_beta)){
    return(NULL)
  }

  error <- as.numeric(estimated_beta) - as.numeric(true_beta)
  data.frame(file = file,
             sim_num = run_meta$sim_num,
             true_mix_num = run_meta$true_mix_num,
             fit_mix_num = run_meta$fit_mix_num,
             model_type = run_meta$model_type,
             sim_scenario = run_meta$sim_scenario,
             days = run_meta$days,
             emission_overlap = run_meta$emission_overlap,
             class = seq_along(error),
             true_beta = as.numeric(true_beta),
             estimated_beta = as.numeric(estimated_beta),
             error = error,
             squared_error = error^2,
             stringsAsFactors = FALSE)
}

parse_one_file <- function(file){
  file_name <- basename(file)
  file_meta <- extract_filename_metadata(file_name)

  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file,envir = load_env)
  if (!"to_save" %in% loaded_names && !exists("to_save",envir = load_env)){
    stop("Loaded file does not contain to_save")
  }
  to_save <- get("to_save",envir = load_env)

  settings <- section_safe(to_save,"settings")
  data_source <- as_character_or_na(setting_or_file(settings,file_meta,
                                                   "data_source","simulation"))
  real_data <- as_logical_or_na(setting_or_file(settings,file_meta,
                                               "real_data",FALSE))
  sim_scenario <- as_number_or_na(setting_or_file(settings,file_meta,
                                                 "sim_scenario",NA_real_))
  if (identical(data_source,"nhanes") || isTRUE(real_data) ||
      is.na(sim_scenario)){
    return(list(skipped = data.frame(file = file,
                                     file_name = file_name,
                                     reason = "not a simulation result",
                                     stringsAsFactors = FALSE)))
  }

  true_params <- section_safe(to_save,"true_params")
  est_params <- section_safe(to_save,"est_params")
  diagnostics <- section_safe(to_save,"diagnostics")

  scenario <- scenario_lookup(sim_scenario)
  file_info <- file.info(file)

  run_meta <- list(
    sim_num = as_number_or_na(setting_or_file(settings,file_meta,"sim_num",
                                             NA_real_)),
    true_mix_num = as_number_or_na(setting_or_file(settings,file_meta,
                                                  "true_mix_num",NA_real_)),
    fit_mix_num = as_number_or_na(setting_or_file(settings,file_meta,
                                                 "fit_mix_num",NA_real_)),
    model_type = as_character_or_na(setting_or_file(settings,file_meta,
                                                   "model_type",NA_character_)),
    data_source = data_source,
    sim_scenario = sim_scenario,
    days = scenario$days,
    num_people = scenario$num_people,
    missing_perc = scenario$missing_perc,
    emission_overlap =
      as_character_or_na(setting_or_file(settings,file_meta,
                                        "emission_overlap","low"))
  )

  bic <- section_safe(to_save,"bic")
  aic <- section_safe(to_save,"aic")
  new_likelihood <- as_number_or_na(param_safe(est_params,"new_likelihood"))
  parameter_count <- as_number_or_na(
    if (!is.null(diagnostics)) diagnostics$bic_parameter_count else NULL
  )
  bic_sample_size <- as_number_or_na(
    if (!is.null(diagnostics)) diagnostics$bic_sample_size else run_meta$num_people
  )

  if ((is.null(aic) || length(aic) == 0) &&
      is.finite(new_likelihood) && is.finite(parameter_count)){
    aic <- (2 * parameter_count) - (2 * new_likelihood)
  }
  if ((is.null(bic) || length(bic) == 0) &&
      is.finite(new_likelihood) && is.finite(parameter_count) &&
      is.finite(bic_sample_size)){
    bic <- parameter_count * log(bic_sample_size) - (2 * new_likelihood)
  }

  training_ibs <- metric_from_nested(
    diagnostics,"training","ibs",diag_value(diagnostics,"ibs",slot = 2)
  )
  training_cindex <- metric_from_nested(
    diagnostics,"training","cindex",diag_value(diagnostics,"cindex",slot = 1)
  )
  test_ibs <- metric_from_nested(diagnostics,"test","ibs",NA_real_)
  test_cindex <- metric_from_nested(diagnostics,"test","cindex",NA_real_)

  survival_coefficients <- extract_survival_coefficients(
    file,run_meta,true_params,est_params
  )
  survival_beta_comparable <- !is.null(survival_coefficients)
  survival_beta_rmse <- if (survival_beta_comparable){
    sqrt(mean(survival_coefficients$squared_error))
  } else {
    NA_real_
  }
  survival_beta_rmse_non_reference <-
    if (survival_beta_comparable && nrow(survival_coefficients) > 1){
      sqrt(mean(survival_coefficients$squared_error[-1]))
    } else {
      NA_real_
    }
  survival_beta_sse <- if (survival_beta_comparable){
    sum(survival_coefficients$squared_error)
  } else {
    NA_real_
  }

  class_entropy <- if (!is.null(diagnostics)) diagnostics$class_entropy else NULL

  run_level <- data.frame(
    file = file,
    file_name = file_name,
    sim_num = run_meta$sim_num,
    true_mix_num = run_meta$true_mix_num,
    fit_mix_num = run_meta$fit_mix_num,
    model_type = run_meta$model_type,
    data_source = run_meta$data_source,
    sim_scenario = run_meta$sim_scenario,
    days = run_meta$days,
    num_people = run_meta$num_people,
    missing_perc = run_meta$missing_perc,
    emission_overlap = run_meta$emission_overlap,
    init_jitter_scale =
      as_number_or_na(setting_or_file(settings,file_meta,
                                     "init_jitter_scale",NA_real_)),
    class_selection_run =
      as_logical_or_na(setting_or_file(settings,file_meta,
                                      "class_selection_run",FALSE)),
    save_reduced_output =
      as_logical_or_na(setting_or_file(settings,file_meta,
                                      "save_reduced_output",FALSE)),
    bic = as_number_or_na(bic),
    aic = as_number_or_na(aic),
    new_likelihood = new_likelihood,
    parameter_count = parameter_count,
    likelihood_type = as_character_or_na(
      if (!is.null(diagnostics)) diagnostics$bic_likelihood else NA_character_
    ),
    training_ibs = as_number_or_na(training_ibs),
    training_cindex = as_number_or_na(training_cindex),
    test_ibs = as_number_or_na(test_ibs),
    test_cindex = as_number_or_na(test_cindex),
    ibs2 = as_number_or_na(diag_value(diagnostics,"ibs2",slot = 4)),
    class_entropy_mean = as_number_or_na(
      list_value(class_entropy,"normalized_mean")
    ),
    class_entropy_median = as_number_or_na(
      list_value(class_entropy,"normalized_median")
    ),
    class_entropy_q25 = as_number_or_na(
      list_value(class_entropy,"normalized_q25")
    ),
    class_entropy_q75 = as_number_or_na(
      list_value(class_entropy,"normalized_q75")
    ),
    class_entropy_mean_max_posterior = as_number_or_na(
      list_value(class_entropy,"mean_max_posterior")
    ),
    survival_beta_comparable = survival_beta_comparable,
    survival_beta_rmse = survival_beta_rmse,
    survival_beta_rmse_non_reference = survival_beta_rmse_non_reference,
    survival_beta_sse = survival_beta_sse,
    file_size = file_info$size,
    file_modified_time = file_info$mtime,
    stringsAsFactors = FALSE
  )

  list(
    run_level = run_level,
    confusion = extract_confusion_table(file,run_meta,diagnostics),
    parameter_counts = extract_parameter_counts(file,run_meta,diagnostics),
    survival_coefficients = survival_coefficients
  )
}

make_class_selection <- function(run_level){
  if (nrow(run_level) == 0){
    return(empty_class_selection())
  }

  group_cols <- c("true_mix_num","model_type","sim_scenario","days",
                  "num_people","missing_perc","emission_overlap","sim_num")
  criteria <- data.frame(
    criterion = c("bic","aic","training_ibs","training_cindex",
                  "test_ibs","test_cindex"),
    direction = c("min","min","min","max","min","max"),
    stringsAsFactors = FALSE
  )

  group_key <- interaction(run_level[,group_cols],drop = TRUE,sep = "\r")
  groups <- split(seq_len(nrow(run_level)),group_key)
  rows <- list()

  for (indices in groups){
    group_data <- run_level[indices,,drop = FALSE]
    for (criterion_ind in seq_len(nrow(criteria))){
      criterion <- criteria$criterion[[criterion_ind]]
      direction <- criteria$direction[[criterion_ind]]
      values <- group_data[[criterion]]
      valid <- is.finite(values)
      if (!any(valid)){
        next
      }
      valid_data <- group_data[valid,,drop = FALSE]
      valid_values <- values[valid]
      selected_local <- if (direction == "min"){
        order(valid_values,valid_data$fit_mix_num)[[1]]
      } else {
        order(-valid_values,valid_data$fit_mix_num)[[1]]
      }
      selected <- valid_data[selected_local,,drop = FALSE]
      rows[[length(rows) + 1]] <- data.frame(
        true_mix_num = selected$true_mix_num,
        model_type = selected$model_type,
        sim_scenario = selected$sim_scenario,
        days = selected$days,
        num_people = selected$num_people,
        missing_perc = selected$missing_perc,
        emission_overlap = selected$emission_overlap,
        sim_num = selected$sim_num,
        criterion = criterion,
        direction = direction,
        selected_fit_mix_num = selected$fit_mix_num,
        selected_value = selected[[criterion]],
        correct_true_mix_num = selected$fit_mix_num == selected$true_mix_num,
        num_fits_available = nrow(group_data),
        selected_file = selected$file,
        stringsAsFactors = FALSE
      )
    }
  }

  bind_rows(rows,empty_class_selection)
}

files <- if (dir.exists(input_dir)){
  list.files(input_dir,pattern = "\\.rda$",full.names = TRUE,
             recursive = TRUE,ignore.case = TRUE)
} else {
  character()
}

run_rows <- list()
confusion_rows <- list()
parameter_count_rows <- list()
survival_coefficient_rows <- list()
error_rows <- list()
skipped_rows <- list()

for (file in files){
  file_name <- basename(file)
  if (grepl("^Inter",file_name)){
    skipped_rows[[length(skipped_rows) + 1]] <- data.frame(
      file = file,file_name = file_name,reason = "interim checkpoint",
      stringsAsFactors = FALSE
    )
    next
  }

  parsed <- tryCatch(
    parse_one_file(file),
    error = function(e) list(error = data.frame(
      file = file,file_name = file_name,reason = conditionMessage(e),
      stringsAsFactors = FALSE
    ))
  )

  if (!is.null(parsed$error)){
    error_rows[[length(error_rows) + 1]] <- parsed$error
    next
  }
  if (!is.null(parsed$skipped)){
    skipped_rows[[length(skipped_rows) + 1]] <- parsed$skipped
    next
  }

  run_rows[[length(run_rows) + 1]] <- parsed$run_level
  confusion_rows[[length(confusion_rows) + 1]] <- parsed$confusion
  parameter_count_rows[[length(parameter_count_rows) + 1]] <-
    parsed$parameter_counts
  survival_coefficient_rows[[length(survival_coefficient_rows) + 1]] <-
    parsed$survival_coefficients
}

run_level <- bind_rows(run_rows,empty_run_level)
confusion_tables <- bind_rows(confusion_rows,empty_confusion)
parameter_counts <- bind_rows(parameter_count_rows,empty_parameter_counts)
survival_coefficients <- bind_rows(survival_coefficient_rows,
                                   empty_survival_coefficients)
class_selection <- make_class_selection(run_level)
parse_errors <- bind_rows(error_rows,empty_file_issue)
skipped_files <- bind_rows(skipped_rows,empty_file_issue)

tables <- list(
  simulation_run_level = run_level,
  simulation_class_selection = class_selection,
  simulation_confusion_tables = confusion_tables,
  simulation_parameter_counts = parameter_counts,
  simulation_survival_coefficients = survival_coefficients,
  parse_errors = parse_errors,
  skipped_files = skipped_files
)

for (table_name in names(tables)){
  saveRDS(tables[[table_name]],
          file.path(output_dir,paste0(table_name,".rds")))
}
saveRDS(tables,file.path(output_dir,"simulation_parse_tables.rds"))

message("Parsed ",nrow(run_level)," simulation result files")
message("Skipped ",nrow(skipped_files)," files; errors in ",
        nrow(parse_errors)," files")
message("Wrote parsed tables to: ",output_dir)
