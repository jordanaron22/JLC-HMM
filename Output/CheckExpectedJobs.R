#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

manifest_file <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Output","expected_jobs.tsv")
}

parsed_results_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Output","parse_sim_results.rds")
}

output_file <- if (length(args) >= 3){
  args[[3]]
} else {
  file.path("Output","expected_job_check.rds")
}

stop_on_missing <- length(args) >= 4 &&
  tolower(args[[4]]) %in% c("true","t","1","yes","y")

key_cols <- c(
  "sim_num",
  "fit_mix_num",
  "model_type",
  "data_source",
  "run_bootstrap",
  "init_jitter_scale",
  "run_leave_one_out_cv",
  "use_hot_start",
  "simulation_days",
  "num_people",
  "true_mix_num",
  "save_reduced_output",
  "class_selection_run",
  "emission_overlap"
)

numeric_cols <- c(
  "sim_num",
  "fit_mix_num",
  "init_jitter_scale",
  "simulation_days",
  "num_people",
  "true_mix_num"
)

text_value <- function(x){
  x <- trimws(as.character(x))
  x[is.na(x)] <- "<NA>"
  tolower(x)
}

numeric_value <- function(x){
  numeric_x <- suppressWarnings(as.numeric(x))
  out <- as.character(numeric_x)
  out[is.na(numeric_x)] <- text_value(x[is.na(numeric_x)])
  out
}

normalize_column <- function(data,column){
  if (column %in% numeric_cols){
    return(numeric_value(data[[column]]))
  }
  text_value(data[[column]])
}

make_key <- function(data,cols){
  normalized <- lapply(cols,function(column){
    normalize_column(data,column)
  })
  do.call(paste,c(normalized,sep = "\r"))
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

add_numeric_summaries <- function(summary_data,data,value_col,prefix,
                                  scenario_cols){
  if (!value_col %in% names(data)){
    return(summary_data)
  }

  summary_input <- data[,scenario_cols,drop = FALSE]
  summary_input$value <- suppressWarnings(as.numeric(data[[value_col]]))
  summary_input <- summary_input[is.finite(summary_input$value),,drop = FALSE]

  if (nrow(summary_input) == 0){
    return(summary_data)
  }

  max_data <- aggregate(
    value ~ fit_mix_num + model_type + simulation_days +
      num_people + true_mix_num + emission_overlap,
    data = summary_input,
    FUN = max_finite
  )
  names(max_data)[names(max_data) == "value"] <- paste0("max_",prefix)

  q95_data <- aggregate(
    value ~ fit_mix_num + model_type + simulation_days +
      num_people + true_mix_num + emission_overlap,
    data = summary_input,
    FUN = q95_finite
  )
  names(q95_data)[names(q95_data) == "value"] <- paste0("q95_",prefix)

  summary_data <- merge(summary_data,max_data,by = scenario_cols,all.x = TRUE)
  merge(summary_data,q95_data,by = scenario_cols,all.x = TRUE)
}

require_columns <- function(data,cols,label){
  missing_cols <- setdiff(cols,names(data))
  if (length(missing_cols) > 0){
    stop(label," is missing columns: ",
         paste(missing_cols,collapse = ", "))
  }
}

if (!file.exists(manifest_file)){
  stop("Manifest file does not exist: ",manifest_file)
}

if (!file.exists(parsed_results_file)){
  stop("Parsed results file does not exist: ",parsed_results_file)
}

manifest <- read.delim(manifest_file,stringsAsFactors = FALSE,
                       check.names = FALSE)
results <- readRDS(parsed_results_file)

require_columns(manifest,key_cols,"Manifest")
require_columns(results,key_cols,"Parsed results")

manifest$expected_key <- make_key(manifest,key_cols)
results$expected_key <- make_key(results,key_cols)

manifest$expected_file_name <- basename(manifest$expected_file)
manifest$expected_file_exists <- file.exists(manifest$expected_file)
manifest$requested_time_hours <- suppressWarnings(
  as.numeric(manifest$requested_time)
)
manifest$requested_mem_gb <- suppressWarnings(
  as.numeric(manifest$requested_mem)
)

result_lookup_cols <- c("expected_key","file","file_name",
                        "diagnostic_convergence_total_iteration_seconds")
result_lookup_cols <- intersect(result_lookup_cols,names(results))
result_lookup <- unique(results[,result_lookup_cols,drop = FALSE])
result_lookup$parsed <- TRUE

checked <- merge(manifest,result_lookup,by = "expected_key",all.x = TRUE)
checked$parsed[is.na(checked$parsed)] <- FALSE

missing_jobs <- checked[!checked$parsed,,drop = FALSE]
completed_jobs <- checked[checked$parsed,,drop = FALSE]

duplicate_expected <- manifest[
  duplicated(manifest$expected_key) |
    duplicated(manifest$expected_key,fromLast = TRUE),
  ,
  drop = FALSE
]

unexpected_results <- results[
  !results$expected_key %in% manifest$expected_key,
  ,
  drop = FALSE
]

summary <- data.frame(
  total_expected = nrow(manifest),
  total_completed = nrow(completed_jobs),
  total_missing = nrow(missing_jobs),
  expected_files_exist = sum(manifest$expected_file_exists),
  duplicate_expected_rows = nrow(duplicate_expected),
  parsed_rows = nrow(results),
  unexpected_parsed_rows = nrow(unexpected_results),
  stringsAsFactors = FALSE
)

scenario_cols <- c("fit_mix_num","model_type","simulation_days",
                   "num_people","true_mix_num","emission_overlap")

missing_by_scenario <- data.frame()
if (nrow(missing_jobs) > 0){
  missing_by_scenario <- aggregate(
    sim_num ~ fit_mix_num + model_type + simulation_days +
      num_people + true_mix_num + emission_overlap,
    data = missing_jobs,
    FUN = length
  )
  names(missing_by_scenario)[names(missing_by_scenario) == "sim_num"] <-
    "n_missing"
}

completed_by_scenario <- data.frame()
if (nrow(completed_jobs) > 0){
  completed_by_scenario <- aggregate(
    sim_num ~ fit_mix_num + model_type + simulation_days +
      num_people + true_mix_num + emission_overlap,
    data = completed_jobs,
    FUN = length
  )
  names(completed_by_scenario)[names(completed_by_scenario) == "sim_num"] <-
    "n_completed"

  completed_by_scenario <- add_numeric_summaries(
    completed_by_scenario,
    completed_jobs,
    "diagnostic_convergence_total_iteration_seconds",
    "em_seconds",
    scenario_cols
  )

  completed_by_scenario <- add_numeric_summaries(
    completed_by_scenario,
    completed_jobs,
    "requested_time_hours",
    "requested_time_hours",
    scenario_cols
  )

  completed_by_scenario <- add_numeric_summaries(
    completed_by_scenario,
    completed_jobs,
    "requested_mem_gb",
    "requested_mem_gb",
    scenario_cols
  )
}

report <- list(
  summary = summary,
  missing_jobs = missing_jobs,
  completed_jobs = completed_jobs,
  missing_by_scenario = missing_by_scenario,
  completed_by_scenario = completed_by_scenario,
  duplicate_expected = duplicate_expected,
  unexpected_results = unexpected_results
)

dir.create(dirname(output_file),recursive = TRUE,showWarnings = FALSE)
saveRDS(report,output_file)

print(summary)
if (nrow(missing_by_scenario) > 0){
  print(missing_by_scenario)
}
message("Saved check report to: ",output_file)

if (stop_on_missing && nrow(missing_jobs) > 0){
  stop("Missing ",nrow(missing_jobs)," expected simulation outputs.")
}
