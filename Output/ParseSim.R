#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Routputs","Routputs")
}

output_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Output","parse_sim_results.rds")
}

CONTROL_PARAMS <- c(
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
  "emission_overlap",
  "time_limit_hours",
  "memory_limit_gb"
)

is_scalar <- function(x){
  is.atomic(x) && length(x) == 1 && is.null(dim(x))
}

clean_name <- function(x){
  x <- gsub("[^A-Za-z0-9_]+","_",x)
  x <- gsub("_+","_",x)
  x <- gsub("^_|_$","",x)
  ifelse(nchar(x) == 0,"value",x)
}

regex_value <- function(pattern,text,default = NA_character_){
  match <- regexec(pattern,text,perl = TRUE)
  found <- regmatches(text,match)[[1]]
  if (length(found) < 2){
    return(default)
  }
  found[[2]]
}

setting_value <- function(settings,param_name,param_index){
  if (!is.null(settings) && !is.null(settings[[param_name]]) &&
      length(settings[[param_name]]) > 0){
    return(settings[[param_name]][[1]])
  }

  command_args <- settings[["command_args"]]
  if (!is.null(command_args) && length(command_args) >= param_index){
    return(command_args[[param_index]])
  }

  NA
}

flatten_diagnostics <- function(x,prefix = "diagnostic"){
  out <- list()
  if (is.null(x)){
    return(out)
  }

  item_names <- names(x)
  for (i in seq_along(x)){
    item_name <- if (!is.null(item_names) && nzchar(item_names[[i]])){
      clean_name(item_names[[i]])
    } else {
      paste0("value",i)
    }
    column_name <- paste(prefix,item_name,sep = "_")
    value <- x[[i]]

    if (is_scalar(value)){
      out[[column_name]] <- value
    } else if (is.list(value) && !inherits(value,"data.frame")){
      out <- c(out,flatten_diagnostics(value,column_name))
    } else if (is.atomic(value) && length(value) > 1 &&
               !is.null(names(value)) && all(nzchar(names(value)))){
      for (j in seq_along(value)){
        out[[paste(column_name,clean_name(names(value)[[j]]),sep = "_")]] <-
          value[[j]]
      }
    } else {
      out[[column_name]] <- value
    }
  }

  out
}

row_to_dataframe <- function(row,all_columns,list_columns){
  out <- vector("list",length(all_columns))
  names(out) <- all_columns

  for (column in all_columns){
    if (column %in% names(row)){
      value <- row[[column]]
    } else {
      value <- NA
    }

    if (column %in% list_columns){
      out[[column]] <- I(list(value))
    } else {
      out[[column]] <- value
    }
  }

  as.data.frame(out,check.names = FALSE,stringsAsFactors = FALSE)
}

rows_to_dataframe <- function(rows){
  if (length(rows) == 0){
    return(data.frame())
  }

  all_columns <- unique(unlist(lapply(rows,names),use.names = FALSE))
  list_columns <- character()
  for (row in rows){
    for (column in names(row)){
      if (!is_scalar(row[[column]])){
        list_columns <- union(list_columns,column)
      }
    }
  }

  row_frames <- lapply(rows,row_to_dataframe,
                       all_columns = all_columns,
                       list_columns = list_columns)
  result <- do.call(rbind,row_frames)
  row.names(result) <- NULL
  result
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
    value ~ true_mix_num + fit_mix_num + model_type +
      simulation_days + num_people + emission_overlap,
    data = summary_input,
    FUN = max_finite
  )
  names(max_data)[names(max_data) == "value"] <- paste0("max_",prefix)

  q95_data <- aggregate(
    value ~ true_mix_num + fit_mix_num + model_type +
      simulation_days + num_people + emission_overlap,
    data = summary_input,
    FUN = q95_finite
  )
  names(q95_data)[names(q95_data) == "value"] <- paste0("q95_",prefix)

  summary_data <- merge(summary_data,max_data,by = scenario_cols,all.x = TRUE)
  merge(summary_data,q95_data,by = scenario_cols,all.x = TRUE)
}

parse_one_file <- function(file){
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file,envir = load_env)
  if (!"to_save" %in% loaded_names || !exists("to_save",envir = load_env)){
    stop("File does not contain object named to_save")
  }

  saved <- get("to_save",envir = load_env)
  settings <- saved[["settings"]]
  diagnostics <- saved[["diagnostics"]]
  file_name <- basename(file)

  row <- list(
    file = normalizePath(file,mustWork = TRUE),
    file_name = file_name,
    sim_num = if (!is.null(settings[["sim_num"]])){
      settings[["sim_num"]][[1]]
    } else {
      suppressWarnings(as.numeric(regex_value("Seed([0-9]+)",file_name)))
    }
  )

  for (i in seq_along(CONTROL_PARAMS)){
    param <- CONTROL_PARAMS[[i]]
    row[[param]] <- setting_value(settings,param,i)
  }

  row[["aic"]] <- if (is_scalar(saved[["aic"]])) saved[["aic"]] else NA
  row[["bic"]] <- if (is_scalar(saved[["bic"]])) saved[["bic"]] else NA
  row <- c(row,flatten_diagnostics(diagnostics))

  row
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

results <- rows_to_dataframe(rows)

dir.create(dirname(output_file),recursive = TRUE,showWarnings = FALSE)
saveRDS(results,output_file)

if (length(errors) > 0){
  error_file <- sub("[.]rds$","_errors.rds",output_file,ignore.case = TRUE)
  saveRDS(do.call(rbind,errors),error_file)
  message("Saved parse errors to: ",error_file)
}

message("Read ",length(files)," files from: ",input_dir)
message("Saved ",nrow(results)," rows and ",ncol(results),
        " columns to: ",output_file)


scenario_cols <- c("true_mix_num","fit_mix_num","model_type",
                   "simulation_days","num_people","emission_overlap")

if (nrow(results) > 0 && all(scenario_cols %in% names(results))){
  scenario_counts <- aggregate(
    sim_num ~ true_mix_num + fit_mix_num + model_type +
      simulation_days + num_people + emission_overlap,
    data = results,
    FUN = function(x) length(unique(x))
  )
  names(scenario_counts)[names(scenario_counts) == "sim_num"] <- "n_seeds"

  scenario_counts <- add_numeric_summaries(
    scenario_counts,
    results,
    "diagnostic_convergence_total_iteration_seconds",
    "em_seconds",
    scenario_cols
  )

  memory_cols <- c(
    diagnostic_memory_peak_gb = "memory_peak_gb",
    max_rss_gb = "max_rss_gb",
    memory_used_gb = "memory_used_gb",
    mem_used_gb = "mem_used_gb"
  )
  for (memory_col in names(memory_cols)){
    scenario_counts <- add_numeric_summaries(
      scenario_counts,
      results,
      memory_col,
      memory_cols[[memory_col]],
      scenario_cols
    )
  }

  limit_cols <- c(
    time_limit_hours = "time_limit_hours",
    memory_limit_gb = "memory_limit_gb"
  )
  for (limit_col in names(limit_cols)){
    if (!limit_col %in% names(results)){
      next
    }

    limit_data <- results[,scenario_cols,drop = FALSE]
    limit_data$value <- suppressWarnings(as.numeric(results[[limit_col]]))
    limit_data <- limit_data[is.finite(limit_data$value),,drop = FALSE]

    if (nrow(limit_data) == 0){
      next
    }

    max_limit <- aggregate(
      value ~ true_mix_num + fit_mix_num + model_type +
        simulation_days + num_people + emission_overlap,
      data = limit_data,
      FUN = max_finite
    )
    names(max_limit)[names(max_limit) == "value"] <-
      paste0("max_",limit_cols[[limit_col]])
    scenario_counts <- merge(scenario_counts,max_limit,
                             by = scenario_cols,all.x = TRUE)
  }
} else {
  scenario_counts <- data.frame()
}

complete_scenario <- scenario_counts$n_seeds == 100
scenario_counts$recommended_time_hours <- ifelse(
  complete_scenario,
  ceiling(pmax(
    scenario_counts$max_em_seconds / 3600 * 1.10,
    scenario_counts$q95_em_seconds / 3600 * 1.10,
    1
  )),
  ceiling(pmax(
    scenario_counts$max_em_seconds / 3600 * 1.25,
    scenario_counts$q95_em_seconds / 3600 * 1.50,
    1
  ))
)

if (all(c("max_memory_peak_gb","q95_memory_peak_gb") %in%
        names(scenario_counts))){
  scenario_counts$recommended_mem_gb <- ifelse(
    complete_scenario,
    ceiling(pmax(
      scenario_counts$max_memory_peak_gb * 1.10,
      scenario_counts$q95_memory_peak_gb * 1.10,
      1
    )),
    ceiling(pmax(
      scenario_counts$max_memory_peak_gb * 1.25,
      scenario_counts$q95_memory_peak_gb * 1.50,
      1
    ))
  )
} else {
  scenario_counts$recommended_mem_gb <- NA_real_
}

scenario_counts_file <- if (grepl("[.]rds$",output_file,ignore.case = TRUE)){
  sub("[.]rds$","_scenario_counts.rds",output_file,ignore.case = TRUE)
} else {
  paste0(output_file,"_scenario_counts.rds")
}

# saveRDS(scenario_counts,scenario_counts_file)
# message("Saved scenario counts to: ",scenario_counts_file)
