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
  "emission_overlap"
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

