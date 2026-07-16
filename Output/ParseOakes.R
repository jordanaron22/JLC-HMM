#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_root <- if (length(args) >= 1) args[[1]] else file.path("Output","Routputs")
output_file <- if (length(args) >= 2) args[[2]] else file.path("Output","parse_oakes_results.rds")

Z_95 <- 1.96

SAVED_SECTION_SLOT <- c(true_params = 1, est_params = 2, settings = 7)
PARAM_SLOT <- c(beta_vec = 7, beta_se = 17)

is_scalar <- function(x){
  is.atomic(x) && length(x) == 1 && is.null(dim(x))
}

scalar_or_na <- function(x){
  if (is.null(x) || length(x) == 0){
    return(NA)
  }
  x[[1]]
}

as_numeric_or_na <- function(x){
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0) NA_real_ else out
}

finite_or_na <- function(x){
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || !is.finite(x[[1]])){
    return(NA_real_)
  }
  x[[1]]
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

q95_finite <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  unname(quantile(x,0.95,names = FALSE,type = 7))
}

regex_value <- function(pattern,text,default = NA_character_){
  match <- regexec(pattern,text,perl = TRUE)
  found <- regmatches(text,match)[[1]]
  if (length(found) < 2){
    return(default)
  }
  found[[2]]
}

load_to_save <- function(file){
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file,envir = load_env)
  if (!"to_save" %in% loaded_names || !exists("to_save",envir = load_env)){
    stop("File does not contain object named to_save")
  }
  get("to_save",envir = load_env)
}

get_section <- function(to_save,name,slot,required = TRUE){
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

get_param <- function(params,name,slot,required = TRUE){
  if (!is.null(params) && !is.null(names(params)) && name %in% names(params)){
    return(params[[name]])
  }
  if (!is.null(params) && length(params) >= slot){
    return(params[[slot]])
  }
  if (required){
    stop("Missing saved parameter: ",name)
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

get_model_type <- function(to_save){
  settings <- get_section(to_save,"settings",SAVED_SECTION_SLOT[["settings"]],
                          required = FALSE)
  model_type <- get_setting(settings,"model_type",NA_character_)
  if (!is.na(model_type)){
    return(as.character(model_type))
  }

  incl_surv <- get_setting(settings,"incl_surv",NA)
  if (!is.na(incl_surv)){
    if (as.numeric(incl_surv) == 2){
      return("joint")
    }
    if (as.numeric(incl_surv) == 0){
      return("two_stage")
    }
  }

  NA_character_
}

classify_file <- function(to_save){
  model_type <- get_model_type(to_save)
  has_oakes <- !is.null(to_save[["oakes"]])

  if (identical(model_type,"joint") && has_oakes){
    return("joint_oakes")
  }
  if (identical(model_type,"joint") && !has_oakes){
    return("joint_jmhmm")
  }
  if (identical(model_type,"two_stage")){
    return("two_stage_jmhmm")
  }

  "unknown"
}

get_fit_mix_num <- function(settings,est_beta,file_name){
  value <- get_setting(settings,"fit_mix_num",NA)
  if (is.na(value)){
    value <- regex_value("FitMix([0-9]+)",file_name)
  }
  if (is.na(value)){
    value <- length(est_beta)
  }
  as.integer(as_numeric_or_na(value))
}

get_oakes_baseline_mode <- function(to_save,file_name){
  mode <- to_save[["oakes"]][["settings"]][["survival_baseline_mode"]]
  mode <- scalar_or_na(mode)
  if (is.na(mode)){
    mode <- regex_value("_oakes_(fixed|profiled)_se[.]rda$",
                        file_name)
  }
  mode
}

get_diagnostics <- function(to_save,file_type){
  if (file_type == "joint_oakes" && !is.null(to_save[["oakes"]])){
    return(to_save[["oakes"]][["diagnostics"]])
  }
  to_save[["diagnostics"]]
}

inventory_row <- function(file,to_save,file_type){
  file_name <- basename(file)
  settings <- get_section(to_save,"settings",SAVED_SECTION_SLOT[["settings"]],
                          required = FALSE)
  diagnostics <- get_diagnostics(to_save,file_type)

  data.frame(
    file = normalizePath(file,mustWork = TRUE),
    file_name = file_name,
    file_type = file_type,
    model_type = get_model_type(to_save),
    sim_num = as_numeric_or_na(get_setting(
      settings,"sim_num",regex_value("Seed([0-9]+)",file_name)
    )),
    simulation_days = as_numeric_or_na(get_setting(
      settings,"simulation_days",regex_value("Days([0-9]+)",file_name)
    )),
    num_people = as_numeric_or_na(get_setting(
      settings,"num_people",regex_value("People([0-9]+)",file_name)
    )),
    true_mix_num = as_numeric_or_na(get_setting(
      settings,"true_mix_num",regex_value("TrueMix([0-9]+)",file_name)
    )),
    fit_mix_num = as_numeric_or_na(get_setting(
      settings,"fit_mix_num",regex_value("FitMix([0-9]+)",file_name)
    )),
    emission_overlap = as.character(get_setting(
      settings,"emission_overlap",
      tolower(regex_value("Overlap([A-Za-z]+)",file_name))
    )),
    oakes_baseline_mode = if (file_type == "joint_oakes"){
      get_oakes_baseline_mode(to_save,file_name)
    } else {
      NA_character_
    },
    runtime_seconds = finite_or_na(diagnostics[["runtime_seconds"]]),
    memory_peak_gb = finite_or_na(diagnostics[["memory_peak_gb"]]),
    time_limit_hours = as_numeric_or_na(get_setting(settings,
                                                    "time_limit_hours",NA)),
    memory_limit_gb = as_numeric_or_na(get_setting(settings,
                                                   "memory_limit_gb",NA)),
    stringsAsFactors = FALSE
  )
}

file_metadata <- function(file,to_save,file_type,est_beta){
  file_name <- basename(file)
  settings <- get_section(to_save,"settings",SAVED_SECTION_SLOT[["settings"]],
                          required = FALSE)
  data.frame(
    file = normalizePath(file,mustWork = TRUE),
    file_name = file_name,
    file_type = file_type,
    model_type = get_model_type(to_save),
    se_method = if (file_type == "joint_oakes"){
      "oakes_schur"
    } else if (file_type == "two_stage_jmhmm"){
      "non_oakes_beta_se"
    } else {
      NA_character_
    },
    sim_num = as_numeric_or_na(get_setting(
      settings,"sim_num",regex_value("Seed([0-9]+)",file_name)
    )),
    simulation_days = as_numeric_or_na(get_setting(
      settings,"simulation_days",regex_value("Days([0-9]+)",file_name)
    )),
    num_people = as_numeric_or_na(get_setting(
      settings,"num_people",regex_value("People([0-9]+)",file_name)
    )),
    true_mix_num = as_numeric_or_na(get_setting(
      settings,"true_mix_num",regex_value("TrueMix([0-9]+)",file_name)
    )),
    fit_mix_num = get_fit_mix_num(settings,est_beta,file_name),
    emission_overlap = as.character(get_setting(
      settings,"emission_overlap",
      tolower(regex_value("Overlap([A-Za-z]+)",file_name))
    )),
    oakes_baseline_mode = if (file_type == "joint_oakes"){
      get_oakes_baseline_mode(to_save,file_name)
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}

get_class_beta_se <- function(to_save,file_type){
  if (file_type == "joint_oakes"){
    se <- to_save[["oakes"]][["survival_schur"]][["se_surv"]]
    return(suppressWarnings(as.numeric(se)))
  }

  est_params <- get_section(to_save,"est_params",SAVED_SECTION_SLOT[["est_params"]])
  se <- get_param(est_params,"beta_se",PARAM_SLOT[["beta_se"]])
  suppressWarnings(as.numeric(se))
}

parse_class_beta_rows <- function(file,to_save,file_type){
  if (file_type == "joint_jmhmm"){
    stop("joint_jmhmm_skipped_use_oakes_file")
  }
  if (!file_type %in% c("joint_oakes","two_stage_jmhmm")){
    stop("unsupported_file_type: ",file_type)
  }

  true_params <- get_section(to_save,"true_params",
                             SAVED_SECTION_SLOT[["true_params"]])
  est_params <- get_section(to_save,"est_params",
                            SAVED_SECTION_SLOT[["est_params"]])

  true_beta <- suppressWarnings(as.numeric(
    get_param(true_params,"beta_vec",PARAM_SLOT[["beta_vec"]])
  ))
  est_beta <- suppressWarnings(as.numeric(
    get_param(est_params,"beta_vec",PARAM_SLOT[["beta_vec"]])
  ))
  se <- get_class_beta_se(to_save,file_type)

  meta <- file_metadata(file,to_save,file_type,est_beta)
  fit_mix_num <- meta$fit_mix_num[[1]]
  if (!is.finite(fit_mix_num) || fit_mix_num < 2){
    stop("fit_mix_num must be at least 2")
  }

  if (length(est_beta) < fit_mix_num){
    stop("Estimated beta_vec is shorter than fit_mix_num")
  }
  if (length(true_beta) < fit_mix_num){
    stop("True beta_vec is shorter than fit_mix_num")
  }
  if (length(se) < fit_mix_num){
    stop("SE vector is shorter than fit_mix_num")
  }

  class_indices <- 2:fit_mix_num
  rows <- meta[rep(1,length(class_indices)),,drop = FALSE]
  rows$class_index <- class_indices
  rows$param_name <- paste0("class_",class_indices)
  rows$estimate <- est_beta[class_indices]
  rows$truth <- true_beta[class_indices]
  rows$se <- se[class_indices]
  rows$bias <- rows$estimate - rows$truth
  rows$ci_lower <- rows$estimate - Z_95 * rows$se
  rows$ci_upper <- rows$estimate + Z_95 * rows$se
  rows$ci_width <- rows$ci_upper - rows$ci_lower
  rows$covered <- rows$truth >= rows$ci_lower & rows$truth <= rows$ci_upper
  rows$valid <- is.finite(rows$estimate) & is.finite(rows$truth) &
    is.finite(rows$se)
  rows
}

summarize_class_beta <- function(rows){
  if (nrow(rows) == 0){
    return(data.frame())
  }

  group_cols <- c("model_type","se_method","simulation_days","num_people",
                  "true_mix_num","fit_mix_num","emission_overlap",
                  "param_name")
  group_key <- interaction(rows[,group_cols,drop = FALSE],
                           drop = TRUE,lex.order = TRUE)
  split_rows <- split(rows,group_key)

  summaries <- lapply(split_rows,function(dat){
    valid <- dat$valid
    valid_dat <- dat[valid,,drop = FALSE]
    out <- as.list(dat[1,group_cols,drop = FALSE])
    out$n <- nrow(dat)
    out$n_valid <- nrow(valid_dat)

    if (nrow(valid_dat) == 0){
      out$coverage <- NA_real_
      out$mean_bias <- NA_real_
      out$median_bias <- NA_real_
      out$rmse <- NA_real_
      out$mean_se <- NA_real_
      out$median_se <- NA_real_
      out$mean_ci_width <- NA_real_
    } else {
      out$coverage <- mean(valid_dat$covered)
      out$mean_bias <- mean(valid_dat$bias)
      out$median_bias <- median(valid_dat$bias)
      out$rmse <- sqrt(mean(valid_dat$bias^2))
      out$mean_se <- mean(valid_dat$se)
      out$median_se <- median(valid_dat$se)
      out$mean_ci_width <- mean(valid_dat$ci_width)
    }

    as.data.frame(out,check.names = FALSE,stringsAsFactors = FALSE)
  })

  result <- do.call(rbind,summaries)
  row.names(result) <- NULL
  result
}

summarize_resources <- function(file_inventory){
  if (nrow(file_inventory) == 0){
    return(data.frame())
  }

  keep <- file_inventory$file_type %in% c("joint_oakes","two_stage_jmhmm")
  data <- file_inventory[keep,,drop = FALSE]
  if (nrow(data) == 0){
    return(data.frame())
  }

  group_cols <- c("file_type","model_type","simulation_days","num_people",
                  "true_mix_num","fit_mix_num","emission_overlap",
                  "oakes_baseline_mode")
  group_key <- interaction(data[,group_cols,drop = FALSE],
                           drop = TRUE,lex.order = TRUE)
  split_data <- split(data,group_key)

  summaries <- lapply(split_data,function(dat){
    out <- as.list(dat[1,group_cols,drop = FALSE])
    out$n_files <- nrow(dat)
    out$n_seeds <- length(unique(dat$sim_num[is.finite(dat$sim_num)]))
    out$mean_runtime_seconds <- mean_finite(dat$runtime_seconds)
    out$q95_runtime_seconds <- q95_finite(dat$runtime_seconds)
    out$max_runtime_seconds <- max_finite(dat$runtime_seconds)
    out$mean_memory_peak_gb <- mean_finite(dat$memory_peak_gb)
    out$q95_memory_peak_gb <- q95_finite(dat$memory_peak_gb)
    out$max_memory_peak_gb <- max_finite(dat$memory_peak_gb)
    out$max_time_limit_hours <- max_finite(dat$time_limit_hours)
    out$max_memory_limit_gb <- max_finite(dat$memory_limit_gb)

    complete_scenario <- is.finite(out$n_seeds) && out$n_seeds >= 98
    if (is.finite(out$max_runtime_seconds) &&
        is.finite(out$q95_runtime_seconds)){
      out$recommended_time_hours <- if (complete_scenario){
        ceiling(pmax(out$max_runtime_seconds / 3600 * 1.10,
                     out$q95_runtime_seconds / 3600 * 1.15,
                     1))
      } else {
        ceiling(pmax(out$max_runtime_seconds / 3600 * 1.25,
                     out$q95_runtime_seconds / 3600 * 1.50,
                     1))
      }
    } else {
      out$recommended_time_hours <- NA_real_
    }

    if (is.finite(out$max_memory_peak_gb) &&
        is.finite(out$q95_memory_peak_gb)){
      out$recommended_mem_gb <- if (complete_scenario){
        ceiling(pmax(out$max_memory_peak_gb * 1.10,
                     out$q95_memory_peak_gb * 1.10,
                     1))
      } else {
        ceiling(pmax(out$max_memory_peak_gb * 1.25,
                     out$q95_memory_peak_gb * 1.50,
                     1))
      }
    } else {
      out$recommended_mem_gb <- NA_real_
    }

    as.data.frame(out,check.names = FALSE,stringsAsFactors = FALSE)
  })

  result <- do.call(rbind,summaries)
  row.names(result) <- NULL
  result[order(result$model_type,result$file_type,result$simulation_days,
               result$emission_overlap,result$oakes_baseline_mode),,
         drop = FALSE]
}

write_outputs <- function(parsed,output_file){
  dir.create(dirname(output_file),recursive = TRUE,showWarnings = FALSE)
  saveRDS(parsed,output_file)
}

if (!dir.exists(input_root)){
  stop("Input directory does not exist: ",input_root)
}

files <- list.files(input_root,pattern = "[.]rda$",full.names = TRUE,
                    recursive = TRUE,ignore.case = TRUE)
files <- files[!grepl("(^|[/\\\\])Inter",files)]

row_list <- list()
inventory_list <- list()
error_list <- list()

for (file in files){
  to_save <- tryCatch(
    load_to_save(file),
    error = function(e){
      error_list[[length(error_list) + 1L]] <<- data.frame(
        file = normalizePath(file,mustWork = TRUE),
        file_name = basename(file),
        file_type = "unknown",
        reason = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )
  if (is.null(to_save)){
    next
  }

  file_type <- classify_file(to_save)
  inventory_list[[length(inventory_list) + 1L]] <-
    inventory_row(file,to_save,file_type)

  parsed <- tryCatch(
    parse_class_beta_rows(file,to_save,file_type),
    error = function(e){
      error_list[[length(error_list) + 1L]] <<- data.frame(
        file = normalizePath(file,mustWork = TRUE),
        file_name = basename(file),
        file_type = file_type,
        reason = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      NULL
    }
  )

  if (!is.null(parsed)){
    row_list[[length(row_list) + 1L]] <- parsed
  }
}

class_beta_long <- if (length(row_list) > 0){
  do.call(rbind,row_list)
} else {
  data.frame()
}
row.names(class_beta_long) <- NULL

class_beta_summary <- summarize_class_beta(class_beta_long)

file_inventory <- if (length(inventory_list) > 0){
  do.call(rbind,inventory_list)
} else {
  data.frame(file = character(),file_name = character(),
             file_type = character())
}
row.names(file_inventory) <- NULL

resource_summary <- summarize_resources(file_inventory)

parse_errors <- if (length(error_list) > 0){
  do.call(rbind,error_list)
} else {
  data.frame(file = character(),file_name = character(),
             file_type = character(),reason = character())
}
row.names(parse_errors) <- NULL

parsed <- list(
  class_beta_long = class_beta_long,
  class_beta_summary = class_beta_summary,
  file_inventory = file_inventory,
  resource_summary = resource_summary,
  parse_errors = parse_errors
)

write_outputs(parsed,output_file)

message("Read ",length(files)," .rda files from: ",input_root)
message("Parsed ",nrow(class_beta_long)," class-beta rows")
message("Saved: ",output_file)
