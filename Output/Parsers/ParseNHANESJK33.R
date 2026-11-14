#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Output","Routputs","RoutputsNHANESJK33")
}

output_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Output","parse_nhanes_jkn.rds")
}

reference_dir <- if (length(args) >= 3){
  args[[3]]
} else {
  "Data"
}

script_arg <- grep("^--file=",commandArgs(FALSE),value = TRUE)
if (length(script_arg) == 1L){
  script_file <- normalizePath(sub("^--file=","",script_arg),mustWork = TRUE)
  project_dir <- normalizePath(file.path(dirname(script_file),"..",".."),
                               mustWork = TRUE)
  local_library <- file.path(project_dir,"Rlib")
  if (dir.exists(local_library)){
    .libPaths(c(local_library,.libPaths()))
  }
}

EXPECTED_REPLICATES <- 1:61
EXPECTED_CLASSES <- 1:5
ESTIMATED_CLASSES <- 2:5
EXPECTED_DESIGN_DF <- 32L
CLASS_NAMES <- paste0("class_",ESTIMATED_CLASSES)
SURVIVAL_COVARIATE_NAMES <- c(
  "age","gender","race","overall_health","education","bmi",
  "diabetes","chd","chf","heart_attack","stroke","alcohol",
  "smoking","physical_function"
)
SURVIVAL_COEF_LENGTHS <- c(1L,2L,6L,5L,5L,4L,2L,2L,2L,2L,2L,5L,3L,2L)
SURVIVAL_PARAMETER_NAMES <- c(
  "age",
  CLASS_NAMES,
  unlist(lapply(seq.int(2L,length(SURVIVAL_COVARIATE_NAMES)),function(i){
    paste0(
      SURVIVAL_COVARIATE_NAMES[[i]],"_",
      seq.int(2L,SURVIVAL_COEF_LENGTHS[[i]])
    )
  }),use.names = FALSE)
)
stopifnot(length(SURVIVAL_PARAMETER_NAMES) == 34L)
REFERENCE_FILES <- c(
  joint = "JMHMMRandInitReducedFitMix5Seed33len96.rda",
  two_stage = "JMHMMNoSurvRandInitReducedFitMix5Seed33len96.rda"
)
FILE_PATTERN <- paste0(
  "^JMHMMJKn(Joint|TwoStage)ReducedFitMix5RefSeed33",
  "Rep([0-9]{3})len96\\.rda$"
)

if (!dir.exists(input_dir)){
  stop("NHANES JKn replicate directory does not exist: ",input_dir)
}

load_to_save <- function(file){
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file,envir = load_env)
  if (!identical(loaded_names,"to_save") ||
      !exists("to_save",envir = load_env,inherits = FALSE)){
    stop("File must contain exactly one object named to_save")
  }
  get("to_save",envir = load_env,inherits = FALSE)
}

assert_true <- function(condition,message){
  if (!isTRUE(condition)){
    stop(message,call. = FALSE)
  }
}

scalar_character <- function(x,name){
  assert_true(length(x) == 1L && !is.na(x),
              paste0(name," must be one nonmissing value"))
  value <- as.character(x[[1]])
  assert_true(nzchar(value),paste0(name," must not be empty"))
  value
}

scalar_integer <- function(x,name){
  assert_true(length(x) == 1L,
              paste0(name," must be one integer value"))
  value <- suppressWarnings(as.integer(x[[1]]))
  assert_true(!is.na(value),paste0(name," must be a finite integer"))
  value
}

scalar_numeric <- function(x,name){
  assert_true(length(x) == 1L,
              paste0(name," must be one numeric value"))
  value <- suppressWarnings(as.numeric(x[[1]]))
  assert_true(is.finite(value),paste0(name," must be finite"))
  value
}

same_numeric <- function(x,y,tolerance = 1e-12){
  length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
    max(abs(x - y)) <= tolerance
}

canonical_psu_table <- function(x){
  required <- c(
    "replicate_id","deleted_stratum","deleted_psu","stratum_psu_count",
    "retained_psu_multiplier","rscale"
  )
  assert_true(is.data.frame(x),"ordered_psu_table must be a data frame")
  assert_true(all(required %in% names(x)),
              "ordered_psu_table is missing required columns")

  out <- data.frame(
    replicate_id = suppressWarnings(as.integer(x$replicate_id)),
    deleted_stratum = as.character(x$deleted_stratum),
    deleted_psu = as.character(x$deleted_psu),
    stratum_psu_count = suppressWarnings(as.integer(x$stratum_psu_count)),
    retained_psu_multiplier = suppressWarnings(
      as.numeric(x$retained_psu_multiplier)
    ),
    rscale = suppressWarnings(as.numeric(x$rscale)),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$replicate_id),,drop = FALSE]
  rownames(out) <- NULL

  assert_true(nrow(out) == length(EXPECTED_REPLICATES),
              "ordered_psu_table must have 61 rows")
  assert_true(identical(out$replicate_id,EXPECTED_REPLICATES),
              "ordered_psu_table must contain replicate IDs 1 through 61")
  assert_true(all(!is.na(out$deleted_stratum) & nzchar(out$deleted_stratum)),
              "ordered_psu_table has invalid strata")
  assert_true(all(!is.na(out$deleted_psu) & nzchar(out$deleted_psu)),
              "ordered_psu_table has invalid PSUs")
  assert_true(all(out$stratum_psu_count %in% c(2L,3L)),
              "ordered_psu_table has unexpected PSU counts")
  assert_true(all(is.finite(out$retained_psu_multiplier) &
                    out$retained_psu_multiplier > 0),
              "ordered_psu_table has invalid PSU multipliers")
  assert_true(all(is.finite(out$rscale) & out$rscale > 0),
              "ordered_psu_table has invalid rscales")
  deleted_units <- paste(out$deleted_stratum,out$deleted_psu,sep = ":")
  assert_true(length(unique(deleted_units)) == length(EXPECTED_REPLICATES),
              "ordered_psu_table does not contain 61 unique deleted PSUs")
  out
}

same_psu_table <- function(x,y){
  isTRUE(all.equal(x,y,tolerance = 1e-12,check.attributes = FALSE))
}

parse_file_name <- function(file){
  file_name <- basename(file)
  match <- regexec(FILE_PATTERN,file_name,perl = TRUE)
  found <- regmatches(file_name,match)[[1]]
  if (length(found) != 3L){
    return(NULL)
  }
  model_type <- if (found[[2]] == "Joint") "joint" else "two_stage"
  data.frame(
    file = normalizePath(file,mustWork = TRUE),
    file_name = file_name,
    model_type = model_type,
    replicate_id = as.integer(found[[3]]),
    stringsAsFactors = FALSE
  )
}

empty_replicate_diagnostic <- function(file,model_type,replicate_id,
                                       error = NA_character_){
  data.frame(
    file = file,
    model_type = model_type,
    replicate_id = as.integer(replicate_id),
    deleted_stratum = NA_character_,
    deleted_psu = NA_character_,
    scale = NA_real_,
    rscale = NA_real_,
    design_df = NA_integer_,
    design_fingerprint = NA_character_,
    converged = FALSE,
    completed_iterations = NA_integer_,
    alignment_score = NA_real_,
    alignment_margin = NA_real_,
    valid = FALSE,
    error = error,
    stringsAsFactors = FALSE
  )
}

parse_replicate <- function(file_info){
  file <- file_info$file[[1]]
  expected_model <- file_info$model_type[[1]]
  expected_replicate <- file_info$replicate_id[[1]]
  to_save <- load_to_save(file)

  assert_true(is.list(to_save),"to_save must be a list")
  settings <- to_save[["settings"]]
  diagnostics <- to_save[["diagnostics"]]
  jackknife <- to_save[["jackknife"]]
  assert_true(is.list(settings),"Missing settings")
  assert_true(is.list(diagnostics),"Missing diagnostics")
  assert_true(is.list(jackknife),"Missing jackknife metadata")

  saved_model <- scalar_character(settings$model_type,"settings$model_type")
  data_source <- scalar_character(settings$data_source,"settings$data_source")
  fit_mix_num <- scalar_integer(settings$fit_mix_num,"settings$fit_mix_num")
  period_len <- scalar_integer(settings$period_len,"settings$period_len")
  resampling_method <- scalar_character(
    settings$resampling_method,"settings$resampling_method"
  )
  settings_replicate <- scalar_integer(
    settings$replicate_id,"settings$replicate_id"
  )
  sim_num <- scalar_integer(settings$sim_num,"settings$sim_num")
  settings_reference_seed <- scalar_integer(
    settings$reference_seed,"settings$reference_seed"
  )

  assert_true(identical(saved_model,expected_model),
              "Saved model type does not match filename")
  assert_true(identical(data_source,"nhanes"),
              "Replicate is not an NHANES fit")
  assert_true(fit_mix_num == 5L,"Replicate is not a five-class fit")
  assert_true(period_len == 96L,"Replicate does not use 96 periods")
  assert_true(identical(resampling_method,"nhanes_jkn"),
              "Replicate does not use nhanes_jkn resampling")
  assert_true(settings_replicate == expected_replicate,
              "settings$replicate_id does not match filename")
  assert_true(sim_num == expected_replicate,
              "settings$sim_num does not match filename")
  assert_true(settings_reference_seed == 33L,
              "settings$reference_seed is not 33")

  method <- scalar_character(jackknife$method,"jackknife$method")
  replicate_id <- scalar_integer(
    jackknife$replicate_id,"jackknife$replicate_id"
  )
  deleted_stratum <- scalar_character(
    jackknife$deleted_stratum,"jackknife$deleted_stratum"
  )
  deleted_psu <- scalar_character(
    jackknife$deleted_psu,"jackknife$deleted_psu"
  )
  scale <- scalar_numeric(jackknife$scale,"jackknife$scale")
  rscale <- scalar_numeric(jackknife$rscale,"jackknife$rscale")
  design_df <- scalar_integer(jackknife$design_df,"jackknife$design_df")
  design_fingerprint <- scalar_character(
    jackknife$design_fingerprint,"jackknife$design_fingerprint"
  )
  reference_seed <- scalar_integer(
    jackknife$reference_seed,"jackknife$reference_seed"
  )
  reference_file <- scalar_character(
    jackknife$reference_file,"jackknife$reference_file"
  )
  alignment_score <- scalar_numeric(
    jackknife$alignment_score,"jackknife$alignment_score"
  )
  alignment_margin <- scalar_numeric(
    jackknife$alignment_margin,"jackknife$alignment_margin"
  )

  assert_true(identical(method,"NHANES_JKn"),
              "Unexpected jackknife method")
  assert_true(replicate_id == expected_replicate,
              "jackknife$replicate_id does not match filename")
  assert_true(scale > 0,"Jackknife scale must be positive")
  assert_true(rscale > 0,"Jackknife rscale must be positive")
  assert_true(design_df == EXPECTED_DESIGN_DF,
              "Jackknife design degrees of freedom must equal 32")
  assert_true(reference_seed == 33L,"Jackknife reference seed is not 33")
  assert_true(basename(reference_file) == REFERENCE_FILES[[expected_model]],
              "Jackknife reference filename is incorrect")
  assert_true(alignment_score >= 0,"Alignment score must be nonnegative")
  assert_true(alignment_margin > 0,
              "Class alignment is tied or has an invalid margin")

  class_mapping <- jackknife$class_mapping
  assert_true(is.data.frame(class_mapping),
              "class_mapping must be a data frame")
  assert_true(all(c("reference_class","replicate_class") %in%
                    names(class_mapping)),
              "class_mapping is missing required columns")
  reference_classes <- suppressWarnings(
    as.integer(class_mapping$reference_class)
  )
  replicate_classes <- suppressWarnings(
    as.integer(class_mapping$replicate_class)
  )
  assert_true(nrow(class_mapping) == length(EXPECTED_CLASSES) &&
                identical(sort(reference_classes),EXPECTED_CLASSES) &&
                identical(sort(replicate_classes),EXPECTED_CLASSES),
              "class_mapping is not a one-to-one five-class mapping")

  aligned_beta <- suppressWarnings(as.numeric(jackknife$aligned_beta_vec))
  assert_true(length(aligned_beta) == length(EXPECTED_CLASSES),
              "aligned_beta_vec must have length five")
  assert_true(all(is.finite(aligned_beta)),
              "aligned_beta_vec contains nonfinite values")
  assert_true(abs(aligned_beta[[1]]) <= 1e-10,
              "The aligned class-1 coefficient is not zero")
  names(aligned_beta) <- paste0("class_",EXPECTED_CLASSES)

  aligned_survival_raw <- jackknife$aligned_survival_vector
  aligned_survival <- suppressWarnings(as.numeric(aligned_survival_raw))
  assert_true(length(aligned_survival) == length(SURVIVAL_PARAMETER_NAMES),
              "aligned_survival_vector must have length 34")
  assert_true(all(is.finite(aligned_survival)),
              "aligned_survival_vector contains nonfinite values")
  assert_true(identical(names(aligned_survival_raw),
                        SURVIVAL_PARAMETER_NAMES),
              "aligned_survival_vector has unexpected parameter names")
  names(aligned_survival) <- SURVIVAL_PARAMETER_NAMES
  assert_true(
    same_numeric(aligned_survival[CLASS_NAMES],
                 aligned_beta[paste0("class_",ESTIMATED_CLASSES)]),
    "Aligned class coefficients differ between saved jackknife vectors"
  )

  psu_table <- canonical_psu_table(jackknife$ordered_psu_table)
  current_psu <- psu_table[
    psu_table$replicate_id == expected_replicate,,drop = FALSE
  ]
  assert_true(current_psu$deleted_stratum[[1]] == deleted_stratum &&
                current_psu$deleted_psu[[1]] == deleted_psu,
              "Deleted PSU does not match ordered_psu_table")
  assert_true(same_numeric(current_psu$rscale[[1]],rscale),
              "Replicate rscale does not match ordered_psu_table")

  convergence <- diagnostics$convergence
  assert_true(is.list(convergence),"Missing convergence diagnostics")
  converged <- isTRUE(convergence$converged)
  completed_iterations <- scalar_integer(
    convergence$completed_iterations,
    "diagnostics$convergence$completed_iterations"
  )
  assert_true(converged,"Replicate model did not converge")

  diagnostic <- data.frame(
    file = file,
    model_type = expected_model,
    replicate_id = expected_replicate,
    deleted_stratum = deleted_stratum,
    deleted_psu = deleted_psu,
    scale = scale,
    rscale = rscale,
    design_df = design_df,
    design_fingerprint = design_fingerprint,
    converged = converged,
    completed_iterations = completed_iterations,
    alignment_score = alignment_score,
    alignment_margin = alignment_margin,
    valid = TRUE,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  list(
    diagnostic = diagnostic,
    aligned_beta = aligned_beta,
    aligned_survival = aligned_survival,
    ordered_psu_table = psu_table
  )
}

validate_reference <- function(model_type,file){
  assert_true(file.exists(file),paste0("Reference file does not exist: ",file))
  to_save <- load_to_save(file)
  settings <- to_save[["settings"]]
  est_params <- to_save[["est_params"]]
  diagnostics <- to_save[["diagnostics"]]
  assert_true(is.list(settings),"Reference is missing settings")
  assert_true(is.list(est_params),"Reference is missing est_params")
  assert_true(is.list(diagnostics),"Reference is missing diagnostics")

  assert_true(
    identical(scalar_character(settings$model_type,"settings$model_type"),
              model_type),
    "Reference model type is incorrect"
  )
  assert_true(
    identical(scalar_character(settings$data_source,"settings$data_source"),
              "nhanes"),
    "Reference data source is not NHANES"
  )
  assert_true(scalar_integer(settings$sim_num,"settings$sim_num") == 33L,
              "Reference sim_num is not Seed 33")
  assert_true(scalar_integer(settings$fit_mix_num,"settings$fit_mix_num") == 5L,
              "Reference is not a five-class fit")
  assert_true(scalar_integer(settings$period_len,"settings$period_len") == 96L,
              "Reference does not use 96 periods")
  assert_true(isTRUE(diagnostics$convergence$converged),
              "Reference model did not converge")

  beta_vec <- suppressWarnings(as.numeric(est_params$beta_vec))
  assert_true(length(beta_vec) == length(EXPECTED_CLASSES),
              "Reference beta_vec must have length five")
  assert_true(all(is.finite(beta_vec)),
              "Reference beta_vec contains nonfinite values")
  assert_true(abs(beta_vec[[1]]) <= 1e-10,
              "Reference class-1 coefficient is not zero")
  beta <- beta_vec[ESTIMATED_CLASSES]
  names(beta) <- CLASS_NAMES

  surv_coef <- est_params$surv_coef
  assert_true(is.list(surv_coef) &&
                length(surv_coef) == length(SURVIVAL_COEF_LENGTHS),
              "Reference surv_coef has an unexpected structure")
  observed_lengths <- vapply(surv_coef,length,integer(1))
  assert_true(identical(observed_lengths,SURVIVAL_COEF_LENGTHS),
              "Reference surv_coef has unexpected category counts")
  surv_coef <- lapply(surv_coef,function(x){
    suppressWarnings(as.numeric(x))
  })
  assert_true(all(is.finite(unlist(surv_coef,use.names = FALSE))),
              "Reference surv_coef contains nonfinite values")
  survival <- c(
    surv_coef[[1]][[1]],
    beta,
    unlist(lapply(surv_coef[-1],function(x) x[-1]),use.names = FALSE)
  )
  assert_true(length(survival) == length(SURVIVAL_PARAMETER_NAMES),
              "Reference survival vector must have length 34")
  names(survival) <- SURVIVAL_PARAMETER_NAMES

  list(
    beta = beta,
    survival = survival,
    diagnostic = data.frame(
      model_type = model_type,
      file = normalizePath(file,mustWork = TRUE),
      valid = TRUE,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  )
}

error_rows <- list()
add_error <- function(stage,message,model_type = NA_character_,
                      replicate_id = NA_integer_,file = NA_character_){
  error_rows[[length(error_rows) + 1L]] <<- data.frame(
    stage = stage,
    model_type = model_type,
    replicate_id = as.integer(replicate_id),
    file = file,
    message = as.character(message),
    stringsAsFactors = FALSE
  )
}

all_rda_files <- list.files(
  input_dir,pattern = "\\.rda$",full.names = TRUE,ignore.case = TRUE
)
file_info_list <- lapply(all_rda_files,parse_file_name)
recognized <- !vapply(file_info_list,is.null,logical(1))
ignored_files <- normalizePath(
  all_rda_files[!recognized],mustWork = TRUE
)
file_info_list <- file_info_list[recognized]
file_info <- if (length(file_info_list) > 0){
  do.call(rbind,file_info_list)
} else {
  data.frame(
    file = character(),file_name = character(),model_type = character(),
    replicate_id = integer(),stringsAsFactors = FALSE
  )
}

for (model_type in names(REFERENCE_FILES)){
  model_ids <- file_info$replicate_id[file_info$model_type == model_type]
  missing_ids <- setdiff(EXPECTED_REPLICATES,model_ids)
  duplicate_ids <- unique(model_ids[duplicated(model_ids)])
  unexpected_ids <- setdiff(model_ids,EXPECTED_REPLICATES)
  if (length(missing_ids) > 0){
    add_error(
      "manifest",
      paste0("Missing replicate IDs: ",paste(missing_ids,collapse = ", ")),
      model_type
    )
  }
  if (length(duplicate_ids) > 0){
    add_error(
      "manifest",
      paste0("Duplicate replicate IDs: ",paste(duplicate_ids,collapse = ", ")),
      model_type
    )
  }
  if (length(unexpected_ids) > 0){
    add_error(
      "manifest",
      paste0("Unexpected replicate IDs: ",
             paste(unexpected_ids,collapse = ", ")),
      model_type
    )
  }
}

parsed_replicates <- vector("list",nrow(file_info))
replicate_diagnostics <- vector("list",nrow(file_info))
for (i in seq_len(nrow(file_info))){
  info <- file_info[i,,drop = FALSE]
  parsed <- tryCatch(
    parse_replicate(info),
    error = function(e){
      add_error(
        "replicate",conditionMessage(e),info$model_type[[1]],
        info$replicate_id[[1]],info$file[[1]]
      )
      NULL
    }
  )
  parsed_replicates[[i]] <- parsed
  replicate_diagnostics[[i]] <- if (is.null(parsed)){
    message_text <- error_rows[[length(error_rows)]]$message[[1]]
    empty_replicate_diagnostic(
      info$file[[1]],info$model_type[[1]],info$replicate_id[[1]],
      message_text
    )
  } else {
    parsed$diagnostic
  }
}

replicate_diagnostics <- if (length(replicate_diagnostics) > 0){
  do.call(rbind,replicate_diagnostics)
} else {
  empty_replicate_diagnostic(
    NA_character_,NA_character_,NA_integer_,NA_character_
  )[0,,drop = FALSE]
}
rownames(replicate_diagnostics) <- NULL

model_data <- list()
for (model_type in names(REFERENCE_FILES)){
  indices <- which(file_info$model_type == model_type)
  model_parsed <- parsed_replicates[indices]
  complete <- length(indices) == length(EXPECTED_REPLICATES) &&
    all(!vapply(model_parsed,is.null,logical(1)))
  if (!complete){
    next
  }
  model_parsed <- model_parsed[
    order(file_info$replicate_id[indices])
  ]
  ids <- vapply(
    model_parsed,function(x) x$diagnostic$replicate_id[[1]],integer(1)
  )
  model_result <- tryCatch({
    assert_true(identical(ids,EXPECTED_REPLICATES),
                "Parsed replicate IDs are not ordered 1 through 61")
    first_table <- model_parsed[[1]]$ordered_psu_table
    assert_true(
      all(vapply(
        model_parsed,
        function(x) same_psu_table(x$ordered_psu_table,first_table),
        logical(1)
      )),
      "Replicates do not contain the same ordered_psu_table"
    )
    scales <- vapply(
      model_parsed,function(x) x$diagnostic$scale[[1]],numeric(1)
    )
    rscales <- vapply(
      model_parsed,function(x) x$diagnostic$rscale[[1]],numeric(1)
    )
    design_dfs <- vapply(
      model_parsed,function(x) x$diagnostic$design_df[[1]],integer(1)
    )
    fingerprints <- vapply(
      model_parsed,function(x) x$diagnostic$design_fingerprint[[1]],
      character(1)
    )
    deleted_units <- vapply(
      model_parsed,
      function(x) paste(
        x$diagnostic$deleted_stratum[[1]],x$diagnostic$deleted_psu[[1]],
        sep = ":"
      ),
      character(1)
    )
    assert_true(length(unique(scales)) == 1L,
                "Replicates do not share one jackknife scale")
    assert_true(all(design_dfs == EXPECTED_DESIGN_DF),
                "Replicates do not all use 32 design degrees of freedom")
    assert_true(length(unique(fingerprints)) == 1L,
                "Replicates have inconsistent design fingerprints")
    assert_true(length(unique(deleted_units)) == length(EXPECTED_REPLICATES),
                "Replicates do not delete 61 unique PSUs")
    assert_true(same_numeric(rscales,first_table$rscale),
                "Ordered replicate rscales do not match the PSU table")

    estimates <- do.call(
      rbind,lapply(model_parsed,function(x){
        x$aligned_beta[ESTIMATED_CLASSES]
      })
    )
    rownames(estimates) <- sprintf("replicate_%03d",EXPECTED_REPLICATES)
    colnames(estimates) <- CLASS_NAMES
    survival_estimates <- do.call(
      rbind,lapply(model_parsed,function(x) x$aligned_survival)
    )
    rownames(survival_estimates) <- sprintf(
      "replicate_%03d",EXPECTED_REPLICATES
    )
    colnames(survival_estimates) <- SURVIVAL_PARAMETER_NAMES
    list(
      estimates = estimates,
      survival_estimates = survival_estimates,
      scale = scales[[1]],
      rscales = rscales,
      design_df = design_dfs[[1]],
      design_fingerprint = fingerprints[[1]],
      ordered_psu_table = first_table
    )
  },error = function(e){
    add_error("design",conditionMessage(e),model_type)
    NULL
  })
  if (!is.null(model_result)){
    model_data[[model_type]] <- model_result
  }
}

if (all(names(REFERENCE_FILES) %in% names(model_data))){
  cross_model_error <- tryCatch({
    joint <- model_data$joint
    two_stage <- model_data$two_stage
    assert_true(joint$design_fingerprint == two_stage$design_fingerprint,
                "Joint and two-stage design fingerprints differ")
    assert_true(same_numeric(joint$scale,two_stage$scale),
                "Joint and two-stage jackknife scales differ")
    assert_true(same_numeric(joint$rscales,two_stage$rscales),
                "Joint and two-stage rscales differ")
    assert_true(joint$design_df == two_stage$design_df,
                "Joint and two-stage design degrees of freedom differ")
    assert_true(same_psu_table(joint$ordered_psu_table,
                               two_stage$ordered_psu_table),
                "Replicate IDs delete different PSUs across model types")
    NULL
  },error = function(e) e)
  if (inherits(cross_model_error,"error")){
    add_error("cross_model",conditionMessage(cross_model_error))
  }
}

reference_results <- list()
reference_diagnostics <- list()
for (model_type in names(REFERENCE_FILES)){
  reference_file <- file.path(reference_dir,REFERENCE_FILES[[model_type]])
  reference <- tryCatch(
    validate_reference(model_type,reference_file),
    error = function(e){
      add_error(
        "reference",conditionMessage(e),model_type,NA_integer_,
        reference_file
      )
      NULL
    }
  )
  if (is.null(reference)){
    reference_diagnostics[[model_type]] <- data.frame(
      model_type = model_type,
      file = reference_file,
      valid = FALSE,
      error = error_rows[[length(error_rows)]]$message[[1]],
      stringsAsFactors = FALSE
    )
  } else {
    reference_results[[model_type]] <- reference
    reference_diagnostics[[model_type]] <- reference$diagnostic
  }
}
reference_diagnostics <- do.call(rbind,reference_diagnostics)
rownames(reference_diagnostics) <- NULL

if (!requireNamespace("survey",quietly = TRUE)){
  add_error(
    "dependency",
    "The survey package is required to calculate NHANES JKn covariance"
  )
}

parse_errors <- function(){
  if (length(error_rows) == 0){
    return(data.frame(
      stage = character(),model_type = character(),
      replicate_id = integer(),file = character(),message = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind,error_rows)
  rownames(out) <- NULL
  out
}

calculate_vcov <- function(thetas,coef,scale,rscales,parameter_names){
  assert_true(ncol(thetas) == length(coef) &&
                length(coef) == length(parameter_names),
              "Covariance inputs have incompatible dimensions")
  vcov_survey <- as.matrix(survey::svrVar(
    thetas = thetas,
    scale = scale,
    rscales = rscales,
    mse = TRUE,
    coef = coef
  ))
  dimnames(vcov_survey) <- list(parameter_names,parameter_names)

  centered <- sweep(thetas,2,coef,"-")
  vcov_manual <- matrix(0,nrow = length(coef),ncol = length(coef))
  for (i in seq_len(nrow(centered))){
    vcov_manual <- vcov_manual +
      scale * rscales[[i]] * tcrossprod(centered[i,])
  }
  dimnames(vcov_manual) <- list(parameter_names,parameter_names)
  assert_true(
    max(abs(vcov_survey - vcov_manual)) <= 1e-10,
    "survey::svrVar does not match the direct JKn quadratic calculation"
  )
  vcov_survey
}

coefficient_summary <- data.frame(
  model_type = character(),class_index = integer(),
  param_name = character(),estimate = numeric(),se = numeric(),
  design_df = integer(),critical_value = numeric(),
  ci_lower = numeric(),ci_upper = numeric(),hr = numeric(),
  hr_ci_lower = numeric(),hr_ci_upper = numeric(),
  stringsAsFactors = FALSE
)
vcov_results <- list()
joint_survival_summary <- data.frame(
  model_type = character(),parameter_index = integer(),
  param_name = character(),estimate = numeric(),se = numeric(),
  design_df = integer(),critical_value = numeric(),
  ci_lower = numeric(),ci_upper = numeric(),hr = numeric(),
  hr_ci_lower = numeric(),hr_ci_upper = numeric(),
  stringsAsFactors = FALSE
)
joint_survival_vcov <- matrix(numeric(),nrow = 0L,ncol = 0L)

inference_ready <- length(error_rows) == 0L &&
  all(names(REFERENCE_FILES) %in% names(model_data)) &&
  all(names(REFERENCE_FILES) %in% names(reference_results))

if (inference_ready){
  inference_error <- tryCatch({
    summary_rows <- list()
    vcov_candidate <- list()
    for (model_type in names(REFERENCE_FILES)){
      model <- model_data[[model_type]]
      beta <- reference_results[[model_type]]$beta
      vcov_jkn <- calculate_vcov(
        model$estimates,beta,model$scale,model$rscales,CLASS_NAMES
      )
      se <- sqrt(diag(vcov_jkn))
      assert_true(all(is.finite(se)),
                  "Calculated standard errors are nonfinite")
      critical <- stats::qt(0.975,df = model$design_df)
      lower <- beta - critical * se
      upper <- beta + critical * se
      summary_rows[[model_type]] <- data.frame(
        model_type = model_type,
        class_index = ESTIMATED_CLASSES,
        param_name = CLASS_NAMES,
        estimate = as.numeric(beta),
        se = as.numeric(se),
        design_df = model$design_df,
        critical_value = critical,
        ci_lower = as.numeric(lower),
        ci_upper = as.numeric(upper),
        hr = exp(as.numeric(beta)),
        hr_ci_lower = exp(as.numeric(lower)),
        hr_ci_upper = exp(as.numeric(upper)),
        stringsAsFactors = FALSE
      )
      vcov_candidate[[model_type]] <- vcov_jkn
    }

    joint_model <- model_data$joint
    joint_reference <- reference_results$joint$survival
    joint_survival_vcov_candidate <- calculate_vcov(
      joint_model$survival_estimates,
      joint_reference,
      joint_model$scale,
      joint_model$rscales,
      SURVIVAL_PARAMETER_NAMES
    )
    assert_true(
      max(abs(
        joint_survival_vcov_candidate[CLASS_NAMES,CLASS_NAMES,drop = FALSE] -
          vcov_candidate$joint
      )) <= 1e-10,
      "Class covariance does not match its full survival covariance block"
    )
    joint_survival_se <- sqrt(diag(joint_survival_vcov_candidate))
    assert_true(all(is.finite(joint_survival_se)),
                "Full joint survival standard errors are nonfinite")
    joint_critical <- stats::qt(0.975,df = joint_model$design_df)
    joint_lower <- joint_reference - joint_critical * joint_survival_se
    joint_upper <- joint_reference + joint_critical * joint_survival_se
    joint_survival_summary_candidate <- data.frame(
      model_type = "joint",
      parameter_index = seq_along(SURVIVAL_PARAMETER_NAMES),
      param_name = SURVIVAL_PARAMETER_NAMES,
      estimate = as.numeric(joint_reference),
      se = as.numeric(joint_survival_se),
      design_df = joint_model$design_df,
      critical_value = joint_critical,
      ci_lower = as.numeric(joint_lower),
      ci_upper = as.numeric(joint_upper),
      hr = exp(as.numeric(joint_reference)),
      hr_ci_lower = exp(as.numeric(joint_lower)),
      hr_ci_upper = exp(as.numeric(joint_upper)),
      stringsAsFactors = FALSE
    )
    list(
      coefficient_summary = do.call(rbind,summary_rows),
      vcov = vcov_candidate,
      joint_survival_summary = joint_survival_summary_candidate,
      joint_survival_vcov = joint_survival_vcov_candidate
    )
  },error = function(e) e)

  if (inherits(inference_error,"error")){
    add_error("covariance",conditionMessage(inference_error))
    inference_ready <- FALSE
  } else {
    coefficient_summary <- inference_error$coefficient_summary
    rownames(coefficient_summary) <- NULL
    vcov_results <- inference_error$vcov
    joint_survival_summary <- inference_error$joint_survival_summary
    rownames(joint_survival_summary) <- NULL
    joint_survival_vcov <- inference_error$joint_survival_vcov
  }
}

reference_estimates <- lapply(reference_results,function(x) x$beta)
replicate_estimates <- lapply(model_data,function(x) x$estimates)
design_metadata <- lapply(model_data,function(x){
  list(
    method = "NHANES_JKn",
    scale = x$scale,
    rscales = x$rscales,
    design_df = x$design_df,
    critical_value = stats::qt(0.975,df = x$design_df),
    design_fingerprint = x$design_fingerprint,
    ordered_psu_table = x$ordered_psu_table
  )
})

out <- list(
  inference_available = isTRUE(inference_ready),
  coefficient_summary = if (inference_ready){
    coefficient_summary
  } else {
    coefficient_summary[0,,drop = FALSE]
  },
  replicate_estimates = replicate_estimates,
  vcov = if (inference_ready) vcov_results else list(),
  reference_estimates = reference_estimates,
  joint_survival = list(
    coefficient_summary = if (inference_ready){
      joint_survival_summary
    } else {
      joint_survival_summary[0,,drop = FALSE]
    },
    replicate_estimates = if ("joint" %in% names(model_data)){
      model_data$joint$survival_estimates
    } else {
      matrix(numeric(),nrow = 0L,ncol = 0L)
    },
    vcov = if (inference_ready){
      joint_survival_vcov
    } else {
      matrix(numeric(),nrow = 0L,ncol = 0L)
    },
    reference_estimates = if ("joint" %in% names(reference_results)){
      reference_results$joint$survival
    } else {
      numeric()
    }
  ),
  replicate_diagnostics = replicate_diagnostics,
  reference_diagnostics = reference_diagnostics,
  design_metadata = design_metadata,
  parse_errors = parse_errors(),
  ignored_files = ignored_files,
  input_dir = normalizePath(input_dir,mustWork = TRUE),
  reference_dir = normalizePath(reference_dir,mustWork = FALSE),
  output_file = normalizePath(output_file,mustWork = FALSE),
  generated_at = Sys.time()
)

dir.create(dirname(output_file),recursive = TRUE,showWarnings = FALSE)
saveRDS(out,output_file)

message("Read ",nrow(file_info)," NHANES JKn replicate files from: ",
        input_dir)
message("Valid joint replicates: ",
        sum(replicate_diagnostics$model_type == "joint" &
              replicate_diagnostics$valid))
message("Valid two-stage replicates: ",
        sum(replicate_diagnostics$model_type == "two_stage" &
              replicate_diagnostics$valid))
message("Saved: ",output_file)

if (!inference_ready){
  stop(
    "NHANES JKn inference was withheld because ",
    nrow(out$parse_errors)," validation error(s) were found. ",
    "Inspect parse_errors in ",output_file,".",
    call. = FALSE
  )
}

print(coefficient_summary)
