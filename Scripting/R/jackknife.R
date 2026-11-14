NHANESJKnReferenceFilename <- function(model_type){
  if (!model_type %in% names(MODEL_TYPE_CODES)){
    stop("Unknown model_type for NHANES JKn reference: ",model_type)
  }
  prefix <- if (model_type == "joint"){
    "JMHMM"
  } else {
    "JMHMMNoSurv"
  }
  paste0(
    prefix,"RandInitReducedFitMix5Seed",NHANES_JKN_REFERENCE_SEED,
    "len",DEFAULT_PERIODS_PER_DAY,".rda"
  )
}

NHANESJKnReferencePath <- function(model_type,data_dir = "Data"){
  file.path(data_dir,NHANESJKnReferenceFilename(model_type))
}

ValidateNHANESJKnReference <- function(to_save,reference_file,model_type){
  expected_name <- NHANESJKnReferenceFilename(model_type)
  if (!file.exists(reference_file)){
    stop("NHANES JKn reference file does not exist: ",reference_file)
  }
  if (basename(reference_file) != expected_name){
    stop(
      "Unexpected NHANES JKn reference filename: ",basename(reference_file),
      "; expected ",expected_name
    )
  }

  validate_saved_results(
    to_save,
    required_sections = c("est_params","settings"),
    source_name = reference_file
  )
  saved_settings <- get_saved_section(to_save,"settings")
  expected_settings <- list(
    fit_mix_num = 5,
    model_type = model_type,
    data_source = DATA_SOURCE[["nhanes"]],
    sim_num = NHANES_JKN_REFERENCE_SEED,
    period_len = DEFAULT_PERIODS_PER_DAY
  )
  for (setting_name in names(expected_settings)){
    saved_value <- saved_settings[[setting_name]]
    if (is.null(saved_value) && setting_name == "period_len"){
      saved_value <- saved_settings$periods_per_day
    }
    if (is.null(saved_value) ||
        !identical(as.character(saved_value),
                   as.character(expected_settings[[setting_name]]))){
      stop(
        "NHANES JKn reference setting mismatch for ",setting_name,
        ": saved ",if (is.null(saved_value)) "<missing>" else saved_value,
        ", expected ",expected_settings[[setting_name]]
      )
    }
  }

  est_params <- get_saved_section(to_save,"est_params")
  beta_vec <- as.numeric(get_saved_param(est_params,"beta_vec"))
  re_prob <- as.matrix(get_saved_param(est_params,"re_prob"))
  if (length(beta_vec) != 5L || any(!is.finite(beta_vec))){
    stop("NHANES JKn reference beta_vec must contain five finite values")
  }
  if (!all(dim(re_prob) == c(NHANES_JKN_EXPECTED_PARTICIPANTS,5L)) ||
      any(!is.finite(re_prob))){
    stop(
      "NHANES JKn reference re_prob must be a finite ",
      NHANES_JKN_EXPECTED_PARTICIPANTS," x 5 matrix"
    )
  }
  if (max(abs(rowSums(re_prob) - 1)) > 1e-6){
    stop("NHANES JKn reference re_prob rows do not sum to one")
  }

  invisible(TRUE)
}

NHANESJKnFingerprint <- function(design_data,ordered_psu_table){
  design_text <- paste(
    paste(
      design_data$SEQN,
      as.character(design_data$SDMVSTRA),
      as.character(design_data$SDMVPSU),
      format(design_data$weight,digits = 17,scientific = TRUE),
      sep = ":"
    ),
    collapse = "|"
  )
  psu_text <- paste(
    paste(
      ordered_psu_table$replicate_id,
      ordered_psu_table$deleted_stratum,
      ordered_psu_table$deleted_psu,
      format(ordered_psu_table$rscale,digits = 17),
      sep = ":"
    ),
    collapse = "|"
  )
  bytes <- as.integer(charToRaw(paste(design_text,psu_text,sep = "||")))
  hash <- 0
  for (byte in bytes){
    hash <- (hash * 131 + byte) %% 2147483647
  }
  paste0("nhanes-jkn-",sprintf("%08x",as.integer(hash)))
}

BuildNHANESJKnDesign <- function(id,sweights_vec,replicate_id){
  if (!requireNamespace("survey",quietly = TRUE)){
    stop("The survey package is required for NHANES JKn runs")
  }
  required_columns <- c("SEQN","SDMVSTRA","SDMVPSU")
  missing_columns <- setdiff(required_columns,names(id))
  if (length(missing_columns) > 0){
    stop(
      "NHANES ID data are missing JKn fields: ",
      paste(missing_columns,collapse = ", ")
    )
  }
  if (length(sweights_vec) != nrow(id) || any(!is.finite(sweights_vec)) ||
      any(sweights_vec <= 0)){
    stop("NHANES JKn base weights must be finite, positive, and match id")
  }
  replicate_id <- as.integer(replicate_id)
  if (length(replicate_id) != 1 || is.na(replicate_id) ||
      replicate_id < 1 || replicate_id > NHANES_JKN_EXPECTED_PSUS){
    stop("NHANES JKn replicate_id must be between 1 and 61")
  }

  design_data <- data.frame(
    SEQN = id$SEQN,
    SDMVSTRA = factor(id$SDMVSTRA),
    SDMVPSU = factor(id$SDMVPSU),
    weight = as.numeric(sweights_vec)
  )
  if (anyNA(design_data) || anyDuplicated(design_data$SEQN)){
    stop("NHANES JKn design fields must be complete with unique SEQN values")
  }

  nhanes_design <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~weight,
    nest = TRUE,
    data = design_data
  )
  jkn_design <- survey::as.svrepdesign(
    nhanes_design,
    type = "JKn",
    mse = TRUE
  )
  replicate_weight_matrix <- stats::weights(
    jkn_design,type = "analysis"
  )
  psu_count <- nrow(unique(design_data[c("SDMVSTRA","SDMVPSU")]))
  design_df <- survey::degf(jkn_design)

  expected_checks <- c(
    participants = nrow(design_data) == NHANES_JKN_EXPECTED_PARTICIPANTS,
    strata = nlevels(design_data$SDMVSTRA) == NHANES_JKN_EXPECTED_STRATA,
    psus = psu_count == NHANES_JKN_EXPECTED_PSUS,
    replicate_columns = ncol(replicate_weight_matrix) ==
      NHANES_JKN_EXPECTED_PSUS,
    design_df = design_df == NHANES_JKN_EXPECTED_DF
  )
  if (!all(expected_checks)){
    stop(
      "NHANES JKn design validation failed: ",
      paste(names(expected_checks)[!expected_checks],collapse = ", ")
    )
  }

  ordered_psu_rows <- lapply(seq_len(ncol(replicate_weight_matrix)),
    function(current_replicate){
      current_weights <- replicate_weight_matrix[,current_replicate]
      zero_pair <- unique(design_data[
        current_weights <= sqrt(.Machine$double.eps),
        c("SDMVSTRA","SDMVPSU"),drop = FALSE
      ])
      if (nrow(zero_pair) != 1L){
        stop(
          "JKn replicate ",current_replicate,
          " does not delete exactly one stratum-PSU pair"
        )
      }
      deleted_stratum <- as.character(zero_pair$SDMVSTRA[[1]])
      deleted_psu <- as.character(zero_pair$SDMVPSU[[1]])
      in_stratum <- as.character(design_data$SDMVSTRA) == deleted_stratum
      stratum_psus <- unique(as.character(
        design_data$SDMVPSU[in_stratum]
      ))
      m <- length(stratum_psus)
      expected_multiplier <- m / (m - 1)
      retained_in_stratum <- in_stratum & current_weights > 0
      outside_stratum <- !in_stratum
      ratio <- current_weights / design_data$weight
      if (any(abs(ratio[retained_in_stratum] - expected_multiplier) > 1e-8) ||
          any(abs(ratio[outside_stratum] - 1) > 1e-8)){
        stop("Unexpected JKn weight multiplier in replicate ",current_replicate)
      }
      data.frame(
        replicate_id = current_replicate,
        deleted_stratum = deleted_stratum,
        deleted_psu = deleted_psu,
        stratum_psu_count = m,
        retained_psu_multiplier = expected_multiplier,
        rscale = as.numeric(jkn_design$rscales[[current_replicate]]),
        stringsAsFactors = FALSE
      )
    }
  )
  ordered_psu_table <- do.call(rbind,ordered_psu_rows)
  expected_rscales <- (ordered_psu_table$stratum_psu_count - 1) /
    ordered_psu_table$stratum_psu_count
  if (any(abs(ordered_psu_table$rscale - expected_rscales) > 1e-12)){
    stop("JKn rscales do not match the deleted stratum PSU counts")
  }

  replicate_weights <- as.numeric(replicate_weight_matrix[,replicate_id])
  keep <- replicate_weights > sqrt(.Machine$double.eps)
  selected_psu <- ordered_psu_table[replicate_id,,drop = FALSE]
  deleted_rows <-
    as.character(design_data$SDMVSTRA) == selected_psu$deleted_stratum &
    as.character(design_data$SDMVPSU) == selected_psu$deleted_psu
  if (!identical(which(!keep),which(deleted_rows))){
    stop("JKn keep vector does not match the selected deleted PSU")
  }

  list(
    replicate_id = replicate_id,
    keep = keep,
    replicate_weights = replicate_weights,
    base_weights = design_data$weight,
    reference_seqn = as.character(design_data$SEQN),
    reference_row_index = which(keep),
    deleted_stratum = selected_psu$deleted_stratum[[1]],
    deleted_psu = selected_psu$deleted_psu[[1]],
    scale = as.numeric(jkn_design$scale),
    rscale = selected_psu$rscale[[1]],
    design_df = as.integer(design_df),
    ordered_psu_table = ordered_psu_table,
    design_fingerprint = NHANESJKnFingerprint(
      design_data,ordered_psu_table
    )
  )
}

AllPermutations <- function(values){
  values <- as.integer(values)
  if (length(values) == 1L){
    return(matrix(values,nrow = 1L))
  }
  do.call(rbind,lapply(seq_along(values),function(i){
    remaining <- values[-i]
    cbind(values[[i]],AllPermutations(remaining))
  }))
}

AlignNHANESJKnClasses <- function(reference_probability,
                                  replicate_probability,
                                  reference_seqn,replicate_seqn,
                                  base_weight){
  reference_probability <- as.matrix(reference_probability)
  replicate_probability <- as.matrix(replicate_probability)
  reference_seqn <- as.character(reference_seqn)
  replicate_seqn <- as.character(replicate_seqn)
  base_weight <- as.numeric(base_weight)
  if (anyDuplicated(reference_seqn) || anyDuplicated(replicate_seqn)){
    stop("SEQN values must be unique for JKn class alignment")
  }
  if (nrow(reference_probability) != length(reference_seqn) ||
      length(base_weight) != length(reference_seqn)){
    stop("Reference probabilities, SEQN, and base weights do not align")
  }
  match_index <- match(replicate_seqn,reference_seqn)
  if (anyNA(match_index)){
    stop("A replicate SEQN was not found in the Seed 159 reference cohort")
  }
  if (nrow(replicate_probability) != length(replicate_seqn)){
    stop("Replicate probabilities do not match replicate SEQN values")
  }
  mix_num <- ncol(reference_probability)
  if (ncol(replicate_probability) != mix_num){
    stop("Reference and replicate probabilities have different class counts")
  }
  reference_retained <- reference_probability[match_index,,drop = FALSE]
  retained_base_weight <- base_weight[match_index]
  overlap <- crossprod(
    reference_retained * retained_base_weight,
    replicate_probability
  )
  permutations <- AllPermutations(seq_len(mix_num))
  scores <- apply(permutations,1L,function(permutation){
    sum(overlap[cbind(seq_len(mix_num),permutation)])
  })
  score_order <- order(scores,decreasing = TRUE)
  selected <- as.integer(permutations[score_order[[1]],])
  alignment_margin <- if (length(scores) > 1L){
    scores[score_order[[1]]] - scores[score_order[[2]]]
  } else {
    Inf
  }
  list(
    permutation = selected,
    class_mapping = data.frame(
      reference_class = seq_len(mix_num),
      replicate_class = selected
    ),
    overlap = overlap,
    alignment_score = unname(scores[score_order[[1]]]),
    alignment_margin = unname(alignment_margin)
  )
}

NHANESSurvivalParameterNames <- function(surv_coef,fit_mix_num,
                                         covariate_names = NULL){
  if (is.null(covariate_names) || length(covariate_names) != length(surv_coef)){
    covariate_names <- c(
      "age",paste0("surv_covar_",seq_len(length(surv_coef) - 1L))
    )
  }
  class_names <- if (fit_mix_num > 1L){
    paste0("class_",2:fit_mix_num)
  } else {
    character(0)
  }
  other_names <- unlist(lapply(seq.int(2L,length(surv_coef)),function(i){
    if (length(surv_coef[[i]]) <= 1L){
      return(character(0))
    }
    paste0(covariate_names[[i]],"_",2:length(surv_coef[[i]]))
  }),use.names = FALSE)
  c(covariate_names[[1]],class_names,other_names)
}

MakeNHANESJKnResult <- function(jackknife_design,reference_to_save,
                                reference_file,replicate_seqn,re_prob,
                                beta_vec,surv_coef,covariate_names){
  reference_est_params <- get_saved_section(reference_to_save,"est_params")
  reference_re_prob <- get_saved_param(reference_est_params,"re_prob")
  alignment <- AlignNHANESJKnClasses(
    reference_probability = reference_re_prob,
    replicate_probability = re_prob,
    reference_seqn = jackknife_design$reference_seqn,
    replicate_seqn = replicate_seqn,
    base_weight = jackknife_design$base_weights
  )
  beta_permuted <- as.numeric(beta_vec[alignment$permutation])
  beta_aligned <- beta_permuted - beta_permuted[[1]]
  aligned_survival_vector <- IntoBetaSurvCoef(
    beta_aligned,surv_coef,length(beta_aligned)
  )
  names(aligned_survival_vector) <- NHANESSurvivalParameterNames(
    surv_coef,length(beta_aligned),covariate_names
  )
  if (length(aligned_survival_vector) != 34L ||
      any(!is.finite(aligned_survival_vector))){
    stop("Aligned NHANES survival vector must contain 34 finite coefficients")
  }

  list(
    method = "NHANES_JKn",
    replicate_id = jackknife_design$replicate_id,
    deleted_stratum = jackknife_design$deleted_stratum,
    deleted_psu = jackknife_design$deleted_psu,
    scale = jackknife_design$scale,
    rscale = jackknife_design$rscale,
    design_df = jackknife_design$design_df,
    ordered_psu_table = jackknife_design$ordered_psu_table,
    design_fingerprint = jackknife_design$design_fingerprint,
    reference_seed = NHANES_JKN_REFERENCE_SEED,
    reference_file = normalizePath(reference_file,mustWork = TRUE),
    class_mapping = alignment$class_mapping,
    overlap = alignment$overlap,
    alignment_score = alignment$alignment_score,
    alignment_margin = alignment$alignment_margin,
    aligned_beta_vec = beta_aligned,
    aligned_survival_vector = aligned_survival_vector
  )
}
