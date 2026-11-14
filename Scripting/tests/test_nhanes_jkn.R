script_arg <- grep("^--file=",commandArgs(FALSE),value = TRUE)
test_dir <- if (length(script_arg) == 1L){
  dirname(normalizePath(sub("^--file=","",script_arg)))
} else {
  file.path("Scripting","tests")
}
script_dir <- normalizePath(file.path(test_dir,".."))
project_dir <- normalizePath(file.path(script_dir,".."))

local_library <- file.path(project_dir,"Rlib")
if (dir.exists(local_library)){
  .libPaths(c(local_library,.libPaths()))
}

for (module_file in c("constants.R","saved_results.R","validation.R",
                      "settings.R","jackknife.R")){
  source(file.path(script_dir,"R",module_file))
}

stopifnot(requireNamespace("survey",quietly = TRUE))

# Argument 16 is parsed only when JMHMM explicitly enables it.
jkn_args <- c(
  "5","joint","nhanes","false","0","false","true","7","5000",
  "5","true","true","low","48","20","nhanes_jkn"
)
jkn_settings <- build_settings(
  cli_args = jkn_args,sim_num = 1,parse_jmhmm_resampling = TRUE
)
stopifnot(
  identical(jkn_settings$resampling_method,"nhanes_jkn"),
  identical(jkn_settings$replicate_id,1L),
  identical(
    build_model_name(jkn_settings),
    "JMHMMJKnJointReducedFitMix5RefSeed33Rep001len96.rda"
  )
)
se_settings <- build_settings(cli_args = jkn_args,sim_num = 1)
stopifnot(identical(se_settings$resampling_method,"none"))

# Build the design from the matched actigraphy/mortality cohort without
# materializing the time-series matrices.
data_dir <- file.path(project_dir,"Data")
wave_g_env <- new.env(parent = emptyenv())
wave_h_env <- new.env(parent = emptyenv())
mortality_env <- new.env(parent = emptyenv())
load(file.path(data_dir,"Wavedata_G.rda"),envir = wave_g_env)
load(file.path(data_dir,"Wavedata_H.rda"),envir = wave_h_env)
load(file.path(data_dir,"NHANES_2011_2012_2013_2014.rda"),
     envir = mortality_env)
id <- rbind(wave_g_env$wave_data_G[[3]],wave_h_env$wave_data_H[[3]])
nhanes_mort_list <- mortality_env$NHANES_mort_list
lmf_data <- rbind(
  nhanes_mort_list[[1]][nhanes_mort_list[[1]]$eligstat == 1,,drop = FALSE],
  nhanes_mort_list[[2]][nhanes_mort_list[[2]]$eligstat == 1,,drop = FALSE]
)
id <- id[id$SEQN %in% lmf_data$seqn,,drop = FALSE]
base_weights <- id$sweights / NHANES_NUM_WAVES
base_weights <- base_weights / mean(base_weights)
jkn <- BuildNHANESJKnDesign(id,base_weights,replicate_id = 1L)
jkn_last <- BuildNHANESJKnDesign(id,base_weights,replicate_id = 61L)
stopifnot(
  length(jkn$keep) == 6483L,
  sum(!jkn$keep) > 0L,
  nrow(jkn$ordered_psu_table) == 61L,
  identical(jkn$design_df,32L),
  identical(jkn$replicate_id,1L),
  jkn$scale == 1,
  jkn$rscale %in% c(0.5,2 / 3),
  length(unique(jkn$ordered_psu_table$replicate_id)) == 61L,
  identical(jkn$design_fingerprint,jkn_last$design_fingerprint),
  identical(jkn$ordered_psu_table,jkn_last$ordered_psu_table),
  !identical(
    c(jkn$deleted_stratum,jkn$deleted_psu),
    c(jkn_last$deleted_stratum,jkn_last$deleted_psu)
  )
)

# Recover a known five-class column permutation from soft probabilities.
set.seed(20260915)
reference_probability <- matrix(
  rexp(500),nrow = 100,ncol = 5
)
reference_probability <- reference_probability /
  rowSums(reference_probability)
replicate_column_order <- c(3L,1L,5L,2L,4L)
replicate_probability <-
  reference_probability[,replicate_column_order,drop = FALSE]
expected_mapping <- match(1:5,replicate_column_order)
alignment <- AlignNHANESJKnClasses(
  reference_probability,
  replicate_probability,
  reference_seqn = seq_len(100),
  replicate_seqn = seq_len(100),
  base_weight = seq(0.5,1.5,length.out = 100)
)
stopifnot(
  identical(alignment$permutation,as.integer(expected_mapping)),
  identical(sort(alignment$class_mapping$replicate_class),1:5),
  alignment$alignment_margin > 0
)
raw_beta <- c(0,0.2,-0.1,0.4,0.7)
aligned_beta <- raw_beta[alignment$permutation]
aligned_beta <- aligned_beta - aligned_beta[[1]]
stopifnot(aligned_beta[[1]] == 0)

# The current NHANES covariate layout yields 34 estimable coefficients.
surv_coef_lengths <- c(1,2,6,5,5,4,2,2,2,2,2,5,3,2)
fake_surv_coef <- lapply(surv_coef_lengths,numeric)
parameter_names <- NHANESSurvivalParameterNames(
  fake_surv_coef,5L,
  c("age","gender","race","overall_health","education","bmi",
    "diabetes","chd","chf","heart_attack","stroke","alcohol",
    "smoking","physical_function")
)
stopifnot(
  length(parameter_names) == 34L,
  identical(parameter_names[1:5],
            c("age","class_2","class_3","class_4","class_5"))
)

# Verify survey::svrVar against its direct MSE-centered quadratic form.
synthetic_thetas <- matrix(
  c(1.1,2.2,0.9,1.8,1.2,2.1),nrow = 3,byrow = TRUE
)
synthetic_coef <- c(1,2)
synthetic_scale <- 0.75
synthetic_rscales <- c(0.5,2 / 3,0.5)
survey_vcov <- survey::svrVar(
  thetas = synthetic_thetas,
  scale = synthetic_scale,
  rscales = synthetic_rscales,
  mse = TRUE,
  coef = synthetic_coef
)
centered <- sweep(synthetic_thetas,2,synthetic_coef,"-")
manual_vcov <- matrix(0,2,2)
for (i in seq_len(nrow(centered))){
  manual_vcov <- manual_vcov + synthetic_scale * synthetic_rscales[[i]] *
    tcrossprod(centered[i,])
}
stopifnot(max(abs(as.matrix(survey_vcov) - manual_vcov)) < 1e-12)

# Strict reference validation is exercised when the fixed files are present.
reference_dir <- Sys.getenv("NHANES_JKN_REFERENCE_DIR",data_dir)
reference_paths <- vapply(
  c("joint","two_stage"),NHANESJKnReferencePath,
  FUN.VALUE = character(1),data_dir = reference_dir
)
if (all(file.exists(reference_paths))){
  for (model_type in names(reference_paths)){
    reference <- load_to_save_from_path(reference_paths[[model_type]])
    ValidateNHANESJKnReference(
      reference,reference_paths[[model_type]],model_type
    )
  }
} else {
  message(
    "Reference-file checks skipped in unit test; run ",
    "validate_nhanes_jkn_preflight.R after placing both Seed 33 files in Data."
  )
}

message("All NHANES JKn unit tests passed.")
