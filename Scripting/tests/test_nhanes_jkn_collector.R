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
                      "jackknife.R","survival.R")){
  source(file.path(script_dir,"R",module_file))
}

reference_dir <- Sys.getenv(
  "NHANES_JKN_REFERENCE_DIR",file.path(project_dir,"Data")
)
reference_paths <- vapply(
  c("joint","two_stage"),NHANESJKnReferencePath,
  FUN.VALUE = character(1),data_dir = reference_dir
)
if (!all(file.exists(reference_paths))){
  message(
    "Collector integration test skipped because both fixed Seed 33 ",
    "reference files are not present."
  )
  quit(save = "no",status = 0)
}

stopifnot(requireNamespace("survey",quietly = TRUE))
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
jkn_design <- BuildNHANESJKnDesign(id,base_weights,1L)

test_root <- tempfile("nhanes-jkn-collector-")
replicate_dir <- file.path(test_root,"replicates")
output_dir <- file.path(test_root,"results")
dir.create(replicate_dir,recursive = TRUE)
on.exit(unlink(test_root,recursive = TRUE,force = TRUE),add = TRUE)

model_labels <- c(joint = "Joint",two_stage = "TwoStage")
covariate_names <- c(
  "age","gender","race","overall_health","education","bmi",
  "diabetes","chd","chf","heart_attack","stroke","alcohol",
  "smoking","physical_function"
)
for (model_type in names(model_labels)){
  reference <- load_to_save_from_path(reference_paths[[model_type]])
  ValidateNHANESJKnReference(
    reference,reference_paths[[model_type]],model_type
  )
  reference_est <- get_saved_section(reference,"est_params")
  beta_vec <- get_saved_param(reference_est,"beta_vec")
  surv_coef <- get_saved_param(reference_est,"surv_coef")
  survival_vector <- IntoBetaSurvCoef(beta_vec,surv_coef,5L)
  names(survival_vector) <- NHANESSurvivalParameterNames(
    surv_coef,5L,covariate_names
  )

  for (replicate_id in 1:61){
    psu_row <- jkn_design$ordered_psu_table[replicate_id,,drop = FALSE]
    perturbation <- sin(replicate_id + seq_along(survival_vector)) / 100
    aligned_survival_vector <- survival_vector + perturbation
    aligned_beta_vec <- c(0,as.numeric(aligned_survival_vector[2:5]))
    to_save <- list(
      est_params = list(),
      diagnostics = list(
        convergence = list(converged = TRUE,completed_iterations = 7L),
        runtime_seconds = 1
      ),
      settings = list(
        fit_mix_num = 5,
        model_type = model_type,
        data_source = "nhanes",
        sim_num = replicate_id,
        period_len = 96,
        resampling_method = "nhanes_jkn"
      ),
      jackknife = list(
        method = "NHANES_JKn",
        replicate_id = replicate_id,
        deleted_stratum = psu_row$deleted_stratum[[1]],
        deleted_psu = psu_row$deleted_psu[[1]],
        scale = jkn_design$scale,
        rscale = psu_row$rscale[[1]],
        design_df = 32L,
        ordered_psu_table = jkn_design$ordered_psu_table,
        design_fingerprint = jkn_design$design_fingerprint,
        reference_seed = 33L,
        reference_file = reference_paths[[model_type]],
        class_mapping = data.frame(
          reference_class = 1:5,replicate_class = 1:5
        ),
        alignment_score = 100,
        alignment_margin = 10,
        aligned_beta_vec = aligned_beta_vec,
        aligned_survival_vector = aligned_survival_vector
      )
    )
    output_name <- sprintf(
      "JMHMMJKn%sReducedFitMix5RefSeed33Rep%03dlen96.rda",
      model_labels[[model_type]],replicate_id
    )
    save(to_save,file = file.path(replicate_dir,output_name))
  }
}

rscript <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)
collector_status <- system2(
  rscript,
  args = shQuote(c(
    file.path(script_dir,"collect_nhanes_jkn.R"),
    replicate_dir,reference_dir,output_dir
  ))
)
stopifnot(collector_status == 0L)

joint_result <- readRDS(file.path(output_dir,"nhanes_jkn_joint.rds"))
two_stage_result <- readRDS(
  file.path(output_dir,"nhanes_jkn_two_stage.rds")
)
stopifnot(
  all(dim(joint_result$replicate_estimates) == c(61L,34L)),
  all(dim(two_stage_result$replicate_estimates) == c(61L,34L)),
  all(is.finite(joint_result$vcov)),
  all(is.finite(two_stage_result$vcov)),
  joint_result$design_df == 32L,
  two_stage_result$design_df == 32L,
  file.exists(file.path(
    output_dir,"nhanes_jkn_joint_vs_two_stage_latent_class_hr.png"
  ))
)

diagnostics <- read.csv(file.path(
  output_dir,"nhanes_jkn_replicate_diagnostics.csv"
))
stopifnot(nrow(diagnostics) == 122L,all(diagnostics$valid))

message("NHANES JKn collector integration test passed.")
