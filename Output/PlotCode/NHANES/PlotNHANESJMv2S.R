#!/usr/bin/env Rscript

script_arg <- grep("^--file=",commandArgs(FALSE),value = TRUE)
if (length(script_arg) == 1L){
  script_file <- normalizePath(sub("^--file=","",script_arg),mustWork = TRUE)
  project_dir <- normalizePath(file.path(dirname(script_file),"..","..",".."),
                               mustWork = TRUE)
  local_library <- file.path(project_dir,"Rlib")
  if (dir.exists(local_library)){
    .libPaths(c(local_library,.libPaths()))
  }
}

library(dplyr)
library(forestploter)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
joint_reference_file <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path(
    "Data",
    "JMHMMLoadInFitMix5Seed159len96_oakes_profiled_se.rda"
  )
}
two_stage_reference_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path(
    "Data",
    "JMHMMNoSurvRandInitReducedFitMix5Seed159len96.rda"
  )
}
output_prefix <- if (length(args) >= 3){
  args[[3]]
} else {
  file.path("Output","Figures","nhanes_jm_vs_two_stage")
}

assert_true <- function(condition,message){
  if (!isTRUE(condition)){
    stop(message,call. = FALSE)
  }
}

read_joint_model_based <- function(reference_file){
  assert_true(file.exists(reference_file),
              paste0("Joint Seed 159 reference does not exist: ",reference_file))
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(reference_file,envir = load_env)
  assert_true(identical(loaded_names,"to_save"),
              "Joint Seed 159 reference must contain only to_save")
  to_save <- load_env$to_save
  assert_true(identical(to_save$settings$model_type,"joint") &&
                identical(to_save$settings$data_source,"nhanes") &&
                as.integer(to_save$settings$sim_num) == 159L,
              "Class-size reference is not the joint NHANES Seed 159 fit")
  re_prob <- to_save$est_params$re_prob
  beta_vec <- to_save$est_params$beta_vec
  assert_true(is.matrix(re_prob) && ncol(re_prob) == 5L,
              "Joint Seed 159 re_prob must be a five-column matrix")
  assert_true(length(beta_vec) == 5L,
              "Joint Seed 159 beta_vec must have length five")

  model_vcov <- to_save$oakes$model_based_survival_vcov
  parameter_map <- to_save$oakes$parameter_map
  survival_map <- parameter_map[parameter_map$block == "survival",,drop = FALSE]
  assert_true(is.matrix(model_vcov) &&
                all(dim(model_vcov) == c(34L,34L)),
              "Oakes model-based survival covariance must be 34 by 34")
  assert_true(nrow(survival_map) == 34L &&
                all(c("param_name","block") %in% names(survival_map)),
              "Oakes survival parameter map is invalid")
  assert_true(all(is.finite(model_vcov)) && all(diag(model_vcov) > 0),
              "Oakes model-based survival covariance is invalid")

  class_names <- paste0("class_",2:5)
  class_indices <- match(class_names,survival_map$param_name)
  assert_true(!anyNA(class_indices),
              "Oakes survival parameter map is missing latent-class effects")
  class_coef <- as.numeric(beta_vec[2:5])
  class_se <- sqrt(diag(model_vcov))[class_indices]
  critical_value <- qnorm(0.975)

  class_hr <- data.frame(
    latent_class = 2:5,
    coef = class_coef,
    se = class_se
  ) |>
    mutate(
      hr = exp(coef),
      ci_low = exp(coef - critical_value * se),
      ci_high = exp(coef + critical_value * se),
      hr_text = sprintf("%.2f (%.2f-%.2f)",hr,ci_low,ci_high)
    )

  list(
    class_hr = class_hr,
    class_size = tabulate(apply(re_prob,1,which.max),nbins = 5L)
  )
}

read_two_stage_naive_cox <- function(reference_file){
  assert_true(file.exists(reference_file),
              paste0("Two-stage Seed 159 reference does not exist: ",
                     reference_file))
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(reference_file,envir = load_env)
  assert_true(identical(loaded_names,"to_save"),
              "Two-stage Seed 159 reference must contain only to_save")
  to_save <- load_env$to_save
  assert_true(identical(to_save$settings$model_type,"two_stage") &&
                identical(to_save$settings$data_source,"nhanes") &&
                as.integer(to_save$settings$sim_num) == 159L,
              "Reference is not the two-stage NHANES Seed 159 fit")

  beta_vec <- to_save$est_params$beta_vec
  beta_se <- to_save$est_params$beta_se
  assert_true(length(beta_vec) == 5L,
              "Two-stage Seed 159 beta_vec must have length five")
  assert_true(length(beta_se) == 34L &&
                all(is.finite(beta_se)) && all(beta_se > 0),
              "Two-stage standard Cox SE vector is invalid")

  class_coef <- as.numeric(beta_vec[2:5])
  class_se <- as.numeric(beta_se[2:5])
  critical_value <- qnorm(0.975)

  data.frame(
    latent_class = 2:5,
    coef = class_coef,
    se = class_se
  ) |>
    mutate(
      hr = exp(coef),
      ci_low = exp(coef - critical_value * se),
      ci_high = exp(coef + critical_value * se),
      hr_text = sprintf("%.2f (%.2f-%.2f)",hr,ci_low,ci_high)
    )
}

joint_model_based <- read_joint_model_based(joint_reference_file)
two_stage_naive_cox <- read_two_stage_naive_cox(two_stage_reference_file)
class_size <- joint_model_based$class_size

jlcm_df <- joint_model_based$class_hr |>
  select(
    latent_class,
    jlcm_hr = hr,
    jlcm_low = ci_low,
    jlcm_high = ci_high,
    jlcm_text = hr_text
  )

two_stage_df <- two_stage_naive_cox |>
  select(
    latent_class,
    two_stage_hr = hr,
    two_stage_low = ci_low,
    two_stage_high = ci_high,
    two_stage_text = hr_text
  )

forest_df <- merge(jlcm_df,two_stage_df,by = "latent_class") |>
  arrange(latent_class) |>
  mutate(
    `Latent Class` = paste("Class",latent_class),
    Size = as.character(class_size[latent_class]),
    ` ` = paste(rep(" ",24),collapse = " "),
    `JLMHMM HR (95% CI)` = jlcm_text,
    `Two-Stage HR (95% CI)` = two_stage_text
  )

plot_table <- forest_df |>
  select(
    `Latent Class`,
    Size,
    ` `,
    `JLMHMM HR (95% CI)`,
    `Two-Stage HR (95% CI)`
  )

# Wrap the two long result-column headers so they do not collide.
names(plot_table)[4:5] <- c(
  "JLMHMM HR\n(95% CI)",
  "Two-Stage HR\n(95% CI)"
)

plot_theme <- forest_theme(
  base_size = 11,
  ci_pch = 15,
  ci_col = c("#440154","#7AD151"),
  ci_fill = c("#440154","#7AD151"),
  ci_lwd = 1.2,
  ci_Theight = 0.15,
  refline_gp = gpar(lty = "dashed",col = "grey35",lwd = 1),
  xaxis_gp = gpar(fontsize = 10),
  xlab_gp = gpar(fontsize = 11),
  legend_name = "Model",
  legend_value = c("JLMHMM","Two-Stage"),
  legend_position = "right"
)

jm_vs_two_stage_plot <- forest(
  plot_table,
  est = list(forest_df$jlcm_hr,forest_df$two_stage_hr),
  lower = list(forest_df$jlcm_low,forest_df$two_stage_low),
  upper = list(forest_df$jlcm_high,forest_df$two_stage_high),
  ci_column = 3,
  ref_line = 1,
  xlim = c(0.5,4.5),
  ticks_at = c(1,2,3,4),
  ticks_digits = 1L,
  xlab = "Hazard Ratio",
  nudge_y = 0.18,
  theme = plot_theme
)

jm_vs_two_stage_plot <- edit_plot(
  jm_vs_two_stage_plot,
  row = c(1,3),
  col = seq_along(plot_table),
  which = "background",
  gp = gpar(fill = "#edf2f0",col = NA)
)

jm_vs_two_stage_plot <- add_border(
  jm_vs_two_stage_plot,
  part = "header",
  where = "bottom",
  gp = gpar(lwd = 1)
)

if (interactive()){
  plot(jm_vs_two_stage_plot,autofit = TRUE)
}

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
output_file <- paste0(output_prefix,".png")
png(
  output_file,
  width = 11,
  height = 5,
  units = "in",
  res = 300
)
plot(jm_vs_two_stage_plot,autofit = TRUE)
invisible(dev.off())

message("Read Seed 159 model-based Oakes results from: ",joint_reference_file)
message("Read Seed 159 naive Cox results from: ",two_stage_reference_file)
message("Saved: ",output_file)
