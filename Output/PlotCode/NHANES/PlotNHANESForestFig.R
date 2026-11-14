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
jkn_file <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Output","parse_nhanes_jkn.rds")
}
model_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Data","JMHMMLoadInFitMix5Seed159len96_oakes_profiled_se.rda")
}
nhanes_file <- file.path("Data","NHANES_2011_2012_2013_2014.rda")
wave_g_file <- file.path("Data","Wavedata_G.rda")
wave_h_file <- file.path("Data","Wavedata_H.rda")
output_prefix <- if (length(args) >= 3){
  args[[3]]
} else {
  file.path("Output","Figures","nhanes_forest_fig")
}

assert_true <- function(condition,message){
  if (!isTRUE(condition)){
    stop(message,call. = FALSE)
  }
}

read_joint_jkn <- function(file){
  assert_true(file.exists(file),
              paste0("NHANES JKn parser output does not exist: ",file))
  parsed <- readRDS(file)
  assert_true(is.list(parsed) && isTRUE(parsed$inference_available),
              "NHANES JKn parser output does not contain valid inference")
  joint <- parsed$joint_survival
  assert_true(is.list(joint),
              "Parser output does not contain the joint_survival block")
  summary <- joint$coefficient_summary
  vcov <- joint$vcov
  required <- c(
    "param_name","estimate","se","design_df","critical_value",
    "ci_lower","ci_upper","hr","hr_ci_lower","hr_ci_upper"
  )
  assert_true(is.data.frame(summary) && nrow(summary) == 34L,
              "Joint survival summary must contain 34 coefficients")
  assert_true(all(required %in% names(summary)),
              "Joint survival summary is missing required columns")
  assert_true(!anyDuplicated(summary$param_name),
              "Joint survival parameter names are duplicated")
  assert_true(is.matrix(vcov) && all(dim(vcov) == c(34L,34L)),
              "Joint survival covariance must be a 34 by 34 matrix")
  assert_true(identical(rownames(vcov),summary$param_name) &&
                identical(colnames(vcov),summary$param_name),
              "Joint survival covariance names do not match the summary")
  assert_true(all(is.finite(as.matrix(vcov))) &&
                all(is.finite(summary$estimate)) &&
                all(is.finite(summary$se)),
              "Joint survival results contain nonfinite values")
  assert_true(all(summary$design_df == 32L) &&
                length(unique(summary$critical_value)) == 1L,
              "Joint survival intervals do not use 32 design degrees of freedom")
  assert_true(max(abs(summary$se - sqrt(diag(vcov)))) <= 1e-10,
              "Joint survival SEs do not match the saved covariance")
  list(
    summary = summary,
    vcov = vcov,
    critical_value = unique(summary$critical_value)[[1]]
  )
}

joint_jkn <- read_joint_jkn(jkn_file)

model_env <- new.env(parent = emptyenv())
loaded_names <- load(model_file,envir = model_env)
assert_true(identical(loaded_names,"to_save"),
            "Joint Seed 159 reference must contain only to_save")
to_save <- model_env$to_save
assert_true(identical(to_save$settings$model_type,"joint") &&
              identical(to_save$settings$data_source,"nhanes") &&
              as.integer(to_save$settings$sim_num) == 159L,
            "Class-size model is not the joint NHANES Seed 159 fit")

est_params <- to_save$est_params
re_prob <- est_params$re_prob
mix_num <- ncol(re_prob)
assert_true(is.matrix(re_prob) && mix_num == 5L,
            "Joint Seed 159 re_prob must be a five-column matrix")

nhanes_env <- new.env()
load(nhanes_file,envir = nhanes_env)
nhanes_mort_list <- nhanes_env$NHANES_mort_list
lmf_data <- bind_rows(
  filter(nhanes_mort_list[[1]],eligstat == 1),
  filter(nhanes_mort_list[[2]],eligstat == 1)
)

wave_g_env <- new.env()
wave_h_env <- new.env()
load(wave_g_file,envir = wave_g_env)
load(wave_h_file,envir = wave_h_env)

id <- bind_rows(
  wave_g_env$wave_data_G[[3]],
  wave_h_env$wave_data_H[[3]]
)

id <- id[id$SEQN %in% lmf_data$seqn,]
lmf_data <- lmf_data[lmf_data$seqn %in% id$SEQN,]
stopifnot(sum(id$SEQN - lmf_data$seqn) == 0)

class_size <- tabulate(
  apply(re_prob,1,which.max),
  nbins = mix_num
)

count_level <- function(variable,level){
  sum(id[[variable]] == level,na.rm = TRUE)
}

joint_summary <- joint_jkn$summary
seed159_coef_values <- c(
  as.numeric(est_params$surv_coef[[1]][[1]]),
  as.numeric(est_params$beta_vec[-1]),
  unlist(lapply(est_params$surv_coef[-1],function(x) as.numeric(x[-1])))
)
assert_true(length(seed159_coef_values) == nrow(joint_summary) &&
              all(is.finite(seed159_coef_values)),
            "Could not reconstruct the 34 Seed 159 survival estimates")
coef_values <- setNames(seed159_coef_values,joint_summary$param_name)

# ## Uncomment for jkn se
# se_values <- setNames(joint_summary$se,joint_summary$param_name)
# critical_value <- joint_jkn$critical_value

## Uncomment for Oakes model-based covariance
model_based_survival_vcov <- to_save$oakes$model_based_survival_vcov
assert_true(is.matrix(model_based_survival_vcov) &&
              all(dim(model_based_survival_vcov) == c(34L,34L)),
            "Oakes model-based survival covariance must be 34 by 34")
assert_true(all(is.finite(model_based_survival_vcov)) &&
              all(diag(model_based_survival_vcov) > 0),
            "Oakes model-based survival covariance is invalid")
dimnames(model_based_survival_vcov) <- list(
  joint_summary$param_name,
  joint_summary$param_name
)
se_values <- sqrt(diag(model_based_survival_vcov))
critical_value <- qnorm(0.975)


linear_contrast <- function(weights){
  assert_true(all(names(weights) %in% names(coef_values)),
              "A requested contrast contains an unknown coefficient")
  contrast <- setNames(numeric(length(coef_values)),names(coef_values))
  contrast[names(weights)] <- weights
  estimate <- sum(contrast * coef_values)
  variance <- drop(
    t(contrast) %*%
      model_based_survival_vcov[names(contrast),names(contrast)] %*%
      contrast
  )
  assert_true(is.finite(variance) && variance >= -1e-12,
              "A calculated contrast variance is invalid")
  c(coef = estimate,se = sqrt(max(variance,0)))
}

bmi_under <- linear_contrast(c(bmi_2 = -1))
bmi_25_30 <- linear_contrast(c(bmi_2 = -1,bmi_3 = 1))
bmi_gt_30 <- linear_contrast(c(bmi_2 = -1,bmi_4 = 1))

coef_values <- c(
  coef_values,
  bmi_under_18_5_ref_18_5_25 = bmi_under[["coef"]],
  bmi_25_30_ref_18_5_25 = bmi_25_30[["coef"]],
  bmi_gt_30_ref_18_5_25 = bmi_gt_30[["coef"]]
)
se_values <- c(
  se_values,
  bmi_under_18_5_ref_18_5_25 = bmi_under[["se"]],
  bmi_25_30_ref_18_5_25 = bmi_25_30[["se"]],
  bmi_gt_30_ref_18_5_25 = bmi_gt_30[["se"]]
)

add_header <- function(label,size = ""){
  data.frame(
    subgroup = label,
    size = as.character(size),
    coef_name = NA_character_,
    is_header = TRUE
  )
}

add_coef <- function(label,size,coef_name,indent = TRUE){
  data.frame(
    subgroup = if (indent) paste0("  ",label) else label,
    size = as.character(size),
    coef_name = coef_name,
    is_header = FALSE
  )
}

forest_df <- bind_rows(
  add_header("Latent Class (Ref: Class 1)",class_size[1]),
  add_coef("Class 2",class_size[2],"class_2"),
  add_coef("Class 3",class_size[3],"class_3"),
  add_coef("Class 4",class_size[4],"class_4"),
  add_coef("Class 5",class_size[5],"class_5"),
  add_coef("Age (1-year Increase)","", "age",indent = FALSE),
  add_coef("Female",count_level("gender",1),"gender_2",indent = FALSE),
  add_header("Race (Ref: Mexican American)",count_level("race",0)),
  add_coef("Other Hispanic",count_level("race",1),"race_2"),
  add_coef("Non-Hispanic White",count_level("race",2),"race_3"),
  add_coef("Non-Hispanic Black",count_level("race",3),"race_4"),
  add_coef("Non-Hispanic Asian",count_level("race",4),"race_5"),
  add_coef("Other Race - Incl Multi-Racial",count_level("race",5),"race_6"),
  add_header("Overall Health (Ref: Excellent)",count_level("overall_health",0)),
  add_coef("Very Good",count_level("overall_health",1),"overall_health_2"),
  add_coef("Good",count_level("overall_health",2),"overall_health_3"),
  add_coef("Fair",count_level("overall_health",3),"overall_health_4"),
  add_coef("Poor",count_level("overall_health",4),"overall_health_5"),
  add_header("Education (Ref: <9th Grade)",count_level("education",0)),
  add_coef("9-11th Grade",count_level("education",1),"education_2"),
  add_coef("High School/GED",count_level("education",2),"education_3"),
  add_coef("Some College/AA Degree",count_level("education",3),"education_4"),
  add_coef("College Graduate or Above",count_level("education",4),"education_5"),
  add_header("BMI (Ref: 18.5-25)",count_level("bmi_disc",1)),
  add_coef("<18.5",count_level("bmi_disc",0),
           "bmi_under_18_5_ref_18_5_25"),
  add_coef("25-30",count_level("bmi_disc",2),
           "bmi_25_30_ref_18_5_25"),
  add_coef(">30",count_level("bmi_disc",3),
           "bmi_gt_30_ref_18_5_25"),
  add_coef("Diabetes",count_level("diabetes",1),"diabetes_2",indent = FALSE),
  add_coef("CHD",count_level("CHD",1),"chd_2",indent = FALSE),
  add_coef("CHF",count_level("CHF",1),"chf_2",indent = FALSE),
  add_coef("Heart Attack",count_level("heart_attack",1),"heart_attack_2",
           indent = FALSE),
  add_coef("Stroke",count_level("stroke",1),"stroke_2",indent = FALSE),
  add_header("Alcohol Use (Ref: Never)",count_level("alcohol",0)),
  add_coef("Former",count_level("alcohol",1),"alcohol_2"),
  add_coef("Moderate",count_level("alcohol",2),"alcohol_3"),
  add_coef("Heavy",count_level("alcohol",3),"alcohol_4"),
  add_coef("Other/Unknown",count_level("alcohol",4),"alcohol_5"),
  add_header("Smoking Status (Ref: Never)",count_level("smoking",0)),
  add_coef("Former",count_level("smoking",1),"smoking_2"),
  add_coef("Current",count_level("smoking",2),"smoking_3"),
  add_coef("Reduced Mobility",count_level("phyfunc",1),"physical_function_2",
           indent = FALSE)
)

required_coef_names <- forest_df$coef_name[!is.na(forest_df$coef_name)]
assert_true(all(required_coef_names %in% names(coef_values)) &&
              all(required_coef_names %in% names(se_values)),
            "The forest plot requests coefficients missing from JKn output")

forest_df <- forest_df |>
  mutate(
    coef = coef_values[coef_name],
    se = se_values[coef_name],
    hr = exp(coef),
    ci_low = exp(coef - critical_value * se),
    ci_high = exp(coef + critical_value * se),
    hr_text = ifelse(
      is.na(hr),
      "",
      sprintf("%.2f (%.2f-%.2f)",hr,ci_low,ci_high)
    ),
    row_number = seq_len(n()),
    shade = row_number %% 2 == 1
  )

plot_max_hr <- 3
plot_ticks <- 0:3

plot_table <- data.frame(
  Subgroup = forest_df$subgroup,
  ` ` = rep(paste(rep(" ",12),collapse = ""),nrow(forest_df)),
  Size = forest_df$size,
  `  ` = rep(paste(rep(" ",24),collapse = " "),nrow(forest_df)),
  `HR (95% CI)` = forest_df$hr_text,
  check.names = FALSE
)

plot_theme <- forest_theme(
  base_size = 10,
  ci_pch = 15,
  ci_col = "black",
  ci_fill = "black",
  ci_lwd = 1,
  ci_Theight = 0.15,
  refline_gp = gpar(lty = "dashed",col = "grey35",lwd = 1),
  xaxis_gp = gpar(fontsize = 9),
  xlab_gp = gpar(fontsize = 10),
  legend_position = "none"
)

forest_plot <- forest(
  plot_table,
  est = forest_df$hr,
  lower = forest_df$ci_low,
  upper = forest_df$ci_high,
  ci_column = 4,
  ref_line = 1,
  xlim = c(0,plot_max_hr),
  ticks_at = plot_ticks,
  ticks_digits = 1L,
  xlab = "Hazard Ratio",
  theme = plot_theme
)

forest_plot <- edit_plot(
  forest_plot,
  row = which(forest_df$shade),
  col = seq_along(plot_table),
  which = "background",
  gp = gpar(fill = "#edf2f0",col = NA)
)

forest_plot <- add_border(
  forest_plot,
  part = "header",
  where = "bottom",
  gp = gpar(lwd = 1)
)

if (interactive()){
  plot(forest_plot,autofit = TRUE)
}
dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
output_file <- paste0(output_prefix,".png")
png(
  output_file,
  width = 8.5,
  height = 12,
  units = "in",
  res = 300
)
plot(forest_plot,autofit = TRUE)
invisible(dev.off())

message("Read joint survival JKn results from: ",jkn_file)
message("Saved: ",output_file)
