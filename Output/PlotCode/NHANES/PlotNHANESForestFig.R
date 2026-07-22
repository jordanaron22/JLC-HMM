#!/usr/bin/env Rscript

library(dplyr)
library(forestploter)
library(grid)

model_file <- file.path("Data","JMHMMFitMix5Seed.rda")
nhanes_file <- file.path("Data","NHANES_2011_2012_2013_2014.rda")
wave_g_file <- file.path("Data","Wavedata_G.rda")
wave_h_file <- file.path("Data","Wavedata_H.rda")
output_prefix <- file.path("Output","Figures","nhanes_forest_fig")

plot_max_hr <- 3

load(model_file)

est_params <- to_save$est_params
beta_vec <- est_params$beta_vec
surv_coef <- est_params$surv_coef
beta_se <- est_params$beta_se
re_prob <- est_params$re_prob
mix_num <- length(beta_vec)

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

coef_values <- c(
  age = surv_coef[[1]],
  setNames(beta_vec[-1],paste0("latent_class",2:mix_num)),
  gender_1 = surv_coef[[2]][2],
  race_1 = surv_coef[[3]][2],
  race_2 = surv_coef[[3]][3],
  race_3 = surv_coef[[3]][4],
  race_4 = surv_coef[[3]][5],
  race_5 = surv_coef[[3]][6],
  overall_health_1 = surv_coef[[4]][2],
  overall_health_2 = surv_coef[[4]][3],
  overall_health_3 = surv_coef[[4]][4],
  overall_health_4 = surv_coef[[4]][5],
  education_1 = surv_coef[[5]][2],
  education_2 = surv_coef[[5]][3],
  education_3 = surv_coef[[5]][4],
  education_4 = surv_coef[[5]][5],
  bmi_1 = surv_coef[[6]][2],
  bmi_2 = surv_coef[[6]][3],
  bmi_3 = surv_coef[[6]][4],
  diabetes_1 = surv_coef[[7]][2],
  chd_1 = surv_coef[[8]][2],
  chf_1 = surv_coef[[9]][2],
  heart_attack_1 = surv_coef[[10]][2],
  stroke_1 = surv_coef[[11]][2],
  alcohol_1 = surv_coef[[12]][2],
  alcohol_2 = surv_coef[[12]][3],
  alcohol_3 = surv_coef[[12]][4],
  alcohol_4 = surv_coef[[12]][5],
  smoking_1 = surv_coef[[13]][2],
  smoking_2 = surv_coef[[13]][3],
  physical_function_1 = surv_coef[[14]][2]
)

coef_values <- c(
  coef_values,
  bmi_under_18_5_ref_18_5_25 = -coef_values[["bmi_1"]],
  bmi_25_30_ref_18_5_25 =
    coef_values[["bmi_2"]] - coef_values[["bmi_1"]],
  bmi_gt_30_ref_18_5_25 =
    coef_values[["bmi_3"]] - coef_values[["bmi_1"]]
)

se_values <- beta_se
se_values <- c(
  se_values,
  bmi_under_18_5_ref_18_5_25 = beta_se[["bmi_1"]],
  # The saved object has marginal SEs, not the full covariance matrix.
  # Use the conservative independent-difference approximation here.
  bmi_25_30_ref_18_5_25 =
    sqrt(beta_se[["bmi_2"]]^2 + beta_se[["bmi_1"]]^2),
  bmi_gt_30_ref_18_5_25 =
    sqrt(beta_se[["bmi_3"]]^2 + beta_se[["bmi_1"]]^2)
)

add_header <- function(label,size = ""){
  data.frame(
    subgroup = label,
    size = as.character(size),
    coef_name = NA_character_,
    is_header = TRUE
  )
}

add_coef <- function(label,size,coef_name){
  data.frame(
    subgroup = paste0("  ",label),
    size = as.character(size),
    coef_name = coef_name,
    is_header = FALSE
  )
}

forest_df <- bind_rows(
  add_header("Latent Class (Ref: Class 1)",class_size[1]),
  add_coef("Class 2",class_size[2],"latent_class2"),
  add_coef("Class 3",class_size[3],"latent_class3"),
  add_coef("Class 4",class_size[4],"latent_class4"),
  add_coef("Class 5",class_size[5],"latent_class5"),
  add_coef("Age (1-year Increase)","", "age"),
  add_coef("Female",count_level("gender",1),"gender_1"),
  add_header("Race (Ref: Mexican American)",count_level("race",0)),
  add_coef("Other Hispanic",count_level("race",1),"race_1"),
  add_coef("Non-Hispanic White",count_level("race",2),"race_2"),
  add_coef("Non-Hispanic Black",count_level("race",3),"race_3"),
  add_coef("Non-Hispanic Asian",count_level("race",4),"race_4"),
  add_coef("Other Race - Incl Multi-Racial",count_level("race",5),"race_5"),
  add_header("Overall Health (Ref: Excellent)",count_level("overall_health",0)),
  add_coef("Very Good",count_level("overall_health",1),"overall_health_1"),
  add_coef("Good",count_level("overall_health",2),"overall_health_2"),
  add_coef("Fair",count_level("overall_health",3),"overall_health_3"),
  add_coef("Poor",count_level("overall_health",4),"overall_health_4"),
  add_header("Education (Ref: <9th Grade)",count_level("education",0)),
  add_coef("9-11th Grade",count_level("education",1),"education_1"),
  add_coef("High School/GED",count_level("education",2),"education_2"),
  add_coef("Some College/AA Degree",count_level("education",3),"education_3"),
  add_coef("College Graduate or Above",count_level("education",4),"education_4"),
  add_header("BMI (Ref: 18.5-25)",count_level("bmi_disc",1)),
  add_coef("<18.5",count_level("bmi_disc",0),
           "bmi_under_18_5_ref_18_5_25"),
  add_coef("25-30",count_level("bmi_disc",2),
           "bmi_25_30_ref_18_5_25"),
  add_coef(">30",count_level("bmi_disc",3),
           "bmi_gt_30_ref_18_5_25"),
  add_coef("Diabetes",count_level("diabetes",1),"diabetes_1"),
  add_coef("CHD",count_level("CHD",1),"chd_1"),
  add_coef("CHF",count_level("CHF",1),"chf_1"),
  add_coef("Heart Attack",count_level("heart_attack",1),"heart_attack_1"),
  add_coef("Stroke",count_level("stroke",1),"stroke_1"),
  add_header("Alcohol Use (Ref: Never)",count_level("alcohol",0)),
  add_coef("Former",count_level("alcohol",1),"alcohol_1"),
  add_coef("Moderate",count_level("alcohol",2),"alcohol_2"),
  add_coef("Heavy",count_level("alcohol",3),"alcohol_3"),
  add_coef("Other/Unknown",count_level("alcohol",4),"alcohol_4"),
  add_header("Smoking Status (Ref: Never)",count_level("smoking",0)),
  add_coef("Former",count_level("smoking",1),"smoking_1"),
  add_coef("Current",count_level("smoking",2),"smoking_2"),
  add_coef("Reduced Mobility",count_level("phyfunc",1),"physical_function_1")
)

forest_df <- forest_df |>
  mutate(
    coef = coef_values[coef_name],
    se = se_values[coef_name],
    hr = exp(coef),
    ci_low = exp(coef - 1.96 * se),
    ci_high = exp(coef + 1.96 * se),
    hr_text = ifelse(
      is.na(hr),
      "",
      sprintf("%.2f (%.2f-%.2f)",hr,ci_low,ci_high)
    ),
    row_number = seq_len(n()),
    shade = row_number %% 2 == 1
  )

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
  ticks_at = c(0,1,2,3),
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

forest_plot <- edit_plot(
  forest_plot,
  row = which(forest_df$is_header),
  col = seq_along(plot_table),
  which = "text",
  gp = gpar(fontface = "bold")
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
png(
  paste0(output_prefix,".png"),
  width = 8.5,
  height = 12,
  units = "in",
  res = 300
)
plot(forest_plot,autofit = TRUE)
dev.off()
