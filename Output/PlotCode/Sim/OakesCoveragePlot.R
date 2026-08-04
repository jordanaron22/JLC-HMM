#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)

input_file <- file.path("Output","parse_oakes_results.rds")
scenario_output_prefix <- file.path(
  "Output","Figures","sim_oakes_scenario_coverage"
)
class_output_prefix <- file.path(
  "Output","Figures","sim_oakes_class_beta_coverage"
)


title_case <- function(x){
  paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
}

make_days_label <- function(days){
  paste0(days," Day",ifelse(days == 1,"","s"))
}

make_model_label <- function(model_type){
  ifelse(model_type == "joint","Joint Oakes","Two-Stage H1")
  if (model_type == "joint"){
    "Joint Oaks"
  } else if (model_type == "two_stage naive"){
    "Two-Stage Naive"
  } else if (model_type == "two_stage MT"){
    "Two-Stage MT"
  } else {
    "Unknown Model Type"
  }
}


x <- readRDS(input_file)

cov_summary <- x$class_beta_summary


wide_cov <- cov_summary %>%
  select(
    simulation_days, num_people, true_mix_num, fit_mix_num,
    emission_overlap, param_name, model_type,
    coverage, mean_bias, rmse, median_se, se_method
  ) %>%
  pivot_wider(
    names_from = c(model_type,se_method),
    values_from = c(coverage, mean_bias, rmse, median_se)
  )




# wide_cov$coverage_diff_joint_minus_two_stage <-
#   wide_cov$coverage.joint - wide_cov$coverage.two_stage

scenario_cov <- wide_cov %>%
  group_by(simulation_days, emission_overlap) %>%
  summarise(
    across(
      c(
        coverage_joint_oakes_schur,
        coverage_joint_joint_naive,
        coverage_two_stage_non_oakes_beta_se,
        coverage_two_stage_two_stage_murphy_topel,
        rmse_joint_oakes_schur,
        rmse_joint_joint_naive,
        rmse_two_stage_non_oakes_beta_se,
        rmse_two_stage_two_stage_murphy_topel,
        median_se_joint_oakes_schur,
        median_se_joint_joint_naive,
        median_se_two_stage_non_oakes_beta_se,
        median_se_two_stage_two_stage_murphy_topel
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )



make_coverage_plot_data <- function(wide_data){
  joint_data <- wide_data[,c("simulation_days","emission_overlap",
                             "coverage_joint_oakes_schur"),drop = FALSE]
  names(joint_data)[names(joint_data) == "coverage_joint_oakes_schur"] <- "coverage"
  joint_data$model_type <- "Joint"
  
  joint_data_naive <- wide_data[,c("simulation_days","emission_overlap",
                             "coverage_joint_joint_naive"),drop = FALSE]
  names(joint_data_naive)[names(joint_data_naive) == "coverage_joint_joint_naive"] <- "coverage"
  joint_data_naive$model_type <- "Joint Naive"

  two_stage_data_naive <- wide_data[,c("simulation_days","emission_overlap",
                                 "coverage_two_stage_non_oakes_beta_se"),drop = FALSE]
  names(two_stage_data_naive)[names(two_stage_data_naive) == "coverage_two_stage_non_oakes_beta_se"] <-
    "coverage"
  two_stage_data_naive$model_type <- "Two-Stage Naive"

  two_stage_data_MT <- wide_data[,c("simulation_days","emission_overlap",
                                 "coverage_two_stage_two_stage_murphy_topel"),drop = FALSE]
  names(two_stage_data_MT)[names(two_stage_data_MT) == "coverage_two_stage_two_stage_murphy_topel"] <-
    "coverage"
  two_stage_data_MT$model_type <- "Two-Stage MT"

  plot_data <- rbind(joint_data,joint_data_naive,two_stage_data_naive,two_stage_data_MT)
  plot_data <- plot_data[is.finite(plot_data$coverage),,drop = FALSE]

  day_values <- sort(unique(plot_data$simulation_days))
  plot_data$days_label <- factor(
    make_days_label(plot_data$simulation_days),
    levels = make_days_label(day_values)
  )

  overlap_values <- sort(unique(plot_data$emission_overlap))
  plot_data$overlap_label <- factor(
    paste0("Overlap: ",title_case(plot_data$emission_overlap)),
    levels = paste0("Overlap: ",title_case(overlap_values))
  )

  plot_data
}

scenario_plot_data <- make_coverage_plot_data(scenario_cov)

scenario_coverage_plot <- ggplot(
  scenario_plot_data,
  aes(x = days_label,y = coverage,fill = model_type)
) +
  geom_col(position = position_dodge(width = 0.8),width = 0.7,
           color = "black") +
  geom_hline(yintercept = 0.95,linetype = "dashed",color = "gray35") +
  facet_grid(rows = vars(overlap_label)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_fill_viridis_d(end = .85,name = "Model") +
  labs(x = "Simulation Days",
       y = "Coverage",
       fill = "Model",
       title = "Average Survival Class-Beta Coverage") +
  theme_bw()

print(scenario_coverage_plot)

make_class_coverage_plot_data <- function(wide_data){
  id_cols <- c("simulation_days","emission_overlap","param_name")

  joint_data <- wide_data[,c(id_cols,"coverage_joint_oakes_schur"),drop = FALSE]
  names(joint_data)[names(joint_data) == "coverage_joint_oakes_schur"] <- "coverage"
  joint_data$model_type <- "Joint"

  joint_data_naive <- wide_data[,c(id_cols,"coverage_joint_joint_naive"),drop = FALSE]
  names(joint_data_naive)[names(joint_data_naive) == "coverage_joint_joint_naive"] <- "coverage"
  joint_data_naive$model_type <- "Joint Naive"

  two_stage_data_naive <- wide_data[,c(id_cols,"coverage_two_stage_non_oakes_beta_se"),drop = FALSE]
  names(two_stage_data_naive)[names(two_stage_data_naive) == "coverage_two_stage_non_oakes_beta_se"] <-
    "coverage"
  two_stage_data_naive$model_type <- "Two-Stage Naive"

  two_stage_data_MT <- wide_data[,c(id_cols,"coverage_two_stage_two_stage_murphy_topel"),drop = FALSE]
  names(two_stage_data_MT)[names(two_stage_data_MT) == "coverage_two_stage_two_stage_murphy_topel"] <-
    "coverage"
  two_stage_data_MT$model_type <- "Two-Stage MT"

  plot_data <- rbind(joint_data,joint_data_naive,two_stage_data_naive, two_stage_data_MT)
  plot_data <- plot_data[is.finite(plot_data$coverage),,drop = FALSE]

  day_values <- sort(unique(plot_data$simulation_days))
  plot_data$days_label <- factor(
    make_days_label(plot_data$simulation_days),
    levels = make_days_label(day_values)
  )

  overlap_values <- sort(unique(plot_data$emission_overlap))
  plot_data$overlap_label <- factor(
    paste0("Overlap: ",title_case(plot_data$emission_overlap)),
    levels = paste0("Overlap: ",title_case(overlap_values))
  )

  plot_data$param_label <- factor(
    gsub("class_","Class ",plot_data$param_name),
    levels = gsub("class_","Class ",sort(unique(plot_data$param_name)))
  )

  plot_data
}

class_plot_data <- make_class_coverage_plot_data(wide_cov)

class_coverage_plot <- ggplot(
  class_plot_data,
  aes(x = param_label,y = coverage,fill = model_type)
) +
  geom_col(position = position_dodge(width = 0.8),width = 0.7,
           color = "black") +
  geom_hline(yintercept = 0.95,linetype = "dashed",color = "gray35") +
  facet_grid(rows = vars(overlap_label),cols = vars(days_label)) +
  coord_cartesian(ylim = c(0,1)) +
  scale_fill_viridis_d(end = .85,name = "Model") +
  labs(x = "Latent Class Survival Coefficient",
       y = "Coverage",
       fill = "Model",
       title = "Survival Class-Beta Coverage by Coefficient") +
  theme_bw()

print(class_coverage_plot)

dir.create(dirname(scenario_output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(
  paste0(scenario_output_prefix,".png"),
  scenario_coverage_plot,
  width = 10,
  height = 6,
  dpi = 300
)
ggsave(
  paste0(class_output_prefix,".png"),
  class_coverage_plot,
  width = 12,
  height = 7,
  dpi = 300
)
