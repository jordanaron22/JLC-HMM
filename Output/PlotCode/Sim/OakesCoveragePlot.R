#!/usr/bin/env Rscript

library(ggplot2)

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
}


x <- readRDS(input_file)

cov_summary <- x$class_beta_summary

wide_cov <- reshape(
  cov_summary[
    , c("simulation_days", "num_people", "true_mix_num", "fit_mix_num",
        "emission_overlap", "param_name", "model_type", "coverage",
        "mean_bias", "rmse", "median_se")
  ],
  idvar = c("simulation_days", "num_people", "true_mix_num", "fit_mix_num",
            "emission_overlap", "param_name"),
  timevar = "model_type",
  direction = "wide"
)

wide_cov$coverage_diff_joint_minus_two_stage <-
  wide_cov$coverage.joint - wide_cov$coverage.two_stage

scenario_cov <- aggregate(
  cbind(
    coverage.joint,
    coverage.two_stage,
    coverage_diff_joint_minus_two_stage,
    rmse.joint,
    rmse.two_stage,
    median_se.joint,
    median_se.two_stage
  ) ~ simulation_days + emission_overlap,
  data = wide_cov,
  FUN = function(x) mean(x,na.rm = TRUE)
)

make_coverage_plot_data <- function(wide_data){
  joint_data <- wide_data[,c("simulation_days","emission_overlap",
                             "coverage.joint"),drop = FALSE]
  names(joint_data)[names(joint_data) == "coverage.joint"] <- "coverage"
  joint_data$model_type <- "joint"

  two_stage_data <- wide_data[,c("simulation_days","emission_overlap",
                                 "coverage.two_stage"),drop = FALSE]
  names(two_stage_data)[names(two_stage_data) == "coverage.two_stage"] <-
    "coverage"
  two_stage_data$model_type <- "two_stage"

  plot_data <- rbind(joint_data,two_stage_data)
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

  plot_data$model_type_label <- factor(
    make_model_label(plot_data$model_type),
    levels = c("Joint Oakes","Two-Stage H1")
  )

  plot_data
}

scenario_plot_data <- make_coverage_plot_data(scenario_cov)

scenario_coverage_plot <- ggplot(
  scenario_plot_data,
  aes(x = days_label,y = coverage,fill = model_type_label)
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

  joint_data <- wide_data[,c(id_cols,"coverage.joint"),drop = FALSE]
  names(joint_data)[names(joint_data) == "coverage.joint"] <- "coverage"
  joint_data$model_type <- "joint"

  two_stage_data <- wide_data[,c(id_cols,"coverage.two_stage"),drop = FALSE]
  names(two_stage_data)[names(two_stage_data) == "coverage.two_stage"] <-
    "coverage"
  two_stage_data$model_type <- "two_stage"

  plot_data <- rbind(joint_data,two_stage_data)
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

  plot_data$model_type_label <- factor(
    make_model_label(plot_data$model_type),
    levels = c("Joint Oakes","Two-Stage H1")
  )

  plot_data
}

class_plot_data <- make_class_coverage_plot_data(wide_cov)

class_coverage_plot <- ggplot(
  class_plot_data,
  aes(x = param_label,y = coverage,fill = model_type_label)
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
