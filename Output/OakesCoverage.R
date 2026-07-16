add_local_lib <- function(path){
  if (dir.exists(path)){
    .libPaths(unique(c(normalizePath(path,winslash = "/",mustWork = TRUE),
                       .libPaths())))
  }
}

for (lib_path in c(".Rlibrary","Rlib")){
  add_local_lib(lib_path)
}

library(ggplot2)

title_case <- function(x){
  paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
}

make_days_label <- function(days){
  paste0(days," Day",ifelse(days == 1,"","s"))
}

make_model_label <- function(model_type){
  ifelse(model_type == "joint","Joint Oakes","Two-Stage H1")
}

x <- readRDS("Output/parse_oakes_results.rds")

cov_summary <- x$class_beta_summary

cov_summary[
  , c("model_type", "se_method", "simulation_days", "num_people",
      "emission_overlap", "param_name", "n_valid", "coverage",
      "mean_bias", "rmse", "median_se")
]

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

wide_cov[
  , c("simulation_days", "emission_overlap", "param_name",
      "coverage.joint", "coverage.two_stage",
      "coverage_diff_joint_minus_two_stage",
      "rmse.joint", "rmse.two_stage",
      "median_se.joint", "median_se.two_stage")
]

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

scenario_cov


############################

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
  geom_col(position = position_dodge(width = 0.8),width = 0.7) +
  geom_hline(yintercept = 0.95,linetype = "dashed",color = "gray35") +
  facet_grid(rows = vars(overlap_label)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Simulation Days",
       y = "Coverage",
       fill = "Model",
       title = "Average Survival Class-Beta Coverage") +
  theme_bw()

scenario_coverage_plot

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
  geom_col(position = position_dodge(width = 0.8),width = 0.7) +
  geom_hline(yintercept = 0.95,linetype = "dashed",color = "gray35") +
  facet_grid(rows = vars(overlap_label),cols = vars(days_label)) +
  coord_cartesian(ylim = c(0,1)) +
  labs(x = "Latent Class Survival Coefficient",
       y = "Coverage",
       fill = "Model",
       title = "Survival Class-Beta Coverage by Coefficient") +
  theme_bw()

class_coverage_plot

# dir.create(file.path("Output","Plots"),recursive = TRUE,showWarnings = FALSE)
# ggsave("oakes_scenario_coverage.png",
#        plot = scenario_coverage_plot,
#        path = file.path("Output","Plots"),
#        width = 10,height = 6,dpi = 300)
# ggsave("oakes_class_beta_coverage.png",
#        plot = class_coverage_plot,
#        path = file.path("Output","Plots"),
#        width = 12,height = 7,dpi = 300)
