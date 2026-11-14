#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)

input_file <- file.path("Output","parse_oakes_results.rds")
simulation_input_file <- file.path("Output","parse_sim_results.rds")
simulation_fit_root <- file.path("Output","Routputs","RoutputsSim")
scenario_output_prefix <- file.path(
  "Output","Figures","sim_oakes_scenario_coverage"
)
class_output_prefix <- file.path(
  "Output","Figures","sim_oakes_class_beta_coverage"
)


title_case <- function(x){
  x <- as.character(x)
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


get_two_stage_naive_summary <- function(simulation_results,target_scenarios,
                                        fit_root,max_seed = 250){
  scenario_cols <- c(
    "simulation_days","num_people","true_mix_num","fit_mix_num",
    "emission_overlap"
  )

  two_stage_runs <- simulation_results %>%
    filter(
      model_type == "two_stage",
      sim_num >= 1,
      sim_num <= max_seed
    ) %>%
    semi_join(target_scenarios,by = scenario_cols)

  fit_files <- list.files(
    fit_root,
    pattern = "[.]rda$",
    recursive = TRUE,
    full.names = TRUE
  )
  duplicate_basenames <- unique(
    basename(fit_files)[duplicated(basename(fit_files))]
  )
  if (length(duplicate_basenames) > 0){
    stop(
      "Duplicate simulation-fit basenames found under ",fit_root,": ",
      paste(utils::head(duplicate_basenames,5),collapse = ", ")
    )
  }
  fit_file_by_basename <- setNames(fit_files,basename(fit_files))

  coefficient_rows <- vector("list",nrow(two_stage_runs))
  missing_files <- character()
  for (run_ind in seq_len(nrow(two_stage_runs))){
    saved_file <- unname(
      fit_file_by_basename[basename(two_stage_runs$file[[run_ind]])]
    )
    if (length(saved_file) == 0 || is.na(saved_file)){
      missing_files <- c(
        missing_files,
        basename(two_stage_runs$file[[run_ind]])
      )
      next
    }

    load_env <- new.env(parent = emptyenv())
    load(saved_file,envir = load_env)
    saved <- load_env[["to_save"]]
    if (is.null(saved)){
      stop("Simulation fit does not contain to_save: ",saved_file)
    }

    true_beta <- suppressWarnings(as.numeric(
      saved[["true_params"]][["beta_vec"]]
    ))
    estimated_beta <- suppressWarnings(as.numeric(
      saved[["est_params"]][["beta_vec"]]
    ))
    naive_se <- suppressWarnings(as.numeric(
      saved[["est_params"]][["beta_se"]]
    ))
    fit_mix_num <- two_stage_runs$fit_mix_num[[run_ind]]
    class_indices <- 2:fit_mix_num

    if (
      length(true_beta) < fit_mix_num ||
      length(estimated_beta) < fit_mix_num ||
      length(naive_se) < fit_mix_num
    ){
      stop("Incomplete class-beta estimates or SEs in ",saved_file)
    }

    run_data <- two_stage_runs[
      rep(run_ind,length(class_indices)),
      c(scenario_cols,"sim_num"),
      drop = FALSE
    ]
    run_data$class_index <- class_indices
    run_data$param_name <- paste0("class_",class_indices)
    run_data$estimate <- estimated_beta[class_indices]
    run_data$truth <- true_beta[class_indices]
    run_data$se <- naive_se[class_indices]
    coefficient_rows[[run_ind]] <- run_data
  }

  if (length(missing_files) > 0){
    warning(
      "Skipped ",length(missing_files),
      " two-stage runs whose simulation-fit file was not found."
    )
  }

  coefficient_data <- bind_rows(coefficient_rows) %>%
    mutate(
      bias = estimate - truth,
      valid = is.finite(estimate) & is.finite(truth) & is.finite(se),
      covered = truth >= estimate - qnorm(.975) * se &
        truth <= estimate + qnorm(.975) * se
    )

  if (nrow(coefficient_data) == 0){
    stop("No two-stage simulation fits were available for naive coverage")
  }

  run_counts <- coefficient_data %>%
    distinct(across(all_of(c(scenario_cols,"sim_num")))) %>%
    count(across(all_of(scenario_cols)),name = "n_runs")
  incomplete_scenarios <- run_counts %>% filter(n_runs < max_seed)
  if (nrow(incomplete_scenarios) > 0){
    warning(
      "Some two-stage scenarios have fewer than ",max_seed,
      " successful saved fits; coverage uses every available fit."
    )
  }

  coefficient_data %>%
    group_by(across(all_of(c(scenario_cols,"param_name")))) %>%
    summarise(
      model_type = "two_stage",
      se_method = "two_stage_naive",
      n = n(),
      n_valid = sum(valid),
      coverage = if (any(valid)) mean(covered[valid]) else NA_real_,
      mean_bias = if (any(valid)) mean(bias[valid]) else NA_real_,
      rmse = if (any(valid)) sqrt(mean(bias[valid]^2)) else NA_real_,
      median_se = if (any(valid)) median(se[valid]) else NA_real_,
      .groups = "drop"
    )
}


x <- readRDS(input_file)

cov_summary <- x$class_beta_summary

# The Oakes results directory contains only a subset of the two-stage fits.
# Recalculate the naive two-stage intervals from all seeds 1--250 in the
# ordinary simulation output, then replace the incomplete parsed rows.
simulation_results <- readRDS(simulation_input_file)
target_scenarios <- cov_summary %>%
  distinct(
    simulation_days,num_people,true_mix_num,fit_mix_num,emission_overlap
  )
two_stage_naive_summary <- get_two_stage_naive_summary(
  simulation_results,
  target_scenarios,
  simulation_fit_root,
  max_seed = 250
)
cov_summary <- cov_summary %>%
  filter(!(model_type == "two_stage" & se_method == "two_stage_naive")) %>%
  bind_rows(two_stage_naive_summary)


wide_cov <- cov_summary %>%
  dplyr::select(
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
        coverage_two_stage_two_stage_naive,
        coverage_two_stage_two_stage_murphy_topel,
        rmse_joint_oakes_schur,
        rmse_joint_joint_naive,
        rmse_two_stage_two_stage_naive,
        rmse_two_stage_two_stage_murphy_topel,
        median_se_joint_oakes_schur,
        median_se_joint_joint_naive,
        median_se_two_stage_two_stage_naive,
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
                                 "coverage_two_stage_two_stage_naive"),drop = FALSE]
  names(two_stage_data_naive)[names(two_stage_data_naive) == "coverage_two_stage_two_stage_naive"] <-
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

  overlap_values <- c("low","mid","high")
  plot_data$emission_overlap <- factor(
    plot_data$emission_overlap,
    levels = overlap_values
  )
  plot_data$overlap_label <- factor(
    paste0("Overlap: ",title_case(plot_data$emission_overlap)),
    levels = paste0("Overlap: ",title_case(overlap_values))
  )

  plot_data
}

scenario_plot_data <- make_coverage_plot_data(scenario_cov)

# library(stringr)
# scenario_plot_data <- scenario_plot_data %>% 
#   filter(model_type == "Joint" | model_type == "Two-Stage Naive") %>% 
#   filter(overlap_label == 'Overlap: Low') %>% 
#   mutate(model_type = str_replace_all(model_type, "Two-Stage Naive", "Two-Stage"))%>% 
#   mutate(model_type = str_replace_all(model_type, "Joint", "JLCM"))

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
       title = "Latent Class Survival Coefficient Coverage") +
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

  two_stage_data_naive <- wide_data[,c(id_cols,"coverage_two_stage_two_stage_naive"),drop = FALSE]
  names(two_stage_data_naive)[names(two_stage_data_naive) == "coverage_two_stage_two_stage_naive"] <-
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
# ggsave(
#   paste0(scenario_output_prefix,".png"),
#   scenario_coverage_plot,
#   width = 10,
#   height = 6,
#   dpi = 300
# )
# ggsave(
#   paste0(class_output_prefix,".png"),
#   class_coverage_plot,
#   width = 12,
#   height = 7,
#   dpi = 300
# )


# scenario_plot_data %>% filter(emission_overlap == "high", simulation_days == 7)

