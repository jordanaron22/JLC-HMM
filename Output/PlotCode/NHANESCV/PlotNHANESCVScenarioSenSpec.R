#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

input_file <- file.path("Output","parse_nhanes_cv_scenarios.rds")
output_prefix <- file.path("Output","Figures","nhanes_cv_scenario_senspec")
model_type_to_plot <- "joint"
scenarios_to_plot <- c("Sleep-Wake Cycle","Null Transition","Standard")

parsed_cv <- readRDS(input_file)

plot_df <- parsed_cv$macro_summary |>
  filter(model_type %in% model_type_to_plot) |>
  filter(testing_data  %in% scenarios_to_plot) |>
  select(
    fit_mix_num,
    testing_data,
    sensitivity = macro_sensitivity,
    specificity = macro_specificity
  ) |>
  pivot_longer(
    cols = c(sensitivity,specificity),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = recode(
      metric,
      sensitivity = "Sensitivity",
      specificity = "Specificity"
    ),
    metric = factor(metric,levels = c("Sensitivity","Specificity")),
    testing_data = factor(testing_data,levels = scenarios_to_plot)
  )

scenario_senspec_plot <- ggplot(
  plot_df,
  aes(x = metric,y = value,fill = testing_data)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.75,
    color = "black"
  ) +
  scale_fill_viridis_d(end = .85) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw() +
  labs(
    title = "Cross-Validated Mixture Sensitivity/Specificity by Testing Data",
    x = NULL,
    y = "Specificity / Sensitivity",
    fill = "Testing Data"
  )

print(scenario_senspec_plot)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(
  paste0(output_prefix,".png"),
  scenario_senspec_plot,
  width = 9,
  height = 6,
  dpi = 300
)
