#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

input_file <- file.path("Output","parse_nhanes_cv_scenarios.rds")
output_prefix <- file.path("Output","Figures","nhanes_cv_light_sensitivity")
model_type_to_plot <- "joint"
scenarios_to_plot <- c("No Light","Standard")
activity_strata_to_plot <- c(
  "Activity Above Individual Median",
  "Activity Below Individual Median",
  "Total"
)

parsed_cv <- readRDS(input_file)

plot_df <- parsed_cv$state_sensitivity_summary |>
  filter(model_type %in% model_type_to_plot) |>
  filter(testing_data %in% scenarios_to_plot) |>
  filter(activity_stratum %in% activity_strata_to_plot) |>
  mutate(
    testing_data = factor(testing_data,levels = scenarios_to_plot),
    activity_stratum = factor(
      activity_stratum,
      levels = activity_strata_to_plot
    ),
    state = factor(state,levels = c("Sleep","Wake"))
  )

light_sensitivity_plot <- ggplot(
  plot_df,
  aes(x = state,y = sensitivity,fill = testing_data)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.75,
    color = "black"
  ) +
  scale_fill_viridis_d(end = .85) +
  coord_cartesian(ylim = c(0,1)) +
  facet_wrap(~ activity_stratum,nrow = 1) +
  theme_bw() +
  labs(
    title = "Cross-Validated Latent State Sensitivity by Testing Data",
    x = NULL,
    y = "Sensitivity",
    fill = "Testing Data"
  )

print(light_sensitivity_plot)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(
  paste0(output_prefix,".png"),
  light_sensitivity_plot,
  width = 9,
  height = 3,
  dpi = 300
)
