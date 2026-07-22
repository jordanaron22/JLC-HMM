#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)

input_file <- file.path("Output","parse_nhanes_cv.rds")
output_prefix <- file.path("Output","Figures","nhanes_cv_model_selection_joint")
model_type_to_plot <- "joint"

parsed_cv <- readRDS(input_file)

model_selection_df <- parsed_cv$cv_model_selection |>
  filter(selection_group == model_type_to_plot) |>
  select(
    fit_mix_num,
    partial_loglik =
      cv_partial_loglik_pooled_weighted,
    pooled_oof_cindex =
      cindex_scenario6_pooled_oof,
    weighted_ibs =
      cv_ibs_pooled_weighted
  )

diagnostic_info <- data.frame(
  diagnostic = c(
    "partial_loglik",
    "pooled_oof_cindex",
    "weighted_ibs"
  ),
  diagnostic_label = c(
    "Partial Weighted Log-Likelihood",
    "Pooled OOF Weighted C-index",
    "Pooled OOF Weighted IBS"
  ),
  best_direction = c("max","max","min")
)

plot_df <- model_selection_df |>
  pivot_longer(
    cols = -fit_mix_num,
    names_to = "diagnostic",
    values_to = "score"
  ) |>
  left_join(diagnostic_info,by = "diagnostic") |>
  filter(is.finite(score)) |>
  mutate(
    diagnostic_label = factor(
      diagnostic_label,
      levels = diagnostic_info$diagnostic_label
    ),
    fit_mix_num_position = as.numeric(factor(fit_mix_num))
  ) |>
  group_by(diagnostic) |>
  mutate(
    score_range = max(score,na.rm = TRUE) - min(score,na.rm = TRUE),
    score_range = if_else(score_range > 0,score_range,1),
    # Bars start near the smallest score so each facet can show differences.
    bar_base = min(score,na.rm = TRUE) - 0.05 * score_range
  ) |>
  ungroup()

best_df <- plot_df |>
  group_by(diagnostic,diagnostic_label,best_direction) |>
  summarise(
    best_score = if_else(
      first(best_direction) == "min",
      min(score,na.rm = TRUE),
      max(score,na.rm = TRUE)
    ),
    .groups = "drop"
  )

model_selection_plot <- ggplot(plot_df) +
  geom_rect(
    aes(
      xmin = fit_mix_num_position - 0.4,
      xmax = fit_mix_num_position + 0.4,
      ymin = bar_base,
      ymax = score
    ),
    fill = "grey45"
  ) +
  geom_hline(
    data = best_df,
    aes(yintercept = best_score),
    color = "red",
    linetype = "dashed"
  ) +
  scale_x_continuous(
    breaks = sort(unique(plot_df$fit_mix_num_position)),
    labels = sort(unique(plot_df$fit_mix_num))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  facet_wrap(~ diagnostic_label,scales = "free_y",ncol = 3) +
  theme_bw() +
  labs(
    x = "Number of latent classes",
    y = "Score"
  )

print(model_selection_plot)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(paste0(output_prefix,".png"),model_selection_plot,
       width = 9,height = 3,dpi = 300)
