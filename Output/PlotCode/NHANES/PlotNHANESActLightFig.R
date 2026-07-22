#!/usr/bin/env Rscript

library(ggplot2)

input_file <- file.path("Data","JMHMMFitMix5Seed.rda")
output_prefix <- file.path("Output","Figures","nhanes_act_light_fig_weekday")
day_type_to_plot <- 1
day_type_label <- "Weekday"
grid_size <- 500
lod_act <- -5.52
lod_light <- -1.56

load(input_file)

est_params <- to_save$est_params
emit_act <- est_params$emit_act
emit_light <- est_params$emit_light
corr_mat <- est_params$corr_mat
mix_num <- dim(emit_act)[3]

bivariate_normal_density <- function(x,y,mu_x,sd_x,mu_y,sd_y,rho){
  x_std <- (x - mu_x) / sd_x
  y_std <- (y - mu_y) / sd_y
  rho_term <- 1 - rho^2

  exp(
    -(x_std^2 - 2 * rho * x_std * y_std + y_std^2) /
      (2 * rho_term)
  ) /
    (2 * pi * sd_x * sd_y * sqrt(rho_term))
}

activity_grid <- seq(
  min(emit_act[,1,,day_type_to_plot]) -
    3.5 * max(emit_act[,2,,day_type_to_plot]),
  max(emit_act[,1,,day_type_to_plot]) +
    3.5 * max(emit_act[,2,,day_type_to_plot]),
  length.out = grid_size
)

light_grid <- seq(
  min(emit_light[,1,,day_type_to_plot]) -
    3.5 * max(emit_light[,2,,day_type_to_plot]),
  max(emit_light[,1,,day_type_to_plot]) +
    3.5 * max(emit_light[,2,,day_type_to_plot]),
  length.out = grid_size
)

grid_df <- expand.grid(
  Activity = activity_grid,
  Light = light_grid
)

plot_rows <- list()
for (class_index in seq_len(mix_num)){
  for (state_index in 1:2){
    density <- bivariate_normal_density(
      x = grid_df$Activity,
      y = grid_df$Light,
      mu_x = emit_act[state_index,1,class_index,day_type_to_plot],
      sd_x = emit_act[state_index,2,class_index,day_type_to_plot],
      mu_y = emit_light[state_index,1,class_index,day_type_to_plot],
      sd_y = emit_light[state_index,2,class_index,day_type_to_plot],
      rho = corr_mat[class_index,state_index,day_type_to_plot]
    )

    plot_rows[[length(plot_rows) + 1L]] <- data.frame(
      Activity = grid_df$Activity,
      Light = grid_df$Light,
      Density = density / max(density),
      state = ifelse(state_index == 1,"Wake","Sleep"),
      latent_class = paste("Class",class_index)
    )
  }
}

plot_df <- do.call(rbind,plot_rows)
plot_df$state <- factor(plot_df$state,levels = c("Wake","Sleep"))
plot_df$latent_class <- factor(
  plot_df$latent_class,
  levels = paste("Class",seq_len(mix_num))
)

act_light_plot <- ggplot(
  plot_df,
  aes(x = Activity,y = Light,z = Density)
) +
  geom_contour(aes(color = after_stat(level)),bins = 150,size = 0.35) +
  geom_hline(yintercept = lod_light,linetype = "dashed",size = 0.5) +
  geom_vline(xintercept = lod_act,linetype = "dashed",size = 0.5) +
  scale_color_viridis_c(limits = c(0,1),end = .95) +
  facet_grid(latent_class ~ state,scales = "free_x") +
  theme_bw() +
  labs(
    title = paste(day_type_label,"Activity and Light by Mixture and Wake/Sleep"),
    x = "Activity (Log MIMS)",
    y = "Light (Log Lux)",
    color = "Density"
  )

print(act_light_plot)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(
  paste0(output_prefix,".png"),
  act_light_plot,
  width = 10,
  height = 6,
  dpi = 300
)
