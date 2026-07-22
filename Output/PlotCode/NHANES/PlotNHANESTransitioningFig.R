#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)

model_file <- file.path("Data","JMHMMFitMix5Seed.rda")
output_prefix <- file.path("Output","Figures","nhanes_transitioning_fig")

load(model_file)

tran_df <- to_save$est_params$tran_df

stopifnot(
  is.data.frame(tran_df),
  all(c("prob","type","time","weekend","mixture") %in% names(tran_df))
)

plot_df <- tran_df |>
  filter(weekend == 1) |>
  mutate(
    transition = ifelse(type == "Falling Asleep","Wake to Sleep","Sleep to Wake"),
    transition = factor(transition,levels = c("Wake to Sleep","Sleep to Wake")),
    latent_class = factor(mixture)
  )

transitioning_plot <- ggplot(
  plot_df,
  aes(x = time,y = prob,color = latent_class)
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ transition,nrow = 1) +
  scale_color_viridis_d(end = .95,name = "Latent Class") +
  scale_x_continuous(
    limits = c(0,25),
    breaks = seq(0,25,5),
    name = "Hour"
  ) +
  scale_y_continuous(
    limits = c(0,.375),
    breaks = seq(0,.3,.1),
    name = "Probability of Transitioning States"
  ) +
  labs(title = "Probability of Transitioning States by Hour and Latent Class") +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(size = 25),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    panel.grid.minor = element_line(color = "grey90"),
    legend.position = "right"
  )

print(transitioning_plot)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(
  paste0(output_prefix,".png"),
  transitioning_plot,
  width = 10.5,
  height = 5,
  dpi = 300
)
