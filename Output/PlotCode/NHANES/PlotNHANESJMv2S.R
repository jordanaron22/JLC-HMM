#!/usr/bin/env Rscript

library(dplyr)
library(forestploter)
library(grid)

jlcm_file <- file.path("Data","JMHMMFitMix5Seed.rda")
two_stage_file <- file.path("Data","JMHMMNoSurvFitMix5Seed.rda")
output_prefix <- file.path("Output","Figures","nhanes_jm_vs_two_stage")

read_class_hr <- function(model_file,model_label){
  load_env <- new.env()
  load(model_file,envir = load_env)
  est_params <- load_env$to_save$est_params

  beta_vec <- est_params$beta_vec
  beta_se <- est_params$beta_se
  mix_num <- length(beta_vec)

  class_coef <- beta_vec[2:mix_num]
  class_se <- if (!is.null(names(beta_se)) &&
                  all(paste0("latent_class",2:mix_num) %in% names(beta_se))){
    beta_se[paste0("latent_class",2:mix_num)]
  } else {
    beta_se[2:mix_num]
  }

  data.frame(
    latent_class = 2:mix_num,
    model = model_label,
    coef = as.numeric(class_coef),
    se = as.numeric(class_se)
  ) |>
    mutate(
      hr = exp(coef),
      ci_low = exp(coef - 1.96 * se),
      ci_high = exp(coef + 1.96 * se),
      hr_text = sprintf("%.2f (%.2f-%.2f)",hr,ci_low,ci_high)
    )
}

load(jlcm_file)
class_size <- tabulate(
  apply(to_save$est_params$re_prob,1,which.max),
  nbins = length(to_save$est_params$beta_vec)
)

jlcm_df <- read_class_hr(jlcm_file,"JLCM") |>
  select(
    latent_class,
    jlcm_hr = hr,
    jlcm_low = ci_low,
    jlcm_high = ci_high,
    jlcm_text = hr_text
  )

two_stage_df <- read_class_hr(two_stage_file,"Two-Stage") |>
  select(
    latent_class,
    two_stage_hr = hr,
    two_stage_low = ci_low,
    two_stage_high = ci_high,
    two_stage_text = hr_text
  )

forest_df <- merge(jlcm_df,two_stage_df,by = "latent_class") |>
  arrange(latent_class) |>
  mutate(
    `Latent Class` = paste("Class",latent_class),
    Size = as.character(class_size[latent_class]),
    ` ` = paste(rep(" ",24),collapse = " "),
    `JLCM HR (95% CI)` = jlcm_text,
    `Two-Stage HR (95% CI)` = two_stage_text
  )

plot_table <- forest_df |>
  select(
    `Latent Class`,
    Size,
    ` `,
    `JLCM HR (95% CI)`,
    `Two-Stage HR (95% CI)`
  )

plot_theme <- forest_theme(
  base_size = 11,
  ci_pch = 15,
  ci_col = c("#440154","#7AD151"),
  ci_fill = c("#440154","#7AD151"),
  ci_lwd = 1.2,
  ci_Theight = 0.15,
  refline_gp = gpar(lty = "dashed",col = "grey35",lwd = 1),
  xaxis_gp = gpar(fontsize = 10),
  xlab_gp = gpar(fontsize = 11),
  legend_name = "Model",
  legend_value = c("JLCM","Two-Stage"),
  legend_position = "right"
)

jm_vs_two_stage_plot <- forest(
  plot_table,
  est = list(forest_df$jlcm_hr,forest_df$two_stage_hr),
  lower = list(forest_df$jlcm_low,forest_df$two_stage_low),
  upper = list(forest_df$jlcm_high,forest_df$two_stage_high),
  ci_column = 3,
  ref_line = 1,
  xlim = c(0.5,3),
  ticks_at = c(1,2,3),
  ticks_digits = 1L,
  xlab = "Hazard Ratio",
  nudge_y = 0.18,
  theme = plot_theme
)

jm_vs_two_stage_plot <- edit_plot(
  jm_vs_two_stage_plot,
  row = c(1,3),
  col = seq_along(plot_table),
  which = "background",
  gp = gpar(fill = "#edf2f0",col = NA)
)

jm_vs_two_stage_plot <- add_border(
  jm_vs_two_stage_plot,
  part = "header",
  where = "bottom",
  gp = gpar(lwd = 1)
)

plot(jm_vs_two_stage_plot,autofit = TRUE)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
png(
  paste0(output_prefix,".png"),
  width = 11,
  height = 3,
  units = "in",
  res = 300
)
plot(jm_vs_two_stage_plot,autofit = TRUE)
dev.off()
