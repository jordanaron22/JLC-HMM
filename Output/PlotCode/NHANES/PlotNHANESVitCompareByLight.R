#!/usr/bin/env Rscript

library(ggplot2)
library(patchwork)

model_file <- file.path("Data","JMHMMFitMix5Seed.rda")
nhanes_file <- file.path("Data","NHANES_2011_2012_2013_2014.rda")
wave_g_file <- file.path("Data","Wavedata_G.rda")
wave_h_file <- file.path("Data","Wavedata_H.rda")
output_prefix <- file.path("Output","Figures","nhanes_vit_compare_by_light")

plot_ind <- 6001
window_start_index <- 65
window_length <- 96

first_day_to_week_ind <- function(first_day){
  weekday <- numeric(96)
  friday <- c(rep(0,68),rep(1,28))
  saturday <- numeric(96) + 1
  sunday <- c(rep(1,68),rep(0,28))

  if (first_day == 1){
    c(sunday,rep(weekday,4),friday,saturday,sunday,weekday)
  } else if (first_day == 2){
    c(rep(weekday,4),friday,saturday,sunday,rep(weekday,2))
  } else if (first_day == 3){
    c(rep(weekday,3),friday,saturday,sunday,rep(weekday,3))
  } else if (first_day == 4){
    c(rep(weekday,2),friday,saturday,sunday,rep(weekday,4))
  } else if (first_day == 5){
    c(weekday,friday,saturday,sunday,rep(weekday,4),friday)
  } else if (first_day == 6){
    c(friday,saturday,sunday,rep(weekday,4),friday,saturday)
  } else if (first_day == 7){
    c(saturday,sunday,rep(weekday,4),friday,saturday,sunday)
  } else {
    stop("PAXDAYWM must be in 1:7.")
  }
}

expit <- function(x){
  exp(x) / (1 + exp(x))
}

params_to_tran <- function(params_tran_array,time,latent_class,weekend_index){
  params_tran <- params_tran_array[,,weekend_index]
  if (dim(params_tran_array)[1] == 1){
    params_tran <- matrix(params_tran,nrow = 1)
  }

  param_matrix <- matrix(
    params_tran[latent_class,],
    nrow = 2,
    ncol = 3,
    byrow = TRUE
  )

  p12 <- param_matrix[1,1] +
    param_matrix[1,2] * cos(2 * pi * time / 96) +
    param_matrix[1,3] * sin(2 * pi * time / 96)

  p21 <- param_matrix[2,1] +
    param_matrix[2,2] * cos(2 * pi * time / 96) +
    param_matrix[2,3] * sin(2 * pi * time / 96)

  tran <- matrix(0,nrow = 2,ncol = 2)
  tran[1,2] <- expit(p12)
  tran[1,1] <- 1 - tran[1,2]
  tran[2,1] <- expit(p21)
  tran[2,2] <- 1 - tran[2,1]
  tran
}

activity_log_density <- function(obs,mu,sigma,lod_act){
  if (is.na(obs)){
    return(0)
  }
  if (obs <= lod_act + 1e-12){
    return(pnorm(lod_act,mean = mu,sd = sigma,log.p = TRUE))
  }
  dnorm(obs,mean = mu,sd = sigma,log = TRUE)
}

viterbi_activity_only <- function(
    ind,
    latent_class,
    act,
    vcovar_mat,
    init,
    params_tran_array,
    emit_act,
    lod_act
){
  day_length <- nrow(act)

  log_class <- matrix(NA,nrow = 2,ncol = day_length)
  for (time in seq_len(day_length)){
    weekend_index <- vcovar_mat[time,ind] + 1
    for (state_index in 1:2){
      log_class[state_index,time] <- activity_log_density(
        obs = act[time,ind],
        mu = emit_act[state_index,1,latent_class,weekend_index],
        sigma = emit_act[state_index,2,latent_class,weekend_index],
        lod_act = lod_act
      )
    }
  }

  viterbi_mat <- matrix(NA,nrow = 2,ncol = day_length)
  viterbi_mat[1,1] <- log(init[latent_class,1]) + log_class[1,1]
  viterbi_mat[2,1] <- log(init[latent_class,2]) + log_class[2,1]

  viterbi_ind_mat <- matrix(NA,nrow = 2,ncol = day_length)

  for (time in 2:day_length){
    tran <- params_to_tran(
      params_tran_array = params_tran_array,
      time = time,
      latent_class = latent_class,
      weekend_index = vcovar_mat[time,ind] + 1
    )

    viterbi_mat[1,time] <- log_class[1,time] +
      max(
        viterbi_mat[1,time - 1] + log(tran[1,1]),
        viterbi_mat[2,time - 1] + log(tran[2,1])
      )
    viterbi_mat[2,time] <- log_class[2,time] +
      max(
        viterbi_mat[1,time - 1] + log(tran[1,2]),
        viterbi_mat[2,time - 1] + log(tran[2,2])
      )

    viterbi_ind_mat[1,time] <- which.max(c(
      viterbi_mat[1,time - 1] + log(tran[1,1]),
      viterbi_mat[2,time - 1] + log(tran[2,1])
    ))
    viterbi_ind_mat[2,time] <- which.max(c(
      viterbi_mat[1,time - 1] + log(tran[1,2]),
      viterbi_mat[2,time - 1] + log(tran[2,2])
    ))
  }

  decoded_state <- c(which.max(viterbi_mat[,day_length]))
  for (time in day_length:2){
    decoded_state <- c(viterbi_ind_mat[decoded_state[1],time],decoded_state)
  }

  decoded_state - 1
}

load(nhanes_file)
lmf_data <- rbind(
  NHANES_mort_list[[1]][NHANES_mort_list[[1]]$eligstat == 1,],
  NHANES_mort_list[[2]][NHANES_mort_list[[2]]$eligstat == 1,]
)

load(wave_g_file)
load(wave_h_file)

act <- t(rbind(wave_data_G[[1]],wave_data_H[[1]])[,-1])
act0 <- act == 0
act <- log(act)
lod_act <- min(act[act != -Inf],na.rm = TRUE) - 1e-5
act[act0] <- lod_act

light <- t(rbind(wave_data_G[[2]],wave_data_H[[2]])[,-1])
light0 <- light == 0
light <- log(light)
lod_light <- min(light[light != -Inf],na.rm = TRUE) - 1e-5
light[light0] <- lod_light

id <- rbind(wave_data_G[[3]],wave_data_H[[3]])

seqn_com_id <- id$SEQN %in% lmf_data$seqn
seqn_com_lmf <- lmf_data$seqn %in% id$SEQN

id <- id[seqn_com_id,]
act <- act[,seqn_com_id]
light <- light[,seqn_com_id]
lmf_data <- lmf_data[seqn_com_lmf,]

stopifnot(sum(id$SEQN - lmf_data$seqn) == 0)

id$age_disc <- ifelse(
  id$age <= 30,
  1,
  ifelse(id$age <= 50,2,ifelse(id$age <= 65,3,4))
)

vcovar_mat <- sapply(as.numeric(id$PAXDAYWM),first_day_to_week_ind)
to_keep_inds <- !is.na(id$age_disc)

id <- id[to_keep_inds,]
act <- act[,to_keep_inds]
light <- light[,to_keep_inds]
vcovar_mat <- vcovar_mat[,to_keep_inds]

load(model_file)

est_params <- to_save$est_params
init <- est_params$init
params_tran_array <- est_params$params_tran_array
emit_act <- est_params$emit_act
re_prob <- est_params$re_prob
decoded_mat <- est_params$decoded_mat
mix_assignment <- apply(re_prob,1,which.max)

stopifnot(
  plot_ind >= 1,
  plot_ind <= ncol(act),
  nrow(decoded_mat) == nrow(act),
  ncol(decoded_mat) == ncol(act),
  window_start_index >= 1,
  window_start_index + window_length - 1 <= nrow(act)
)

latent_class <- mix_assignment[plot_ind]
full_decode <- decoded_mat[,plot_ind,latent_class]
activity_only_decode <- viterbi_activity_only(
  ind = plot_ind,
  latent_class = latent_class,
  act = act,
  vcovar_mat = vcovar_mat,
  init = init,
  params_tran_array = params_tran_array,
  emit_act = emit_act,
  lod_act = lod_act
)

window_index <- window_start_index + seq_len(window_length) - 1
time_start <- as.POSIXct("2024-01-01 16:00:00",tz = "UTC")
time_vec <- time_start + (seq_len(window_length) - 1) * 15 * 60

plot_data <- function(measure_name,prediction_name,observed,decode,other_decode){
  data.frame(
    time = time_vec,
    observed = observed[window_index],
    decode = decode[window_index],
    differs = decode[window_index] != other_decode[window_index],
    measure = measure_name,
    prediction = prediction_name
  )
}

plot_rows <- rbind(
  plot_data(
    "Activity",
    "Activity and Light",
    act[,plot_ind],
    full_decode,
    activity_only_decode
  ),
  plot_data(
    "Activity",
    "Activity Only",
    act[,plot_ind],
    activity_only_decode,
    full_decode
  ),
  plot_data(
    "Light",
    "Activity and Light",
    light[,plot_ind],
    full_decode,
    activity_only_decode
  ),
  plot_data(
    "Light",
    "Activity Only",
    light[,plot_ind],
    activity_only_decode,
    full_decode
  )
)

plot_rows$state <- ifelse(plot_rows$decode == 0,"Wake","Sleep")
plot_rows$state_label <- ifelse(
  plot_rows$differs,
  paste0(plot_rows$state," (Differ)"),
  plot_rows$state
)
plot_rows$state_label <- factor(
  plot_rows$state_label,
  levels = c("Wake","Sleep","Wake (Differ)","Sleep (Differ)")
)

make_panel <- function(measure_name,prediction_name,tag_title){
  panel_df <- plot_rows[
    plot_rows$measure == measure_name &
      plot_rows$prediction == prediction_name,
  ]

  if (measure_name == "Activity"){
    sleep_y <- -6
    wake_y <- 4
    right_axis_name <- "Activity (Log MIMS)"
    right_axis_breaks <- seq(-6,4,2)
  } else {
    sleep_y <- -2
    wake_y <- 8
    right_axis_name <- "Light (Log Lux)"
    right_axis_breaks <- seq(-2,8,2)
  }

  panel_df$state_y <- ifelse(panel_df$decode == 0,wake_y,sleep_y)

  ggplot(panel_df,aes(x = time)) +
    geom_line(aes(y = state_y),color = "#1f7a8c",linewidth = 0.8) +
    geom_point(
      aes(y = observed,color = state_label,shape = state_label),
      size = 1.6
    ) +
    scale_color_manual(
      values = c(
        "Wake" = "#8bd63a",
        "Sleep" = "#440154",
        "Wake (Differ)" = "#ff7f50",
        "Sleep (Differ)" = "#c51b7d"
      ),
      name = "Predicted State",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "Wake" = 16,
        "Sleep" = 16,
        "Wake (Differ)" = 8,
        "Sleep (Differ)" = 8
      ),
      name = "Predicted State",
      drop = FALSE
    ) +
    scale_x_datetime(
      breaks = seq(time_start,time_start + 24 * 60 * 60,by = 4 * 60 * 60),
      limits = c(time_start,time_start + 24 * 60 * 60),
      date_labels = "%H:%M",
      expand = c(0,0),
      name = "Time"
    ) +
    scale_y_continuous(
      breaks = c(sleep_y,wake_y),
      labels = c("Sleep","Wake"),
      limits = c(sleep_y - 0.5,wake_y + 0.5),
      name = "Predicted State",
      sec.axis = sec_axis(
        ~ .,
        name = right_axis_name,
        breaks = right_axis_breaks
      )
    ) +
    labs(
      title = paste("Predicted State and Observed",measure_name,"by Time"),
      subtitle = prediction_name,
      tag = tag_title
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(size = 13),
      plot.tag = element_text(size = 14,face = "bold"),
      axis.text.x = element_text(angle = 45,hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_line(color = "grey90")
    )
}

vit_compare_plot <-
  (
    make_panel("Activity","Activity and Light","A") +
      make_panel("Activity","Activity Only","B")
  ) /
  (
    make_panel("Light","Activity and Light","C") +
      make_panel("Light","Activity Only","D")
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(vit_compare_plot)

dir.create(dirname(output_prefix),recursive = TRUE,showWarnings = FALSE)
ggsave(
  paste0(output_prefix,".png"),
  vit_compare_plot,
  width = 12,
  height = 8,
  dpi = 300
)
