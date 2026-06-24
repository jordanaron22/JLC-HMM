input_file <- file.path("Output","parse_sim_results.rds")


group_key <- function(data,cols){
  apply(data[,cols,drop = FALSE],1,function(row){
    row[is.na(row)] <- "<NA>"
    paste(row,collapse = "\r")
  })
}



dat <- readRDS(input_file)


criteria <- data.frame(
  criterion = c("bic","aic",
                "test_ibs","test_cindex"),
  value_col = c("bic","aic","diagnostic_test_ibs",
                "diagnostic_test_cindex"),
  direction = c("min","min","min","max"),
  stringsAsFactors = FALSE
)

seed_cols <- c("true_mix_num","model_type","simulation_days","num_people",
               "emission_overlap","sim_num")

class_levels <- sort(unique(dat$fit_mix_num[is.finite(dat$fit_mix_num)]))

selected_rows <- list()
for (criterion_ind in seq_len(nrow(criteria))){
  criterion <- criteria$criterion[[criterion_ind]]
  value_col <- criteria$value_col[[criterion_ind]]
  direction <- criteria$direction[[criterion_ind]]

  seed_key <- group_key(dat,seed_cols)
  groups <- split(seq_len(nrow(dat)),seed_key)

  for (indices in groups){
    group_data <- dat[indices,,drop = FALSE]
    values <- group_data[[value_col]]
    valid <- is.finite(values)
    if (!any(valid)){
      next
    }

    group_data <- group_data[valid,,drop = FALSE]
    values <- values[valid]
    selected_ind <- if (direction == "min"){
      order(values,group_data$fit_mix_num)[[1]]
    } else {
      order(-values,group_data$fit_mix_num)[[1]]
    }
    selected <- group_data[selected_ind,,drop = FALSE]

    selected_rows[[length(selected_rows) + 1]] <- data.frame(
      true_mix_num = selected$true_mix_num,
      model_type = selected$model_type,
      days = selected$simulation_days,
      num_people = selected$num_people,
      emission_overlap = selected$emission_overlap,
      sim_num = selected$sim_num,
      criterion = criterion,
      selected_fit_mix_num = selected$fit_mix_num,
      stringsAsFactors = FALSE
    )
  }
}

class_selection <- do.call(rbind,selected_rows)
if (is.null(class_selection) || nrow(class_selection) == 0){
  stop("No class selections could be calculated.")
}

scenario_cols <- c("true_mix_num","model_type","days","num_people",
                   "emission_overlap","criterion")

scenario_key <- group_key(class_selection,scenario_cols)
scenario_groups <- split(seq_len(nrow(class_selection)),scenario_key)

plot_rows <- list()
for (indices in scenario_groups){
  scenario_data <- class_selection[indices,,drop = FALSE]
  n_seeds <- length(unique(scenario_data$sim_num))

  for (class_level in class_levels){
    n <- sum(scenario_data$selected_fit_mix_num == class_level)
    row <- scenario_data[1,scenario_cols,drop = FALSE]
    row$selected_fit_mix_num <- class_level
    row$n_seeds <- n_seeds
    row$n <- n
    row$percent <- 100 * n / n_seeds
    plot_rows[[length(plot_rows) + 1]] <- row
  }
}

plot_data <- do.call(rbind,plot_rows)


ggplot(
  plot_data |> dplyr::filter(model_type == "joint"),
  aes(x = selected_fit_mix_num, y = n,
               fill = factor(model_type))
) +
  geom_col(
    position = position_dodge2()
  ) +
  facet_grid(rows = vars(criterion),
                      cols = vars(days, emission_overlap)) +
  labs(x = "Selected latent classes",
                y = "Selections (n)",
                fill = "Model")

ggsave("class_selection_percentages.png",  path = "Output/Plots",width = 10, height = 6, dpi = 300)
