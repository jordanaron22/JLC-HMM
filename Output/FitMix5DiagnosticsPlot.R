library(ggplot2)
input_file <- file.path("Output","parse_sim_results.rds")


title_case <- function(x){
  paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
}


dat <- readRDS(input_file)
dat <- dat[dat$fit_mix_num == 5,]

diagnostics <- data.frame(
  diagnostic = c("test_cindex","test_ibs"),
  value_col = c("diagnostic_test_cindex","diagnostic_test_ibs"),
  diagnostic_label = c("Test C-index","Test IBS"),
  stringsAsFactors = FALSE
)

plot_rows <- list()
for (i in seq_len(nrow(diagnostics))){
  value_col <- diagnostics$value_col[[i]]

  if (!value_col %in% names(dat)){
    next
  }

  plot_rows[[length(plot_rows) + 1]] <- data.frame(
    true_mix_num = dat$true_mix_num,
    fit_mix_num = dat$fit_mix_num,
    model_type = dat$model_type,
    days = dat$simulation_days,
    num_people = dat$num_people,
    emission_overlap = dat$emission_overlap,
    sim_num = dat$sim_num,
    diagnostic = diagnostics$diagnostic[[i]],
    diagnostic_label = diagnostics$diagnostic_label[[i]],
    value = dat[[value_col]],
    stringsAsFactors = FALSE
  )
}

diagnostic_data <- do.call(rbind,plot_rows)
diagnostic_data <- diagnostic_data[is.finite(diagnostic_data$value),]

scenario_cols <- c("true_mix_num","fit_mix_num","model_type","days",
                   "num_people","emission_overlap","diagnostic",
                   "diagnostic_label")

plot_data <- aggregate(
  value ~ true_mix_num + fit_mix_num + model_type + days +
    num_people + emission_overlap + diagnostic + diagnostic_label,
  data = diagnostic_data,
  FUN = mean
)

n_seeds <- aggregate(
  sim_num ~ true_mix_num + fit_mix_num + model_type + days +
    num_people + emission_overlap + diagnostic + diagnostic_label,
  data = diagnostic_data,
  FUN = function(x) length(unique(x))
)
names(n_seeds)[names(n_seeds) == "sim_num"] <- "n_seeds"

plot_data <- merge(plot_data,n_seeds,by = scenario_cols)

day_values <- sort(unique(plot_data$days))
plot_data$days_label <- factor(
  paste0(plot_data$days," Day",ifelse(plot_data$days == 1,"","s")),
  levels = paste0(day_values," Day",ifelse(day_values == 1,"","s"))
)

plot_data$model_type_label <- factor(
  ifelse(plot_data$model_type == "joint","JLCM","Two-Stage"),
  levels = c("JLCM","Two-Stage")
)

overlap_values <- sort(unique(plot_data$emission_overlap))
plot_data$overlap_label <- factor(
  paste0("Overlap: ",title_case(plot_data$emission_overlap)),
  levels = paste0("Overlap: ",title_case(overlap_values))
)

people_values <- sort(unique(plot_data$num_people))
plot_data$people_label <- factor(
  paste0("n = ",format(plot_data$num_people,big.mark = ",",
                       scientific = FALSE,trim = TRUE)),
  levels = paste0("n = ",format(people_values,big.mark = ",",
                                scientific = FALSE,trim = TRUE))
)

diagnostic_levels <- diagnostics$diagnostic_label[
  diagnostics$diagnostic %in% unique(plot_data$diagnostic)
]
plot_data$diagnostic_label <- factor(
  plot_data$diagnostic_label,
  levels = diagnostic_levels
)

plot_data <- plot_data[order(plot_data$diagnostic,
                             plot_data$emission_overlap,
                             plot_data$days,
                             plot_data$model_type),]
row.names(plot_data) <- NULL

ggplot(
  plot_data,
  aes(x = days_label,y = value,fill = model_type_label)
) +
  geom_col(position = position_dodge()) +
  facet_grid(rows = vars(diagnostic_label),
                      cols = vars(overlap_label),
                      scales = "free_y") +
  labs(x = "Days",y = "Mean diagnostic",fill = "Model")


ggsave("SurvivalDiagnosticsSim.png", path = "Output/Plots", width = 10, height = 6, dpi = 300)
