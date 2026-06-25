input_file <- file.path("Output","parse_sim_results.rds")


title_case <- function(x){
  paste0(toupper(substr(x,1,1)),substr(x,2,nchar(x)))
}


all_permutations <- function(x){
  if (length(x) == 1){
    return(matrix(x,nrow = 1))
  }

  do.call(rbind,lapply(seq_along(x),function(i){
    cbind(x[[i]],all_permutations(x[-i]))
  }))
}


best_class_match <- function(confusion_table){
  tab <- as.matrix(confusion_table)
  class_num <- nrow(tab)
  permutations <- all_permutations(seq_len(class_num))
  scores <- apply(permutations,1,function(perm){
    sum(tab[cbind(seq_len(class_num),perm)])
  })
  best_ind <- which.max(scores)
  matched_n <- scores[[best_ind]]
  total_n <- sum(tab)

  list(
    fitted_class_for_true_class = permutations[best_ind,],
    matched_n = matched_n,
    total_n = total_n,
    matched_percent = 100 * matched_n / total_n
  )
}


get_beta_data <- function(file){
  load_env <- new.env(parent = emptyenv())
  load(file,envir = load_env)
  saved <- load_env[["to_save"]]

  true_beta <- as.numeric(saved[["true_params"]][["beta_vec"]])
  fit_beta <- as.numeric(saved[["est_params"]][["beta_vec"]])
  beta_se <- as.numeric(saved[["est_params"]][["beta_se"]])
  class_match <- best_class_match(saved[["diagnostics"]][["confusion_table"]])

  beta_se_for_beta_vec <- rep(NA_real_,length(true_beta))
  beta_se_for_beta_vec[-1] <- beta_se[2:length(true_beta)]

  true_beta <- true_beta - true_beta[[1]]
  fit_beta <- fit_beta[class_match$fitted_class_for_true_class]
  beta_se_for_beta_vec <- beta_se_for_beta_vec[
    class_match$fitted_class_for_true_class
  ]
  reference_beta_se <- beta_se_for_beta_vec[[1]]
  fit_beta <- fit_beta - fit_beta[[1]]

  data.frame(
    beta_num = seq_along(true_beta),
    true_beta = true_beta,
    fit_beta = fit_beta,
    beta_se = beta_se_for_beta_vec,
    reference_beta_se = reference_beta_se,
    matched_fitted_class = class_match$fitted_class_for_true_class,
    matched_n = class_match$matched_n,
    total_n = class_match$total_n,
    matched_percent = class_match$matched_percent,
    stringsAsFactors = FALSE
  )
}


dat <- readRDS(input_file)
dat <- dat[dat$fit_mix_num == 5 & dat$true_mix_num == 5,]

run_cols <- c("file","true_mix_num","fit_mix_num","model_type",
              "simulation_days","num_people","emission_overlap","sim_num")

coef_rows <- list()
for (i in seq_len(nrow(dat))){
  beta_data <- get_beta_data(dat$file[[i]])
  beta_data <- beta_data[beta_data$beta_num > 1,]

  run_data <- dat[i,run_cols,drop = FALSE]
  run_data <- run_data[rep(1,nrow(beta_data)),]

  coef_rows[[length(coef_rows) + 1]] <- cbind(run_data,beta_data)
}

coef_data <- do.call(rbind,coef_rows)
row.names(coef_data) <- NULL
coef_data$error <- coef_data$fit_beta - coef_data$true_beta
coef_data$squared_error <- coef_data$error^2

scenario_cols <- c("true_mix_num","fit_mix_num","model_type",
                   "simulation_days","num_people","emission_overlap")

match_data <- unique(coef_data[,c(scenario_cols,"sim_num",
                                  "matched_percent")])

large_beta_se <- coef_data[is.finite(coef_data$beta_se) &
                             coef_data$beta_se > 100 |
                             is.finite(coef_data$reference_beta_se) &
                             coef_data$reference_beta_se > 100,]

coef_data <- coef_data[!(is.finite(coef_data$beta_se) &
                          coef_data$beta_se > 100 |
                          is.finite(coef_data$reference_beta_se) &
                          coef_data$reference_beta_se > 100),]

plot_data <- aggregate(
  squared_error ~ true_mix_num + fit_mix_num + model_type +
    simulation_days + num_people + emission_overlap,
  data = coef_data,
  FUN = function(x) sqrt(mean(x,na.rm = TRUE))
)
names(plot_data)[names(plot_data) == "squared_error"] <- "rmse"

n_seeds <- aggregate(
  sim_num ~ true_mix_num + fit_mix_num + model_type +
    simulation_days + num_people + emission_overlap,
  data = coef_data,
  FUN = function(x) length(unique(x))
)
names(n_seeds)[names(n_seeds) == "sim_num"] <- "n_seeds"

n_coef <- aggregate(
  squared_error ~ true_mix_num + fit_mix_num + model_type +
    simulation_days + num_people + emission_overlap,
  data = coef_data,
  FUN = length
)
names(n_coef)[names(n_coef) == "squared_error"] <- "n"

mean_matched_percent <- aggregate(
  matched_percent ~ true_mix_num + fit_mix_num + model_type +
    simulation_days + num_people + emission_overlap,
  data = match_data,
  FUN = mean
)
names(mean_matched_percent)[names(mean_matched_percent) ==
                              "matched_percent"] <- "mean_matched_percent"

plot_data <- merge(plot_data,n_seeds,by = scenario_cols)
plot_data <- merge(plot_data,n_coef,by = scenario_cols)
plot_data <- merge(plot_data,mean_matched_percent,by = scenario_cols)

plot_data$days <- plot_data$simulation_days

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

plot_data <- plot_data[,c("true_mix_num","fit_mix_num","model_type","days",
                          "num_people","emission_overlap","rmse","n_seeds",
                          "n","mean_matched_percent",
                          "days_label","model_type_label",
                          "overlap_label","people_label")]
row.names(plot_data) <- NULL


ggplot(plot_data, aes(x = factor(days), y = rmse, fill = factor(model_type))) + 
  geom_col(position = position_dodge()) + 
  facet_grid(rows = vars(emission_overlap),scales = "free_y") 

ggsave("survival_rmse.png",  path = "Output/Plots",width = 10, height = 6, dpi = 300)
