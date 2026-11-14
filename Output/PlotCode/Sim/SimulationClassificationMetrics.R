#!/usr/bin/env Rscript

# Calculate hard latent-class classification metrics for the five-class
# simulation study. The input contains one true-class-by-fitted-class
# confusion matrix for each fitted simulation replicate.
#
# Class labels are arbitrary. Within each replicate, fitted labels are aligned
# to true labels using the permutation that maximizes total classification
# agreement. Metrics are then calculated as follows:
#
#   accuracy          = total correctly classified / total participants
#   sensitivity_g     = correctly classified in true class g / true class g
#   macro sensitivity = unweighted mean of the five sensitivities
#   PPV_g             = correctly classified in fitted class g / fitted class g
#   macro PPV         = unweighted mean of the five PPVs
#
# Scenario summaries average replicate-level metrics, so every successful
# simulation replicate receives equal weight. These are in-sample, modal-class
# diagnostics. They complement entropy but do not fully describe uncertainty
# in the posterior class probabilities used by the two-stage Cox model.

input_file <- file.path("Output","parse_sim_results.rds")
output_rds <- file.path("Output","simulation_classification_metrics.rds")
output_csv <- file.path("Output","simulation_classification_metrics.csv")

true_mix_num_target <- 5L
fit_mix_num_target <- 5L
follow_up_values <- c(1,3,7)
overlap_values <- c("low","mid","high")
model_values <- c("joint","two_stage")
class_count <- 5L

separation_labels <- c(
  low = "High",
  mid = "Middle",
  high = "Low"
)

model_labels <- c(
  joint = "JLCM",
  two_stage = "Two-stage"
)


stop_with <- function(...){
  stop(paste0(...),call. = FALSE)
}


all_permutations <- function(x){
  if (length(x) == 1L){
    return(matrix(x,nrow = 1L))
  }

  do.call(rbind,lapply(seq_along(x),function(i){
    cbind(x[[i]],all_permutations(x[-i]))
  }))
}


best_class_match <- function(confusion_table){
  tab <- as.matrix(confusion_table)
  if (nrow(tab) != ncol(tab)){
    stop_with("Confusion matrix must be square")
  }

  class_num <- nrow(tab)
  permutations <- all_permutations(seq_len(class_num))
  scores <- apply(permutations,1L,function(permutation){
    sum(tab[cbind(seq_len(class_num),permutation)])
  })
  best_index <- which.max(scores)

  list(
    fitted_class_for_true_class =
      as.integer(permutations[best_index,,drop = TRUE]),
    matched_n = as.numeric(scores[[best_index]]),
    total_n = as.numeric(sum(tab))
  )
}


validate_confusion_matrix <- function(confusion_table,expected_total,context){
  if (is.null(confusion_table) || length(confusion_table) == 0L){
    stop_with("Missing confusion matrix for ",context)
  }

  tab <- as.matrix(confusion_table)
  storage.mode(tab) <- "numeric"

  if (!identical(dim(tab),c(class_count,class_count))){
    stop_with(
      "Expected a ",class_count," x ",class_count,
      " confusion matrix for ",context,"; found ",
      paste(dim(tab),collapse = " x ")
    )
  }
  if (any(!is.finite(tab))){
    stop_with("Confusion matrix contains non-finite counts for ",context)
  }
  if (any(tab < 0)){
    stop_with("Confusion matrix contains negative counts for ",context)
  }
  if (any(abs(tab - round(tab)) > sqrt(.Machine$double.eps))){
    stop_with("Confusion matrix contains non-integer counts for ",context)
  }

  axis_names <- names(dimnames(tab))
  if (!is.null(axis_names) &&
      !identical(axis_names,c("true_class","fitted_class"))){
    stop_with(
      "Confusion matrix axes are not true_class by fitted_class for ",context
    )
  }

  total_count <- sum(tab)
  if (!is.finite(expected_total) || length(expected_total) != 1L){
    stop_with("Invalid num_people value for ",context)
  }
  if (abs(total_count - expected_total) > sqrt(.Machine$double.eps)){
    stop_with(
      "Confusion matrix total does not match num_people for ",context,
      ": ",total_count," versus ",expected_total
    )
  }
  if (any(rowSums(tab) <= 0)){
    stop_with("At least one true class is empty for ",context)
  }
  if (any(colSums(tab) <= 0)){
    stop_with("At least one fitted class is empty for ",context)
  }

  tab
}


classification_metrics <- function(confusion_table){
  class_match <- best_class_match(confusion_table)
  aligned_table <- confusion_table[
    ,class_match$fitted_class_for_true_class,
    drop = FALSE
  ]

  true_positive <- diag(aligned_table)
  true_count <- rowSums(aligned_table)
  predicted_count <- colSums(aligned_table)
  sensitivity <- true_positive / true_count
  positive_predictive_value <- true_positive / predicted_count
  accuracy <- sum(true_positive) / sum(aligned_table)
  matched_accuracy <- class_match$matched_n / class_match$total_n

  if (!isTRUE(all.equal(accuracy,matched_accuracy,tolerance = 1e-12))){
    stop_with("Accuracy does not equal the optimal label-matching percentage")
  }

  list(
    aligned_table = aligned_table,
    fitted_class_for_true_class =
      class_match$fitted_class_for_true_class,
    matched_n = class_match$matched_n,
    total_n = class_match$total_n,
    accuracy = accuracy,
    sensitivity = sensitivity,
    macro_sensitivity = mean(sensitivity),
    positive_predictive_value = positive_predictive_value,
    macro_positive_predictive_value = mean(positive_predictive_value),
    true_positive = true_positive,
    true_count = true_count,
    predicted_count = predicted_count
  )
}


check_label_invariance <- function(){
  test_table <- matrix(
    c(
      80, 4, 3, 2, 1,
       3,75, 4, 2, 1,
       2, 4,70, 3, 1,
       1, 2, 3,68, 2,
       1, 1, 2, 3,65
    ),
    nrow = class_count,
    byrow = TRUE
  )
  column_permutation <- c(3,5,1,4,2)

  original <- classification_metrics(test_table)
  permuted <- classification_metrics(test_table[,column_permutation])
  metric_names <- c(
    "accuracy","macro_sensitivity","macro_positive_predictive_value"
  )
  original_values <- unlist(original[metric_names],use.names = FALSE)
  permuted_values <- unlist(permuted[metric_names],use.names = FALSE)

  if (!isTRUE(all.equal(original_values,permuted_values,tolerance = 1e-12))){
    stop_with("Label-invariance validation failed")
  }

  invisible(TRUE)
}


metric_summary <- function(x,prefix){
  x <- as.numeric(x)
  if (length(x) == 0L || any(!is.finite(x))){
    stop_with("Cannot summarize non-finite or empty ",prefix," values")
  }

  metric_sd <- if (length(x) > 1L) stats::sd(x) else NA_real_
  result <- c(
    mean(x),
    metric_sd,
    metric_sd / sqrt(length(x))
  )
  names(result) <- paste0(c("mean_","sd_","mcse_"),prefix)
  result
}


group_rows <- function(data,group_row,group_columns){
  keep <- rep(TRUE,nrow(data))
  for (column in group_columns){
    keep <- keep & data[[column]] == group_row[[column]]
  }
  data[keep,,drop = FALSE]
}


make_scenario_summary <- function(replication_metrics){
  group_columns <- c(
    "simulation_days","emission_overlap","separation_label",
    "model_type","model_label"
  )
  groups <- unique(replication_metrics[,group_columns,drop = FALSE])
  groups <- groups[
    order(
      match(groups$emission_overlap,overlap_values),
      groups$simulation_days,
      match(groups$model_type,model_values)
    ),
    ,drop = FALSE
  ]

  rows <- vector("list",nrow(groups))
  for (group_index in seq_len(nrow(groups))){
    current_group <- groups[group_index,,drop = FALSE]
    current <- group_rows(replication_metrics,current_group,group_columns)
    summary_values <- c(
      metric_summary(current$accuracy,"accuracy"),
      metric_summary(current$macro_sensitivity,"macro_sensitivity"),
      metric_summary(
        current$macro_positive_predictive_value,
        "macro_positive_predictive_value"
      )
    )

    rows[[group_index]] <- cbind(
      current_group,
      data.frame(
        n_runs = nrow(current),
        as.list(summary_values),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }

  result <- do.call(rbind,rows)
  row.names(result) <- NULL
  result
}


make_class_scenario_summary <- function(replication_class_metrics){
  group_columns <- c(
    "simulation_days","emission_overlap","separation_label",
    "model_type","model_label","latent_class"
  )
  groups <- unique(replication_class_metrics[,group_columns,drop = FALSE])
  groups <- groups[
    order(
      match(groups$emission_overlap,overlap_values),
      groups$simulation_days,
      match(groups$model_type,model_values),
      groups$latent_class
    ),
    ,drop = FALSE
  ]

  rows <- vector("list",nrow(groups))
  for (group_index in seq_len(nrow(groups))){
    current_group <- groups[group_index,,drop = FALSE]
    current <- group_rows(
      replication_class_metrics,current_group,group_columns
    )
    summary_values <- c(
      metric_summary(current$sensitivity,"sensitivity"),
      metric_summary(
        current$positive_predictive_value,
        "positive_predictive_value"
      )
    )

    rows[[group_index]] <- cbind(
      current_group,
      data.frame(
        n_runs = nrow(current),
        as.list(summary_values),
        mean_true_count = mean(current$true_count),
        mean_predicted_count = mean(current$predicted_count),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    )
  }

  result <- do.call(rbind,rows)
  row.names(result) <- NULL
  result
}


make_wide_summary <- function(scenario_summary){
  id_columns <- c(
    "simulation_days","emission_overlap","separation_label"
  )
  measure_columns <- c(
    "n_runs",
    "mean_accuracy","sd_accuracy","mcse_accuracy",
    "mean_macro_sensitivity","sd_macro_sensitivity",
    "mcse_macro_sensitivity",
    "mean_macro_positive_predictive_value",
    "sd_macro_positive_predictive_value",
    "mcse_macro_positive_predictive_value"
  )

  joint <- scenario_summary[
    scenario_summary$model_type == "joint",
    c(id_columns,measure_columns),
    drop = FALSE
  ]
  names(joint)[match(measure_columns,names(joint))] <-
    paste0(measure_columns,"_jlcm")

  two_stage <- scenario_summary[
    scenario_summary$model_type == "two_stage",
    c(id_columns,measure_columns),
    drop = FALSE
  ]
  names(two_stage)[match(measure_columns,names(two_stage))] <-
    paste0(measure_columns,"_two_stage")

  result <- merge(joint,two_stage,by = id_columns,all = TRUE,sort = FALSE)
  result <- result[
    order(
      match(result$emission_overlap,overlap_values),
      result$simulation_days
    ),
    ,drop = FALSE
  ]
  row.names(result) <- NULL
  result
}


validate_scenarios <- function(simulation_data){
  scenario_columns <- c(
    "simulation_days","emission_overlap","model_type"
  )
  expected <- expand.grid(
    simulation_days = follow_up_values,
    emission_overlap = overlap_values,
    model_type = model_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  actual <- unique(simulation_data[,scenario_columns,drop = FALSE])

  make_key <- function(x){
    paste(x$simulation_days,x$emission_overlap,x$model_type,sep = "|")
  }
  expected_keys <- make_key(expected)
  actual_keys <- make_key(actual)

  missing_keys <- setdiff(expected_keys,actual_keys)
  unexpected_keys <- setdiff(actual_keys,expected_keys)
  if (length(missing_keys) > 0L){
    stop_with("Missing model-scenarios: ",paste(missing_keys,collapse = ", "))
  }
  if (length(unexpected_keys) > 0L){
    stop_with(
      "Unexpected model-scenarios: ",paste(unexpected_keys,collapse = ", ")
    )
  }
  if (nrow(actual) != 18L){
    stop_with("Expected 18 model-scenarios; found ",nrow(actual))
  }

  scenario_count <- nrow(unique(
    simulation_data[,c("simulation_days","emission_overlap"),drop = FALSE]
  ))
  if (scenario_count != 9L){
    stop_with("Expected nine simulation scenarios; found ",scenario_count)
  }

  invisible(TRUE)
}


if (!file.exists(input_file)){
  stop_with("Input file does not exist: ",input_file)
}

check_label_invariance()

parsed_results <- readRDS(input_file)
required_columns <- c(
  "file","sim_num","true_mix_num","fit_mix_num","model_type",
  "simulation_days","num_people","emission_overlap",
  "diagnostic_confusion_table"
)
missing_columns <- setdiff(required_columns,names(parsed_results))
if (length(missing_columns) > 0L){
  stop_with(
    "Parsed simulation results are missing required columns: ",
    paste(missing_columns,collapse = ", ")
  )
}

keep <- parsed_results$true_mix_num == true_mix_num_target &
  parsed_results$fit_mix_num == fit_mix_num_target &
  parsed_results$model_type %in% model_values &
  parsed_results$simulation_days %in% follow_up_values &
  parsed_results$emission_overlap %in% overlap_values
simulation_data <- parsed_results[keep,required_columns,drop = FALSE]

if (nrow(simulation_data) == 0L){
  stop_with("No five-class simulation runs matched the requested scenarios")
}
if (length(unique(simulation_data$num_people)) != 1L){
  stop_with("Expected one num_people value across the five-class simulations")
}

run_key <- paste(
  simulation_data$simulation_days,
  simulation_data$emission_overlap,
  simulation_data$model_type,
  simulation_data$sim_num,
  sep = "|"
)
if (anyDuplicated(run_key)){
  stop_with("Duplicate simulation run keys were found in the parsed results")
}

validate_scenarios(simulation_data)

replication_rows <- vector("list",nrow(simulation_data))
class_rows <- vector("list",nrow(simulation_data) * class_count)
class_row_index <- 0L

for (run_index in seq_len(nrow(simulation_data))){
  run <- simulation_data[run_index,,drop = FALSE]
  context <- paste0(
    "days=",run$simulation_days,
    ", overlap=",run$emission_overlap,
    ", model=",run$model_type,
    ", seed=",run$sim_num
  )
  confusion_table <- validate_confusion_matrix(
    run$diagnostic_confusion_table[[1]],
    expected_total = run$num_people,
    context = context
  )
  metrics <- classification_metrics(confusion_table)

  replication_rows[[run_index]] <- data.frame(
    file = run$file,
    sim_num = run$sim_num,
    true_mix_num = run$true_mix_num,
    fit_mix_num = run$fit_mix_num,
    model_type = run$model_type,
    model_label = unname(model_labels[[run$model_type]]),
    simulation_days = run$simulation_days,
    num_people = run$num_people,
    emission_overlap = run$emission_overlap,
    separation_label = unname(separation_labels[[run$emission_overlap]]),
    accuracy = metrics$accuracy,
    macro_sensitivity = metrics$macro_sensitivity,
    macro_positive_predictive_value =
      metrics$macro_positive_predictive_value,
    matched_n = metrics$matched_n,
    total_n = metrics$total_n,
    fitted_class_permutation = paste(
      metrics$fitted_class_for_true_class,
      collapse = ","
    ),
    stringsAsFactors = FALSE
  )

  for (latent_class in seq_len(class_count)){
    class_row_index <- class_row_index + 1L
    true_positive <- metrics$true_positive[[latent_class]]
    true_count <- metrics$true_count[[latent_class]]
    predicted_count <- metrics$predicted_count[[latent_class]]

    class_rows[[class_row_index]] <- data.frame(
      file = run$file,
      sim_num = run$sim_num,
      model_type = run$model_type,
      model_label = unname(model_labels[[run$model_type]]),
      simulation_days = run$simulation_days,
      emission_overlap = run$emission_overlap,
      separation_label = unname(separation_labels[[run$emission_overlap]]),
      latent_class = latent_class,
      matched_fitted_class =
        metrics$fitted_class_for_true_class[[latent_class]],
      true_positive = true_positive,
      false_negative = true_count - true_positive,
      false_positive = predicted_count - true_positive,
      true_count = true_count,
      predicted_count = predicted_count,
      total_count = metrics$total_n,
      sensitivity = metrics$sensitivity[[latent_class]],
      positive_predictive_value =
        metrics$positive_predictive_value[[latent_class]],
      stringsAsFactors = FALSE
    )
  }
}

replication_metrics <- do.call(rbind,replication_rows)
row.names(replication_metrics) <- NULL
replication_class_metrics <- do.call(rbind,class_rows)
row.names(replication_class_metrics) <- NULL

metric_columns <- c(
  "accuracy","macro_sensitivity","macro_positive_predictive_value"
)
class_metric_columns <- c("sensitivity","positive_predictive_value")
if (any(!is.finite(as.matrix(replication_metrics[,metric_columns])))){
  stop_with("Replication-level metrics contain non-finite values")
}
if (any(!is.finite(as.matrix(
  replication_class_metrics[,class_metric_columns]
)))){
  stop_with("Class-level metrics contain non-finite values")
}
if (any(as.matrix(replication_metrics[,metric_columns]) < 0 | 
        as.matrix(replication_metrics[,metric_columns]) > 1)){
  stop_with("Replication-level metrics fall outside [0, 1]")
}
if (any(as.matrix(replication_class_metrics[,class_metric_columns]) < 0 |
        as.matrix(replication_class_metrics[,class_metric_columns]) > 1)){
  stop_with("Class-level metrics fall outside [0, 1]")
}

scenario_summary <- make_scenario_summary(replication_metrics)
scenario_class_summary <- make_class_scenario_summary(
  replication_class_metrics
)
scenario_summary_wide <- make_wide_summary(scenario_summary)

if (nrow(scenario_summary) != 18L){
  stop_with("Expected 18 rows in the long scenario summary")
}
if (nrow(scenario_summary_wide) != 9L){
  stop_with("Expected nine rows in the wide scenario summary")
}

if (any(scenario_summary$n_runs < 245L | scenario_summary$n_runs > 250L)){
  warning(
    "At least one scenario has a run count outside the expected 245--250 range",
    call. = FALSE
  )
}

output <- list(
  metadata = list(
    input_file = normalizePath(input_file,mustWork = TRUE),
    metric_scale = "proportion",
    true_mix_num = true_mix_num_target,
    fit_mix_num = fit_mix_num_target,
    class_count = class_count,
    follow_up_values = follow_up_values,
    emission_overlap_values = overlap_values,
    model_values = model_values,
    label_alignment = "per-replicate permutation maximizing accuracy",
    scenario_aggregation = "unweighted mean of replicate-level metrics"
  ),
  replication_metrics = replication_metrics,
  replication_class_metrics = replication_class_metrics,
  scenario_class_summary = scenario_class_summary,
  scenario_summary = scenario_summary,
  scenario_summary_wide = scenario_summary_wide
)

dir.create(dirname(output_rds),recursive = TRUE,showWarnings = FALSE)
saveRDS(output,output_rds)
utils::write.csv(scenario_summary_wide,output_csv,row.names = FALSE,na = "")

display_summary <- scenario_summary[,c(
  "separation_label","simulation_days","model_label","n_runs",
  "mean_accuracy","mean_macro_sensitivity",
  "mean_macro_positive_predictive_value"
)]
names(display_summary) <- c(
  "Separation","Days","Model","Runs",
  "Accuracy","Macro sensitivity","Macro PPV"
)
for (column in c("Accuracy","Macro sensitivity","Macro PPV")){
  display_summary[[column]] <- sprintf("%.1f%%",100 * display_summary[[column]])
}

cat("\nFive-class simulation classification metrics\n")
print(display_summary,row.names = FALSE,right = FALSE)
cat("\nSaved detailed results to: ",output_rds,"\n",sep = "")
cat("Saved nine-scenario summary to: ",output_csv,"\n",sep = "")
