#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1){
  args[[1]]
} else {
  file.path("Output","Routputs","RoutputsNHANESCVscenarios")
}

output_file <- if (length(args) >= 2){
  args[[2]]
} else {
  file.path("Output","parse_nhanes_cv_scenarios.rds")
}


support_plot_file <- if (length(args) >= 3){
  args[[3]]
} else {
  file.path("SubmissionDraft2025","Support",
            "SensSpecbyDataForPaperScenarios.png")
}

light_support_plot_file <- if (length(args) >= 4){
  args[[4]]
} else {
  file.path("SubmissionDraft2025","Support","LightSensitivityScenarios.png")
}

light_plot_file <- file.path("Output","Plots","LightSensitivity.png")

SAVED_SECTION_SLOT <- c(true_params = 1, est_params = 2, bic = 3,
                        leave_out = 4, simulated_hmm = 5,
                        diagnostics = 6, settings = 7)

LEAVE_OUT_SLOT <- c(leave_out_inds = 1, conf_mat_list = 2,
                    cindex_new_list = 3, ibs_new_list = 4,
                    senspec_list = 5, ibs2_new_list = 6,
                    senspec_mix_list = 7)

LEAVE_OUT_LABELS <- c(
  "1" = "Sleep-Wake Cycle",
  "2" = "No Light",
  "3" = "No Activity",
  "4" = "Null Transition",
  "5" = "Activity Only",
  "6" = "Standard"
)

PLOT_PALETTE <- c(
  "Sleep-Wake Cycle" = "#440154",
  "No Light" = "#414487",
  "No Activity" = "#2A788E",
  "Null Transition" = "#22A884",
  "Activity Only" = "#7AD151",
  "Standard" = "#FDE725"
)

LIGHT_PAPER_PALETTE <- c(
  "No Light" = "#440154",
  "Standard" = "#95D840"
)

STATE_STRATA_LABELS <- c(
  "1" = "Activity Below Individual Median",
  "2" = "Activity Above Individual Median",
  "3" = "Total"
)

STATE_STRATA_PLOT_ORDER <- c(
  "Activity Above Individual Median",
  "Activity Below Individual Median",
  "Total"
)

STATE_LABELS <- c("Sleep","Wake")

regex_value <- function(pattern,text,default = NA_character_){
  match <- regexec(pattern,text,perl = TRUE)
  found <- regmatches(text,match)[[1]]
  if (length(found) < 2){
    return(default)
  }
  found[[2]]
}

scalar_or_na <- function(x){
  if (is.null(x) || length(x) == 0){
    return(NA)
  }
  x[[1]]
}

numeric_or_na <- function(x){
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0){
    return(NA_real_)
  }
  out[[1]]
}

integer_or_na <- function(x){
  out <- suppressWarnings(as.integer(x))
  if (length(out) == 0){
    return(NA_integer_)
  }
  out[[1]]
}

mean_or_na <- function(x){
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0){
    return(NA_real_)
  }
  mean(x)
}

rbind_fill <- function(data_list){
  data_list <- Filter(function(x) !is.null(x) && nrow(x) > 0,data_list)
  if (length(data_list) == 0){
    return(data.frame())
  }

  all_names <- unique(unlist(lapply(data_list,names),use.names = FALSE))
  aligned <- lapply(data_list,function(current){
    missing_names <- setdiff(all_names,names(current))
    for (missing_name in missing_names){
      current[[missing_name]] <- NA
    }
    current[,all_names,drop = FALSE]
  })

  out <- do.call(rbind,aligned)
  row.names(out) <- NULL
  out
}

load_to_save <- function(file){
  if (grepl("[.]rds$",file,ignore.case = TRUE)){
    return(readRDS(file))
  }

  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file,envir = load_env)
  if ("to_save" %in% loaded_names && exists("to_save",envir = load_env)){
    return(get("to_save",envir = load_env))
  }
  if (length(loaded_names) == 1){
    return(get(loaded_names[[1]],envir = load_env))
  }
  stop("File does not contain object named to_save")
}

get_section <- function(to_save,name,slot,required = FALSE){
  if (!is.null(to_save) && !is.null(names(to_save)) && name %in% names(to_save)){
    return(to_save[[name]])
  }
  if (!is.null(to_save) && length(to_save) >= slot){
    return(to_save[[slot]])
  }
  if (required){
    stop("Missing saved section: ",name)
  }
  NULL
}

get_setting <- function(settings,name,default = NA){
  if (!is.null(settings) && !is.null(settings[[name]]) &&
      length(settings[[name]]) > 0){
    return(settings[[name]][[1]])
  }
  default
}

get_leave_out_field <- function(leave_out,name,slot,default = NULL){
  if (!is.null(leave_out) &&
      !is.null(names(leave_out)) &&
      name %in% names(leave_out)){
    return(leave_out[[name]])
  }
  if (!is.null(leave_out) && length(leave_out) >= slot){
    return(leave_out[[slot]])
  }
  default
}

model_type_from_file <- function(file_name){
  if (grepl("NoSurv",file_name,fixed = TRUE)){
    return("two_stage")
  }
  "joint"
}

normalize_confusion_matrix <- function(conf_mat,fit_mix_num = NA_integer_){
  conf_mat <- as.matrix(conf_mat)
  storage.mode(conf_mat) <- "numeric"

  if (length(conf_mat) == 0 || nrow(conf_mat) == 0 || ncol(conf_mat) == 0){
    stop("Empty confusion matrix")
  }

  row_labels <- rownames(conf_mat)
  col_labels <- colnames(conf_mat)

  if (is.null(row_labels) || any(!nzchar(row_labels))){
    row_labels <- as.character(seq_len(nrow(conf_mat)))
  }
  if (is.null(col_labels) || any(!nzchar(col_labels))){
    col_labels <- as.character(seq_len(ncol(conf_mat)))
  }

  class_labels <- union(row_labels,col_labels)
  if (is.finite(fit_mix_num) && fit_mix_num > 0){
    class_labels <- as.character(seq_len(fit_mix_num))
  }

  full_mat <- matrix(
    0,
    nrow = length(class_labels),
    ncol = length(class_labels),
    dimnames = list(class_labels,class_labels)
  )

  row_match <- match(row_labels,class_labels)
  col_match <- match(col_labels,class_labels)
  if (any(is.na(row_match)) || any(is.na(col_match))){
    stop("Confusion matrix contains class labels outside fit_mix_num")
  }

  full_mat[row_match,col_match] <- conf_mat
  full_mat
}

confusion_matrix_to_long <- function(conf_mat,file,file_name,fit_mix_num,
                                     model_type,cv_fold_id,leave_out_type,
                                     leave_out_n){
  conf_mat <- normalize_confusion_matrix(conf_mat,fit_mix_num)

  if (!is.na(leave_out_n) && sum(conf_mat) != leave_out_n){
    warning(
      "Confusion matrix count does not match held-out count in ",
      file_name,
      " scenario ",
      leave_out_type,
      ": sum(confusion)=",
      sum(conf_mat),
      ", length(leave_out_inds)=",
      leave_out_n
    )
  }

  cells <- expand.grid(
    mix_assignment_pred = rownames(conf_mat),
    mix_assignment_true_ind = colnames(conf_mat),
    stringsAsFactors = FALSE
  )
  cells$count <- as.numeric(conf_mat)
  cells$file <- normalizePath(file,mustWork = TRUE)
  cells$file_name <- file_name
  cells$fit_mix_num <- fit_mix_num
  cells$model_type <- model_type
  cells$cv_fold_id <- cv_fold_id
  cells$leave_out_type <- leave_out_type
  cells$testing_data <- unname(LEAVE_OUT_LABELS[as.character(leave_out_type)])
  cells$heldout_n <- leave_out_n

  cells
}

state_from_label <- function(x){
  out <- rep(NA_character_,length(x))
  out[grepl("Sleep",x,ignore.case = TRUE)] <- "Sleep"
  out[grepl("Wake",x,ignore.case = TRUE)] <- "Wake"
  out
}

state_confusion_matrix_to_long <- function(conf_mat,file,file_name,fit_mix_num,
                                           model_type,cv_fold_id,
                                           leave_out_type,stratum_index){
  conf_mat <- as.matrix(conf_mat)
  storage.mode(conf_mat) <- "numeric"

  if (nrow(conf_mat) != 2 || ncol(conf_mat) != 2){
    stop(
      "Expected a 2x2 latent-state confusion matrix, got ",
      paste(dim(conf_mat),collapse = "x")
    )
  }

  predicted_state <- state_from_label(rownames(conf_mat))
  true_state <- state_from_label(colnames(conf_mat))

  if (any(is.na(predicted_state)) || any(is.na(true_state))){
    stop("Could not identify Wake/Sleep labels in latent-state matrix")
  }

  full_mat <- matrix(
    0,
    nrow = length(STATE_LABELS),
    ncol = length(STATE_LABELS),
    dimnames = list(STATE_LABELS,STATE_LABELS)
  )
  full_mat[predicted_state,true_state] <- conf_mat

  cells <- expand.grid(
    predicted_state = rownames(full_mat),
    true_state = colnames(full_mat),
    stringsAsFactors = FALSE
  )
  cells$count <- as.numeric(full_mat)
  cells$file <- normalizePath(file,mustWork = TRUE)
  cells$file_name <- file_name
  cells$fit_mix_num <- fit_mix_num
  cells$model_type <- model_type
  cells$cv_fold_id <- cv_fold_id
  cells$leave_out_type <- leave_out_type
  cells$testing_data <- unname(LEAVE_OUT_LABELS[as.character(leave_out_type)])
  cells$activity_stratum_index <- stratum_index
  cells$activity_stratum <- unname(
    STATE_STRATA_LABELS[as.character(stratum_index)]
  )

  cells
}

state_senspec_to_long <- function(senspec_list,file,file_name,fit_mix_num,
                                  model_type,cv_fold_id){
  if (is.null(senspec_list) || length(senspec_list) == 0){
    return(data.frame())
  }

  rows <- list()
  for (leave_out_type in seq_along(senspec_list)){
    if (!as.character(leave_out_type) %in% names(LEAVE_OUT_LABELS)){
      next
    }
    scenario_senspec <- senspec_list[[leave_out_type]]
    if (is.null(scenario_senspec) || length(scenario_senspec) == 0){
      next
    }

    stratum_indices <- seq_along(scenario_senspec)
    stratum_indices <- stratum_indices[
      as.character(stratum_indices) %in% names(STATE_STRATA_LABELS)
    ]
    for (stratum_index in stratum_indices){
      conf_mat <- scenario_senspec[[stratum_index]]
      if (is.null(conf_mat) || length(conf_mat) == 0){
        next
      }
      rows[[length(rows) + 1L]] <- state_confusion_matrix_to_long(
        conf_mat = conf_mat,
        file = file,
        file_name = file_name,
        fit_mix_num = fit_mix_num,
        model_type = model_type,
        cv_fold_id = cv_fold_id,
        leave_out_type = leave_out_type,
        stratum_index = stratum_index
      )
    }
  }

  rbind_fill(rows)
}

calc_class_metrics <- function(conf_mat){
  conf_mat <- as.matrix(conf_mat)
  storage.mode(conf_mat) <- "numeric"

  total_count <- sum(conf_mat)
  if (!is.finite(total_count) || total_count <= 0){
    stop("Confusion matrix has no observations")
  }

  true_count <- colSums(conf_mat)
  predicted_count <- rowSums(conf_mat)
  true_positive <- diag(conf_mat)
  false_negative <- true_count - true_positive
  false_positive <- predicted_count - true_positive
  true_negative <- total_count - true_positive - false_negative - false_positive

  sensitivity_denominator <- true_positive + false_negative
  specificity_denominator <- true_negative + false_positive
  precision_denominator <- true_positive + false_positive

  data.frame(
    latent_class = colnames(conf_mat),
    true_positive = true_positive,
    false_negative = false_negative,
    false_positive = false_positive,
    true_negative = true_negative,
    true_count = true_count,
    predicted_count = predicted_count,
    total_count = total_count,
    sensitivity = ifelse(
      sensitivity_denominator > 0,
      true_positive / sensitivity_denominator,
      NA_real_
    ),
    specificity = ifelse(
      specificity_denominator > 0,
      true_negative / specificity_denominator,
      NA_real_
    ),
    positive_predictive_value = ifelse(
      precision_denominator > 0,
      true_positive / precision_denominator,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
}

build_group_metrics <- function(confusion_long,
                                group_cols = c(
                                  "fit_mix_num",
                                  "model_type",
                                  "leave_out_type",
                                  "testing_data"
                                )){
  if (nrow(confusion_long) == 0){
    return(list(class_summary = data.frame(),macro_summary = data.frame()))
  }

  groups <- unique(confusion_long[,group_cols,drop = FALSE])

  class_rows <- list()
  macro_rows <- list()

  for (group_index in seq_len(nrow(groups))){
    current_group <- groups[group_index,,drop = FALSE]
    keep <- rep(TRUE,nrow(confusion_long))
    for (group_col in group_cols){
      keep <- keep & confusion_long[[group_col]] == current_group[[group_col]]
    }
    current <- confusion_long[keep,,drop = FALSE]

    predicted_levels <- sort(unique(current$mix_assignment_pred))
    true_levels <- sort(unique(current$mix_assignment_true_ind))
    class_levels <- union(predicted_levels,true_levels)
    conf_mat <- matrix(
      0,
      nrow = length(class_levels),
      ncol = length(class_levels),
      dimnames = list(class_levels,class_levels)
    )
    for (row_index in seq_len(nrow(current))){
      conf_mat[
        current$mix_assignment_pred[[row_index]],
        current$mix_assignment_true_ind[[row_index]]
      ] <- conf_mat[
        current$mix_assignment_pred[[row_index]],
        current$mix_assignment_true_ind[[row_index]]
      ] + current$count[[row_index]]
    }

    class_metrics <- calc_class_metrics(conf_mat)
    class_metrics <- cbind(
      current_group[rep(1,nrow(class_metrics)),,drop = FALSE],
      class_metrics
    )
    class_rows[[length(class_rows) + 1L]] <- class_metrics

    macro_rows[[length(macro_rows) + 1L]] <- data.frame(
      current_group,
      macro_sensitivity = mean_or_na(class_metrics$sensitivity),
      macro_specificity = mean_or_na(class_metrics$specificity),
      macro_positive_predictive_value =
        mean_or_na(class_metrics$positive_predictive_value),
      accuracy = sum(class_metrics$true_positive) /
        unique(class_metrics$total_count),
      total_count = unique(class_metrics$total_count),
      class_count = nrow(class_metrics),
      fold_count = length(unique(current$cv_fold_id)),
      stringsAsFactors = FALSE
    )
  }

  class_summary <- rbind_fill(class_rows)
  macro_summary <- rbind_fill(macro_rows)

  order_cols <- order(
    macro_summary$model_type,
    macro_summary$fit_mix_num,
    macro_summary$leave_out_type
  )
  macro_summary <- macro_summary[order_cols,,drop = FALSE]
  row.names(macro_summary) <- NULL

  order_class <- order(
    class_summary$model_type,
    class_summary$fit_mix_num,
    class_summary$leave_out_type,
    suppressWarnings(as.numeric(class_summary$latent_class))
  )
  class_summary <- class_summary[order_class,,drop = FALSE]
  row.names(class_summary) <- NULL

  list(class_summary = class_summary,macro_summary = macro_summary)
}

summarize_state_sensitivity <- function(state_confusion_long,
                                        group_cols = c(
                                          "fit_mix_num",
                                          "model_type",
                                          "leave_out_type",
                                          "testing_data"
                                        )){
  if (nrow(state_confusion_long) == 0){
    return(data.frame())
  }

  context_cols <- c(
    group_cols,
    "activity_stratum_index",
    "activity_stratum"
  )
  context_cols <- unique(context_cols)
  groups <- unique(state_confusion_long[,context_cols,drop = FALSE])

  rows <- list()
  for (group_index in seq_len(nrow(groups))){
    current_group <- groups[group_index,,drop = FALSE]
    keep <- rep(TRUE,nrow(state_confusion_long))
    for (group_col in context_cols){
      keep <- keep & state_confusion_long[[group_col]] ==
        current_group[[group_col]]
    }
    current <- state_confusion_long[keep,,drop = FALSE]

    total_count <- sum(current$count)
    for (state in STATE_LABELS){
      denominator <- sum(
        current$count[current$true_state == state],
        na.rm = TRUE
      )
      correct <- sum(
        current$count[
          current$true_state == state &
            current$predicted_state == state
        ],
        na.rm = TRUE
      )

      rows[[length(rows) + 1L]] <- data.frame(
        current_group,
        state = state,
        correct_count = correct,
        true_count = denominator,
        total_count = total_count,
        sensitivity = ifelse(denominator > 0,correct / denominator,NA_real_),
        fold_count = length(unique(current$cv_fold_id)),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- rbind_fill(rows)
  if (nrow(out) == 0){
    return(out)
  }
  out$activity_stratum <- factor(
    out$activity_stratum,
    levels = STATE_STRATA_PLOT_ORDER
  )
  out$state <- factor(out$state,levels = STATE_LABELS)

  out <- out[order(
    out$model_type,
    out$fit_mix_num,
    out$leave_out_type,
    out$activity_stratum,
    out$state
  ),,drop = FALSE]
  row.names(out) <- NULL
  out
}

parse_one_file <- function(file){
  file_name <- basename(file)
  to_save <- load_to_save(file)
  settings <- get_section(
    to_save,
    "settings",
    SAVED_SECTION_SLOT[["settings"]],
    required = FALSE
  )
  leave_out <- get_section(
    to_save,
    "leave_out",
    SAVED_SECTION_SLOT[["leave_out"]],
    required = TRUE
  )

  fit_mix_num <- integer_or_na(
    get_setting(
      settings,
      "fit_mix_num",
      regex_value("FitMix([0-9]+)",file_name)
    )
  )
  if (!is.finite(fit_mix_num)){
    stop("Cannot determine fit_mix_num from settings or filename")
  }

  model_type <- as.character(
    get_setting(settings,"model_type",model_type_from_file(file_name))
  )
  cv_fold_id <- integer_or_na(
    get_leave_out_field(
      leave_out,
      "cv_fold_id",
      slot = length(leave_out) + 1L,
      default = get_setting(
        settings,
        "cv_fold_id",
        get_setting(settings,"sim_num",regex_value("Seed([0-9]+)",file_name))
      )
    )
  )

  leave_out_inds <- get_leave_out_field(
    leave_out,
    "leave_out_inds",
    LEAVE_OUT_SLOT[["leave_out_inds"]],
    default = NULL
  )
  leave_out_n <- if (is.null(leave_out_inds)) NA_integer_ else length(leave_out_inds)

  conf_mat_list <- get_leave_out_field(
    leave_out,
    "conf_mat_list",
    LEAVE_OUT_SLOT[["conf_mat_list"]],
    default = NULL
  )
  if (is.null(conf_mat_list) || length(conf_mat_list) == 0){
    stop("Missing leave_out$conf_mat_list")
  }

  scenario_indices <- seq_along(LEAVE_OUT_LABELS)
  available_scenarios <- scenario_indices[
    scenario_indices <= length(conf_mat_list) &
      vapply(
        conf_mat_list[scenario_indices],
        function(x) !is.null(x) && length(x) > 0,
        logical(1)
      )
  ]
  missing_scenarios <- setdiff(scenario_indices,available_scenarios)
  if (length(missing_scenarios) > 0){
    warning(
      "Missing leave-out scenarios in ",
      file_name,
      ": ",
      paste(missing_scenarios,collapse = ", ")
    )
  }

  rows <- list()
  for (leave_out_type in available_scenarios){
    conf_mat <- conf_mat_list[[leave_out_type]]
    rows[[length(rows) + 1L]] <- confusion_matrix_to_long(
      conf_mat = conf_mat,
      file = file,
      file_name = file_name,
      fit_mix_num = fit_mix_num,
      model_type = model_type,
      cv_fold_id = cv_fold_id,
      leave_out_type = leave_out_type,
      leave_out_n = leave_out_n
    )
  }

  senspec_list <- get_leave_out_field(
    leave_out,
    "senspec_list",
    LEAVE_OUT_SLOT[["senspec_list"]],
    default = NULL
  )

  list(
    confusion_long = rbind_fill(rows),
    state_confusion_long = state_senspec_to_long(
      senspec_list = senspec_list,
      file = file,
      file_name = file_name,
      fit_mix_num = fit_mix_num,
      model_type = model_type,
      cv_fold_id = cv_fold_id
    )
  )
}

make_plot_data <- function(macro_summary){
  sensitivity_rows <- macro_summary
  sensitivity_rows$metric <- "Sensitivity"
  sensitivity_rows$value <- sensitivity_rows$macro_sensitivity

  specificity_rows <- macro_summary
  specificity_rows$metric <- "Specificity"
  specificity_rows$value <- specificity_rows$macro_specificity

  plot_data <- rbind(
    sensitivity_rows,
    specificity_rows
  )
  plot_data$metric <- factor(
    plot_data$metric,
    levels = c("Sensitivity","Specificity")
  )
  plot_data$testing_data <- factor(
    plot_data$testing_data,
    levels = unname(LEAVE_OUT_LABELS)
  )
  plot_data$model_label <- paste0(
    ifelse(plot_data$model_type == "joint","Joint","Two-stage"),
    ", ",
    plot_data$fit_mix_num,
    " LC"
  )
  plot_data
}

make_state_plot_data <- function(state_summary,testing_data_filter = NULL){
  if (!is.null(testing_data_filter)){
    state_summary <- state_summary[
      as.character(state_summary$testing_data) %in% testing_data_filter,
      ,drop = FALSE
    ]
  }
  state_summary$activity_stratum <- factor(
    as.character(state_summary$activity_stratum),
    levels = STATE_STRATA_PLOT_ORDER
  )
  state_summary$state <- factor(
    as.character(state_summary$state),
    levels = STATE_LABELS
  )
  state_summary$testing_data <- factor(
    as.character(state_summary$testing_data),
    levels = unname(LEAVE_OUT_LABELS)
  )
  state_summary$model_label <- paste0(
    ifelse(state_summary$model_type == "joint","Joint","Two-stage"),
    ", ",
    state_summary$fit_mix_num,
    " LC"
  )
  state_summary
}

save_base_state_sensitivity_plot <- function(filename,plot_data,device,
                                             title,palette = PLOT_PALETTE){
  if (device == "png"){
    grDevices::png(filename,width = 3600,height = 2100,res = 300)
  } else if (device == "pdf"){
    grDevices::pdf(filename,width = 12,height = 7)
  } else {
    stop("Unsupported plot device: ",device)
  }
  on.exit(grDevices::dev.off(),add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par),add = TRUE)

  plot_data <- plot_data[is.finite(plot_data$sensitivity),,drop = FALSE]
  model_labels <- unique(plot_data$model_label)
  strata_labels <- STATE_STRATA_PLOT_ORDER[
    STATE_STRATA_PLOT_ORDER %in% as.character(unique(plot_data$activity_stratum))
  ]
  testing_labels <- unname(LEAVE_OUT_LABELS)[
    unname(LEAVE_OUT_LABELS) %in%
      as.character(unique(plot_data$testing_data))
  ]

  if (length(model_labels) == 0 ||
      length(strata_labels) == 0 ||
      length(testing_labels) == 0){
    warning("No latent-state sensitivity rows available to plot: ",filename)
    return(invisible(NULL))
  }

  panel_count <- length(model_labels) * length(strata_labels)
  layout_matrix <- matrix(
    seq_len(panel_count),
    nrow = length(model_labels),
    ncol = length(strata_labels),
    byrow = TRUE
  )
  layout_matrix <- rbind(
    layout_matrix,
    rep(panel_count + 1L,length(strata_labels))
  )

  graphics::layout(
    layout_matrix,
    heights = c(rep(4,length(model_labels)),1.0)
  )

  graphics::par(mar = c(5,4.5,4,1),oma = c(0,0,2.25,0))

  for (model_label in model_labels){
    for (stratum_label in strata_labels){
      current <- plot_data[
        plot_data$model_label == model_label &
          as.character(plot_data$activity_stratum) == stratum_label,
        ,drop = FALSE
      ]

      value_mat <- matrix(
        NA_real_,
        nrow = length(testing_labels),
        ncol = length(STATE_LABELS),
        dimnames = list(testing_labels,STATE_LABELS)
      )
      for (row_index in seq_len(nrow(current))){
        value_mat[
          as.character(current$testing_data[[row_index]]),
          as.character(current$state[[row_index]])
        ] <- current$sensitivity[[row_index]]
      }

      panel_title <- if (length(model_labels) > 1){
        paste(stratum_label,model_label,sep = "\n")
      } else {
        stratum_label
      }

      graphics::barplot(
        value_mat,
        beside = TRUE,
        ylim = c(0,1.05),
        col = palette[rownames(value_mat)],
        border = "black",
        ylab = "Sensitivity",
        main = panel_title,
        las = 1,
        cex.names = 1.05,
        cex.axis = 0.95,
        cex.lab = 1.2
      )
    }
  }

  graphics::par(mar = c(0,0,0,0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = testing_labels,
    fill = palette[testing_labels],
    border = "black",
    title = "Testing Data",
    horiz = length(testing_labels) <= 3,
    bty = "n",
    cex = 1.0
  )

  graphics::mtext(title,outer = TRUE,cex = 1.45,line = 0)
  invisible(NULL)
}

save_state_sensitivity_plot <- function(plot_data,file_png,file_pdf = NULL,
                                        title =
                                          paste(
                                            "Cross-Validated Latent State",
                                            "Sensitivity by Testing Data"
                                          ),
                                        palette = PLOT_PALETTE){
  plot_data <- plot_data[is.finite(plot_data$sensitivity),,drop = FALSE]
  if (nrow(plot_data) == 0){
    warning("No latent-state sensitivity rows available to plot")
    return(invisible(NULL))
  }

  if (!is.null(file_png)){
    plot_dir <- dirname(file_png)
    if (!dir.exists(plot_dir)){
      dir.create(plot_dir,recursive = TRUE,showWarnings = FALSE)
    }
  }
  if (!is.null(file_pdf)){
    plot_dir <- dirname(file_pdf)
    if (!dir.exists(plot_dir)){
      dir.create(plot_dir,recursive = TRUE,showWarnings = FALSE)
    }
  }

  if (requireNamespace("ggplot2",quietly = TRUE)){
    state_plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = state,y = sensitivity,fill = testing_data)
    ) +
      ggplot2::geom_col(
        position = ggplot2::position_dodge(width = 0.78),
        width = 0.72,
        color = "black",
        linewidth = 0.25
      ) +
      ggplot2::facet_grid(
        ggplot2::vars(model_label),
        ggplot2::vars(activity_stratum)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0,1.05),
        breaks = seq(0,1,by = 0.25),
        name = "Sensitivity"
      ) +
      ggplot2::scale_fill_manual(
        values = palette,
        drop = FALSE,
        name = "Testing Data"
      ) +
      ggplot2::labs(x = NULL,title = title) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 20),
        axis.title.y = ggplot2::element_text(size = 16),
        axis.text = ggplot2::element_text(size = 14),
        legend.position = "bottom",
        legend.title = ggplot2::element_text(size = 15),
        legend.text = ggplot2::element_text(size = 13)
      )

    ggplot2::ggsave(
      filename = file_png,
      plot = state_plot,
      width = 12,
      height = 7,
      dpi = 300
    )
    if (!is.null(file_pdf)){
      ggplot2::ggsave(
        filename = file_pdf,
        plot = state_plot,
        width = 12,
        height = 7
      )
    }
    return(invisible(state_plot))
  }

  save_base_state_sensitivity_plot(file_png,plot_data,"png",title,palette)
  if (!is.null(file_pdf)){
    save_base_state_sensitivity_plot(file_pdf,plot_data,"pdf",title,palette)
  }

  invisible(NULL)
}

files <- list.files(input_dir,pattern = "[.](rda|rds)$",
                    full.names = TRUE,ignore.case = TRUE)
if (length(files) == 0){
  stop("No .rda or .rds files found in: ",input_dir)
}

message("Parsing ",length(files)," files from ",normalizePath(input_dir))

rows <- list()
state_rows <- list()
errors <- list()
for (file in files){
  parsed <- tryCatch(
    parse_one_file(file),
    error = function(e){
      errors[[length(errors) + 1L]] <<- data.frame(
        file = normalizePath(file,mustWork = TRUE),
        file_name = basename(file),
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      list(confusion_long = data.frame(),state_confusion_long = data.frame())
    }
  )

  if (!is.null(parsed$confusion_long) && nrow(parsed$confusion_long) > 0){
    rows[[length(rows) + 1L]] <- parsed$confusion_long
  }
  if (!is.null(parsed$state_confusion_long) &&
      nrow(parsed$state_confusion_long) > 0){
    state_rows[[length(state_rows) + 1L]] <- parsed$state_confusion_long
  }
}

confusion_long <- rbind_fill(rows)
state_confusion_long <- rbind_fill(state_rows)
error_data <- rbind_fill(errors)

if (nrow(confusion_long) == 0){
  stop("No confusion matrices could be parsed")
}

fold_metrics <- build_group_metrics(
  confusion_long[,c(
    "file_name",
    "fit_mix_num",
    "model_type",
    "cv_fold_id",
    "leave_out_type",
    "testing_data",
    "mix_assignment_pred",
    "mix_assignment_true_ind",
    "count"
  ),drop = FALSE],
  group_cols = c(
    "file_name",
    "fit_mix_num",
    "model_type",
    "cv_fold_id",
    "leave_out_type",
    "testing_data"
  )
)
names(fold_metrics$class_summary)[
  names(fold_metrics$class_summary) == "total_count"
] <- "fold_total_count"
names(fold_metrics$macro_summary)[
  names(fold_metrics$macro_summary) == "total_count"
] <- "fold_total_count"

pooled_metrics <- build_group_metrics(
  confusion_long[,c(
    "fit_mix_num",
    "model_type",
    "leave_out_type",
    "testing_data",
    "mix_assignment_pred",
    "mix_assignment_true_ind",
    "count",
    "cv_fold_id"
  ),drop = FALSE]
)

state_sensitivity_by_fold <- summarize_state_sensitivity(
  state_confusion_long,
  group_cols = c(
    "file_name",
    "fit_mix_num",
    "model_type",
    "cv_fold_id",
    "leave_out_type",
    "testing_data"
  )
)
state_sensitivity_summary <- summarize_state_sensitivity(
  state_confusion_long,
  group_cols = c(
    "fit_mix_num",
    "model_type",
    "leave_out_type",
    "testing_data"
  )
)

results <- list(
  input_dir = normalizePath(input_dir,mustWork = TRUE),
  output_prefix = output_file,
  support_plot_file = support_plot_file,
  light_support_plot_file = light_support_plot_file,
  light_plot_file = light_plot_file,
  file_count = length(files),
  parsed_file_count = length(unique(confusion_long$file_name)),
  error_count = nrow(error_data),
  errors = error_data,
  confusion_long = confusion_long,
  fold_class_summary = fold_metrics$class_summary,
  fold_macro_summary = fold_metrics$macro_summary,
  class_summary = pooled_metrics$class_summary,
  macro_summary = pooled_metrics$macro_summary,
  state_confusion_long = state_confusion_long,
  state_sensitivity_by_fold = state_sensitivity_by_fold,
  state_sensitivity_summary = state_sensitivity_summary
)

dir.create(dirname(output_file),recursive = TRUE,showWarnings = FALSE)
saveRDS(results,output_file)

