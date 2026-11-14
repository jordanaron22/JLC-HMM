#!/usr/bin/env Rscript

# Reproduce the NHANES LC4/LC5/LC6 neighboring-class tables directly from
# the final joint-model fits. Run with:
#
#   Rscript Output/PlotCode/NHANES/CalculateNHANESNeighboringClassTables.R
#
# Classes are ordered by increasing estimated mortality coefficient. Each
# participant is assigned to the class with the largest posterior probability.
# The saved fits do not contain participant IDs, so cross-model comparisons
# rely on the common participant row order used when the models were fitted.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_file <- normalizePath(
    sub("^--file=", "", script_arg), mustWork = TRUE
  )
  project_dir <- normalizePath(
    file.path(dirname(script_file), "..", "..", ".."), mustWork = TRUE
  )
} else {
  project_dir <- normalizePath(".", mustWork = TRUE)
}

data_dir <- file.path(project_dir, "Data")
table_dir <- file.path(project_dir, "Output", "PlotCode", "NHANES")
class_table_file <- file.path(table_dir, "NHANESClassPAHRSensitivityTable.tex")
stability_table_file <- file.path(
  table_dir, "NHANESJointClassStabilityTables.tex"
)

model_files <- c(
  `4` = file.path(data_dir, "JMHMMFitMix4Seed.rda"),
  `5` = file.path(data_dir, "JMHMMFitMix5Seed.rda"),
  `6` = file.path(data_dir, "JMHMMFitMix6Seed.rda")
)

EXPECTED_PARTICIPANTS <- 6483L
PROBABILITY_TOLERANCE <- 1e-6
SECONDARY_LC6_THRESHOLD <- 0.10

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

load_to_save <- function(file) {
  assert_true(file.exists(file), paste0("Model file not found: ", file))
  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(file, envir = load_env)
  assert_true(
    "to_save" %in% loaded_names,
    paste0("Expected an object named to_save in ", basename(file))
  )
  load_env$to_save
}

load_model <- function(class_count, file) {
  fit <- load_to_save(file)
  settings <- fit$settings
  estimates <- fit$est_params

  assert_true(is.list(settings), paste0("Missing settings in ", basename(file)))
  assert_true(
    is.list(estimates), paste0("Missing est_params in ", basename(file))
  )
  assert_true(
    identical(tolower(as.character(settings$model_type)), "joint"),
    paste0(basename(file), " is not a joint-model fit")
  )
  assert_true(
    as.integer(settings$fit_mix_num) == class_count,
    paste0(basename(file), " does not contain the expected LC", class_count,
           " fit")
  )

  posterior <- as.matrix(estimates$re_prob)
  storage.mode(posterior) <- "numeric"
  assert_true(
    identical(dim(posterior), c(EXPECTED_PARTICIPANTS, class_count)),
    paste0("Unexpected re_prob dimensions in ", basename(file))
  )
  assert_true(
    all(is.finite(posterior)) && all(posterior >= 0) && all(posterior <= 1),
    paste0("Invalid posterior probabilities in ", basename(file))
  )
  assert_true(
    max(abs(rowSums(posterior) - 1)) <= PROBABILITY_TOLERANCE,
    paste0("Posterior rows do not sum to one in ", basename(file))
  )

  beta <- as.numeric(estimates$beta_vec)
  assert_true(
    length(beta) == class_count && all(is.finite(beta)),
    paste0("Invalid beta_vec in ", basename(file))
  )

  activity <- estimates$emit_act
  assert_true(
    is.array(activity) && identical(dim(activity), c(2L, 2L, class_count, 2L)) &&
      all(is.finite(activity)),
    paste0("Invalid emit_act array in ", basename(file))
  )

  # Apply the intended increasing-mortality-risk class order to all quantities.
  class_order <- order(beta)
  beta <- beta[class_order] - min(beta)
  posterior <- posterior[, class_order, drop = FALSE]
  activity <- activity[, , class_order, , drop = FALSE]
  assignment <- max.col(posterior, ties.method = "first")

  list(
    classes = class_count,
    file = file,
    original_class_order = class_order,
    beta = beta,
    hazard_ratio = exp(beta),
    posterior = posterior,
    assignment = assignment,
    activity = activity
  )
}

# Base-R adjusted Rand index for two hard class assignments.
adjusted_rand_index <- function(first, second) {
  assert_true(length(first) == length(second), "ARI inputs have different lengths")
  contingency <- table(first, second)
  choose_two <- function(x) sum(x * (x - 1) / 2)

  total <- sum(contingency)
  total_pairs <- total * (total - 1) / 2
  if (total_pairs == 0) {
    return(1)
  }

  observed <- choose_two(contingency)
  row_pairs <- choose_two(rowSums(contingency))
  column_pairs <- choose_two(colSums(contingency))
  expected <- row_pairs * column_pairs / total_pairs
  denominator <- (row_pairs + column_pairs) / 2 - expected

  if (abs(denominator) <= .Machine$double.eps) {
    return(if (abs(observed - expected) <= .Machine$double.eps) 1 else 0)
  }
  (observed - expected) / denominator
}

class_summary <- function(model) {
  counts <- tabulate(model$assignment, nbins = model$classes)
  data.frame(
    model = model$classes,
    class = seq_len(model$classes),
    count = counts,
    proportion = counts / length(model$assignment),
    hazard_ratio = model$hazard_ratio,
    weekday_wake_pa = model$activity[1L, 1L, , 1L],
    weekday_sleep_pa = model$activity[2L, 1L, , 1L],
    weekend_wake_pa = model$activity[1L, 1L, , 2L],
    weekend_sleep_pa = model$activity[2L, 1L, , 2L],
    row.names = NULL
  )
}

cross_classification <- function(reference_model, neighboring_model) {
  counts <- table(
    factor(reference_model$assignment, levels = seq_len(reference_model$classes)),
    factor(neighboring_model$assignment,
           levels = seq_len(neighboring_model$classes))
  )
  row_proportions <- counts / rowSums(counts)

  assert_true(
    identical(as.integer(rowSums(counts)),
              tabulate(reference_model$assignment,
                       nbins = reference_model$classes)),
    "Cross-classification rows do not reproduce the LC5 class totals"
  )
  assert_true(
    sum(counts) == EXPECTED_PARTICIPANTS,
    "Cross-classification does not contain all participants"
  )

  list(
    counts = counts,
    row_proportions = row_proportions,
    ari = adjusted_rand_index(
      reference_model$assignment, neighboring_model$assignment
    )
  )
}

format_count <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_percent <- function(x) {
  paste0(sprintf("%.1f", 100 * x), "\\%")
}

format_count_percent <- function(count, proportion) {
  paste0(format_count(count), " (", format_percent(proportion), ")")
}

format_decimal <- function(x) sprintf("%.3f", x)

format_pa <- function(x) {
  value <- format_decimal(x)
  if (x < 0) paste0("$", value, "$") else value
}

make_class_table <- function(summary_data) {
  lines <- c(
    "% Requires \\usepackage{booktabs} and \\usepackage{graphicx}.",
    "",
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0(
      "\\caption{Class size, mortality risk, and physical activity profiles ",
      "across the four-, five-, and six-class joint models}"
    ),
    "\\label{tab:nhanes-class-pa-hr-sensitivity}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{ccrrrrrr}",
    "\\toprule",
    paste0(
      "Model ($G$) & Class & $n$ (\\%) & HR & ",
      "\\multicolumn{2}{c}{Weekday PA mean} & ",
      "\\multicolumn{2}{c}{Weekend PA mean} \\\\"
    ),
    "\\cmidrule(lr){5-6} \\cmidrule(lr){7-8}",
    "& & & & Wake & Sleep & Wake & Sleep \\\\",
    "\\midrule"
  )

  for (class_count in 4:6) {
    rows <- summary_data[summary_data$model == class_count, , drop = FALSE]
    for (row_index in seq_len(nrow(rows))) {
      row <- rows[row_index, ]
      cells <- c(
        row$model,
        row$class,
        format_count_percent(row$count, row$proportion),
        format_decimal(row$hazard_ratio),
        format_pa(row$weekday_wake_pa),
        format_pa(row$weekday_sleep_pa),
        format_pa(row$weekend_wake_pa),
        format_pa(row$weekend_sleep_pa)
      )
      lines <- c(lines, paste0(paste(cells, collapse = " & "), " \\\\"))
    }
    if (class_count < 6L) {
      lines <- c(lines, "\\addlinespace")
    }
  }

  c(
    lines,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\begin{minipage}{0.92\\linewidth}",
    "\\footnotesize",
    paste0(
      "\\textit{Note.} Classes within each model are ordered by increasing ",
      "estimated"
    ),
    "mortality hazard. HR denotes the class-specific mortality hazard ratio relative",
    paste0(
      "to Class 1 of the same model. Class sizes use maximum-posterior ",
      "assignments."
    ),
    "Wake and sleep PA means are estimated class-specific weekday and weekend",
    paste0(
      "means on the log-MIMS scale; consequently, a less negative sleep PA mean ",
      "indicates"
    ),
    "greater sleep-period activity.",
    "\\end{minipage}",
    "\\end{table}"
  )
}

make_cross_table <- function(comparison, neighboring_classes) {
  model_word <- if (neighboring_classes == 4L) "four" else "six"
  counts <- comparison$counts
  percentages <- comparison$row_proportions
  last_column <- neighboring_classes + 1L
  column_spec <- paste0("l", paste(rep("r", neighboring_classes + 1L),
                                   collapse = ""))

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0(
      "\\caption{Cross-classification of the five-class and ", model_word,
      "-class joint models}"
    ),
    paste0("\\label{tab:nhanes-five-to-", model_word, "}"),
    "\\resizebox{\\linewidth}{!}{%",
    paste0("\\begin{tabular}{", column_spec, "}"),
    "\\toprule",
    paste0(
      "& \\multicolumn{", neighboring_classes, "}{c}{",
      tools::toTitleCase(model_word), "-class model} & \\\\"
    ),
    paste0("\\cmidrule(lr){2-", last_column, "}"),
    "Five-class model",
    paste0(
      "  & ",
      paste(paste("Class", seq_len(neighboring_classes)), collapse = " & "),
      " & Total \\\\"
    ),
    "\\midrule"
  )

  for (row_index in seq_len(nrow(counts))) {
    cells <- vapply(seq_len(ncol(counts)), function(column_index) {
      format_count_percent(
        counts[row_index, column_index],
        percentages[row_index, column_index]
      )
    }, character(1))
    lines <- c(
      lines,
      paste0(
        "Class ", row_index, " & ", paste(cells, collapse = " & "),
        " & ", format_count(sum(counts[row_index, ])), " \\\\"
      )
    )
  }

  lines <- c(
    lines,
    "\\midrule",
    paste0(
      "Total & ", paste(format_count(colSums(counts)), collapse = " & "),
      " & ", format_count(sum(counts)), " \\\\"
    ),
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\begin{minipage}{\\linewidth}",
    "\\footnotesize",
    paste0(
      "\\textit{Note.} Values are numbers of participants, with percentages ",
      "calculated"
    ),
    "within each row of the five-class model. Classes in both models are ordered by",
    paste0(
      "increasing estimated mortality hazard. The adjusted Rand index is ",
      sprintf("%.3f", comparison$ari), "."
    ),
    "\\end{minipage}",
    "\\end{table}"
  )
  lines
}

correspondence_entries <- function(row_proportions, model, include_secondary) {
  dominant <- which.max(row_proportions)
  selected <- dominant
  if (include_secondary) {
    selected <- which(row_proportions >= SECONDARY_LC6_THRESHOLD)
    selected <- unique(c(dominant, selected))
    selected <- selected[order(row_proportions[selected], decreasing = TRUE)]
  }

  vapply(selected, function(class_number) {
    paste0(
      "Class ", class_number, ": ", format_percent(row_proportions[class_number]),
      ", HR ", format_decimal(model$hazard_ratio[class_number])
    )
  }, character(1))
}

correspondence_text <- function(row_proportions, model, include_secondary) {
  entries <- correspondence_entries(row_proportions, model, include_secondary)
  if (length(entries) == 1L) {
    return(entries)
  }
  paste0(
    "\\begin{tabular}[t]{@{}l@{}}",
    paste(entries, collapse = " \\\\\n      "),
    "\\end{tabular}"
  )
}

make_anchored_table <- function(models, summary_data, comparison_4,
                                comparison_6) {
  lc5_summary <- summary_data[summary_data$model == 5L, , drop = FALSE]
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    paste0(
      "\\caption{Five-class-anchored comparison of latent classes and ",
      "mortality hazard ratios}"
    ),
    "\\label{tab:nhanes-five-class-anchored}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrll}",
    "\\toprule",
    paste0(
      "Five-class class & $n$ (\\%) & HR & Four-class correspondence & ",
      "Six-class correspondence \\\\"
    ),
    "\\midrule"
  )

  for (class_number in 1:5) {
    row <- lc5_summary[class_number, ]
    lc4_text <- correspondence_text(
      comparison_4$row_proportions[class_number, ], models[["4"]], FALSE
    )
    lc6_text <- correspondence_text(
      comparison_6$row_proportions[class_number, ], models[["6"]], TRUE
    )
    lines <- c(
      lines,
      paste0(
        "Class ", class_number, " & ",
        format_count_percent(row$count, row$proportion), " & ",
        format_decimal(row$hazard_ratio)
      ),
      paste0("  & ", lc4_text),
      paste0("  & ", lc6_text, " \\\\ ")
    )
  }

  c(
    lines,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    "\\begin{minipage}{\\linewidth}",
    "\\footnotesize",
    "\\textit{Note.} HR denotes the class-specific mortality hazard ratio relative",
    "to Class 1 of the same model. Classes are ordered by increasing estimated",
    "mortality hazard, so Class 1 is the reference and lowest-risk class in every",
    "model. Correspondence percentages are row percentages conditional on membership",
    "in the indicated five-class class. The dominant correspondence is reported,",
    paste0(
      "along with any secondary six-class correspondence containing at least 10\\% ",
      "of"
    ),
    "the five-class class.",
    "\\end{minipage}",
    "\\end{table}"
  )
}

models <- lapply(names(model_files), function(class_count) {
  load_model(as.integer(class_count), model_files[[class_count]])
})
names(models) <- names(model_files)

assert_true(
  all(vapply(models, function(model) {
    length(model$assignment) == EXPECTED_PARTICIPANTS
  }, logical(1))),
  "The fitted models do not share the expected participant count"
)

summary_data <- do.call(rbind, lapply(models, class_summary))
row.names(summary_data) <- NULL
comparison_4 <- cross_classification(models[["5"]], models[["4"]])
comparison_6 <- cross_classification(models[["5"]], models[["6"]])

# Regression checks for the final fits used when this transparent script was made.
assert_true(
  abs(comparison_4$ari - 0.507) < 0.001,
  paste0("Unexpected LC5-vs-LC4 ARI: ", format(comparison_4$ari, digits = 8))
)
assert_true(
  abs(comparison_6$ari - 0.774) < 0.001,
  paste0("Unexpected LC5-vs-LC6 ARI: ", format(comparison_6$ari, digits = 8))
)

cat("NHANES neighboring-class analysis\n")
cat("===================================\n")
cat("Hard assignments use each participant's largest posterior probability.\n")
cat("Participant IDs are not saved; comparisons use the common saved row order.\n\n")

for (class_count in names(models)) {
  model <- models[[class_count]]
  cat("LC", class_count, " original classes in increasing-risk order: ",
      paste(model$original_class_order, collapse = ", "), "\n", sep = "")
  cat("Ordered mortality coefficients: ",
      paste(format_decimal(model$beta), collapse = ", "), "\n", sep = "")
  print(summary_data[summary_data$model == as.integer(class_count), ],
        row.names = FALSE)
  cat("\n")
}

cat("LC5 by LC4 counts:\n")
print(comparison_4$counts)
cat("LC5 by LC4 row percentages:\n")
print(round(100 * comparison_4$row_proportions, 1))
cat("Adjusted Rand index:", format_decimal(comparison_4$ari), "\n\n")

cat("LC5 by LC6 counts:\n")
print(comparison_6$counts)
cat("LC5 by LC6 row percentages:\n")
print(round(100 * comparison_6$row_proportions, 1))
cat("Adjusted Rand index:", format_decimal(comparison_6$ari), "\n\n")

cat("LC5-anchored correspondences:\n")
for (class_number in 1:5) {
  lc4_entries <- correspondence_entries(
    comparison_4$row_proportions[class_number, ], models[["4"]], FALSE
  )
  lc6_entries <- correspondence_entries(
    comparison_6$row_proportions[class_number, ], models[["6"]], TRUE
  )
  lc4_console <- gsub("\\%", "%", paste(lc4_entries, collapse = "; "),
                      fixed = TRUE)
  lc6_console <- gsub("\\%", "%", paste(lc6_entries, collapse = "; "),
                      fixed = TRUE)
  cat(
    "LC5 Class ", class_number, ": LC4 ",
    lc4_console,
    "; LC6 ",
    lc6_console,
    "\n", sep = ""
  )
}

writeLines(make_class_table(summary_data), class_table_file, useBytes = TRUE)
stability_lines <- c(
  "% Requires \\usepackage{booktabs} and \\usepackage{graphicx}.",
  "% Classes in every model are ordered by increasing estimated mortality hazard.",
  "% Each cross-classification cell is: hard-assignment count (row percentage).",
  "",
  make_cross_table(comparison_4, 4L),
  "",
  make_cross_table(comparison_6, 6L),
  "",
  make_anchored_table(models, summary_data, comparison_4, comparison_6)
)
writeLines(stability_lines, stability_table_file, useBytes = TRUE)

cat("\nWrote:\n")
cat(" - ", class_table_file, "\n", sep = "")
cat(" - ", stability_table_file, "\n", sep = "")
