#!/usr/bin/env Rscript

# Summarize transition curves and predicted sleep duration for the final
# five-class NHANES joint model. This script prints results only; it does not
# refit the model or write any output files.
#
# Run from any directory with:
#   Rscript Output/PlotCode/NHANES/SummarizeNHANESLC5Interpretation.R

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
fit_file <- file.path(data_dir, "JMHMMFitMix5Seed.rda")
mortality_file <- file.path(
  data_dir, "NHANES_2011_2012_2013_2014.rda"
)
wave_g_file <- file.path(data_dir, "Wavedata_G.rda")
wave_h_file <- file.path(data_dir, "Wavedata_H.rda")

CLASS_COUNT <- 5L
PARTICIPANT_COUNT <- 6483L
PERIODS_PER_DAY <- 96L
OBSERVATION_DAYS <- 9L
COMPLETE_DAY_COUNT <- 8L
FIRST_COMPLETE_EPOCH <- 69L
MINIMUM_OBSERVED_EPOCHS <- PERIODS_PER_DAY * 6L
MAXIMUM_MISSING_EPOCHS_PER_DAY <- PERIODS_PER_DAY / 2L
EXPECTED_COMPLETE_PARTICIPANTS <- 6179L

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

load_environment <- function(file) {
  assert_true(file.exists(file), paste0("File not found: ", file))
  environment <- new.env(parent = emptyenv())
  load(file, envir = environment)
  environment
}

median_or_na <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(NA_real_)
  }
  median(values)
}

# Return the first value whose normalized cumulative weight reaches 0.5.
weighted_median <- function(values, weights) {
  keep <- is.finite(values) & is.finite(weights) & weights > 0
  values <- values[keep]
  weights <- weights[keep]
  assert_true(length(values) > 0L, "Weighted median has no usable values")

  order_index <- order(values)
  values <- values[order_index]
  weights <- weights[order_index] / sum(weights)
  values[which(cumsum(weights) >= 0.5)[1L]]
}

# Reproduce the weekday/weekend calendar used when the NHANES model was fit.
# Codes are 0 for weekday and 1 for weekend. Weekend begins Friday at 5 p.m.
# and ends Sunday at 5 p.m., so epochs 69:864 contain eight complete days.
first_day_to_weekend <- function(first_day) {
  weekday <- numeric(PERIODS_PER_DAY)
  friday <- c(numeric(68L), rep(1, 28L))
  saturday <- rep(1, PERIODS_PER_DAY)
  sunday <- c(rep(1, 68L), numeric(28L))

  calendar <- switch(
    as.character(as.integer(first_day)),
    `1` = c(sunday, rep(weekday, 4L), friday, saturday, sunday, weekday),
    `2` = c(rep(weekday, 4L), friday, saturday, sunday, rep(weekday, 2L)),
    `3` = c(rep(weekday, 3L), friday, saturday, sunday, rep(weekday, 3L)),
    `4` = c(rep(weekday, 2L), friday, saturday, sunday, rep(weekday, 4L)),
    `5` = c(weekday, friday, saturday, sunday, rep(weekday, 4L), friday),
    `6` = c(friday, saturday, sunday, rep(weekday, 4L), friday, saturday),
    `7` = c(saturday, sunday, rep(weekday, 4L), friday, saturday, sunday),
    NULL
  )

  assert_true(
    !is.null(calendar) && length(calendar) == PERIODS_PER_DAY * OBSERVATION_DAYS,
    paste0("Invalid first observation day: ", first_day)
  )
  calendar
}

# -----------------------------------------------------------------------------
# Load and validate the fitted model.
# -----------------------------------------------------------------------------

fit_environment <- load_environment(fit_file)
assert_true(
  exists("to_save", envir = fit_environment, inherits = FALSE),
  "The fitted model file does not contain to_save"
)
fit <- fit_environment$to_save
settings <- fit$settings
estimates <- fit$est_params

assert_true(is.list(settings), "The fitted model has no settings list")
assert_true(is.list(estimates), "The fitted model has no est_params list")
assert_true(
  identical(tolower(as.character(settings$model_type)), "joint") &&
    as.integer(settings$fit_mix_num) == CLASS_COUNT,
  "The fitted model is not the expected five-class joint model"
)
assert_true(
  as.integer(settings$period_len) == PERIODS_PER_DAY,
  "The fitted model does not use 96 observations per day"
)

transition_data <- estimates$tran_df
posterior <- as.matrix(estimates$re_prob)
decoded <- estimates$decoded_mat
beta <- as.numeric(estimates$beta_vec)

assert_true(
  is.data.frame(transition_data) &&
    all(c("prob", "type", "time", "weekend", "mixture") %in%
          names(transition_data)),
  "tran_df is missing required columns"
)
assert_true(
  identical(dim(posterior), c(PARTICIPANT_COUNT, CLASS_COUNT)),
  "re_prob does not have dimensions 6483 x 5"
)
assert_true(
  is.array(decoded) &&
    identical(
      dim(decoded),
      c(PERIODS_PER_DAY * OBSERVATION_DAYS,
        PARTICIPANT_COUNT, CLASS_COUNT)
    ),
  "decoded_mat does not have dimensions 864 x 6483 x 5"
)
assert_true(
  all(is.finite(posterior)) && all(posterior >= 0) && all(posterior <= 1) &&
    max(abs(rowSums(posterior) - 1)) <= 1e-6,
  "re_prob contains invalid posterior probabilities"
)
assert_true(
  all(is.na(decoded) | decoded %in% c(0, 1)),
  "decoded_mat contains values other than wake=0 and sleep=1"
)
assert_true(
  length(beta) == CLASS_COUNT && all(is.finite(beta)),
  "beta_vec is not a finite five-element vector"
)

# Classes in the paper are ordered from lowest to highest mortality coefficient.
class_order <- order(beta)
original_to_ordered <- integer(CLASS_COUNT)
original_to_ordered[class_order] <- seq_len(CLASS_COUNT)
posterior <- posterior[, class_order, drop = FALSE]
decoded <- decoded[, , class_order, drop = FALSE]
transition_data$latent_class <- original_to_ordered[
  as.integer(transition_data$mixture)
]

# -----------------------------------------------------------------------------
# Transition peaks and troughs.
# -----------------------------------------------------------------------------

transition_data$day_type <- ifelse(
  transition_data$weekend == 1, "Weekday",
  ifelse(transition_data$weekend == 2, "Weekend", NA_character_)
)
transition_data$transition <- ifelse(
  transition_data$type == "Falling Asleep", "Wake to Sleep",
  ifelse(transition_data$type == "Waking", "Sleep to Wake", NA_character_)
)

assert_true(
  !anyNA(transition_data$day_type) && !anyNA(transition_data$transition) &&
    all(is.finite(transition_data$prob)) &&
    all(is.finite(transition_data$time)) &&
    all(transition_data$latent_class %in% seq_len(CLASS_COUNT)),
  "tran_df contains invalid day, transition, class, time, or probability values"
)

transition_rows <- list()
row_index <- 1L
for (day_type in c("Weekday", "Weekend")) {
  for (transition in c("Wake to Sleep", "Sleep to Wake")) {
    for (latent_class in seq_len(CLASS_COUNT)) {
      group <- transition_data[
        transition_data$day_type == day_type &
          transition_data$transition == transition &
          transition_data$latent_class == latent_class,
        , drop = FALSE
      ]
      group <- group[order(group$time), , drop = FALSE]
      assert_true(
        nrow(group) == PERIODS_PER_DAY,
        paste0("Expected 96 transition points for ", day_type, ", ",
               transition, ", Class ", latent_class)
      )

      peak_probability <- max(group$prob)
      trough_probability <- min(group$prob)
      peak_row <- which(group$prob == peak_probability)[1L]
      trough_row <- which(group$prob == trough_probability)[1L]

      transition_rows[[row_index]] <- data.frame(
        class = latent_class,
        transition = transition,
        day_type = day_type,
        peak_time = group$time[peak_row],
        peak_probability = peak_probability,
        trough_time = group$time[trough_row],
        trough_probability = trough_probability,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
}

transition_summary <- do.call(rbind, transition_rows)
row.names(transition_summary) <- NULL
assert_true(nrow(transition_summary) == 20L, "Transition summary must have 20 rows")

# -----------------------------------------------------------------------------
# Reconstruct participant order, missingness, and weekday/weekend calendars.
# -----------------------------------------------------------------------------

mortality_environment <- load_environment(mortality_file)
wave_g_environment <- load_environment(wave_g_file)
wave_h_environment <- load_environment(wave_h_file)

assert_true(
  exists("NHANES_mort_list", envir = mortality_environment, inherits = FALSE),
  "Mortality file does not contain NHANES_mort_list"
)
assert_true(
  exists("wave_data_G", envir = wave_g_environment, inherits = FALSE) &&
    exists("wave_data_H", envir = wave_h_environment, inherits = FALSE),
  "Wave files do not contain wave_data_G and wave_data_H"
)

mortality_list <- mortality_environment$NHANES_mort_list
mortality_data <- rbind(
  mortality_list[[1]][mortality_list[[1]]$eligstat == 1, , drop = FALSE],
  mortality_list[[2]][mortality_list[[2]]$eligstat == 1, , drop = FALSE]
)
wave_g <- wave_g_environment$wave_data_G
wave_h <- wave_h_environment$wave_data_H

activity <- t(rbind(wave_g[[1]], wave_h[[1]])[, -1])
light <- t(rbind(wave_g[[2]], wave_h[[2]])[, -1])
participant_data <- rbind(wave_g[[3]], wave_h[[3]])

included <- participant_data$SEQN %in% mortality_data$seqn
mortality_included <- mortality_data$seqn %in% participant_data$SEQN
participant_data <- participant_data[included, , drop = FALSE]
activity <- activity[, included, drop = FALSE]
light <- light[, included, drop = FALSE]
mortality_data <- mortality_data[mortality_included, , drop = FALSE]

assert_true(
  identical(
    match(participant_data$SEQN, mortality_data$seqn),
    seq_len(nrow(mortality_data))
  ),
  "Participant order does not match the mortality-data order"
)
assert_true(
  nrow(participant_data) == PARTICIPANT_COUNT &&
    identical(dim(activity), c(PERIODS_PER_DAY * OBSERVATION_DAYS,
                               PARTICIPANT_COUNT)) &&
    identical(dim(light), dim(activity)),
  "Raw-data reconstruction does not match the 6,483 fitted participants"
)

first_day <- as.numeric(participant_data$PAXDAYWM)
weekend_calendar <- vapply(
  first_day, first_day_to_weekend,
  numeric(PERIODS_PER_DAY * OBSERVATION_DAYS)
)

complete_participant <-
  colSums(!is.na(activity)) > MINIMUM_OBSERVED_EPOCHS &
  colSums(!is.na(light)) > MINIMUM_OBSERVED_EPOCHS
complete_index <- which(complete_participant)
assert_true(
  length(complete_index) == EXPECTED_COMPLETE_PARTICIPANTS,
  paste0("Expected 6,179 complete participants; found ",
         length(complete_index))
)

day_start <- FIRST_COMPLETE_EPOCH +
  PERIODS_PER_DAY * (seq_len(COMPLETE_DAY_COUNT) - 1L)
day_indices <- lapply(day_start, function(start) {
  start:(start + PERIODS_PER_DAY - 1L)
})

valid_day <- t(vapply(day_indices, function(index) {
  colSums(is.na(activity[index, complete_index, drop = FALSE])) <=
    MAXIMUM_MISSING_EPOCHS_PER_DAY
}, logical(length(complete_index))))

weekend_day <- t(vapply(day_indices, function(index) {
  weekend_count <- colSums(
    weekend_calendar[index, complete_index, drop = FALSE]
  )
  assert_true(
    all(weekend_count %in% c(0L, PERIODS_PER_DAY)),
    "A complete 5 p.m.-to-5 p.m. window contains mixed day types"
  )
  weekend_count == PERIODS_PER_DAY
}, logical(length(complete_index))))

# -----------------------------------------------------------------------------
# Participant-first, posterior-probability-weighted sleep summaries.
# -----------------------------------------------------------------------------

sleep_rows <- vector("list", CLASS_COUNT)
for (latent_class in seq_len(CLASS_COUNT)) {
  class_decode <- decoded[, complete_index, latent_class, drop = FALSE]
  dim(class_decode) <- c(nrow(decoded), length(complete_index))

  sleep_hours <- t(vapply(day_indices, function(index) {
    colSums(class_decode[index, , drop = FALSE] == 1) /
      (PERIODS_PER_DAY / 24)
  }, numeric(length(complete_index))))
  sleep_hours[!valid_day] <- NA_real_

  participant_overall <- apply(sleep_hours, 2L, median_or_na)
  participant_weekday <- vapply(seq_along(complete_index), function(person) {
    median_or_na(sleep_hours[!weekend_day[, person], person])
  }, numeric(1))
  participant_weekend <- vapply(seq_along(complete_index), function(person) {
    median_or_na(sleep_hours[weekend_day[, person], person])
  }, numeric(1))

  class_weights <- posterior[complete_index, latent_class]
  available_overall <- is.finite(participant_overall)
  available_weekday <- is.finite(participant_weekday)
  available_weekend <- is.finite(participant_weekend)

  sleep_rows[[latent_class]] <- data.frame(
    class = latent_class,
    overall_median_hours = weighted_median(
      participant_overall, class_weights
    ),
    weekday_median_hours = weighted_median(
      participant_weekday, class_weights
    ),
    weekend_median_hours = weighted_median(
      participant_weekend, class_weights
    ),
    overall_available_n = sum(available_overall),
    overall_posterior_weight = sum(class_weights[available_overall]),
    weekday_available_n = sum(available_weekday),
    weekday_posterior_weight = sum(class_weights[available_weekday]),
    weekend_available_n = sum(available_weekend),
    weekend_posterior_weight = sum(class_weights[available_weekend]),
    stringsAsFactors = FALSE
  )
}

sleep_summary <- do.call(rbind, sleep_rows)
row.names(sleep_summary) <- NULL

expected_sleep_medians <- data.frame(
  class = seq_len(CLASS_COUNT),
  overall_median_hours = c(8.625, 12.000, 8.250, 8.625, 9.250),
  weekday_median_hours = c(8.500, 11.750, 8.250, 8.625, 8.750),
  weekend_median_hours = c(9.250, 12.750, 8.750, 9.000, 12.500)
)
median_columns <- c(
  "overall_median_hours", "weekday_median_hours", "weekend_median_hours"
)
assert_true(
  isTRUE(all.equal(
    sleep_summary[, median_columns],
    expected_sleep_medians[, median_columns],
    tolerance = 1e-10,
    check.attributes = FALSE
  )),
  "Sleep-duration medians do not match the expected current-fit values"
)

options(width = max(260L, getOption("width")))

cat("LC5 transition extrema\n")
cat("========================\n")
cat("Times use the model's 0.25-hour transition grid.\n")
print(transition_summary, row.names = FALSE, digits = 6)

cat("\nPosterior-weighted predicted sleep duration\n")
cat("=============================================\n")
cat(
  "Participant-first medians; eight complete 5 p.m.-to-5 p.m. windows; ",
  length(complete_index), " eligible participants.\n", sep = ""
)
print(sleep_summary, row.names = FALSE, digits = 7)
