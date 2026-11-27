#By default, this script displays the saved summary results and opens any available
#analysis figures in an interactive R session. Run from the project root.
#Only M2_summary_2021_2023.rds is required in results/2021_2023_models/ for this mode.
#Set run_models <- TRUE to load downloaded model fits or fit any missing models
#and recompute the results. This requires the analysis data and dependencies.
run_models <- FALSE

library(data.table)
results_dir <- "results/2021_2023_models"
figures_dir <- "Figures"
summary_file <- file.path(results_dir, "M2_summary_2021_2023.rds")

if (!run_models) {
  saved_summary <- readRDS(summary_file)
} else {
  library(glmmTMB)
  library(ggplot2)
  if (!requireNamespace("statmod", quietly = TRUE)) stop("Install statmod for normal quadrature.")
  source("R/rg_inference.R")
  source("R/poststrat_helpers.R")
  
  #Social isolation 2021-2023####
  years <- 2021:2023
  input_file <- "Intermediate Data/NWSG_data_2021_2023.RData"
  n_workers <- 10L
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  inputs <- new.env(parent = emptyenv())
  objects <- load(input_file, envir = inputs)
  tables <- objects[vapply(objects, function(x) {
    all(c("alone", "TUCASEID", "year_normalized") %in% names(inputs[[x]]))
  }, logical(1))]
  stopifnot(length(tables) == 1L)
  episodes <- as.data.table(inputs[[tables]])
  year <- as.integer(as.character(episodes$year_normalized))
  if ("TUYEAR" %in% names(episodes) && !all(is.na(episodes$TUYEAR))) {
    year <- as.integer(as.character(episodes$TUYEAR))
  } else if (all(na.omit(year) %in% 0:(length(years) - 1L))) {
    year <- year + min(years)
  } else if (max(year, na.rm = TRUE) < 100) {
    year <- year + 2003L
  }
  rows <- which(year %in% years)
  NWSG_data <- episodes[rows, .(
    alone, year_normalized, TUCASEID, TUACTDUR24, TUFNWGTP,
    age_category, HouseholdSize_category, education_category, hours_worked_category,
    PESEX, race_category, hispanic_category, marital_status, income_category,
    category, subcategory, TRCODEP, TUMONTH, TUDIARYDAY, State
  )]
  
  NWSG_data[, year_normalized := year[rows] - min(years)]
  rm(inputs, episodes)
  NWSG_data[, alone := as.numeric(as.character(alone))]
  NWSG_data[, PESEX := factor(PESEX)]
  NWSG_data <- droplevels(NWSG_data[complete.cases(NWSG_data) & alone %in% 0:1 &
                                      TUFNWGTP > 0 & TUACTDUR24 > 0])
  stopifnot(setequal(NWSG_data$year_normalized, years - min(years)))
  NWSG_data[, scaled_weight := TUFNWGTP * TUACTDUR24 / mean(TUFNWGTP * TUACTDUR24)]
  
  NWSG_data[, subj := TUCASEID]
  
  #Episode duration by year####
  frequencies <- NWSG_data[, .(weight = sum(TUFNWGTP)),
                           by = .(year = 2021L + year_normalized, TUACTDUR24)]
  upper <- log10(max(frequencies$TUACTDUR24))
  curves <- frequencies[, {
    x <- log10(TUACTDUR24)
    w <- weight / sum(weight)
    #Reflect at one minute, the shortest eligible duration
    curve <- density(c(x, -x), weights = c(w, w) / 2, bw = 0.09,
                     from = 0, to = upper, n = 512)
    .(minutes = 10^curve$x, density = 2 * curve$y)
  }, by = year]
  
  duration_plot <- ggplot(curves, aes(minutes, density, colour = factor(year))) +
    geom_line(linewidth = 0.9) +
    scale_x_log10(breaks = c(1, 5, 15, 30, 60, 120, 480, 960)) +
    scale_colour_manual(values = c("2021" = "#0072B2", "2022" = "#D55E00", "2023" = "#009E73")) +
    labs(x = "Episode duration (minutes, log scale)", y = "Density (per log10 minute)",
         colour = "Year") +
    theme_classic() +
    theme(legend.position = "bottom")
  
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figures_dir, "episode_duration_density_2021_2023.pdf"),
         duration_plot, width = 7, height = 4.5)
  print(duration_plot)
  
  #Rank activities by survey-weighted episode frequency####
  
  ranking <- NWSG_data[, .(frequency = sum(TUFNWGTP)), by = TRCODEP]
  ranking[, code := sprintf("%06d", as.integer(as.character(TRCODEP)))]
  ranking[, activity := unname(atus_activity_labels[code])]
  ranking[is.na(activity), activity := paste("Activity", code)]
  setorder(ranking, -frequency, code)
  ranking[, share := frequency / sum(frequency)]
  selected <- head(ranking, 20L)
  stopifnot(nrow(selected) == 20L, !anyDuplicated(selected$TRCODEP))
  activities <- as.character(selected$TRCODEP)
  
  #M1 weighted adjusted####
  fit_file <- file.path(results_dir, "M1_adjusted.RData")
  
  if (file.exists(fit_file)) {
    saved <- load(fit_file)
    if (!"time_weighted" %in% saved) time_weighted <- c(elapsed = NA_real_)
  } else {
    time_weighted <- system.time({
      M1_adjusted <- glmmTMB(
        alone ~ year_normalized + age_category + HouseholdSize_category + education_category +
          hours_worked_category + PESEX + race_category + hispanic_category +
          marital_status + income_category +
          (1 | TUCASEID) +
          (1 + year_normalized || category) +
          (1 + year_normalized || subcategory) +
          (1 + year_normalized || TRCODEP) +
          (1 + year_normalized || TUMONTH) +
          (1 + year_normalized || TUDIARYDAY) +
          (1 + year_normalized || State),
        data = NWSG_data, family = binomial(link = "logit"), weights = scaled_weight,
        control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
      )
    })
    save(M1_adjusted, time_weighted, file = fit_file)
  }
  
  
  #M1 weighted baseline####
  fit_file <- file.path(results_dir, "M1_unadjusted.RData")
  
  if (file.exists(fit_file)) {
    saved <- load(fit_file)
    if (!"time_weighted_baseline" %in% saved) time_weighted_baseline <- c(elapsed = NA_real_)
  } else {
    time_weighted_baseline <- system.time({
      M1_unadjusted <- glmmTMB(alone ~ year_normalized + (1 | TUCASEID) +
                                 (1 + year_normalized || category) + (1 + year_normalized || subcategory) +
                                 (1 + year_normalized || TRCODEP) + (1 + year_normalized || TUMONTH) +
                                 (1 + year_normalized || TUDIARYDAY),
                               data = NWSG_data, family = binomial("logit"), weights = scaled_weight,
                               control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000)))
    })
    save(M1_unadjusted, time_weighted_baseline, file = fit_file)
  }
  
  stopifnot(setequal(all.vars(formula(M1_unadjusted)),
                     c("alone", "year_normalized", "TUCASEID", "category", "subcategory",
                       "TRCODEP", "TUMONTH", "TUDIARYDAY")))
  
  #M2####
  #Compute the basis once so all prediction grids use the fitted coordinates
  duration_basis <- poly(log(NWSG_data$TUACTDUR24), degree = 3) * sqrt(nrow(NWSG_data))
  NWSG_data[, duration_linear := duration_basis[, 1]]
  NWSG_data[, duration_quadratic := duration_basis[, 2]]
  NWSG_data[, duration_cubic := duration_basis[, 3]]
  rm(duration_basis)
  NWSG_data[, duration_weight := TUACTDUR24 / mean(TUACTDUR24)]
  
  fit_file <- file.path(results_dir, "M2.RData")
  if (file.exists(fit_file)) {
    saved <- load(fit_file)
    stopifnot("M2" %in% saved, inherits(M2, "glmmTMB"), "duration_cubic" %in% all.vars(formula(M2)))
    if (!"time_M2" %in% saved) time_M2 <- c(elapsed = NA_real_)
  } else {
    time_M2 <- system.time({
      M2 <- glmmTMB(
        alone ~ year_normalized * (age_category + HouseholdSize_category + education_category +
                                     hours_worked_category + PESEX + race_category + hispanic_category +
                                     marital_status + income_category) + duration_linear + duration_quadratic + duration_cubic +
          (1 | TUCASEID) +
          (1 + year_normalized || category) +
          (1 + year_normalized || subcategory) +
          (1 + year_normalized || TRCODEP) +
          (1 + year_normalized || TUMONTH) +
          (1 + year_normalized || TUDIARYDAY) +
          (1 + year_normalized || State),
        data = NWSG_data, family = binomial(link = "logit"), weights = duration_weight,
        control = glmmTMBControl(optCtrl = list(iter.max = 2000, eval.max = 2000))
      )
    })
    save(M2, time_M2, file = fit_file)
  }
  
  
  #M1 predictions with fixed and observed activity shares####
  span <- diff(range(years))
  m1 <- rbind(M1_adjusted = atus_trends(M1_adjusted, NWSG_data, span = span),
              M1_unadjusted = atus_trends(M1_unadjusted, NWSG_data, span = span))
  
  #Historical M1 analytical intervals for the fixed year coefficient####
  year_base <- summary(M1_unadjusted)$coefficients$cond["year_normalized", ]
  year_adj <- summary(M1_adjusted)$coefficients$cond["year_normalized", ]
  m1_results <- data.table(
    term = c("M1_unadjusted", "M1_adjusted", "M1_difference"),
    estimate = c(year_base[1], year_adj[1], year_adj[1] - year_base[1]),
    se = c(year_base[2], year_adj[2], sqrt(year_base[2]^2 + year_adj[2]^2))
  )
  #The historical difference assumes zero covariance between fits
  m1_results[, `:=`(lower = estimate - 1.96 * se, upper = estimate + 1.96 * se)]
  m1_results[, ci_width := upper - lower]
  
  #Full-sample fitting times####
  fit_times <- data.table(
    model = c("M1 baseline", "M1 adjusted", "M2"),
    fit_seconds = c(time_weighted_baseline[["elapsed"]], time_weighted[["elapsed"]],
                    time_M2[["elapsed"]])
  )
  fit_times[, fit_minutes := fit_seconds / 60]
  
  
  model <- M2
  frame <- model.frame(model)
  stopifnot(inherits(model, "glmmTMB"), nrow(frame) == nrow(NWSG_data),
            identical(as.character(frame$TUCASEID), as.character(NWSG_data$TUCASEID)),
            identical(as.character(frame$TRCODEP), as.character(NWSG_data$TRCODEP)),
            all(frame$year_normalized == NWSG_data$year_normalized),
            all(model.response(frame) == NWSG_data$alone),
            isTRUE(all.equal(frame$duration_linear, NWSG_data$duration_linear)),
            isTRUE(all.equal(frame$duration_quadratic, NWSG_data$duration_quadratic)),
            isTRUE(all.equal(frame$duration_cubic, NWSG_data$duration_cubic)),
            isTRUE(all.equal(as.numeric(model.weights(frame)), NWSG_data$duration_weight)))
  rm(frame)
  
  #Accuracy on observed diaries####
  diary <- NWSG_data[, .(year = min(years) + year_normalized, TUCASEID, TUACTDUR24, TUFNWGTP, alone)]
  diary[, probability := as.numeric(predict(model, type = "response"))]
  daily <- diary[, .(weight = first(TUFNWGTP), observed = sum(TUACTDUR24 * alone),
                     predicted = sum(TUACTDUR24 * probability)), by = .(year, TUCASEID)]
  
  accuracy <- daily[, {
    obs <- weighted.mean(observed, weight)
    pred <- weighted.mean(predicted, weight)
    mse <- weighted.mean((predicted - observed)^2, weight)
    .(respondents = .N, observed_minutes = obs, predicted_minutes = pred,
      absolute_error_minutes = abs(pred - obs),
      R2 = 1 - mse / weighted.mean((observed - obs)^2, weight), RMSE_minutes = sqrt(mse))
  }, by = year][order(year)]
  
  
  #Observed and predicted time alone by episode-duration decile####
  #Split the survey weight at tied durations across adjacent pooled deciles
  pooled <- diary[, .(weight = sum(TUFNWGTP)), by = TUACTDUR24][order(TUACTDUR24)]
  pooled[, upper := cumsum(weight) / sum(weight)]
  pooled[, lower := shift(upper, fill = 0)]
  shares <- rbindlist(lapply(1:10, function(k) {
    pooled[upper > (k - 1) / 10 & lower < k / 10,
           .(TUACTDUR24, decile = k,
             fraction = (pmin(upper, k / 10) - pmax(lower, (k - 1) / 10)) / (upper - lower))]
  }))
  
  totals <- diary[, .(episodes = sum(TUFNWGTP),
                      minutes = sum(TUFNWGTP * TUACTDUR24),
                      observed = sum(TUFNWGTP * TUACTDUR24 * alone),
                      predicted = sum(TUFNWGTP * TUACTDUR24 * probability)),
                  by = TUACTDUR24]
  parts <- merge(totals, shares, by = "TUACTDUR24", allow.cartesian = TRUE)
  decile_check <- parts[, .(episodes = sum(fraction * episodes),
                            min_minutes = min(TUACTDUR24), max_minutes = max(TUACTDUR24)),
                        by = decile][order(decile)]
  decile_check[, percent := 100 * episodes / sum(episodes)]
  stopifnot(identical(decile_check$decile, 1:10), all(abs(decile_check$percent - 10) < 1e-7))
  
  calibration <- parts[, .(Observed = 100 * sum(fraction * observed) / sum(fraction * minutes),
                           M2 = 100 * sum(fraction * predicted) / sum(fraction * minutes)),
                       by = decile][order(decile)]
  points <- melt(calibration, id.vars = "decile",
                 variable.name = "series", value.name = "percent")
  
  decile_plot <- ggplot(points, aes(decile, percent, colour = series)) +
    geom_line(linewidth = 0.7) +
    geom_point(aes(shape = series), size = 2) +
    scale_x_continuous(breaks = 1:10) +
    scale_colour_manual(values = c(Observed = "black", M2 = "#0072B2")) +
    labs(title = "Time spent alone by episode duration, 2021–2023",
         x = "Episode duration decile (1 = shortest, 10 = longest)",
         y = "Eligible time spent alone (%)", colour = NULL, shape = NULL) +
    theme_classic() +
    theme(legend.position = "bottom")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figures_dir, "observed_predicted_alone_duration_deciles_pooled_2021_2023.pdf"),
         decile_plot, width = 6, height = 4.5)
  print(decile_plot)
  
  #Predicted daily minutes for Tests 1 to 3####
  minute_trends <- function(model, data, span) {
    people <- copy(data[!duplicated(TUCASEID)])
    people[, survey_year := year_normalized]
    people[, `:=`(duration_linear = 0, duration_quadratic = 0, duration_cubic = 0)]
    
    #Use one common eligible-day length and pooled durations within each activity
    episodes <- data[, .(TUCASEID, year_normalized, category, subcategory, TRCODEP,
                         TUACTDUR24, TUFNWGTP, duration_linear, duration_quadratic, duration_cubic)]
    episodes[, eligible := sum(TUACTDUR24), by = TUCASEID]
    diaries <- episodes[!duplicated(TUCASEID)]
    day_minutes <- weighted.mean(diaries$eligible, diaries$TUFNWGTP)
    episodes[, share := TUFNWGTP * TUACTDUR24 / eligible]
    
    mix <- episodes[, .(pooled = sum(share),
                        first = sum(share[year_normalized == 0]),
                        last = sum(share[year_normalized == span])),
                    by = .(category, subcategory, TRCODEP, duration_linear, duration_quadratic, duration_cubic)]
    mix[, within_activity := pooled / sum(pooled), by = TRCODEP]
    mix[, first := within_activity * sum(first), by = TRCODEP]
    mix[, last := within_activity * sum(last), by = TRCODEP]
    
    #Average over respondent heterogeneity using normal quadrature
    sd_person <- unname(attr(VarCorr(model)$cond$TUCASEID, "stddev")[1])
    normal <- statmod::gauss.quad.prob(40, dist = "normal")
    beta <- fixef(model)$cond
    duration_effect <- mix$duration_linear * beta[["duration_linear"]] +
      mix$duration_quadratic * beta[["duration_quadratic"]] +
      mix$duration_cubic * beta[["duration_cubic"]]
    
    endpoint <- function(t, activity_weights) {
      people[, year_normalized := t]
      mix[, year_normalized := t]
      person_eta <- as.numeric(predict(model, newdata = people, type = "link",
                                       re.form = NA)) + re_eta(model, people, "State")
      activity_eta <- duration_effect + re_eta(model, mix, "category") +
        re_eta(model, mix, "subcategory") + re_eta(model, mix, "TRCODEP")
      
      #Interpolate the population mean probability to avoid a person by episode grid
      shifts <- outer(activity_eta, sd_person * normal$nodes, "+")
      knots <- seq(min(shifts), max(shifts) + 0.02, by = 0.02)
      probability <- function(rows) {
        mean_at_knots <- vapply(knots, function(x) {
          weighted.mean(plogis(person_eta[rows] + x), people$TUFNWGTP[rows])
        }, numeric(1))
        interpolated <- approx(knots, mean_at_knots, xout = as.vector(shifts), rule = 2)$y
        drop(matrix(interpolated, nrow = nrow(mix)) %*% normal$weights)
      }
      
      pooled <- probability(seq_len(nrow(people)))
      by_year <- probability(which(people$survey_year == t))
      day_minutes * c(fixed = weighted.mean(pooled, mix$pooled),
                      population = weighted.mean(by_year, mix$pooled),
                      activities = weighted.mean(pooled, activity_weights))
    }
    
    first <- endpoint(0, mix$first)
    last <- endpoint(span, mix$last)
    change <- (last - first) / span
    c(setNames(first, paste0("minutes_2021_", names(first))),
      setNames(last, paste0("minutes_2023_", names(last))),
      minutes_trend_fixed_composition = unname(change["fixed"]),
      minutes_trend_varying_population = unname(change["population"]),
      minutes_population_composition_difference = unname(change["population"] - change["fixed"]),
      minutes_trend_varying_activities = unname(change["activities"]),
      minutes_activity_composition_difference = unname(change["activities"] - change["fixed"]))
  }
  
  #Subgroup estimates####
  atus_postfit <- function(models, data, activities, span) {
    model <- models[["M2"]]
    terms <- c("trend_fixed_composition", "trend_varying_population", "trend_varying_activities",
               "population_composition_difference", "activity_composition_difference",
               paste0("trend_activity_", activities))
    overall <- terms[1:5]
    minute_terms <- c(paste0("minutes_2021_", c("fixed", "population", "activities")),
                      paste0("minutes_2023_", c("fixed", "population", "activities")),
                      paste0("minutes_", terms[c(1, 2, 4, 3, 5)]))
    missing_minutes <- setNames(rep(NA_real_, length(minute_terms)), minute_terms)
    minutes <- missing_minutes
    missing <- setNames(rep(NA_real_, length(terms)), terms)
    estimates <- missing
    
    if (isTRUE(is.finite(model$fit$objective))) {
      #Hold duration fixed in both endpoint predictions, including the population contrast
      grid <- copy(data)
      grid[, `:=`(duration_linear = 0, duration_quadratic = 0, duration_cubic = 0)]
      estimates <- tryCatch(atus_trends(model, grid, activities, span),
                            error = function(e) { warning(conditionMessage(e)); missing })
      minutes <- tryCatch(minute_trends(model, data, span),
                          error = function(e) { warning(conditionMessage(e)); missing_minutes })
    }
    #Empirical pooled within-year SD of observed minutes alone, one row per diary
    daily <- data[, .(observed = sum(TUACTDUR24 * alone), weight = first(TUFNWGTP)),
                  by = .(year_normalized, TUCASEID)]
    sums <- daily[, .(ss = sum(weight * (observed - weighted.mean(observed, weight))^2),
                      denom = sum(weight) - sum(weight^2) / sum(weight)), by = year_normalized]
    pooled_sd <- sqrt(sum(sums$ss) / sum(sums$denom))
    period_change <- unname(minutes[paste0("minutes_", overall)]) * span
    d <- setNames(period_change / pooled_sd, paste0("d_", overall))
    estimates <- c(estimates, minutes, d)
    estimates[!is.finite(estimates)] <- NA_real_
    estimates
  }
  
  
  #Fit M2 adjusted within each group####
  fit_group <- function(data, calls) {
    data <- droplevels(data)
    data[, duration_weight := TUACTDUR24 / mean(TUACTDUR24)]
    messages <- character()
    error <- ""
    
    seconds <- system.time({
      model <- tryCatch(
        withCallingHandlers(eval(calls$M2),
                            warning = function(w) messages <<- c(messages, conditionMessage(w))),
        error = function(e) {
          error <<- conditionMessage(e)
          NULL
        }
      )
    })[["elapsed"]]
    
    details <- data.table(model = "M2", fit_seconds = seconds, error = error,
                          warnings = paste(unique(messages), collapse = " | "),
                          converged = isTRUE(model$fit$convergence == 0L),
                          pdHess = isTRUE(model$sdr$pdHess))
    
    list(models = list(M2 = model), details = details)
  }
  
  #Fit and save the subgroup models####
  s <- 10L
  B <- 2L
  seed <- 2023L
  set.seed(seed)
  units <- unique(NWSG_data[, .(TUCASEID, year_normalized)])
  stopifnot(!anyDuplicated(units$TUCASEID))
  assignments <- lapply(seq_len(B), function(b) {
    ids <- copy(units)
    ids[, group := rg_partition(TUCASEID, s, list(year_normalized))]
    ids
  })
  
  calls <- list(M2 = getCall(M2))
  calls$M2$data <- quote(data)
  partition_files <- file.path(results_dir, sprintf("M2_subgroups_%02d.rds", seq_len(B)))
  input_md5 <- unname(tools::md5sum(input_file))
  
  for (b in seq_len(B)) {
    assignment <- assignments[[b]]
    file <- partition_files[b]
    if (file.exists(file)) {
      saved <- readRDS(file)
      old_call <- saved$calls[["M2"]]
      new_call <- calls$M2
      #Compare model arguments without formula environments or source references
      call_spec <- function(x) {
        x$start <- NULL
        vapply(as.list(x)[-1L], function(arg) paste(deparse(arg), collapse = " "), character(1))
      }
      stopifnot(identical(saved$input_md5, input_md5), identical(call_spec(old_call), call_spec(new_call)),
                identical(as.data.frame(saved$assignment), as.data.frame(assignment)),
                length(saved$fits) == s)
    } else {
      group <- assignment$group[match(NWSG_data$TUCASEID, assignment$TUCASEID)]
      groups <- split(NWSG_data, factor(group, levels = seq_len(s)))
      elapsed <- system.time({
        cl <- parallel::makeCluster(min(n_workers, s))
        fits <- tryCatch({
          parallel::clusterEvalQ(cl, { library(glmmTMB); library(data.table); setDTthreads(1L); NULL })
          parallel::parLapplyLB(cl, groups, fit_group, calls = calls)
        }, finally = parallel::stopCluster(cl))
      })[["elapsed"]]
      saved <- list(assignment = assignment, fits = fits, calls = calls,
                    input_md5 = input_md5, fit_seconds = elapsed)
      tmp <- tempfile(tmpdir = results_dir)
      saveRDS(saved, tmp, compress = FALSE)
      stopifnot(file.rename(tmp, file))
      rm(fits, groups)
    }
    
    rm(saved)
    gc()
  }
  
  #Post-fit inference from the same saved models####
  full_models <- list(M2 = model)
  full <- atus_postfit(full_models, NWSG_data, activities, span)
  main_terms <- c("trend_fixed_composition", "trend_varying_population",
                  "population_composition_difference", "trend_varying_activities",
                  "activity_composition_difference")
  stopifnot(all(is.finite(full[c(main_terms, paste0("d_", main_terms))])))
  rg <- rg_inference(NWSG_data, atus_postfit, unit = "TUCASEID", theta_full = full,
                     full_models = full_models, partition_files = partition_files,
                     keep_groups = FALSE, estimate_args = list(activities = activities, span = span))
  
  results <- data.table(term = names(full), estimate = unname(full), se = NA_real_,
                        lower = NA_real_, upper = NA_real_, p = NA_real_)
  if (NROW(rg)) {
    results[rg, on = "term", `:=`(se = i.se, lower = i.ci_lo, upper = i.ci_hi, p = i.p)]
  }
  
  #Tests 1 to 3 together. Slopes in annual log-odds and d for the full period####
  labels <- c("Change with fixed composition", "Change with year-specific population",
              "Difference in change: population minus fixed", "Change with year-specific activities",
              "Difference in change: activities minus fixed")
  log_odds <- results[match(main_terms, term)]
  standardized <- results[match(paste0("d_", main_terms), term)]
  trends <- data.table(test = c(1L, 2L, 2L, 3L, 3L), term = labels,
                       annual_log_odds = log_odds$estimate, se = log_odds$se,
                       lower = log_odds$lower, upper = log_odds$upper, p = log_odds$p)
  
  
  #Daily minutes. Changes and paired differences over the full two years####
  minute_rows <- results[match(paste0("minutes_", main_terms), term)]
  minutes <- data.table(test = trends$test, term = labels,
                        change_2021_2023 = span * minute_rows$estimate,
                        lower = span * minute_rows$lower, upper = span * minute_rows$upper)
  minutes[, `:=`(d = standardized$estimate, d_lower = standardized$lower, d_upper = standardized$upper)]
  
  #Activity slopes and pointwise 95% confidence intervals####
  selected[, term := paste0("trend_activity_", TRCODEP)]
  plot_data <- results[selected, on = "term"]
  plot_data[, direction := factor(sign(estimate), levels = c(-1, 0, 1),
                                  labels = c("Negative", "Zero", "Positive"))]
  plot_data[, activity := factor(activity, levels = unique(rev(activity)))]
  
  activity_plot <- ggplot(plot_data, aes(activity, estimate, colour = direction)) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, na.rm = TRUE) +
    geom_point(size = 2.6) +
    coord_flip() +
    scale_x_discrete(labels = function(x) vapply(x, function(label)
      paste(strwrap(label, width = 38), collapse = "\n"), character(1))) +
    scale_colour_manual(values = c(Negative = "#0072B2", Zero = "grey40", Positive = "#CC79A7")) +
    labs(title = paste("Time spent alone by activity,", min(years), "–", max(years)),
         x = NULL, y = "Annual change in log-odds of being alone (95% CI)", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank())
  
  if (interactive() || dev.cur() > 1L) print(activity_plot)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(figures_dir, "activity_slopes_2021_2023.pdf"), activity_plot,
         width = 11, height = max(6, 0.45 * nrow(plot_data) + 2), units = "in")
  
  
  
  
  #Compare the same post-stratified estimands across models####
  comparison_terms <- c("trend_fixed_composition", "trend_varying_population",
                        "trend_varying_activities", "population_composition_difference",
                        "activity_composition_difference")
  model_comparison <- data.table(term = comparison_terms,
                                 M1_unadjusted = unname(m1["M1_unadjusted", comparison_terms]),
                                 M1_adjusted = unname(m1["M1_adjusted", comparison_terms]),
                                 M2 = unname(full[comparison_terms]))
  model_comparison[, M1_adjusted_minus_unadjusted := M1_adjusted - M1_unadjusted]
  model_comparison[, M2_minus_M1_adjusted := M2 - M1_adjusted]
  
  #Save the compact post-fit results####
  saved_summary <- list(accuracy = accuracy, trends = trends, minutes = minutes, estimates = results,
                        m1_fixed_year = m1_results, model_comparison = model_comparison,
                        calibration = calibration)
  saveRDS(saved_summary, summary_file)
}

#Final tables####
cat("\nM1 fixed year coefficients and adjusted minus unadjusted (historical analytical intervals)\n")
print(saved_summary$m1_fixed_year, digits = 4)
cat("\nPost-stratified annual log-odds: M1 and M2 point estimates and differences\n")
print(saved_summary$model_comparison, digits = 4)
cat("\nM2 accuracy by year: observed and predicted daily minutes\n")
print(saved_summary$accuracy, digits = 5)
cat("\nTests 1-3: annual slopes and contrasts in log-odds\n")
print(saved_summary$trends, digits = 4)
cat("\nTests 1-3: changes in daily minutes from 2021 to 2023 and Cohen's d (95% CIs)\n")
print(saved_summary$minutes, digits = 4)

cat("\nTest 4: activity-specific annual slopes\n")
print(saved_summary$estimates[startsWith(as.character(term), "trend_activity_")], digits = 4)
cat("\nObserved and predicted time alone by episode-duration decile\n")
print(saved_summary$calibration, digits = 4)

# Open the existing PDFs when viewing the saved results.
if (!run_models) {
  figure_files <- file.path(figures_dir, c(
    "episode_duration_density_2021_2023.pdf",
    "observed_predicted_alone_duration_deciles_pooled_2021_2023.pdf",
    "activity_slopes_2021_2023.pdf"
  ))
  for (figure in figure_files[file.exists(figure_files)]) {
    if (interactive()) {
      if (.Platform$OS.type == "windows") shell.exec(normalizePath(figure))
      else utils::browseURL(figure)
    } else {
      message("Saved figure: ", figure)
    }
  }
}



