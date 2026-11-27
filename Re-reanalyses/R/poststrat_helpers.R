#Activity contributions with fixed and observed time shares####
activity_mix <- function(intercept, slope, first, last, pooled, span) {
  stopifnot(span > 0)
  fixed <- weighted.mean(slope, pooled)
  start <- weighted.mean(intercept, first)
  end <- weighted.mean(intercept + span * slope, last)

  c(fixed = fixed, observed = (end - start) / span)
}

#Post-stratified ATUS trends on the log-odds scale####
atus_trends <- function(model, data, activities = character(), span) {
  rhs <- formula(delete.response(terms(lme4::nobars(formula(model)))))
  beta <- fixef(model)$cond
  effects <- ranef(model)$cond

  #Predict the same participants at the first and last year####
  people <- subj_grid(data, unique(c(all.vars(rhs), "State")), wcol = "TUFNWGTP")
  first <- people$year_normalized == 0
  last <- people$year_normalized == span

  person_prediction <- function(year) {
    grid <- copy(people)
    grid[, year_normalized := year]
    X <- model.matrix(rhs, grid)[, names(beta), drop = FALSE]
    eta <- drop(X %*% beta)

    if (!is.null(effects$State)) {
      state <- effects$State[match(grid$State, rownames(effects$State)), , drop = FALSE]
      eta <- eta + state[["(Intercept)"]] + year * state[["year_normalized"]]
    }
    eta
  }

  start <- person_prediction(0)
  end <- person_prediction(span)
  population_fixed <- weighted.mean(end - start, people$w_subj) / span
  population_observed <- (weighted.mean(end[last], people$w_subj[last]) -
                          weighted.mean(start[first], people$w_subj[first])) / span

  #Split each participant's survey weight across the activities they recorded####
  diary <- subj_grid(data, c("category", "subcategory", "year_normalized"),
                     wcol = "TUFNWGTP", act = "TRCODEP", durcol = "TUACTDUR24")

  act <- diary[, .(pooled = sum(w_ps),
                   first = sum(w_ps[year_normalized == 0]),
                   last = sum(w_ps[year_normalized == span])),
               by = .(TRCODEP, subcategory, category)]
  stopifnot(!anyDuplicated(act$TRCODEP))

  act[, `:=`(intercept = 0, slope = 0)]
  for (group in c("category", "subcategory", "TRCODEP")) {
    re <- effects[[group]]
    stopifnot(all(c("(Intercept)", "year_normalized") %in% names(re)))
    index <- match(as.character(act[[group]]), rownames(re))
    act[, intercept := intercept + re[["(Intercept)"]][index]]
    act[, slope := slope + re[["year_normalized"]][index]]
  }

  mix <- with(act, activity_mix(intercept, slope, first, last, pooled, span))
  fixed <- population_fixed + mix[["fixed"]]
  varying_population <- population_observed + mix[["fixed"]]
  varying_activities <- population_fixed + mix[["observed"]]
  selected <- match(activities, as.character(act$TRCODEP))
  activity_names <- if (length(activities)) paste0("trend_activity_", activities)

  c(trend_fixed_composition = fixed,
    trend_varying_population = varying_population,
    trend_varying_activities = varying_activities,
    population_composition_difference = varying_population - fixed,
    activity_composition_difference = varying_activities - fixed,
    setNames(population_fixed + act$slope[selected], activity_names))
}

#Subject-level grid####
subj_grid <- function(dat, cols, wcol = "w_subj", act = NULL, durcol = "dur") {
  if (is.null(act)) {
    return(dat[, c(lapply(.SD, `[`, 1), list(w_subj = get(wcol)[1])), by = subj, .SDcols = cols])
  }
  keep <- setdiff(cols, act)
  g <- dat[, c(lapply(.SD, `[`, 1),
               list(w_subj = get(wcol)[1], .dur = sum(get(durcol)))),
           by = c("subj", act), .SDcols = keep]
  g[, w_ps := w_subj * .dur / sum(.dur), by = subj]
  g[, .dur := NULL]
  g[]
}

#Activity-level random-effect contribution####
re_eta <- function(fit, grid, group = "act") {
  z <- numeric(nrow(grid))
  if (is.null(group) || !group %in% names(grid)) return(z)
  re <- tryCatch(ranef(fit)$cond[[group]], error = function(e) NULL)
  if (is.null(re) || !nrow(re)) return(z)
  i <- match(as.character(grid[[group]]), rownames(re))
  for (v in names(re)) {
    b <- re[[v]][i]
    b[is.na(b)] <- 0
    z <- z + if (v == "(Intercept)") b else b * as.numeric(grid[[v]])
  }
  z
}

#Population-average contrast####
pop_contrast <- function(fit, rhs, sub, set1, set0 = NULL, scale = 1, by = NULL, se = TRUE,
                         wcol = "w_subj", re_group = NULL) {
  #grid_1, grid_0: n rows x same columns as sub, only set1/set0 overridden
  grid_1 <- copy(sub)
  for (nm in names(set1)) set(grid_1, j = nm, value = set1[[nm]])
  grid_0 <- copy(sub)
  if (!is.null(set0)) for (nm in names(set0)) set(grid_0, j = nm, value = set0[[nm]])
  
  #diff_mm: n x p, one row per grid row, one column per fixed-effect coefficient
  diff_mm <- (model.matrix(rhs, grid_1) - model.matrix(rhs, grid_0)) / scale
  b_hat <- fixef(fit)$cond #length p
  #d_re: length n, the same contrast in the re_group random effects, zero if absent
  d_re <- if (is.null(re_group)) numeric(nrow(sub)) else
    (re_eta(fit, grid_1, re_group) - re_eta(fit, grid_0, re_group)) / scale
  #wv: length n, the post-stratification weights used for every average below
  wv <- sub[[wcol]]
  
  if (is.null(by)) {
    #pw: length n, sums to 1
    pw <- wv / sum(wv)
    #contrast: 1 x p, the weighted average of diff_mm's rows
    contrast <- matrix(colSums(pw * diff_mm), nrow = 1,
                       dimnames = list("pop", colnames(diff_mm)))
    #re_add: length 1, the weighted average of the random-effect contrast
    re_add <- sum(pw * d_re)
  } else {
    lv <- levels(sub[[by]])
    #contrast: k x p, one weighted average per level of by
    contrast <- t(vapply(lv, function(a) {
      keep <- sub[[by]] == a
      pw <- wv[keep] / sum(wv[keep])
      colSums(pw * diff_mm[keep, , drop = FALSE])
    }, numeric(ncol(diff_mm))))
    re_add <- vapply(lv, function(a) {
      keep <- sub[[by]] == a
      pw <- wv[keep] / sum(wv[keep])
      sum(pw * d_re[keep])
    }, numeric(1))
    #prepend the unweighted average across levels as the "pop" row: (k+1) x p
    contrast <- rbind(pop = colMeans(contrast), contrast)
    rownames(contrast) <- c("pop", paste0(by, lv))
    re_add <- c(pop = mean(re_add), re_add)
  }
  
  contrast <- contrast[, names(b_hat), drop = FALSE]
  #est: length k+1 (or 1), one number per row of contrast
  est <- as.vector(contrast %*% b_hat) + unname(re_add)
  
  if (!se) {
    out <- data.table(target = rownames(contrast), est = est,
                      se = NA_real_, lo = NA_real_, hi = NA_real_)
    attr(out, "contrast") <- contrast
    return(out[])
  }
  
  #cov_est: (k+1) x (k+1) (or 1x1), the covariance between all rows' estimates
  v_hat <- vcov(fit)$cond #p x p
  cov_est <- contrast %*% v_hat %*% t(contrast)
  se_v <- sqrt(diag(cov_est))
  out <- data.table(target = rownames(contrast), est = est, se = se_v,
                    lo = est - 1.96 * se_v, hi = est + 1.96 * se_v)
  attr(out, "contrast") <- contrast
  attr(out, "vcov") <- cov_est
  out[]
}

#Slope of a continuous covariate####
pop_slope <- function(fit, rhs, sub, var = "year", by = NULL, h = 1e-4, se = TRUE,
                      wcol = "w_subj", re_group = NULL) {
  pop_contrast(fit, rhs, sub, set1 = setNames(list(sub[[var]] + h), var),
               scale = h, by = by, se = se, wcol = wcol, re_group = re_group)
}

#Contrast between two levels of a factor####
#Example: pop_diff(fit, ~ sex, sub, "sex", "m", "f") builds two grids from sub,
#one with sex forced to "m" for every subject and one with sex forced to "f",
#leaving year (and anything else) at each subject's own observed value. The
#model predicts on both grids, takes grid_1 - grid_0 for each subject, and
#averages those differences with the survey weights: the population male vs
#female contrast, holding the composition of the sample fixed.
pop_diff <- function(fit, rhs, sub, var, lev1, lev0, by = NULL, se = TRUE,
                     wcol = "w_subj", re_group = NULL) {
  lv <- levels(sub[[var]])
  pop_contrast(fit, rhs, sub,
               set1 = setNames(list(factor(lev1, levels = lv)), var),
               set0 = setNames(list(factor(lev0, levels = lv)), var), by = by, se = se,
               wcol = wcol, re_group = re_group)
}

#Linear combination of rows####
row_combo <- function(x, weights) {
  cov_est <- attr(x, "vcov")
  est <- sum(weights * x$est)
  se <- sqrt(drop(t(weights) %*% cov_est %*% weights))
  c(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se)
}

#Extra term from the drifting activity composition####
marg_act_term <- function(fit, act_share_first, act_share_last, n_year, include_slope = TRUE) {
  re <- ranef(fit)$cond$act
  j <- as.integer(rownames(re))
  drift <- (sum(act_share_last[j] * re[["(Intercept)"]]) -
              sum(act_share_first[j] * re[["(Intercept)"]])) / (n_year - 1)
  if (!include_slope) return(drift)
  drift + sum(act_share_last[j] * re[["year"]])
}

#Marginal year slope over a drifting fixed-effect composition####
pop_marg_slope <- function(fit, rhs, sub, by, share_first, share_last,
                           var = "year", year_first, year_last) {
  b_hat <- fixef(fit)$cond
  lv <- levels(sub[[by]])
  
  #level_at: weighted-average predicted level per by-level, at a fixed year value
  level_at <- function(year_val) {
    g <- copy(sub); set(g, j = var, value = year_val)
    mm <- model.matrix(rhs, g)[, names(b_hat), drop = FALSE]
    pred <- as.vector(mm %*% b_hat)
    vapply(lv, function(a) {
      keep <- sub[[by]] == a
      weighted.mean(pred[keep], sub$w_subj[keep])
    }, numeric(1))
  }
  
  level_first <- level_at(year_first)
  level_last <- level_at(year_last)
  (sum(share_last * level_last) - sum(share_first * level_first)) / (year_last - year_first)
}

#Post-stratified annual slopes on the logit scale####
activity_slopes <- function(model, data, activities = character()) {
  stopifnot(model$fit$convergence == 0L, isTRUE(model$sdr$pdHess))
  rhs <- formula(delete.response(terms(lme4::nobars(formula(model)))))
  people <- subj_grid(data, all.vars(rhs), wcol = "TUFNWGTP")
  fixed <- pop_slope(model, rhs, people, var = "year_normalized", h = 1, se = FALSE)$est
  grid <- subj_grid(data, c("category", "subcategory"), wcol = "TUFNWGTP",
                    act = "TRCODEP", durcol = "TUACTDUR24")
  codes <- unique(grid[, .(TRCODEP, category, subcategory)])
  stopifnot(!anyDuplicated(codes$TRCODEP))
  codes[, slope := as.numeric(fixed)]
  for (group in c("category", "subcategory", "TRCODEP")) {
    codes[, year_normalized := 0]
    zero <- re_eta(model, codes, group)
    codes[, year_normalized := 1]
    codes[, slope := slope + re_eta(model, codes, group) - zero]
  }
  overall <- weighted.mean(codes$slope[match(grid$TRCODEP, codes$TRCODEP)], grid$w_ps)
  if (!length(activities)) return(c(overall = overall))
  by_activity <- codes$slope[match(activities, as.character(codes$TRCODEP))]
  c(overall = overall, setNames(by_activity, paste0("activity_", activities)))
}

#ATUS activity labels from the 2003-2023 multi-year coding lexicon####
atus_activity_labels <- c(
  "010101" = "Sleeping",
  "010102" = "Sleeplessness",
  "010199" = "Sleeping, n.e.c.*",
  "010201" = "Washing, dressing and grooming oneself",
  "010299" = "Grooming, n.e.c.*",
  "010301" = "Health-related self care",
  "010399" = "Self care, n.e.c.*",
  "010401" = "Personal/Private activities",
  "010499" = "Personal activities, n.e.c.*",
  "010501" = "Personal emergencies",
  "010599" = "Personal care emergencies, n.e.c.*",
  "019999" = "Personal care, n.e.c.*",
  "020101" = "Interior cleaning",
  "020102" = "Laundry",
  "020103" = "Sewing, repairing, & maintaining textiles",
  "020104" = "Storing interior hh items, inc. food",
  "020199" = "Housework, n.e.c.*",
  "020201" = "Food and drink preparation",
  "020202" = "Food presentation",
  "020203" = "Kitchen and food clean-up",
  "020299" = "Food & drink prep, presentation, & clean-up, n.e.c.*",
  "020301" = "Interior arrangement, decoration, & repairs",
  "020302" = "Building and repairing furniture",
  "020303" = "Heating and cooling",
  "020399" = "Interior maintenance, repair, & decoration, n.e.c.*",
  "020401" = "Exterior cleaning",
  "020402" = "Exterior repair, improvements, & decoration",
  "020499" = "Exterior maintenance, repair & decoration, n.e.c.*",
  "020501" = "Lawn, garden, and houseplant care",
  "020502" = "Ponds, pools, and hot tubs",
  "020599" = "Lawn and garden, n.e.c.*",
  "020681" = "Care for animals and pets (not veterinary care)",
  "020699" = "Pet and animal care, n.e.c.*",
  "020701" = "Vehicle repair and maintenance (by self)",
  "020799" = "Vehicles, n.e.c.*",
  "020801" = "Appliance, tool, and toy set-up, repair, & maintenance (by self)",
  "020899" = "Appliances and tools, n.e.c.*",
  "020901" = "Financial management",
  "020902" = "Household & personal organization and planning",
  "020903" = "HH & personal mail & messages (except e-mail)",
  "020904" = "HH & personal e-mail and messages",
  "020905" = "Home security",
  "020999" = "Household management, n.e.c.*",
  "029999" = "Household activities, n.e.c.*",
  "030101" = "Physical care for hh children",
  "030102" = "Reading to/with hh children",
  "030103" = "Playing with hh children, not sports",
  "030104" = "Arts and crafts with hh children",
  "030105" = "Playing sports with hh children",
  "030186" = "Talking with/listening to hh children",
  "030108" = "Organization & planning for hh children",
  "030109" = "Looking after hh children (as a primary activity)",
  "030110" = "Attending hh children's events",
  "030111" = "Waiting for/with hh children",
  "030112" = "Picking up/dropping off hh children",
  "030199" = "Caring for & helping hh children, n.e.c.*",
  "030201" = "Homework (hh children)",
  "030202" = "Meetings and school conferences (hh children)",
  "030203" = "Home schooling of hh children",
  "030204" = "Waiting associated with hh children's education",
  "030299" = "Activities related to hh child's education, n.e.c.*",
  "030301" = "Providing medical care to hh children",
  "030302" = "Obtaining medical care for hh children",
  "030303" = "Waiting associated with hh children's health",
  "030399" = "Activities related to hh child's health, n.e.c.*",
  "030401" = "Physical care for hh adults",
  "030402" = "Looking after hh adult (as a primary activity)",
  "030403" = "Providing medical care to hh adult",
  "030404" = "Obtaining medical and care services for hh adult",
  "030405" = "Waiting associated with caring for household adults",
  "030499" = "Caring for household adults, n.e.c.*",
  "030501" = "Helping hh adults",
  "030502" = "Organization & planning for hh adults",
  "030503" = "Picking up/dropping off hh adult",
  "030504" = "Waiting associated with helping hh adults",
  "030599" = "Helping household adults, n.e.c.*",
  "039999" = "Caring for & helping hh members, n.e.c.*",
  "040101" = "Physical care for nonhh children",
  "040102" = "Reading to/with nonhh children",
  "040103" = "Playing with nonhh children, not sports",
  "040104" = "Arts and crafts with nonhh children",
  "040105" = "Playing sports with nonhh children",
  "040186" = "Talking with/listening to nonhh children",
  "040108" = "Organization & planning for nonhh children",
  "040109" = "Looking after nonhh children (as primary activity)",
  "040110" = "Attending nonhh children's events",
  "040111" = "Waiting for/with nonhh children",
  "040112" = "Dropping off/picking up nonhh children",
  "040199" = "Caring for and helping nonhh children, n.e.c.*",
  "040201" = "Homework (nonhh children)",
  "040202" = "Meetings and school conferences (nonhh children)",
  "040203" = "Home schooling of nonhh children",
  "040204" = "Waiting associated with nonhh children's education",
  "040299" = "Activities related to nonhh child's educ., n.e.c.*",
  "040301" = "Providing medical care to nonhh children",
  "040302" = "Obtaining medical care for nonhh children",
  "040303" = "Waiting associated with nonhh children's health",
  "040399" = "Activities related to nonhh child's health, n.e.c.*",
  "040401" = "Physical care for nonhh adults",
  "040402" = "Looking after nonhh adult (as a primary activity)",
  "040403" = "Providing medical care to nonhh adult",
  "040404" = "Obtaining medical and care services for nonhh adult",
  "040405" = "Waiting associated with caring for nonhh adults",
  "040499" = "Caring for nonhh adults, n.e.c.*",
  "040501" = "Housework, cooking, & shopping assistance for nonhh adults",
  "040502" = "House & lawn maintenance & repair assistance for nonhh adults",
  "040503" = "Animal & pet care assistance for nonhh adults",
  "040504" = "Vehicle & appliance maintenance/repair assistance for nonhh adults",
  "040505" = "Financial management assistance for nonhh adults",
  "040506" = "Household management & paperwork assistance for nonhh adults",
  "040507" = "Picking up/dropping off nonhh adult",
  "040508" = "Waiting associated with helping nonhh adults",
  "040599" = "Helping nonhh adults, n.e.c.*",
  "049999" = "Caring for & helping nonhh members, n.e.c.*",
  "050101" = "Work, main job",
  "050102" = "Work, other job(s)",
  "050103" = "Security procedures related to work",
  "050189" = "Working, n.e.c.*",
  "050201" = "Socializing, relaxing, and leisure as part of job",
  "050202" = "Eating and drinking as part of job",
  "050203" = "Sports and exercise as part of job",
  "050204" = "Security procedures as part of job",
  "050289" = "Work-related activities, n.e.c.*",
  "050301" = "Income-generating hobbies, crafts, and food",
  "050302" = "Income-generating performances",
  "050303" = "Income-generating services",
  "050304" = "Income-generating rental property activities",
  "050389" = "Other income-generating activities, n.e.c.*",
  "050481" = "Job search activities",
  "050403" = "Job interviewing",
  "050404" = "Waiting associated with job search or interview",
  "050405" = "Security procedures rel. to job search/interviewing",
  "050499" = "Job search and Interviewing, n.e.c.*",
  "059999" = "Work and work-related activities, n.e.c.*",
  "060101" = "Taking class for degree, certification, or licensure",
  "060102" = "Taking class for personal interest",
  "060103" = "Waiting associated with taking classes",
  "060104" = "Security procedures rel. to taking classes",
  "060199" = "Taking class, n.e.c.*",
  "060201" = "Extracurricular club activities",
  "060202" = "Extracurricular music & performance activities",
  "060203" = "Extracurricular student government activities",
  "060289" = "Education-related extracurricular activities, n.e.c.*",
  "060301" = "Research/homework for class for degree, certification, or licensure",
  "060302" = "Research/homework for class for pers. interest",
  "060303" = "Waiting associated with research/homework",
  "060399" = "Research/homework n.e.c.*",
  "060401" = "Administrative activities: class for degree, certification, or licensure",
  "060402" = "Administrative activities: class for personal interest",
  "060403" = "Waiting associated w/admin. activities (education)",
  "060499" = "Administrative for education, n.e.c.*",
  "069999" = "Education, n.e.c.*",
  "070101" = "Grocery shopping",
  "070102" = "Purchasing gas",
  "070103" = "Purchasing food (not groceries)",
  "070104" = "Shopping, except groceries, food and gas",
  "070105" = "Waiting associated with shopping",
  "070199" = "Shopping, n.e.c.*",
  "070201" = "Comparison shopping",
  "070299" = "Researching purchases, n.e.c.*",
  "070301" = "Security procedures rel. to consumer purchases",
  "070399" = "Security procedures rel. to consumer purchases, n.e.c.*",
  "079999" = "Consumer purchases, n.e.c.*",
  "080101" = "Using paid childcare services",
  "080102" = "Waiting associated w/purchasing childcare svcs",
  "080199" = "Using paid childcare services, n.e.c.*",
  "080201" = "Banking",
  "080202" = "Using other financial services",
  "080203" = "Waiting associated w/banking/financial services",
  "080299" = "Using financial services and banking, n.e.c.*",
  "080301" = "Using legal services",
  "080302" = "Waiting associated with legal services",
  "080399" = "Using legal services, n.e.c.*",
  "080401" = "Using health and care services outside the home",
  "080402" = "Using in-home health and care services",
  "080403" = "Waiting associated with medical services",
  "080499" = "Using medical services, n.e.c.*",
  "080501" = "Using personal care services",
  "080502" = "Waiting associated w/personal care services",
  "080599" = "Using personal care services, n.e.c.*",
  "080601" = "Activities rel. to purchasing/selling real estate",
  "080602" = "Waiting associated w/purchasing/selling real estate",
  "080699" = "Using real estate services, n.e.c.*",
  "080701" = "Using veterinary services",
  "080702" = "Waiting associated with veterinary services",
  "080799" = "Using veterinary services, n.e.c.*",
  "080801" = "Security procedures rel. to professional/personal svcs.",
  "080899" = "Security procedures rel. to professional/personal svcs n.e.c.*",
  "089999" = "Professional and personal services, n.e.c.*",
  "090101" = "Using interior cleaning services",
  "090102" = "Using meal preparation services",
  "090103" = "Using clothing repair and cleaning services",
  "090104" = "Waiting associated with using household services",
  "090199" = "Using household services, n.e.c.*",
  "090201" = "Using home maint/repair/décor/construction svcs",
  "090202" = "Waiting associated w/ home main/repair/décor/constr",
  "090299" = "Using home maint/repair/décor/constr services, n.e.c.*",
  "090301" = "Using pet services",
  "090302" = "Waiting associated with pet services",
  "090399" = "Using pet services, n.e.c.*",
  "090401" = "Using lawn and garden services",
  "090402" = "Waiting associated with using lawn & garden services",
  "090499" = "Using lawn and garden services, n.e.c.*",
  "090501" = "Using vehicle maintenance or repair services",
  "090502" = "Waiting associated with vehicle main. or repair svcs",
  "090599" = "Using vehicle maint. & repair svcs, n.e.c.*",
  "099999" = "Using household services, n.e.c.*",
  "100101" = "Using police and fire services",
  "100102" = "Using social services",
  "100103" = "Obtaining licenses & paying fines, fees, taxes",
  "100199" = "Using government services, n.e.c.*",
  "100201" = "Civic obligations & participation",
  "100299" = "Civic obligations & participation, n.e.c.*",
  "100381" = "Waiting associated with using government services",
  "100383" = "Waiting associated w/civic obligations & participation",
  "100399" = "Waiting assoc. w/govt svcs or civic obligations, n.e.c.*",
  "100401" = "Security procedures rel. to govt svcs/civic obligations",
  "100499" = "Security procedures rel. to govt svcs/civic obligations, n.e.c.*",
  "109999" = "Government services, n.e.c.*",
  "110101" = "Eating and drinking",
  "110199" = "Eating and drinking, n.e.c.*",
  "110281" = "Waiting associated w/eating & drinking",
  "110289" = "Waiting associated with eating & drinking, n.e.c.*",
  "119999" = "Eating and drinking, n.e.c.*",
  "120101" = "Socializing and communicating with others",
  "120199" = "Socializing and communicating, n.e.c.*",
  "120201" = "Attending or hosting parties/receptions/ceremonies",
  "120202" = "Attending meetings for personal interest (not volunteering)",
  "120299" = "Attending/hosting social events, n.e.c.*",
  "120301" = "Relaxing, thinking",
  "120302" = "Tobacco and drug use",
  "120303" = "Television and movies (not religious)",
  "120304" = "Television (religious)",
  "120305" = "Listening to the radio",
  "120306" = "Listening to/playing music (not radio)",
  "120307" = "Playing games",
  "120308" = "Computer use for leisure (exc. Games)",
  "120309" = "Arts and crafts as a hobby",
  "120310" = "Collecting as a hobby",
  "120311" = "Hobbies, except arts & crafts and collecting",
  "120312" = "Reading for personal interest",
  "120313" = "Writing for personal interest",
  "120399" = "Relaxing and leisure, n.e.c.*",
  "120401" = "Attending performing arts",
  "120402" = "Attending museums",
  "120403" = "Attending movies/film",
  "120404" = "Attending gambling establishments",
  "120405" = "Security procedures rel. to arts & entertainment",
  "120499" = "Arts and entertainment, n.e.c.*",
  "120501" = "Waiting assoc. w/socializing & communicating",
  "120502" = "Waiting assoc. w/attending/hosting social events",
  "120503" = "Waiting associated with relaxing/leisure",
  "120504" = "Waiting associated with arts & entertainment",
  "120599" = "Waiting associated with socializing, n.e.c.*",
  "129999" = "Socializing, relaxing, and leisure, n.e.c.*",
  "130101" = "Doing aerobics",
  "130102" = "Playing baseball",
  "130103" = "Playing basketball",
  "130104" = "Biking",
  "130105" = "Playing billiards",
  "130106" = "Boating",
  "130107" = "Bowling",
  "130108" = "Climbing, spelunking, caving",
  "130109" = "Dancing",
  "130110" = "Participating in equestrian sports",
  "130111" = "Fencing",
  "130112" = "Fishing",
  "130113" = "Playing football",
  "130114" = "Golfing",
  "130115" = "Doing gymnastics",
  "130116" = "Hiking",
  "130117" = "Playing hockey",
  "130118" = "Hunting",
  "130119" = "Participating in martial arts",
  "130120" = "Playing racquet sports",
  "130121" = "Participating in rodeo competitions",
  "130122" = "Rollerblading",
  "130123" = "Playing rugby",
  "130124" = "Running",
  "130125" = "Skiing, ice skating, snowboarding",
  "130126" = "Playing soccer",
  "130127" = "Softball",
  "130128" = "Using cardiovascular equipment",
  "130129" = "Vehicle touring/racing",
  "130130" = "Playing volleyball",
  "130131" = "Walking",
  "130132" = "Participating in water sports",
  "130133" = "Weightlifting/strength training",
  "130134" = "Working out, unspecified",
  "130135" = "Wrestling",
  "130136" = "Doing yoga",
  "130199" = "Playing sports n.e.c.*",
  "130201" = "Watching aerobics",
  "130202" = "Watching baseball",
  "130203" = "Watching basketball",
  "130204" = "Watching biking",
  "130205" = "Watching billiards",
  "130206" = "Watching boating",
  "130207" = "Watching bowling",
  "130208" = "Watching climbing, spelunking, caving",
  "130209" = "Watching dancing",
  "130210" = "Watching equestrian sports",
  "130211" = "Watching fencing",
  "130212" = "Watching fishing",
  "130213" = "Watching football",
  "130214" = "Watching golfing",
  "130215" = "Watching gymnastics",
  "130216" = "Watching hockey",
  "130217" = "Watching martial arts",
  "130218" = "Watching racquet sports",
  "130219" = "Watching rodeo competitions",
  "130220" = "Watching rollerblading",
  "130221" = "Watching rugby",
  "130222" = "Watching running",
  "130223" = "Watching skiing, ice skating, snowboarding",
  "130224" = "Watching soccer",
  "130225" = "Watching softball",
  "130226" = "Watching vehicle touring/racing",
  "130227" = "Watching volleyball",
  "130228" = "Watching walking",
  "130229" = "Watching water sports",
  "130230" = "Watching weightlifting/strength training",
  "130231" = "Watching people working out, unspecified",
  "130232" = "Watching wrestling",
  "130299" = "Attending sporting events, n.e.c.*",
  "130301" = "Waiting related to playing sports or exercising",
  "130302" = "Waiting related to attending sporting events",
  "130399" = "Waiting associated with sports, exercise, & recreation, n.e.c.*",
  "130401" = "Security related to playing sports or exercising",
  "130402" = "Security related to attending sporting events",
  "130499" = "Security related to sports, exercise, & recreation, n.e.c.*",
  "139999" = "Sports, exercise, & recreation, n.e.c.*",
  "140101" = "Attending religious services",
  "140102" = "Participation in religious practices",
  "140103" = "Waiting associated w/religious & spiritual activities",
  "140104" = "Security procedures rel. to religious & spiritual activities",
  "140105" = "Religious education activities",
  "149999" = "Religious and spiritual activities, n.e.c.*",
  "150101" = "Computer use",
  "150102" = "Organizing and preparing",
  "150103" = "Reading",
  "150104" = "Telephone calls (except hotline counseling)",
  "150105" = "Writing",
  "150106" = "Fundraising",
  "150199" = "Administrative & support activities, n.e.c.*",
  "150201" = "Food preparation, presentation, clean-up",
  "150202" = "Collecting & delivering clothing & other goods",
  "150203" = "Providing care",
  "150204" = "Teaching, leading, counseling, mentoring",
  "150299" = "Social service & care activities, n.e.c.*",
  "150301" = "Building houses, wildlife sites, & other structures",
  "150302" = "Indoor & outdoor maintenance, repair, & clean-up",
  "150399" = "Indoor & outdoor maintenance, building & clean-up activities, n.e.c.*",
  "150401" = "Performing",
  "150402" = "Serving at volunteer events & cultural activities",
  "150499" = "Participating in performance & cultural activities, n.e.c.*",
  "150501" = "Attending meetings, conferences, & training",
  "150599" = "Attending meetings, conferences, & training, n.e.c.*",
  "150601" = "Public health activities",
  "150602" = "Public safety activities",
  "150699" = "Public health & safety activities, n.e.c.*",
  "159989" = "Volunteer activities, n.e.c.*",
  "160101" = "Telephone calls to/from family members",
  "160102" = "Telephone calls to/from friends, neighbors, or acquaintances",
  "160103" = "Telephone calls to/from education services providers",
  "160104" = "Telephone calls to/from salespeople",
  "160105" = "Telephone calls to/from professional or personal care svcs providers",
  "160106" = "Telephone calls to/from household services providers",
  "160107" = "Telephone calls to/from paid child or adult care providers",
  "160108" = "Telephone calls to/from government officials",
  "169989" = "Telephone calls, n.e.c.*",
  "180101" = "Travel related to personal care",
  "180199" = "Travel related to personal care, n.e.c.*",
  "180280" = "Travel related to household activities",
  "180381" = "Travel related to caring for and helping hh children",
  "180382" = "Travel related to caring for and helping hh adults",
  "180399" = "Travel rel. to caring for & helping hh members, n.e.c.*",
  "180481" = "Travel related to caring for and helping nonhh children",
  "180482" = "Travel related to caring for and helping nonhh adults",
  "180499" = "Travel rel. to caring for & helping nonhh members, n.e.c.*",
  "180501" = "Travel related to working",
  "180502" = "Travel related to work-related activities",
  "180589" = "Travel related to work, n.e.c.*",
  "180601" = "Travel related to taking class",
  "180682" = "Travel related to education (except taking class)",
  "180699" = "Travel related to education, n.e.c.*",
  "180701" = "Travel related to grocery shopping",
  "180782" = "Travel related to shopping (except grocery shopping)",
  "180801" = "Travel related to using childcare services",
  "180802" = "Travel related to using financial services and banking",
  "180803" = "Travel related to using legal services",
  "180804" = "Travel related to using medical services",
  "180805" = "Travel related to using personal care services",
  "180806" = "Travel related to using real estate services",
  "180807" = "Travel related to using veterinary services",
  "180899" = "Travel rel. to using prof. & personal care services, n.e.c.*",
  "180901" = "Travel related to using household services",
  "180902" = "Travel related to using home main./repair/décor./construction svcs",
  "180903" = "Travel related to using pet services (not vet)",
  "180904" = "Travel related to using lawn and garden services",
  "180905" = "Travel related to using vehicle maintenance & repair services",
  "180999" = "Travel related to using household services, n.e.c.*",
  "181081" = "Travel related to using government services",
  "181002" = "Travel related to civic obligations & participation",
  "181099" = "Travel rel. to govt svcs & civic obligations, n.e.c.*",
  "181101" = "Travel related to eating and drinking",
  "181199" = "Travel related to eating and drinking, n.e.c.*",
  "181201" = "Travel related to socializing and communicating",
  "181202" = "Travel related to attending or hosting social events",
  "181283" = "Travel related to relaxing and leisure",
  "181204" = "Travel related to arts and entertainment",
  "181299" = "Travel rel. to socializing, relaxing, & leisure, n.e.c.*",
  "181301" = "Travel related to participating in sports/exercise/recreation",
  "181302" = "Travel related to attending sporting/recreational events",
  "181399" = "Travel related to sports, exercise, & recreation, n.e.c.*",
  "181401" = "Travel related to religious/spiritual practices",
  "181499" = "Travel rel. to religious/spiritual activities, n.e.c.*",
  "181501" = "Travel related to volunteering",
  "181599" = "Travel related to volunteer activities, n.e.c.*",
  "181601" = "Travel related to phone calls",
  "181699" = "Travel rel. to phone calls, n.e.c.*",
  "181801" = "Security procedures related to traveling",
  "181899" = "Security procedures related to traveling, n.e.c.*",
  "189999" = "Traveling, n.e.c.*",
  "500101" = "Insufficient detail in verbatim",
  "500103" = "Missing travel or destination",
  "500104" = "Recorded simultaneous activities incorrectly code was used in 2003-2012 only",
  "500105" = "Respondent refused to provide information/\"none of your business\"",
  "500106" = "Gap/can't remember",
  "500107" = "Unable to code activity at 1st tier",
  "509989" = "Data codes, n.e.c.*"
)
