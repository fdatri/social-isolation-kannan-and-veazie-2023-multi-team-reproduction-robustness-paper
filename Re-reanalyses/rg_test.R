#Setup####
Sys.setenv(OMP_NUM_THREADS = "1")
library(glmmTMB)
library(data.table)
library(pbapply)
library(ggplot2)
library(parallel)

#Keep this script in the project folder, next to R/
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
} else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
  if (nzchar(script_path)) setwd(dirname(script_path))
}

source("R/poststrat_helpers.R")
source("R/rg_inference.R")

results_dir <- "function tests results"
dir.create(results_dir, showWarnings = FALSE)
pboptions(type = "timer", min_time = 30)
cores <- detectCores(logical = FALSE)
n_workers <- if (is.na(cores)) 1L else max(1L, cores - 3L)

#Parallelize replications. Each replication fits its random groups sequentially
run_simulation <- function(n_sims, run_replication, objects) {
  cl <- makeCluster(min(n_workers, n_sims))
  on.exit(stopCluster(cl))
  clusterEvalQ(cl, {
    library(glmmTMB)
    library(data.table)
    source("R/poststrat_helpers.R")
    source("R/rg_inference.R")
    NULL
  })
  clusterExport(cl, objects, envir = parent.frame())
  rbindlist(pblapply(seq_len(n_sims), run_replication, cl = cl))
}

#Test 1: sex only as demographic variable. Weights in the likelihood vs post-stratification####
#Sampling is not representative of sex, and the year slope interacts with sex.
#The setup is deliberately simple: each subject contributes six activities, but these
#play no role here, since activities do not enter the DGP. Being alone depends only
#on the subject, on sex, on year, and on the sex by year interaction.
#Also sex distribution is confounded with the year in this example
#(men concentrated in early years, women in late ones)
#Both additive models return a biased point estimate for the year slope conditional
#on the average subject, whether or not weights are used in the likelihood. The model
#that includes the interaction recovers it once the weights are applied afterwards,
#as a post-stratification step.
#The strong bias is mainly due to the strong unmodelled interaction.


n_subj  <- 30000
n_obs   <- 6 #activities per subject
p_pop_f <- 0.7 #population share of females
p_smp_f <- 0.5 #sampling share of females
b0      <- -0.2 #intercept: log-odds of being alone, female, year 0, average subject
b_m     <- 0.30 #male vs female, at year 0
b_yf    <- 0.090 #year slope for females
b_ym    <- -0.030 #year slope for males
sd_u    <- 1.0 #between-subject sd of the random intercept

#estimand: year slope averaged over the population sex composition, for the average subject
truth_t1_year_slope <- p_pop_f * b_yf + (1 - p_pop_f) * b_ym

result_file <- file.path(results_dir, "single_draw_test1.rds")
if (file.exists(result_file)) {
  draw_t1 <- readRDS(result_file)
} else {
  #DGP####
  #Sex is sampled 50/50 while the population is 70/30, so females are under-sampled
  #and the survey weight is the ratio of the two shares. Year is skewed in opposite
  #directions by sex: men are concentrated in the early years, women in the late ones,
  #so sex and year are strongly confounded.
  set.seed(11)
  male <- rbinom(n_subj, 1, 1 - p_smp_f)
  year <- ifelse(male == 1, 8 * rbeta(n_subj, 2, 5), 8 * rbeta(n_subj, 5, 2))
  w    <- ifelse(male == 1, (1 - p_pop_f) / (1 - p_smp_f), p_pop_f / p_smp_f)
  u    <- rnorm(n_subj, 0, sd_u)

  dat_t1 <- data.table(
    subj = factor(rep(1:n_subj, each = n_obs)),
    male = rep(male, each = n_obs),
    year = rep(year, each = n_obs),
    w = rep(w, each = n_obs),
    u = rep(u, each = n_obs)
  )
  dat_t1[, sex := factor(male, 0:1, c("f", "m"))]
  dat_t1[, log_odds := b0 + b_m * male + b_yf * year + (b_ym - b_yf) * male * year + u]
  dat_t1[, alone := rbinom(.N, 1, plogis(log_odds))]
  dat_t1[, w := w * .N / sum(w)]

  #Fit####
  #the same additive model with and without weights, then the model with the interaction
  fit_t1_additive_weighted <- glmmTMB(alone ~ year + sex + (1 | subj), data = dat_t1,
                                      family = binomial("logit"), weights = w)
  fit_t1_additive <- glmmTMB(alone ~ year + sex + (1 | subj), data = dat_t1,
                             family = binomial("logit"))
  fit_t1_interaction <- glmmTMB(alone ~ year * sex + (1 | subj), data = dat_t1,
                                family = binomial("logit"))

  draw_t1 <- list(
    dat_t1 = dat_t1,
    fit_t1_additive_weighted = fit_t1_additive_weighted,
    fit_t1_additive = fit_t1_additive,
    fit_t1_interaction = fit_t1_interaction
  )
  saveRDS(draw_t1, result_file)
}

print(summary(draw_t1$fit_t1_additive_weighted))
print(summary(draw_t1$fit_t1_additive))
print(summary(draw_t1$fit_t1_interaction))

#Population slope####
#one row per subject, carrying the survey weight, is all post-stratification needs
sub_t1 <- subj_grid(draw_t1$dat_t1, c("sex", "year"), "w")

print(truth_t1_year_slope)
print(pop_slope(draw_t1$fit_t1_additive_weighted, ~ year + sex, sub_t1))
print(pop_slope(draw_t1$fit_t1_additive, ~ year + sex, sub_t1))
print(pop_slope(draw_t1$fit_t1_interaction, ~ year * sex, sub_t1))

#Only the model that includes the interaction, combined with post-stratification,
#recovers an estimate of the year slope compatible with the ground truth


#Test 2: sex and age, two omitted interactions####
#Slightly more complicated simulation: sampling is now non-representative on six sex by age cells,
#and the year slope interacts with both sex and age. Year is again skewed by age, so
#the three age groups are observed over different parts of the window. Being alone
#depends on the subject, on sex, on age, on year and on the two interactions with year.
#Both additive models miss the year slope conditional on the average subject. The model
#with the interactions recovers it once the weights are applied afterwards, as a
#post-stratification step

#b0, sd_u and n_obs are inherited from Test 1

n_subj  <- 30000
n_obs   <- 6 #activities per subject
sd_u    <- 1.0 #between-subject sd of the random intercept
b0      <- -0.2 #intercept: log-odds of being alone, female, age group 1, year 0, average subject
b_sex   <- 0.30 #male vs female, at year 0
b_age2  <- 0.15 #age group 2 vs 1, at year 0
b_age3  <- 0.55 #age group 3 vs 1, at year 0
b_y     <- 0.090 #year slope in the reference cell: female, age group 1
d_m     <- -0.120 #shift of the year slope for males
d_a2    <-  0.020 #shift of the year slope for age group 2
d_a3    <-  0.050 #shift of the year slope for age group 3

p_pop <- c(f1 = .7/3, f2 = .7/3, f3 = .7/3, m1 = .1, m2 = .1, m3 = .1) #population cell shares
p_smp <- c(f1 = .10, f2 = .15, f3 = .25, m1 = .25, m2 = .15, m3 = .10) #sampling cell shares


#estimand: year slope averaged over the population sex and age composition, conditional
#on the average subject
truth_t2_year_slope <- b_y + d_m * 0.3 + d_a2 * (1/3) + d_a3 * (1/3)
print(truth_t2_year_slope)

result_file <- file.path(results_dir, "single_draw_test2.rds")
if (file.exists(result_file)) {
  draw_t2 <- readRDS(result_file)
} else {
  #DGP####
  #cells are drawn with the sampling shares and weighted back to the population ones,
  #so the weights span roughly a factor of six; year is skewed by age group
  set.seed(11)
  cell    <- sample(names(p_smp), n_subj, TRUE, prob = p_smp)
  male    <- as.integer(substr(cell, 1, 1) == "m")
  age_num <- as.integer(substr(cell, 2, 2))
  shp     <- c(2, 5, 8)
  year    <- 8 * rbeta(n_subj, shp[age_num], rev(shp)[age_num])
  w       <- p_pop[cell] / p_smp[cell]
  u       <- rnorm(n_subj, 0, sd_u)

  dat_t2 <- data.table(
    subj = factor(rep(1:n_subj, each = n_obs)),
    male = rep(male, each = n_obs),
    age_num = rep(age_num, each = n_obs),
    year = rep(year, each = n_obs),
    w = rep(w, each = n_obs),
    u = rep(u, each = n_obs)
  )
  dat_t2[, sex := factor(male, 0:1, c("f", "m"))]
  dat_t2[, age := factor(age_num)]
  dat_t2[, log_odds := b0 + b_sex * male + b_age2 * (age_num == 2) + b_age3 * (age_num == 3) +
           (b_y + d_m * male + d_a2 * (age_num == 2) + d_a3 * (age_num == 3)) * year + u]
  dat_t2[, alone := rbinom(.N, 1, plogis(log_odds))]
  dat_t2[, w := w * .N / sum(w)]

  #Fit####
  fit_t2_additive_weighted <- glmmTMB(alone ~ year + sex + age + (1 | subj), data = dat_t2,
                                      family = binomial("logit"), weights = w)
  fit_t2_additive <- glmmTMB(alone ~ year + sex + age + (1 | subj), data = dat_t2,
                             family = binomial("logit"))
  fit_t2_interaction <- glmmTMB(alone ~ year * sex + year * age + (1 | subj), data = dat_t2,
                                family = binomial("logit"))

  draw_t2 <- list(
    dat_t2 = dat_t2,
    fit_t2_additive_weighted = fit_t2_additive_weighted,
    fit_t2_additive = fit_t2_additive,
    fit_t2_interaction = fit_t2_interaction
  )
  saveRDS(draw_t2, result_file)
}

print(summary(draw_t2$fit_t2_additive_weighted))
print(summary(draw_t2$fit_t2_additive))
print(summary(draw_t2$fit_t2_interaction))

#Population slope####
sub_t2 <- subj_grid(draw_t2$dat_t2, c("sex", "age", "year"), "w")

print(truth_t2_year_slope)
print(pop_slope(draw_t2$fit_t2_additive_weighted, ~ year + sex + age, sub_t2))
print(pop_slope(draw_t2$fit_t2_additive, ~ year + sex + age, sub_t2))
print(pop_slope(draw_t2$fit_t2_interaction, ~ year * sex + year * age, sub_t2))


#Test 2 simulation: additive weighted model vs post-stratified interaction model####
#Same DGP as Test 2, with the sampling/year
#confounding sharpened (shp) and the interactions strengthened (d_m, d_a2, d_a3)
n_sims_t2 <- 200
n_subj_t2 <- 5000; n_obs_t2 <- 6; sd_u_t2 <- 1.0
b0_t2 <- -0.2; b_sex_t2 <- 0.30; b_age2_t2 <- 0.15; b_age3_t2 <- 0.35
b_y_t2  <- 0.090
d_m_t2  <- -0.200
d_a2_t2 <- -0.040
d_a3_t2 <- -0.100
shp_t2  <- c(1.5, 4, 8)  #was c(2, 3.5, 5), sharper age/year confounding
p_pop_t2 <- c(f1 = .7/3, f2 = .7/3, f3 = .7/3, m1 = .1, m2 = .1, m3 = .1)
p_smp_t2 <- c(f1 = .10, f2 = .15, f3 = .25, m1 = .25, m2 = .15, m3 = .10)
truth_t2_sim <- b_y_t2 + d_m_t2 * 0.3 + d_a2_t2 * (1 / 3) + d_a3_t2 * (1 / 3)

est_fun_t2 <- function(d) {
  f_add_w <- suppressWarnings(glmmTMB(alone ~ year + sex + age + (1 | subj), data = d,
                                      family = binomial("logit"), weights = w))
  f_int   <- suppressWarnings(glmmTMB(alone ~ year * sex + year * age + (1 | subj), data = d,
                                      family = binomial("logit")))
  sub <- subj_grid(d, c("sex", "age", "year"), "w")
  c(additive_weighted = pop_slope(f_add_w, ~ year + sex + age, sub)$est,
    post_stratified_interaction = pop_slope(f_int, ~ year * sex + year * age, sub)$est)
}

run_rep_t2 <- function(k) {
  set.seed(11 + k)
  cell    <- sample(names(p_smp_t2), n_subj_t2, TRUE, prob = p_smp_t2)
  male    <- as.integer(substr(cell, 1, 1) == "m")
  age_num <- as.integer(substr(cell, 2, 2))
  year    <- 8 * rbeta(n_subj_t2, shp_t2[age_num], rev(shp_t2)[age_num])
  w       <- p_pop_t2[cell] / p_smp_t2[cell]
  u       <- rnorm(n_subj_t2, 0, sd_u_t2)
  d <- data.table(subj = factor(rep(1:n_subj_t2, each = n_obs_t2)),
                  male = rep(male, each = n_obs_t2), age_num = rep(age_num, each = n_obs_t2),
                  year = rep(year, each = n_obs_t2), w = rep(w, each = n_obs_t2))
  d[, sex := factor(male, 0:1, c("f", "m"))]
  d[, age := factor(age_num)]
  d[, log_odds := b0_t2 + b_sex_t2 * male + b_age2_t2 * (age_num == 2) + b_age3_t2 * (age_num == 3) +
      (b_y_t2 + d_m_t2 * male + d_a2_t2 * (age_num == 2) + d_a3_t2 * (age_num == 3)) * year +
      rep(u, each = n_obs_t2)]
  d[, alone := rbinom(.N, 1, plogis(log_odds))]
  d[, w := w * .N / sum(w)]

  f_add_w <- suppressWarnings(glmmTMB(alone ~ year + sex + age + (1 | subj), data = d,
                                      family = binomial("logit"), weights = w))
  f_int   <- suppressWarnings(glmmTMB(alone ~ year * sex + year * age + (1 | subj), data = d,
                                      family = binomial("logit")))
  sub <- subj_grid(d, c("sex", "age", "year"), "w")
  ps_add <- pop_slope(f_add_w, ~ year + sex + age, sub)
  ps_int <- pop_slope(f_int, ~ year * sex + year * age, sub)

  #point estimate is the full-data one for both se types, only the variance differs
  theta_full <- c(additive_weighted = ps_add$est, post_stratified_interaction = ps_int$est)
  rg <- rg_inference(d, est_fun_t2, unit = "subj", s = 10, theta_full = theta_full, seed = k)
  rg <- rg[match(names(theta_full), rg$term)]

  rbind(
    data.table(rep = k, method = names(theta_full), se_type = "delta",
               est = theta_full, se = c(ps_add$se, ps_int$se),
               lo = c(ps_add$lo, ps_int$lo), hi = c(ps_add$hi, ps_int$hi)),
    data.table(rep = k, method = names(theta_full), se_type = "rg_subj",
               est = theta_full, se = rg$se, lo = rg$ci_lo, hi = rg$ci_hi)
  )[, cover := truth_t2_sim >= lo & truth_t2_sim <= hi][]
}

result_file <- file.path(results_dir, "sim_test2.rds")
if (file.exists(result_file)) {
  out_t2 <- readRDS(result_file)
} else {
  objects <- c("run_rep_t2", "est_fun_t2", "n_subj_t2", "n_obs_t2", "sd_u_t2", "b0_t2", "b_sex_t2",
               "b_age2_t2", "b_age3_t2", "b_y_t2", "d_m_t2", "d_a2_t2", "d_a3_t2",
               "shp_t2", "p_pop_t2", "p_smp_t2", "truth_t2_sim")
  out_t2 <- run_simulation(n_sims_t2, run_rep_t2, objects)
  saveRDS(out_t2, result_file)
}

print(out_t2[, .(truth = truth_t2_sim, bias = mean(est) - truth_t2_sim,
           rmse = sqrt(mean((est - truth_t2_sim)^2)),
           coverage = mean(cover)), by = .(method, se_type)])


#Test 3: PSU clustering, model-based SE vs random-groups####
#Same DGP as Test 1, plus a design feature the model does not contain: subjects are
#clustered in primary sampling units, and PSUs differ both in their level and in their
#year trend. The model with the interaction is fitted and the estimand recovered
#by post-stratification, so the point estimate should be fine. What fails, compared to
#previous examples, is the SE taken from the model's Hessian (via vcov(fit)), which
#assumes subjects are independent beyond the variables already accounted for by the model.
#Assign whole PSUs to random groups so subjects in the same PSU stay together.
#Compare coverage against the model-based interval over 200 replications.
n_sims <- 200 #replications
n_psu  <- 50 #primary sampling units
n_subj <- 4000 #subjects, nested in PSUs
n_obs  <- 4 #activities per subject
sd_p   <- 0.40 #PSU intercept
sd_py  <- 0.04 #PSU year slope

#estimand: as in Test 1, the year slope averaged over the population sex composition
truth_t3_year_slope <- p_pop_f * b_yf + (1 - p_pop_f) * b_ym
print(truth_t3_year_slope)

#DGP####
#each PSU has its own level, its own year slope and its own stretch of the window,
#so the clustering loads on the trend and not only on the intercept
gen_psu <- function() {
  psu   <- rep(1:n_psu, each = n_subj / n_psu)
  p_int <- rnorm(n_psu, 0, sd_p)
  p_slp <- rnorm(n_psu, 0, sd_py)
  p_yr  <- runif(n_psu, 2, 6)
  male  <- rbinom(n_subj, 1, 1 - p_smp_f)
  year  <- pmin(8, pmax(0, p_yr[psu] + rnorm(n_subj, 0, 1)))
  w     <- ifelse(male == 1, (1 - p_pop_f) / (1 - p_smp_f), p_pop_f / p_smp_f)
  u     <- rnorm(n_subj, 0, sd_u)

  d <- data.table(
    subj = factor(rep(1:n_subj, each = n_obs)),
    psu = factor(rep(psu, each = n_obs)),
    male = rep(male, each = n_obs),
    year = rep(year, each = n_obs),
    w = rep(w, each = n_obs),
    u = rep(u, each = n_obs)
  )
  d[, sex := factor(male, 0:1, c("f", "m"))]
  q <- as.integer(d$psu)
  d[, log_odds := b0 + b_m * male +
      (b_yf + (b_ym - b_yf) * male + p_slp[q]) * year + p_int[q] + u]
  d[, alone := rbinom(.N, 1, plogis(log_odds))]
  d[, w := w * .N / sum(w)][]
}

est_fun <- function(d) {
  f <- suppressWarnings(glmmTMB(alone ~ year * sex + (1 | subj), data = d,
                                family = binomial("logit")))
  sub <- subj_grid(d, c("sex", "year"), "w")
  c(pop_slope = pop_slope(f, ~ year * sex, sub)$est)
}

#Simulation####
run_rep <- function(k) {
  set.seed(11 + k)
  d <- gen_psu()
  f <- suppressWarnings(glmmTMB(alone ~ year * sex + (1 | subj), data = d,
                                family = binomial("logit")))
  ps <- pop_slope(f, ~ year * sex, subj_grid(d, c("sex", "year"), "w"))
  rg <- rg_inference(d, est_fun, unit = "psu", s = 10,
                     theta_full = c(pop_slope = ps$est), seed = k)
  data.table(
    rep = k, method = c("delta", "rg_psu"),
    est = c(ps$est, rg$est_full),
    se = c(ps$se, rg$se),
    lo = c(ps$lo, rg$ci_lo),
    hi = c(ps$hi, rg$ci_hi)
  )[, cover := truth_t3_year_slope >= lo & truth_t3_year_slope <= hi][]
}


result_file <- file.path(results_dir, "sim_test3.rds")
if (file.exists(result_file)) {
  out_psu <- readRDS(result_file)
} else {
  objects <- c("run_rep", "gen_psu", "est_fun", "n_sims", "n_psu", "n_subj", "n_obs",
               "sd_p", "sd_py", "sd_u", "b0", "b_m", "b_yf", "b_ym", "p_pop_f", "p_smp_f",
               "truth_t3_year_slope")
  out_psu <- run_simulation(n_sims, run_rep, objects)
  saveRDS(out_psu, result_file)
}

out_psu[, emp_sd := sd(est), by = method]
print(out_psu[, .(mean_est = mean(est), bias = mean(est) - truth_t3_year_slope, coverage = mean(cover),
            mc_se = sqrt(mean(cover) * (1 - mean(cover)) / .N), mean_se = mean(se),
            emp_sd = emp_sd[1], ratio = mean(se) / emp_sd[1],
            se_cv = sd(se) / mean(se),
            ci_width = mean(hi - lo)), by = method])


#Test 3b: repeated partitions####
#Compare the variability of the estimated SE for B=4 and B=1.
n_sims_b <- 200
run_rep_b <- function(k) {
  set.seed(700 + k)
  d <- gen_psu()
  f <- suppressWarnings(glmmTMB(alone ~ year * sex + (1 | subj), data = d,
                                family = binomial("logit")))
  ps <- pop_slope(f, ~ year * sex, subj_grid(d, c("sex", "year"), "w"))
  rg <- rg_inference(d, est_fun, unit = "psu", s = 10, B = 4,
                     theta_full = c(pop_slope = ps$est), seed = k)
  data.table(
    rep = k, method = c("delta", "rg_psu"),
    est = c(ps$est, rg$est_full),
    se = c(ps$se, rg$se),
    lo = c(ps$lo, rg$ci_lo),
    hi = c(ps$hi, rg$ci_hi)
  )[, cover := truth_t3_year_slope >= lo & truth_t3_year_slope <= hi][]
}

result_file <- file.path(results_dir, "sim_test3b.rds")
if (file.exists(result_file)) {
  out_psu_b <- readRDS(result_file)
} else {
  objects <- c("run_rep_b", "gen_psu", "est_fun", "n_sims_b", "n_psu", "n_subj", "n_obs",
               "sd_p", "sd_py", "sd_u", "b0", "b_m", "b_yf", "b_ym", "p_pop_f", "p_smp_f",
               "truth_t3_year_slope")
  out_psu_b <- run_simulation(n_sims_b, run_rep_b, objects)
  saveRDS(out_psu_b, result_file)
}

out_psu_b[, emp_sd := sd(est), by = method]
print(out_psu_b[, .(mean_est = mean(est), bias = mean(est) - truth_t3_year_slope, coverage = mean(cover),
              mc_se = sqrt(mean(cover) * (1 - mean(cover)) / .N), mean_se = mean(se),
              emp_sd = emp_sd[1], ratio = mean(se) / emp_sd[1],
              se_cv = sd(se) / mean(se),
              ci_width = mean(hi - lo)), by = method])

#se_cv drops from 0.236 (B=1) to 0.155 (B=4)


#Test 4: drifting activity mix, single draw from the DGP####
#Two compositions move over time while conditional behaviour does not:
#the population ages, and time shifts towards activities done alone.
#Only the youngest group has a genuine trend, and the age gap is wider in men.
n_year         <- 10 #years covered by the window
n_per_year     <- 2500 #subjects per year, single draw
n_sims_t4      <- 200 #replications
n_act          <- 20 #activities in the coding scheme
eps_per        <- 6 #activities per subject

b_intercept <- -0.30 #log-odds of being alone, female, age 1, year 0, average activity and subject
b_age       <- c(0, 0.20, 0.45) #oldest group is the most alone
b_sex_male  <- 0.10 #male vs female, at year 0
b_age_male  <- c(0, 0.05, 0.15) #the age gap is wider in men
b_year_age  <- c(0.04, 0, 0) #trend only in the youngest group
sd_subj     <- 0.80 #between-subject sd of the random intercept
sd_alone    <- 0.60 #spread of activity alone propensity
sd_slope    <- 0.10 #spread of activity year slopes

delta_age <- 0.015 #yearly drift of population age shares
kappa_act <- 0.050 #yearly drift of activity mix towards high-alone activities
p_smp_age <- c(1.5, 1.0, 0.6) #the young are over-sampled
p_smp_sex <- c(f = 1.2, m = 0.8) #sampling factors by sex
p_pop_male <- 0.5 #population share of men, constant over the window

#Activities####
#fixed features of the activity system, drawn once
set.seed(11)
acts <- data.table(act_num = 1:n_act, alone = rnorm(n_act, 0, sd_alone))
acts[, z := as.numeric(scale(alone))]
acts[, dur := 45 * exp(0.5 * z)] #activities done alone last longer
acts[, freq := exp(rnorm(n_act, 0, 0.8))][, freq := freq / sum(freq)]
acts[, share := freq * dur][, share := share / sum(share)]
acts[, slope := rnorm(n_act, 0, sd_slope)]
acts[, slope := slope - sum(share * slope)] #average activity has zero trend
#note that share is the share of total time, freq * dur, while the estimator
#averages the within-subject time shares, which are not the same weight. So the
#slopes do not average to exactly zero under the estimator, and the ground truth
#below computes that residual instead of assuming it away. The DGP is unchanged
print(acts)

#activity frequencies drift towards the high-alone ones
freq_year <- lapply(0:(n_year - 1), function(t) {
  p <- acts$freq * exp(kappa_act * acts$z * t); p / sum(p)
})

#A priori ground truth for the year trend marginalised over the activities####
#Obtained by Monte Carlo over the duration mechanism alone, before any data
#are generated and with no model involved.
n_mc <- 200000
act_share_expected <- vapply(seq_len(n_year), function(t) {
  tot <- numeric(n_act)
  for (i in seq_len(n_mc)) {
    a <- sample.int(n_act, eps_per, prob = freq_year[[t]])
    d <- rgamma(eps_per, shape = 4, scale = acts$dur[a] / 4)
    for (j in seq_along(a)) tot[a[j]] <- tot[a[j]] + d[j] / sum(d)
  }
  tot / sum(tot)
}, numeric(n_act))

#Known activity contributions under each set of shares####
span <- n_year - 1
first <- act_share_expected[, 1]
last <- act_share_expected[, n_year]
year_slope_act_term_prior <- sum(rowMeans(act_share_expected) * acts$slope)
observed <- (sum(last * (acts$alone + span * acts$slope)) - sum(first * acts$alone)) / span
year_slope_marg_act_term_prior <- observed - year_slope_act_term_prior

#year slope, conditional on the average activity and the average subject
truth_t4_year_slope <- c(pop = mean(b_year_age), age1 = b_year_age[1],
                         age2 = b_year_age[2], age3 = b_year_age[3])
#same, conditional on the activity specific effects instead and averaged over
#the time spent in each activity. This is what pop_slope returns on the subject
#by activity grid and is the primary Test 4 estimand
truth_t4_year_slope_act <- truth_t4_year_slope + year_slope_act_term_prior
#same, also marginalised over the drift of the activity composition
truth_t4_year_slope_marg_act <- truth_t4_year_slope_act + year_slope_marg_act_term_prior

#age gap 3 - 1, averaged over years and over sex with population shares
age_gap_truth <- function(male) {
  (b_age[3] - b_age[1]) + (b_age_male[3] - b_age_male[1]) * male +
    (b_year_age[3] - b_year_age[1]) * mean(0:(n_year - 1))
}
truth_t4_age_gap <- c(pop = age_gap_truth(p_pop_male),
                      sexf = age_gap_truth(0), sexm = age_gap_truth(1))
#how much the age gap widens in men
truth_t4_age_gap_sex_diff <- b_age_male[3] - b_age_male[1]

print(truth_t4_year_slope)
print(truth_t4_year_slope_act)
print(truth_t4_year_slope_marg_act)
print(truth_t4_age_gap)
print(truth_t4_age_gap_sex_diff)
print(year_slope_act_term_prior)
print(year_slope_marg_act_term_prior)

#DGP####
gen_t4 <- function(n_per_year) {
  year   <- rep(0:(n_year - 1), each = n_per_year)
  n_subj <- length(year)

  #population age shares start equal and drift towards the oldest group
  p_pop_age <- cbind(1/3 - delta_age * year, 1/3, 1/3 + delta_age * year)

  #sampling over-represents the young, with factors constant across years
  p_smp <- sweep(p_pop_age, 2, p_smp_age, `*`)
  p_smp <- p_smp / rowSums(p_smp)

  age_num <- apply(p_smp, 1, function(p) sample.int(3, 1, prob = p))
  male    <- rbinom(n_subj, 1, p_smp_sex["m"] / sum(p_smp_sex))
  re_subj <- rnorm(n_subj, 0, sd_subj)

  #survey weight normalised within year
  w_subj <- p_pop_age[cbind(seq_len(n_subj), age_num)] /
    (p_smp[cbind(seq_len(n_subj), age_num)] *
       ifelse(male == 1, p_smp_sex["m"], p_smp_sex["f"]))
  w_subj <- ave(w_subj, year, FUN = function(x) x / mean(x))

  act_draw <- vapply(seq_len(n_subj), function(i)
    sample.int(n_act, eps_per, prob = freq_year[[year[i] + 1]]), integer(eps_per))

  d <- data.table(
    subj    = factor(rep(seq_len(n_subj), each = eps_per)),
    year    = rep(year, each = eps_per),
    age_num = rep(age_num, each = eps_per),
    male    = rep(male, each = eps_per),
    w_subj  = rep(w_subj, each = eps_per),
    re_subj = rep(re_subj, each = eps_per),
    act_num = as.vector(act_draw)
  )
  d[, age := factor(age_num)]
  d[, sex := factor(male, 0:1, c("f", "m"))]
  d[, act := factor(act_num, levels = 1:n_act)]

  #episode durations, rescaled so that a diary day sums to a plausible budget
  d[, dur := rgamma(.N, shape = 4, scale = acts$dur[act_num] / 4)]
  d[, budget := rep(runif(n_subj, 600, 1080), each = eps_per)]
  d[, dur := dur / sum(dur) * budget[1], by = subj]
  d[, budget := NULL]

  d[, log_odds := b_intercept + b_age[age_num] + b_sex_male * male +
      b_age_male[age_num] * male + acts$alone[act_num] + re_subj +
      (b_year_age[age_num] + acts$slope[act_num]) * year]
  d[, alone := rbinom(.N, 1, plogis(log_odds))]

  #row weight, used only for post-stratification: the subject's survey weight
  #split across the activities in proportion to the time spent in each, so that a
  #subject with a longer diary does not count for more than one with a shorter one
  d[, w_row := w_subj * dur / sum(dur), by = subj]
  d[, w_row := w_row * .N / sum(w_row)][]
}


#Fixed and observed activity shares in the simulation####
t4_mix <- function(fit, grid) {
  shares <- grid[, .(first = sum(w_ps[year == 0]), last = sum(w_ps[year == n_year - 1]),
                     pooled = sum(w_ps)), by = act]
  re <- ranef(fit)$cond$act
  index <- match(as.character(shares$act), rownames(re))

  activity_mix(re[["(Intercept)"]][index], re[["year"]][index],
               shares$first, shares$last, shares$pooled, n_year - 1)
}

result_file <- file.path(results_dir, "single_draw_test4_endpoints.rds")
if (file.exists(result_file)) {
  draw_t4 <- readRDS(result_file)
} else {
  #Single draw####
  #generated once, fitted once, everything (data, fits, results) saved together
  set.seed(21)
  dat_t4 <- gen_t4(n_per_year)

  comp_t4 <- dat_t4[, .(share_young = mean(age_num == 1),
                        share_old = mean(age_num == 3),
                        pop_young = sum(w_subj * (age_num == 1)) / sum(w_subj),
                        pop_middle = sum(w_subj * (age_num == 2)) / sum(w_subj),
                        pop_old = sum(w_subj * (age_num == 3)) / sum(w_subj),
                        time_high_alone = sum(dur * (acts$z[act_num] > 0)) / sum(dur)),
                    by = year][order(year)]

  share_t4 <- dat_t4[, .(time = sum(w_row)), by = .(year, act_num)]
  share_t4[, prop := time / sum(time), by = year]
  share_t4 <- merge(share_t4, acts[, .(act_num, alone)], by = "act_num")

  fit_t4_additive_weighted <- glmmTMB(alone ~ year + age + sex + (1 | subj) + (1 + year || act),
                                      data = dat_t4, family = binomial("logit"), weights = w_row)
  fit_t4_interaction_with_act <- glmmTMB(alone ~ year * age * sex + (1 | subj) + (1 + year || act),
                                         data = dat_t4, family = binomial("logit"))
  fit_t4_interaction_no_act <- glmmTMB(alone ~ year * age * sex + (1 | subj),
                                       data = dat_t4, family = binomial("logit"))

  #one row per subject and activity, weighted by the share of the subject's time
  sub_t4 <- subj_grid(dat_t4, c("age", "sex", "year"), act = "act", durcol = "dur")

  slope_add_w  <- pop_slope(fit_t4_additive_weighted, ~ year + age + sex, sub_t4, by = "age",
                            wcol = "w_ps", re_group = "act")
  #no activity random effect in this fit, so re_group contributes exactly zero
  slope_int_noact <- pop_slope(fit_t4_interaction_no_act, ~ year * age * sex, sub_t4, by = "age",
                               wcol = "w_ps", re_group = "act")

  #Average demographics within age and apply the same activity mix to every age####
  mix <- t4_mix(fit_t4_interaction_with_act, sub_t4)
  person <- pop_slope(fit_t4_interaction_with_act, ~ year * age * sex, sub_t4,
                      by = "age", wcol = "w_ps")
  slope_int_act <- copy(person)[, est := est + mix[["fixed"]]]
  slope_marg_act <- copy(person)[, est := est + mix[["observed"]]]
  for (tab in list(slope_int_act, slope_marg_act)) {
    tab[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se)]
  }
  year_slope_marg_act_term_hat <- mix[["observed"]] - mix[["fixed"]]

  age_gap_with_act <- pop_diff(fit_t4_interaction_with_act, ~ year * age * sex, sub_t4,
                               "age", "3", "1", by = "sex", wcol = "w_ps")
  age_gap_no_act <- pop_diff(fit_t4_interaction_no_act, ~ year * age * sex, sub_t4,
                             "age", "3", "1", by = "sex", wcol = "w_ps")
  age_gap_add_w <- pop_diff(fit_t4_additive_weighted, ~ year + age + sex, sub_t4,
                            "age", "3", "1", by = "sex", wcol = "w_ps")
  age_gap_sex_diff <- row_combo(age_gap_with_act, c(0, -1, 1))

  draw_t4 <- list(
    dat_t4 = dat_t4, comp_t4 = comp_t4, share_t4 = share_t4,
    fit_t4_additive_weighted = fit_t4_additive_weighted,
    fit_t4_interaction_with_act = fit_t4_interaction_with_act,
    fit_t4_interaction_no_act = fit_t4_interaction_no_act,
    slope_add_w = slope_add_w, slope_int_act = slope_int_act,
    slope_int_noact = slope_int_noact, slope_marg_act = slope_marg_act,
    year_slope_marg_act_term_hat = year_slope_marg_act_term_hat,
    age_gap_with_act = age_gap_with_act, age_gap_no_act = age_gap_no_act,
    age_gap_add_w = age_gap_add_w, age_gap_sex_diff = age_gap_sex_diff
  )
  saveRDS(draw_t4, result_file)
}

print(draw_t4$comp_t4)
print(ggplot(melt(draw_t4$comp_t4, id.vars = "year", measure.vars = c("pop_young", "pop_middle", "pop_old")),
       aes(year, value, color = variable)) +
  geom_line(linewidth = 1) +
  labs(title = "Composition drift", x = "year", y = "share") +
  theme_classic(base_size = 14) + theme(panel.grid = element_blank(), axis.ticks = element_line(color = "black")))

print(ggplot(draw_t4$share_t4, aes(x = year, y = prop, group = act_num, color = alone)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  scale_color_gradient2(low = "#2166ac", mid = "grey80", high = "#b2182b", midpoint = 0,
                        name = "alone\npropensity") +
  labs(title = "Activity time share drift", x = "year", y = "share of total time") +
  theme_classic(base_size = 14) + theme(panel.grid = element_blank(), axis.ticks = element_line(color = "black")))

print(summary(draw_t4$fit_t4_additive_weighted))
print(summary(draw_t4$fit_t4_interaction_with_act))
print(summary(draw_t4$fit_t4_interaction_no_act))

#Year slope, conditional on the activity specific effects, time weighted####
print(truth_t4_year_slope_act)
print(draw_t4$slope_add_w)
print(draw_t4$slope_int_act)
print(draw_t4$slope_int_noact)



#Year slope, marginalised over the activities####
print(truth_t4_year_slope_marg_act)
print(draw_t4$slope_marg_act)

#Shift due to the drifting activity mix alone, holding behaviour fixed:
#prior from the true DGP parameters, hat from the fitted model's activity BLUPs
print(c(prior = year_slope_marg_act_term_prior, hat = draw_t4$year_slope_marg_act_term_hat))


#Age gap 3-1####
print(truth_t4_age_gap)
print(draw_t4$age_gap_with_act)
print(draw_t4$age_gap_no_act)
print(draw_t4$age_gap_add_w)


print(truth_t4_age_gap_sex_diff)
print(draw_t4$age_gap_sex_diff)


#Test 4 simulation: analytical SE vs random-groups####
truth_t4_all <- c(
  slope_pop       = truth_t4_year_slope_act[["pop"]],
  slope_age1      = truth_t4_year_slope_act[["age1"]],
  slope_age2      = truth_t4_year_slope_act[["age2"]],
  slope_age3      = truth_t4_year_slope_act[["age3"]],
  slope_marg_pop  = truth_t4_year_slope_marg_act[["pop"]],
  slope_marg_age1 = truth_t4_year_slope_marg_act[["age1"]],
  slope_marg_age2 = truth_t4_year_slope_marg_act[["age2"]],
  slope_marg_age3 = truth_t4_year_slope_marg_act[["age3"]],
  age_gap         = truth_t4_age_gap[["pop"]]
)

fit_t4 <- function(d) {
  d <- copy(d)[, subj := droplevels(subj)]
  suppressWarnings(glmmTMB(alone ~ year * age * sex + (1 | subj) + (1 + year || act),
                           data = d, family = binomial("logit")))
}

t4_estimates <- function(fit, d, se = FALSE) {
  grid <- subj_grid(d, c("age", "sex", "year"), act = "act", durcol = "dur")
  person <- pop_slope(fit, ~ year * age * sex, grid, by = "age", wcol = "w_ps", se = se)
  mix <- t4_mix(fit, grid)
  gap <- pop_diff(fit, ~ year * age * sex, grid, "age", "3", "1", wcol = "w_ps", se = se)

  estimates <- c(setNames(person$est + mix[["fixed"]], paste0("slope_", person$target)),
                 setNames(person$est + mix[["observed"]], paste0("slope_marg_", person$target)),
                 age_gap = gap$est)
  #Analytical SE includes only the fixed-effect uncertainty
  list(est = estimates, se = setNames(c(person$se, person$se, gap$se), names(estimates)))
}

est_fun_t4 <- function(d) {
  t4_estimates(fit_t4(d), d)$est
}

#Simulation for test 4####
run_rep_t4 <- function(k) {
  set.seed(3000 + k)
  d <- gen_t4(2500)
  full <- t4_estimates(fit_t4(d), d, se = TRUE)
  est <- full$est
  se <- full$se

  rg <- rg_inference(d, est_fun_t4, unit = "subj", s = 25,
                     theta_full = est, seed = k)
  rg <- rg[match(names(est), rg$term)]

  rbind(
    data.table(rep = k, method = "delta", target = names(est), est = est, se = se,
               lo = est - 1.96 * se, hi = est + 1.96 * se),
    data.table(rep = k, method = "rg_subj", target = names(est), est = est, se = rg$se,
               lo = rg$ci_lo, hi = rg$ci_hi)
  )[, cover := truth_t4_all[target] >= lo & truth_t4_all[target] <= hi][]
}

out_t4 <- readRDS(file.path(results_dir, "sim_test4.rds"))

out_t4[, emp_sd := sd(est), by = .(method, target)]
print(out_t4[, .(truth = truth_t4_all[target][1], mean_est = mean(est),
           bias = mean(est) - truth_t4_all[target][1], coverage = mean(cover),
           mc_se = sqrt(mean(cover) * (1 - mean(cover)) / .N), mean_se = mean(se),
           emp_sd = emp_sd[1], ratio = mean(se) / emp_sd[1],
           se_cv = sd(se) / mean(se)),
       by = .(target, method)][order(target, method)])

