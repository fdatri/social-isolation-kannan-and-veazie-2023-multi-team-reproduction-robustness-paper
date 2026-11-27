## Rationale of the re-reanalysis

The aim of this re-reanalysis of the ATUS (America Time Use Survey) is to address potential bias in the estimated annual trend and to improve how we quantified its uncertainty, relative to the mixed-effects logistic regression approach we used in one of our original reanalysis (https://osf.io/p8er4).

In that reanalysis, we used `glmmTMB` to fit separate weighted mixed-effects logistic regression models for six measures of social connectedness (being alone, being with friends, being with household family members, being with non-household family members, being with others, and companionship). We adjusted for age, sex, race, Hispanic origin, education, income, marital status, household size, and hours worked. We weighted each diary episode by its duration multiplied by the respondent’s ATUS survey weight and rescaled the resulting weights to have an average of one.

This approach however has two main limitations. First, `glmmTMB` uses the episode weights as multipliers of the log-likelihood contributions. In this way the models treat larger weights as additional statistical information, even though they reflect episode duration and survey representation. Each episode still provides a single observation, so the model-based standard errors may misstate the uncertainty. Second, even with survey weights, our models assumed a common annual trend across sociodemographic groups. If trends differed substantially across groups, the estimated slope could be biased relative to the population-weighted average of the group-specific slopes (see Test 2 in `rg_test.R`). In our original subgroup analyses, however, the estimated differences in slopes were modest and smaller in magnitude than the overall trend, suggesting that any bias from assuming a common slope was likely small.

To develop the final model for the 2003–2019 re-reanalysis, we first used the 2021–2023 data to perfect model calibration, compare alternative model specifications, and assess the computational feasibility of refitting the model within random subgroups.

## Procedure

To address these limitations, we make two changes to our original analysis to improve validity and see if our conclusions hold.

1. **We use post-stratification to estimate the population-average annual slope.** We fit the model without survey weights and include interactions between year and each sociodemographic predictor. We weight each diary episode by its duration, rescaled to have a mean of one, so that longer episodes contribute to the fit in proportion to the time they represent. We also include linear, quadratic, and cubic terms for log episode duration because, even for the same activity, the probability of being alone can vary with how long the episode lasts. This is especially useful for calibrating the model across episode durations when the aim is to predict daily minutes for each measure of social connectedness. We then use the respondent survey weights to average the estimated annual slopes over the sociodemographic composition of the population, pooled across years. For the overall trend, we also average the activity-specific slope contributions according to the share of each respondent’s eligible diary time spent in each activity. The resulting estimate describes the average annual change in log-odds while allowing that change to vary across groups and activities.

2. **We use random groups to estimate uncertainty.** This approach draws on the method of interpenetrating subsamples introduced by Mahalanobis (1946) and is a well-established method for estimating uncertainty in surveys (see also Wolter, 2007, Chapter 2). We randomly divide respondents into different groups balanced by survey year, keeping all diary episodes from each respondent together. Within each group, we refit the model and repeat the post-stratification. We estimate the variance of the full-sample slope by dividing the variance of the group-specific estimates by the number of valid groups. We repeat the partition twice and average the two variance estimates. This is computationally feasible for our mixed-effects models because each refit uses just portion of the respondents. A full bootstrap would require many refits on samples comparable in size to the original dataset.

## Implementation and testing

### Functions provided by the two scripts

`rg_inference.R` contains functions for random-group inference. Sampling units can be assigned to random groups, and estimates obtained from those groups can be used to quantify uncertainty in a statistic or a paired contrast. `poststrat_helpers.R` contains functions for constructing post-stratification weights and for combining and averaging model-based estimates over a specified reference distribution.

### Application to the ATUS analysis

In the ATUS analysis, respondents are randomly assigned to groups balanced by survey year. All diary episodes from a given respondent are kept together. Within each group, the models are refitted and the same post-stratification procedure used for the full sample is applied.

Annual slopes are estimated on the log-odds scale. For each activity, the annual slope is averaged over the same pooled sociodemographic distribution using respondent survey weights. This estimate includes the fitted year-slope contributions for the activity’s category, subcategory, and activity code. For the overall slope, the activity contributions are first averaged within each respondent using the share of eligible diary time spent in each activity. The resulting estimates are then averaged across respondents using their survey weights.

Uncertainty is estimated separately for each slope or contrast. For each valid partition, the sample variance of the valid group estimates is divided by their number. If `valid_group_estimates` contains the estimates for one partition, its contribution to the variance estimate is

```r
var(valid_group_estimates) / length(valid_group_estimates)
```

If `partition_variances` contains one such value per valid partition, the estimated variance and standard error of the full-sample estimate are

```r
variance_rg <- mean(partition_variances)
se_rg <- sqrt(variance_rg)
```

The standard error is used to construct a 95% *t* confidence interval around the full-sample estimate. For each paired comparison, the fixed-composition trend is subtracted from the corresponding trend with year-specific population composition or activity shares within each random group.



## Purpose of the analysis

As in the original reanalysis, the aim is to estimate the annual trend in each measure of social connectedness after accounting for changes in the sociodemographic composition of the population. If a mechanism operating independently of these changes accounts for a substantial part of the observed trend, a substantial trend would be expected to remain after adjustment.

In the original reanalysis, for all outcomes except time spent with friends, there was no evidence of an annual change after adjustment for sociodemographic variables and conditional on the average activity. The estimated confidence intervals indicated that any remaining change was at most small. The attenuation of the trends after adjustment was consistent with the hypothesis that a substantial part of the observed change reflected shifts in sociodemographic composition or in unmeasured factors that covaried with it.
For social isolation, there was also no evidence that the estimated annual trend differed when the activity distribution was allowed to vary across years. This finding indicated that shifts in how respondents allocated their time across activities had, at most, a small impact on the observed trend in social isolation.


## Tests

To assess the robustness of the findings from the original reanalysis, the following tests will be conducted using post-stratified estimates and random-group inference. Annual slopes will be expressed as changes in log-odds per year.

1. **Adjusted annual trend for each outcome.** For each measure of social connectedness, a mixed-effects logistic regression model will be fitted with `glmmTMB`, including the sociodemographic variables used in the original reanalysis and their interactions with year. The overall annual slope will be estimated by averaging the slopes for individual activities according to the share of eligible diary time each respondent spent on them. Respondents from all survey years will contribute to this average according to their survey weights. This is equivalent to predicting log-odds for each respondent in both the first and last years, keeping their characteristics, state of residence, and activity shares fixed at their observed values. The annual slope is the difference between the last and first years’ survey-weighted mean predictions divided by the number of years elapsed. The adjusted slope and its 95% confidence interval will be reported, and a two-sided *t* test will be used to assess whether the slope differs from zero.

2. **Difference between trends with fixed and year-specific population composition.** For each outcome, the slope from Test 1, which holds the joint distribution of sociodemographic characteristics and state of residence fixed at its pooled value, will be compared with the slope obtained when this distribution is allowed to vary by year. For this second slope, the mean prediction for the first year will use only respondents surveyed in the first year, and the mean prediction for the last year will use only respondents surveyed in the last year. This differs from the original reanalysis, which compared slopes from two separately fitted models. Here, both slopes will be derived from the same adjusted model using respondent survey weights and the same pooled activity shares. The difference (`trend_varying_population - trend_fixed_composition`) will be calculated within each random group to obtain a paired standard error and a 95% confidence interval. A positive difference indicates that changes in population composition are associated with a more positive annual trend.

3. **Contribution of changes in the activity distribution to the trends.** Two slopes from the same adjusted social connectedness measure model will be compared. The first is the slope from Test 1, which keeps each respondent’s observed activity shares fixed in both the first and last years. The second slope will also be calculated from the difference between two predictions for each respondent. For the first prediction, year will be set to the first year and the activity mix to that observed in the first year. For the second prediction, year will be set to the last year and the activity mix to that observed in the last year. Each year’s mix will be calculated by averaging the activity shares of respondents surveyed in that year using their survey weights. For each respondent, sociodemographic characteristics and state of residence will be held fixed. The first prediction will be subtracted from the second on the log-odds scale. These differences will then be averaged across respondents from all survey years using their survey weights and divided by the number of years elapsed to obtain the annual slope.This keeps the sociodemographic and geographic composition fixed while allowing the activity mix to change. The difference between the two slopes (`trend_varying_activities - trend_fixed_composition`), its 95% confidence interval, and the two-sided *p*-value will be reported. A positive difference would indicate that changes in how respondents allocated their time across activities contributed to a more positive trend in the social connectedness measure.

4. **Activity-specific trends in social isolation.** For each of the 20 most frequently reported activities in 2003–2019, the adjusted annual slope will be estimated using the same pooled sociodemographic distribution. The slopes will be plotted with pointwise 95% confidence intervals as in the report of the original reanalysis.

Effect sizes will be interpreted on the log-odds scale, considering both the point estimate and the limits of its 95% confidence interval. To put the magnitude of the trend in context, the estimated change over the study period may also be compared with differences in log-odds associated with sociodemographic characteristics, such as living alone versus living with others or belonging to different age groups.

For each measure of social connectedness, we may also report the estimates and 95% confidence intervals for the first three tests in minutes. This will depend on how well the model reproduces observed daily minutes, in terms of the difference between survey-weighted observed and predicted means in each year (e.g., an absolute error/bias of less than 2 minutes) and the survey-weighted *R squared* across respondents.

For Tests 1-3, we will convert the mean predicted probability at each endpoint to daily minutes by multiplying it by the same survey-weighted mean number of eligible diary minutes, pooled across study years. Eligible diary time is thus held fixed at its pooled value to account for differences in analyzable time across years.

## Main modelling assumptions

The model allows the annual slope to vary by sociodemographic characteristic and by activity. The main assumption is that omitted interactions involving year and two or more other variables have little influence on the slopes and contrasts of interest. The very small differences in annual slopes observed in the original reanalysis provide some support for this assumption. Including all such interactions would exponentially increase the number of parameters, and some combinations of characteristics and activities would have too few observations for reliable estimation, possibly leading to convergence problems.

The interpretation of the adjusted trends also assumes that the model includes the sociodemographic characteristics most relevant to changes in each measure of social connectedness. This assumption received some indirect support from the original reanalysis, in which the model fitted to the 2003–2019 data showed strong predictive performance on the 2021–2023 data.

Under these assumptions, a small and precisely estimated adjusted trend would weaken explanations that predict a substantial average increase in reported social isolation among people in comparable circumstances performing similar activities. It would, however, remain possible that social processes that covaried with those demographic changes also contributed to the overall increase in isolation.

## References

Mahalanobis, P. C. (1946). Recent experiments in statistical sampling in the Indian Statistical Institute. *Journal of the Royal Statistical Society, 109*(4), 325–378.

Wolter, K. M. (2007). The method of random groups. In *Introduction to Variance Estimation* (2nd ed., chap. 2). Springer. <https://doi.org/10.1007/978-0-387-35099-8_2>

