## Rationale of the re-reanalysis

The aim of this re-reanalysis of the ATUS (America Times Usage Survey) is to address potential bias in the estimated annual trend and to improve how we quantified its uncertainty, relative to the mixed-effects logistic regression approach we used in one of our original reanalysis (https://osf.io/p8er4).

In that reanalysis, we used `glmmTMB` to fit separate weighted mixed-effects logistic regression models for six measures of social connectedness (being alone, being with friends, being with household family members, being with non-household family members, being with others, and companionship). We adjusted for age, sex, race, Hispanic origin, education, income, marital status, household size, and hours worked. We weighted each diary episode by its duration multiplied by the respondent’s ATUS survey weight and rescaled the resulting weights to have an average of one.

This approach however has two main limitations. First, `glmmTMB` uses the episode weights as multipliers of the log-likelihood contributions. In this way the models treat larger weights as additional statistical information, even though they reflect episode duration and survey representation. Each episode still provides a single observation, so the model-based standard errors may misstate the uncertainty. Second, even with survey weights, our models assumed a common annual trend across sociodemographic groups. If trends differed substantially across groups, the estimated slope could be biased relative to the population-weighted average of the group-specific slopes (see Test 2 in `rg_test.R`). In our original subgroup analyses, however, the estimated differences in slopes were modest and smaller in magnitude than the overall trend, suggesting that any bias from assuming a common slope was likely small.

## Procedure

To address these limitations, we make two changes to our original analysis to improve validity and see if our conclusions hold.

1. **We use post-stratification to estimate the population-average annual slope.** We fit the model without survey weights and include interactions between year and each sociodemographic predictor. We then use the respondent survey weights to average the estimated annual slopes over the sociodemographic composition of the population, pooled across years. For the overall trend, we also average the activity-specific slope contributions according to the share of each respondent’s eligible diary time spent in each activity. The resulting estimate describes the average annual change in log-odds while allowing that change to vary across groups and activities.

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

The standard error is used to construct a 95% *t* confidence interval around the full-sample estimate. For the paired comparison, the baseline slope is subtracted from the adjusted slope within each group. The same variance calculation is then applied to these within-group differences.

## Purpose of the analysis

As in the original reanalysis, the aim is to estimate the annual trend in each measure of social connectedness after accounting for changes in the sociodemographic composition of the population. If a mechanism operating independently of these changes accounts for a substantial part of the observed trend, a substantial trend would be expected to remain after adjustment.

In the original reanalysis, for all outcomes except time spent with friends, there was no evidence of an annual change after adjustment for sociodemographic variables and conditional on the average activity. The estimated confidence intervals indicated that any remaining change was at most small. The attenuation of the trends after adjustment was consistent with the hypothesis that a substantial part of the observed change reflected shifts in sociodemographic composition or in unmeasured factors that covaried with it.
For social isolation, there was also no evidence that the estimated annual trend differed when the activity distribution was allowed to vary across years. This finding indicated that shifts in how respondents allocated their time across activities had, at most, a small impact on the observed trend in social isolation.


## Tests

To assess the robustness of the findings from the original reanalysis, the following tests will be conducted using post-stratified estimates and random-group inference. Annual slopes will be expressed as changes in log-odds per year.

1. **Adjusted annual trend for each outcome.** For each measure of social connectedness, a mixed-effects logistic regression model will be fitted with `glmmTMB`, including the sociodemographic variables used in the original reanalysis and their interactions with year. The overall annual slope will be estimated by averaging the slopes for individual activities according to how much time each respondent spent on them. Respondents from all survey years will contribute to this average according to their survey weights. The adjusted slope and its 95% confidence interval will be reported, and a two-sided *t* test will be used to assess whether the slope differs from zero.
2. **Difference between trends before and after sociodemographic and geographic adjustment.** For each outcome, the adjusted slope from Test 1 will be compared with the slope from a baseline model that accounts for activity, diary month, day of the week, and repeated observations within respondents. This baseline differs slightly from the one used in the original reanalysis. Accounting for these factors in both models helps isolate how much the estimated trend changes after adjustment for sociodemographic characteristics and state of residence. Both slopes will be averaged using the same activity shares and respondent survey weights. The difference (`adjusted - baseline`) will be calculated within each random group to obtain a paired standard error and a 95% confidence interval.
3. **Contribution of changes in the activity distribution to the social isolation trend.** Two slopes from the same social isolation adjusted model will be compared. The first holds the share of eligible diary time spent in each activity fixed at its pooled value. The second allows these shares to vary by year. Both use the same sociodemographic reference distribution. The difference, its 95% confidence interval, and the two-sided *p*-value will be reported. A positive difference would indicate that changes in how respondents allocated their time across activities contributed to a more positive trend in social isolation.
4. **Activity-specific trends in social isolation.** For each of the 20 most frequently reported activities in 2003–2019, the adjusted annual slope will be estimated using the same pooled sociodemographic distribution. The slopes will be plotted with pointwise 95% confidence intervals as in the report of the original reanalysis.

Effect sizes will be interpreted on the log-odds scale, considering both the point estimate and the limits of its 95% confidence interval. To put the magnitude of the trend in context, the estimated change over the study period may also be compared with differences in log-odds associated with sociodemographic characteristics, such as living alone versus living with others or belonging to different age groups.

## Main modelling assumptions

The model allows the annual slope to vary by sociodemographic characteristic and by activity. The main assumption is that omitted interactions involving year and two or more other variables have little influence on the slopes and contrasts of interest. The very small differences in annual slopes observed in the original reanalysis provide some support for this assumption. Including all such interactions would exponentially increase the number of parameters, and some combinations of characteristics and activities would have too few observations for reliable estimation, possibly leading to convergence problems.

The interpretation of the adjusted trends also assumes that the model includes the sociodemographic characteristics most relevant to changes in each measure of social connectedness. This assumption received some indirect support from the original reanalysis, in which the model fitted to the 2003–2019 data showed strong predictive performance on the 2021–2023 data.

Under these assumptions, a small and precisely estimated adjusted trend would weaken explanations that predict a substantial average increase in reported social isolation among people in comparable circumstances performing similar activities. It would, however, remain possible that social processes that covaried with those demographic changes also contributed to the overall increase in isolation.

