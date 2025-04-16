# Clinical parameter analysis 
# 1. Paired Wilcoxon test for comparing two time points within the same group
wilcox.test(group_time1, group_time2, paired = TRUE, alternative = "two.sided")

# 2. Unpaired Wilcoxon (Mann–Whitney U) test for comparing two groups at the same time point
wilcox.test(group1_data, group2_data, paired = FALSE, alternative = "two.sided")

###### Statistical analysis for continuous CGM data, alpha diversity, and ordination coordinates (PCoA/PCA/aPCoA)

## 1. Compare each day of experimental phase with baseline (normal diet) using paired Wilcoxon test
wilcox.test(day_data, baseline_mean, paired = TRUE, alternative = alternative_set)

## 2. Within-group comparison between two time points (paired Wilcoxon test)
wilcox.test(timepoint1_mean, timepoint2_mean, paired = TRUE, alternative = alternative_set)

## 3. Between-group comparison at the same time point using ANCOVA (baseline as covariate)
# Main-effects ANCOVA model: tests group differences adjusting for baseline
formula_main <- paste(response_variable, " ~ Nor_mean + group")
fit_main <- aov(as.formula(formula_main), data = data_set)
summary(fit_main)

# Interaction ANCOVA model: tests whether the effect of baseline differs by group
formula_interaction <- paste(response_variable, " ~ Nor_mean * group")
fit_interaction <- aov(as.formula(formula_interaction), data = data_set)
summary(fit_interaction)

## 4. Between-group comparison during the same experimental phase using ANCOVA (same as part 3)

## 5. Baseline comparison between groups using unpaired Wilcoxon test
wilcox.test(baseline_T, baseline_N, paired = FALSE, alternative = "two.sided")

