# 1. Analyze abundance dynamics of each CAG during the intervention period using linear mixed-effects models
# 'each_CAG_abundance' contains log10-transformed abundance for one CAG across samples
# The mean value during the baseline (normal diet) period is used as reference

# Fit linear mixed-effects model: day as fixed effect, subject as random effect
model1 <- lmer(log_CAG ~ day + (1 | subject), data = each_CAG_abundance)

# Extract fixed-effect estimates
model_summary <- summary(model1)
model_results <- data.frame(model_summary$coefficients, check.names = FALSE) %>%
  tibble::rownames_to_column("item") %>%
  mutate(adjusted_slope = ifelse(
    item != "(Intercept)",
    Estimate + model_results[model_results$item == "(Intercept)", "Estimate"],
    Estimate
  ))


# 2. Local Similarity Analysis
# LSA was performed in a Linux environment to evaluate similarity between temporal trends
# of CAGs and metabolites based on their LMM-derived slopes.


# 3. Correlation between fecal metabolites and CGM indicators
# 'common_df' contains matched fecal metabolite concentrations and CGM indicators
# for each subject and each time point

# Perform repeated measures correlation between a metabolite and CGM feature
rmcorr(subject = common_df$subject,
       measure1 = common_df[, meta],      # fecal metabolite
       measure2 = common_df[, CGM_item],  # CGM index
       dataset = common_df)