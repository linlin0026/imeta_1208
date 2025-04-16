# 1.  Prepare per-subject abundance data

# Convert abundance matrix to binary presence/absence (1 = present, 0 = absent)
presence_matrix <- individual_relative_abundance
presence_matrix[presence_matrix > 0] <- 1

# Count the number of subjects in which each bin is present
frequency_table <- data.frame(
  bin_id = rownames(presence_matrix),
  count = rowSums(presence_matrix)
)

# Filter bins that are present in at least 45% of subjects
filtered_bins <- frequency_table[frequency_table$count >= 0.45 * ncol(individual_relative_abundance), ]

# 2.  Compute correlations between bins

# FastSpar was executed on a Linux system

# 3. Select stable correlations (p < 0.05) shared by at least 6 subjects in one group

# 'edge_df' contains all possible pairwise bin combinations (rows)
# Each column represents whether the corresponding edge is significant in a subject
# The value is the correlation coefficient if significant, otherwise 0

subject_num <- 6
edge_df %>%
  filter(rowSums(. > 0) >= subject_num | rowSums(. < 0) >= subject_num)


# 4. Community Assembly Group (CAG) detection

# Method based on:
# Wu G, Zhao N, Zhang C, Lam YY, Zhao L. Guild-based analysis for understanding gut microbiome in human health and diseases. 
# Genome Med. 2021;13(1):22. doi:10.1186/s13073-021-00840-y
        
        


# 5. Based on CAG classification, calculate the total abundance of each CAG in each sample

# Note: CAG abundance is computed separately for each  group 
# This step aggregates bin-level relative abundances into CAG-level values
# resulting in a 'CAG_abundance' table

# 6. Analyze the abundance dynamics of each CAG during the intervention period

# 'each_CAG_abundance' contains log10-transformed abundance values of a single CAG for each sample
# The mean abundance during the normal diet period is used as the baseline

# Fit a linear mixed-effects model with subject as a random effect
 lmer(log_CAG ~ day + (1 | subject), data = each_CAG_abundance)

 
# 7. Focus on CAG abundance changes between two time points: day 13 and day 27
# 'cag_abundance' contains CAG-level relative abundances for each sample
# 'metadata' includes corresponding sample metadata with a 'day' column (e.g., day 13 or 27)

linda.obj <- linda(cag_abundance,
                  metadata,
                  formula = '~ day + (1 | subject)',
                  alpha = 0.05,               # significance threshold
                  type = "proportion",        # input data is relative abundance
                  adaptive = TRUE,            # use adaptive shrinkage
                  p.adj.method = "BH",        # Benjamini-Hochberg correction
                  lib.cut = 0,                # no sample filtering by library size
                  prev.cut = 0,               # no taxa filtering by prevalence
                  winsor.quan = NULL,         # no winsorization
                  n.cores = 1)                # number of CPU cores for parallel computation


# 8. Identify temporal patterns of CAGs based on their abundance trajectories over time
# Input: 'CAG_abundance' table (rows = CAGs, columns = time points or samples ordered by time)

# Set DTW parameters
dtw_method <- "sakoechiba"
window.size <- 7

# Define DTW distance function
dtw_distance <- function(x, y) {
 dtw::dtw(x, y,
          keep.internals = FALSE,
          window.type = dtw_method,
          window.size = window.size)$distance
}

# Create a DTW-based distance matrix between CAGs
dtw_dist_matrix <- proxy::dist(CAG_abundance, method = dtw_distance)

# Perform hierarchical clustering using Ward's method
hc <- hclust(dtw_dist_matrix, method = "ward.D2")

# 9. Compare the proportions of functional gene categories (e.g., CAZy or KO) 
# across CAG clusters classified as increasing, decreasing, or unchanged over time

# Example: Compare "Increasing" vs "Decreasing" clusters using Wilcoxon test
wilcox.test(cluster_increasing$gene_ratio,
            cluster_decreasing$gene_ratio,
            paired = FALSE)

# Within-clusters comparison: Test variation among clusters within one trend group
# Example: Within the "Increasing" group, test differences among CAGs
kruskal.test(gene_ratio ~ cluster, data = cluster_increasing)

# Post hoc pairwise comparisons using Dunn's test (with BH correction)
PMCMRplus::kwAllPairsDunnTest(gene_ratio ~ cluster,
                              data = cluster_increasing,
                              p.adjust.method = "BH")
