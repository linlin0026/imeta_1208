# Serum Metabolomics Analysis Workflow (based on MetaboAnalystR)

# 1. Data filtering was performed using MetaboAnalyst (https://www.metaboanalyst.ca)

# 2. PERMANOVA to test group-level differences in metabolomic profiles
# 'scaled_data' contains standardized metabolite concentration data
# 'metadata' includes sample information; select one grouping variable for comparison

distance_matrix <- vegan::vegdist(scaled_data, method = "euclidean")
group <- metadata[[i]]  # select the i-th metadata column as grouping variable

permanova_result <- vegan::adonis2(distance_matrix ~ group, permutations = 999)


# 3. Principal Component Analysis (PCA) for dimensionality reduction
# Note: If 'scaled_data' is already scaled, set center = FALSE, scale. = FALSE

pca_result <- prcomp(scaled_data,
                     center = FALSE,
                     scale. = FALSE,
                     rank. = 10)

# 4. Differential metabolite analysis (pre- vs post-intervention)
# Input: 'Raw_data_for_LocalR.csv' (raw serum metabolite concentrations)
library(MetaboAnalystR)

### local
mSet<-InitDataObjects("conc", "stat", TRUE)
mSet<-Read.TextData(mSet, paste0("Raw_data_for_LocalR.csv"), "rowp", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)

# Normalization: SumNorm + LogNorm + Pareto scaling
mSet <- Normalization(mSet, "SumNorm", "LogNorm", "ParetoNorm", ratio = FALSE, ratioNum = 20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# Fold-change analysis (0 = log2, 1 = threshold = 2-fold change)
mSet <- FC.Anal(mSet, 0, 1, TRUE)
mSet <- PlotFC(mSet, "fc_1_", "png", 72, width = NA)

# Non-parametric paired t-tests (Wilcoxon signed-rank test)
mSet <- Ttests.Anal(mSet, paired = TRUE, thresh = 0.05,
                    equal.var = TRUE, nonlog = TRUE, nonpar = TRUE)
mSet <- PlotTT(mSet, "tt_1_", "png", 72, width = NA)


####oplsda
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", "png", 72, width=NA, 1,2,0.95,0,0, "na")
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "png", 72, width=NA);
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "png", 72, width=NA, "vip", "tscore", 15,FALSE)
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", "png", 72, width=NA)
mSet<-OPLSDA.Permut(mSet, 100)
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_0_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)

### enrichment
# 'combind_Df' refers to the list of metabolites that were significantly upregulated after intervention,
# filtered based on predefined differential metabolite criteria.

mSet<-InitDataObjects("conc", "msetora", FALSE) #paired false
cmpd.vec<- combind_Df$KEGG
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "kegg");
mSet<-CreateMappingResultTable(mSet)
mSet<-CreateMappingResultTable(mSet)
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-SetCurrentMsetLib(mSet, "kegg_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 300, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "pdf", 300, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 300, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "pdf", 300, width=NA)
# mSet<-SaveTransformedData(mSet)

### pathway
mSet<-InitDataObjects("conc", "pathora", FALSE)
cmpd.vec<-combind_Df$KEGG
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "kegg");
mSet<-CreateMappingResultTable(mSet)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")
mSet<-PlotPathSummary(mSet, T, "path_view_0_", "pdf", 300, width=20, NA, NA ) #show.grid TRUE
mSet<-PlotPathSummary(mSet, T, "path_view_0_", "png", 300, width=NA, NA, NA )
