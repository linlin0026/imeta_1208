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

# 4. Differential metabolite analysis between experimental phases per subject
# Input file: 'Raw_data_for_LocalR.csv' (unscaled fecal metabolite concentrations)
library(MetaboAnalystR)
# Initialize MetaboAnalystR object (unpaired design, concentration data)
mSet <- InitDataObjects("conc", "stat", FALSE)
mSet<-Read.TextData(mSet, paste0("Raw_data_for_LocalR.csv"), "rowu", "disc")### 注意这里不配对  rowu而不是rowp
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-PreparePrenormData(mSet)

# Apply normalization: Sum normalization, log transformation, and Pareto scaling
mSet <- Normalization(mSet, "SumNorm", "LogNorm", "ParetoNorm", ratio = FALSE, ratioNum = 20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)

# Perform unpaired non-parametric statistical test 
mSet<-Ttests.Anal(mSet, T, 0.05, FALSE, TRUE, "raw", TRUE)

# OPLS-DA analysis for multivariate discrimination between groups
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
