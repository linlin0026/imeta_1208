###Index Calculation

### Calculate alpha diversity indices from relative abundance data
# Input: relative_abundance (rows = samples, columns = taxa)

# Shannon diversity index
shannon <- vegan::diversity(relative_abundance, index = "shannon")

# Observed richness (OTU/ASV count per sample)
richness <- vegan::specnumber(relative_abundance)

# Pielou's evenness
evenness <- shannon / log(richness)

### Compute beta diversity distance matrix (Bray–Curtis dissimilarity)
distance_matrix <- vegan::vegdist(relative_abundance, method = "bray")

### Perform aPCoA (adjusted Principal Coordinates Analysis)
# subject: subject ID
# metadata: data frame containing sample metadata
# group: main grouping variable (e.g., Low dose phase)
# color: color vector for plotting

aPCoA(distance_matrix ~ subject,
      metadata = metadata,
      maincov = group,
      col = color,
      cex = 1,
      drawCenter = FALSE,
      drawEllipse = FALSE)


# Perform PERMANOVA to test group differences in beta diversity
# relative_abundance: abundance table
# group: grouping variable (e.g., Overwieght, T2D)
# subject: subject ID

vegan::adonis2(relative_abundance ~ group,
               data = metadata,
               method = "bray",
               by = "margin",
               strata = metadata$subject)

