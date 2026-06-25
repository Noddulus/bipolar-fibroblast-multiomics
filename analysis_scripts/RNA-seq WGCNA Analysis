# ==============================================================================
# Script: RNA-seq WGCNA Analysis for Temporal Fibroblast Data
# Description: Weighted Gene Co-expression Network Analysis (WGCNA) pipeline 
#              to cluster temporally expressed genes across 6 fibroblast lines.
# ==============================================================================

# 1. Setup and Package Loading -------------------------------------------------

# Set working directory to the location of the script/data
# Users should change this to their local directory
setwd("./") 

# Load required libraries
# Note: Ensure WGCNA and edgeR are installed via BiocManager prior to running
library(WGCNA)
library(edgeR)

# Enable multi-threading for WGCNA to speed up correlation calculations
enableWGCNAThreads()

# 2. Data Loading and Preprocessing --------------------------------------------

# Load pre-filtered and normalized CPM data
OrderedCPM <- readRDS("CPMFilteredData.rds")
# Assuming CubicGenes and RNAlCPM are available in the environment from previous 
# Spline modeling steps. If not, they should be loaded here.
TimeGenes <- rownames(CubicGenes)
TimeCPM <- RNAlCPM[TimeGenes, 1:78]

# Define a generalized plotting function for temporal data
# Note: Column indices (e.g., c(1, 12, 13, 2:11)) are ordered this way to correct 
# for alphanumeric sorting in the raw data, reordering them to 0h, 4h, 8h...48h.
plot_timecourse <- function(dataset, feature_id, y_label, title) {
  
  if(feature_id %in% rownames(dataset)) {
    data_vec <- dataset[feature_id, ]
  } else if (feature_id %in% colnames(dataset)) {
    data_vec <- dataset[, feature_id]
  } else {
    stop("Feature not found in dataset.")
  }
  
  y_min <- min(data_vec[1:78], na.rm = TRUE)
  y_max <- max(data_vec[1:78], na.rm = TRUE)
  
  plot(as.numeric(data_vec[c(1, 12, 13, 2:11)]), type = "l", ylim = c(y_min, y_max), 
       main = title, ylab = y_label, xlab = "Time (Hours)", axes = FALSE, col = 2, lty = 2)
  
  axis(side = 1, at = 1:13, labels = seq(0, 48, by = 4))
  axis(2)
  box()
  
  # Define plotting indices and colors for the 6 individuals
  individuals <- list(
    ind1 = list(idx = c(1, 12, 13, 2:11), col = 2),
    ind2 = list(idx = c(14, 25, 26, 15:24), col = 3),
    ind3 = list(idx = c(27, 38, 39, 28:37), col = 4),
    ind4 = list(idx = c(40, 51, 52, 41:50), col = 7),
    ind5 = list(idx = c(53, 64, 65, 54:63), col = 5),
    ind6 = list(idx = c(66, 77, 78, 67:76), col = 6)
  )
  
  for (ind in individuals) {
    lines(as.numeric(data_vec[ind$idx]), col = ind$col, lty = 2)
    points(as.numeric(data_vec[ind$idx]), col = ind$col, pch = 16)
  }
}

# Verify data structure with a known circadian clock gene
plot_timecourse(RNAlCPM, "TEF", "lCPM values", "Temporal lCPM for TEF")


# 3. WGCNA Network Construction ------------------------------------------------

datExpr <- as.data.frame(t(TimeCPM))

# Check for outliers
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample Clustering to Detect Outliers", 
     sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# Plot Scale Free Topology Model Fit
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, col = "red")
abline(h = 0.90, col = "red")

# Plot Mean Connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")

# Construct Adjacency and TOM (Topological Overlap Matrix)
softPower <- 12 # Selected based on scale-free topology fit
adjacency <- adjacency(datExpr, power = softPower, type = "signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# 4. Module Identification -----------------------------------------------------

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)

# Calculate module eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

# Cluster module eigengenes to merge highly correlated modules
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# Merge close modules
MEDissThres <- 0.25
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Plot gene dendrogram with initial and merged module colors
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Verify Eigengene structures
par(mfrow = c(1,1))
plot_timecourse(mergedMEs, "MEblack", "Eigengene Value", "Black Module Eigengene Trajectory")

# Save Eigengenes for downstream modeling
save(MEs, mergedMEs, file = "TimelCPMWGCNAEigengnes.RData")

# 5. Exporting Results for Gene Ontology ---------------------------------------

# Recalculate MEs with merged color labels
MEs0 <- moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs <- orderMEs(MEs0)

# Calculate Module Membership
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
names(geneModuleMembership) <- paste0("MM", modNames)
names(MMPvalue) <- paste0("p.MM", modNames)

# Compile Gene Information
geneInfo0 <- data.frame(geneSymbol = colnames(datExpr),
                        moduleColor = mergedColors,
                        geneModuleMembership, 
                        MMPvalue)

write.csv(geneInfo0, file = "geneInfolCPMWGCNACubicGenes.csv", row.names = FALSE)

# 6. Automated Output for Downstream Tools (MetaScape) -------------------------

# Create a clean directory for module outputs
dir.create("Module_Gene_Lists", showWarnings = FALSE)

# Automatically write a .txt file for every unique module identified
modules <- unique(geneInfo0$moduleColor)
for (mod in modules) {
  mod_genes <- geneInfo0$geneSymbol[geneInfo0$moduleColor == mod]
  file_name <- paste0("Module_Gene_Lists/RNAseqWGCNA_", tools::toTitleCase(mod), "Genes.txt")
  write.table(mod_genes, file = file_name, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Aggregate genes by module into a single format for MetaScape multi-list input
metascape_df <- aggregate(geneSymbol ~ moduleColor, data = geneInfo0, 
                          FUN = function(x) paste(x, collapse = ","))

# Rename columns and write to file
colnames(metascape_df) <- c("Module", "Gene_List")
write.table(metascape_df, file = "RNAseqWGCNA_ForMetaScape.txt", 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
