# ==============================================================================
# Script: Figure Generation for Fibroblast Circadian Multi-Omics Paper
# Description: Generates heatmaps, time-series plots with spline fits, 
#              sLDSR enrichment bar plots, and WGCNA Eigengene trajectories.
# ==============================================================================

setwd("./") 

# Load Required Libraries
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(reshape2)
library(tidyverse)
library(gridExtra)

# ==============================================================================
# 1. Circadian Gene Heatmap (Figure)
# ==============================================================================

# Load fitted model summary data
SMESummary <- readRDS("SMEsummaryData.rds")

# Target Circadian Genes
CircadianGenes <- c("NR1D2","PER3","PER2","TRIM5","CDC25B","ARNTL","TEF","DBP",
                    "CRY1","NPAS2")

SMESummaryCircadian <- SMESummary[which(SMESummary$Gene %in% CircadianGenes),]

# Extract and Z-Score the fitted data (Columns 3:15 contain the 13 timepoints)
ZScoredCircadianGenes <- SMESummaryCircadian[, 3:15]
rownames(ZScoredCircadianGenes) <- SMESummaryCircadian$Gene

for (i in 1:nrow(ZScoredCircadianGenes)){
  row_data <- unlist(ZScoredCircadianGenes[i, ])
  ZScoredCircadianGenes[i,] <- (row_data - mean(row_data)) / sd(row_data)    
}

# Manually cluster/order genes to match the publication figure structure
order_idx <- c(2, 6, 9, 7, 8, 4, 10, 3, 1, 5)
ClusterZScoredCircadianGenes <- ZScoredCircadianGenes[order_idx, ]

# Expand grid for ggplot2 heatmap
Timepoints <- seq(0, 48, by = 4)
x <- as.character(Timepoints)
y <- rownames(ClusterZScoredCircadianGenes)
data <- expand.grid(X = x, Y = y)

# Populate grid with Z-scored values
data$Z <- as.numeric(t(ClusterZScoredCircadianGenes))

# Generate Heatmap
p_heatmap <- ggplot(data, aes(X, Y, fill = Z)) + 
  geom_tile() +
  scale_fill_viridis(discrete = FALSE, option = "plasma") +
  theme_ipsum() +
  labs(x = "Timepoints (Hrs)", y = "", fill = "Z-Score") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(p_heatmap)


# ==============================================================================
# 2. Time-Series Expression Plots with Spline Trendlines (Figure)
# ==============================================================================

# Load RNA-seq logCPM normalized data
RNAlCPM <- readRDS("lCPMFilteredData.rds")

# Extract the data for the selected genes
CircadianGeneslCPM <- RNAlCPM[rownames(RNAlCPM) %in% CircadianGenes, ]

# Reorder columns chronologically across the 6 individuals
TimeOrder <- c(1,12,13,2,3,4,5,6,7,8,9,10,11, 14,25,26,15,16,17,18,19,20,21,22,23,24,
               27,38,39,28,29,30,31,32,33,34,35,36,37, 40,51,52,41,42,43,44,45,46,47,48,49,50,
               53,64,65,54,55,56,57,58,59,60,61,62,63, 66,77,78,67,68,69,70,71,72,73,74,75,76)

OderCircadianGeneslCPM <- CircadianGeneslCPM[, TimeOrder]

# Reshape raw data to long format for plotting points
long_data <- melt(as.matrix(OderCircadianGeneslCPM), varnames = c("Gene", "Sample"), value.name = "Expression")
long_data$Time <- rep(seq(0, 48, by = 4), times = 6 * nrow(OderCircadianGeneslCPM))

# Reshape fitted spline data to long format for plotting trendlines
trend_data <- SMESummaryCircadian[, 3:15]
rownames(trend_data) <- SMESummaryCircadian$Gene
trend_data$Gene <- rownames(trend_data)

long_trend_data <- melt(trend_data, id.vars = "Gene", variable.name = "TimeCol", value.name = "Trend")
long_trend_data$Time <- rep(seq(0, 48, by = 4), each = nrow(trend_data))

# Generate Time-Series Plot
p_timeseries <- ggplot() +
  geom_point(data = long_data, aes(x = Time, y = Expression), size = 1, color = "black", alpha = 0.5) +
  geom_line(data = long_trend_data, aes(x = Time, y = Trend), color = "red", size = 1) + 
  facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(x = "Time (Hrs)", y = "Expression (lCPM)") +
  scale_x_continuous(breaks = seq(0, 48, by = 4), limits = c(0, 48)) + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"))

print(p_timeseries)


# ==============================================================================
# 3. Stratified LD Score Regression (sLDSR) Enrichment (Figure)
# ==============================================================================

# Load sLDSR results
PTSD <- read.table("HG19LDTimePeaksPTSD.results", header = TRUE)

# Color palettes for annotations
color_palette <- c(
  c("#53338a","#8b18cd", "#d7a3f5"), # Decreasing
  c("#395921","#30cd18", "#aaf5a3"), # Increasing
  c("#344b8d","#2959db","#7392e7","#c2d1f9")  # Baseline
)

# Plot PTSD Enrichment with Standard Error Bars
selected_annotations <- c(1:7, 30, 35, 64)

base_ptsd <- barplot(PTSD$Enrichment[selected_annotations],
                     main = "PTSD Partitioned Heritability Enrichment",
                     ylab = "Enrichment",
                     col = color_palette,
                     ylim = c(-40, 40))

# Add error bars
arrows(x0 = base_ptsd,                            
       y0 = PTSD$Enrichment[selected_annotations] + PTSD$Enrichment_std_error[selected_annotations],
       y1 = PTSD$Enrichment[selected_annotations] - PTSD$Enrichment_std_error[selected_annotations],
       angle = 90, code = 3, length = 0.1)

# Legend
legend("bottomleft", 
       legend = c("Decreasing Time Peaks", "Decreasing + 1kb", "Decreasing + 10kb", 
                  "Increasing Time Peaks", "Increasing + 1kb", "Increasing + 10kb",
                  "Baseline","H3K1me1","H3K9ac","Conserved Regions"),
       fill = color_palette, cex = 0.7)


# ==============================================================================
# 4. WGCNA Eigengene Trajectories (Figure)
# ==============================================================================

load("TimelCPMWGCNAEigengnes.RData")

# Extract Cell Line and Timepoint from rownames
split_names <- strsplit(rownames(MEs), " ")
MEs$CellLine <- sapply(split_names, `[`, 1) 
MEs$TimePoint <- as.numeric(sapply(split_names, `[`, 3)) 

# Reshape data
df_long <- MEs %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Column", values_to = "Value") %>%
  filter(Column %in% c("MEblue","MEyellow","MEred","MEblack","MEpink","MEbrown"))

# Custom titles including gene counts
titles <- c(MEblue = "Blue (396)", MEyellow ="Yellow (297)", MEred = "Red (216)",
            MEblack = "Black (207)", MEpink = "Pink (161)", MEbrown ="Brown (381)")

# Generate Eigengene Trajectory Plot
p_eigengene <- ggplot(df_long, aes(x = TimePoint, y = Value, group = CellLine, color = CellLine)) +
  geom_line() + 
  geom_point(shape = 16) + 
  scale_color_manual(values = c("red", "green", "blue", "#FFDB58", "cyan", "purple")) + 
  facet_wrap(~ Column, ncol = 2, scales = "free_y", labeller = as_labeller(titles)) + 
  labs(x = "Time Point (Hours)", y = "Eigengene Value", color = "Cell Line") +
  scale_x_continuous(breaks = seq(0, 48, by = 4), limits = c(0, 48)) + 
  theme_minimal() +
  theme(strip.text = element_text(size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.ticks = element_line(color = "black"))

print(p_eigengene)
