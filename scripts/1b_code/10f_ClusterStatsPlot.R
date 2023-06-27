# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: 10f_ClusterStatsPlot.R
# Goal: To plot a stacked bar plot and a heatmap for clustering vs. some condition of interest
# Usage: Rscript 10f_ClusterStatsPlot.R {sampleID} {comparisonID} {cluster_col} {annotation_col} {seuratObj} {outdir}
# WHERE:
# sampleID and comparisonID are useful names used for saving the data
# cluster_col and annotation_col are the column names for the clustering and some annotation of interest
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  cat("ERROR: 6 arguments expected\n")
  cat("example: Rscript 10f_ClusterStatsPlot.R {sampleID} {comparisonID} {cluster_col} {annotation_col} {seuratObj} {outdir}\n")
  quit()
}

sampleID <- args[1]
comparisonID <- args[2]
cluster_col <- args[3]
annotation_col <- args[4]
df <- readRDS(args[5])
outdir <- args[6]

# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "pheatmap", "ggplot2")
## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
  }
}
# ## function to load bioconductor packages
# biocondpkgTest <- function(x) {
#   if (!require(x,character.only = TRUE)) {
#     source("http://www.bioconductor.org/biocLite.R")
#     biocLite(x)
#     if(!require(x,character.only = TRUE)) stop (paste0(x, "bioconductor package not found"))
#   }
# }
# ## load packages
# for (b_pkg in bioc_packages) {
#   biocondpkgTest(b_pkg)
# }
for (r_pkg in r_packages) {
  baseRpkgTest(r_pkg)
}

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

# ------------------------------------------------------------------
# MAKE PLOTS
# ------------------------------------------------------------------

# Get the data
plot <- DimPlot(df, group.by = cluster_col, split.by = annotation_col)

# Plot 1: Stacked bar
data <- as.data.frame(table(plot$data[,3:4]))
colnames(data) <- c("cluster", "condition", "count")
# Stacked + percent
ggplot(data, aes(fill=condition, y=count, x=cluster)) + 
  geom_bar(position="fill", stat="identity")
ggsave(paste0(outdir, sampleID, "_", comparisonID, "_barplot.pdf"))

# Plot 2: Heatmap
mat <- as.matrix(table(plot$data[,3:4]))
pheatmap(mat, scale = "row", filename = paste0(outdir, sampleID, "_", comparisonID, "_heatmap.pdf"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()