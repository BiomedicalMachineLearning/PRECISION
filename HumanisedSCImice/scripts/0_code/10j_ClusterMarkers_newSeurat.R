# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: 10j_ClusterMarkers.R
# Goal: To run marker predictions for a given Seurat object
# Usage: Rscript 10j_ClusterMarkers.R {sampleID} {cluster_col} {seuratObj} {outdir}
# WHERE:
# sampleID and comparisonID are useful names used for saving the data
# cluster_col and annotation_col are the column names for the clustering and some annotation of interest
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 10j_ClusterMarkers.R {sampleID} {cluster_col} {seuratObj} {outdir}\n")
  quit()
}

sampleID <- args[1]
cluster_col <- args[2]
df <- readRDS(args[3])
outdir <- args[4]

# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "dplyr")
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
# SET IDENTS
# ------------------------------------------------------------------

DefaultAssay(df) <- "RNA"
Idents(df) <- cluster_col

# ------------------------------------------------------------------
# FIND ALL MARKERS
# ------------------------------------------------------------------

df.markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
df.markers <- df.markers %>% filter(p_val_adj <= 0.05)
top10.markers <- df.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 10)

write.table(df.markers, file = paste0(outdir, sampleID, "_allmarkers.txt"), quote = FALSE, col.names = NA, sep = "\t")
write.table(top10.markers, file = paste0(outdir, sampleID, "_top10markers.txt"), quote = FALSE, col.names = NA, sep = "\t")

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()