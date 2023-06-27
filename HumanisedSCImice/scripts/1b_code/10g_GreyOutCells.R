# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: 10g_GreyOutCells.R
# Goal: To plot a dimplot of the clustering UMAP but with the treatment cells greyed out
# Usage: Rscript 10g_GreyOutCells.R {sampleID} {clusterCol} {seuratObj} {out_figs} {out_objects}
# WHERE:
# clusterCol = name of the cluster annotation column you want to use
# out_figs = place to save figures/other data
# out_objects = place to save RDS objects
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat("ERROR: 5 arguments expected\n")
  cat("example: Rscript 10g_GreyOutCells.R {sampleID} {clusterCol} {seuratObj} {out_figs} {out_objects}\n")
  quit()
}

sampleID <- args[1]
clusterCol <- args[2]
df <- readRDS(args[3])
outdir <- args[4]
outdir_obj <- args[5]


# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "ggplot2", "dplyr")
## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
  }
}
## function to load bioconductor packages
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

if (endsWith(outdir_obj, "/") == FALSE) {
  outdir_obj <- paste0(outdir_obj, "/", sep="")
}

# ------------------------------------------------------------------
# PERFORM ANALYSIS
# ------------------------------------------------------------------

# make a new metadata column combining state (treatment vs control) and cluster ID
## Get the annotations ready
annotations <- df[[c("state", clusterCol)]]
colnames(annotations) <- c("state", "clustering")
nclust <- length(unique(annotations$clustering))
annotations$cells <- rownames(annotations)
annotations <- dplyr::mutate(annotations, cluststate = paste0(annotations$state, "_", annotations$clustering))
rownames(annotations) <- annotations$cells
# make another column where all "treatment" cells are just called "treat"
annotations$cluststate_simple <- annotations$cluststate
annotations$cluststate_simple <- sub("^treat.*", "treat", annotations$cluststate_simple)
annotations$cluststate_simple <- sub("^cntl_", "", annotations$cluststate_simple)
# make into factor so the labels are in the right order
mylevels <- c("treat", as.character(0:nclust-1))
annotations$cluststate_simple <- factor(annotations$cluststate_simple, levels = mylevels)
annotations <- annotations[ , c("cluststate", "cluststate_simple"), drop = FALSE]
## Add to Seurat object
df <- AddMetaData(df, annotations)

# Plot cells
mypalette <- c(rep("grey", 1), rev(DiscretePalette(nclust, palette = "alphabet")))
mypalette_2tone  <- c(rep("grey40", 1), rep("coral1", nclust))

## Plot with colours for each cluster
plot1 <- DimPlot(df, group.by = "cluststate_simple", cols = mypalette)
plot1$data <- plot1$data[order(plot1$data$cluststate_simple),]
ggsave(plot = plot1, paste0(outdir, sampleID, "_greyedOutTreatment.pdf"))

## Plot treatment vs. control
plot2 <- DimPlot(df, group.by = "cluststate_simple", cols = mypalette_2tone)
plot2$data <- plot2$data[order(plot2$data$cluststate_simple),]
ggsave(plot = plot2, paste0(outdir, sampleID, "_treatVScntl.pdf"))

## Plot treatment vs. control, split
DimPlot(df, group.by = "cluststate_simple", split.by = "state", cols = mypalette_2tone)
ggsave(paste0(outdir, sampleID, "_treatVScntl_split.pdf"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

saveRDS(df, file = paste0(outdir_obj, sampleID, "Clustered_ann_Seurat.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()