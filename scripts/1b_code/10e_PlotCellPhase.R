# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 26th June 2020
# Title: 10e_PlotCellphase.R
# Goal: To make plots and tables showing the breakdown of cell phase across clusters/identities
# Usage: Rscript 10e_PlotCellphase.R {sampleID} {clusteringColName} {seuratObject} {outdir}
# Where "clusteringColName" is the column of metadata you want to group by (i.e. best clustering result)
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 10e_PlotCellphase.R {sampleID} {clusteringColName} {seuratObject} {outdir}\n")
  quit()
}

sampleID <- args[1]
column <- args[2]
df <- readRDS(args[3])
outdir <- args[4]


# load packages
bioc_packages <- c()
r_packages <- c("Seurat")
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
# DATA TABLE
# ------------------------------------------------------------------

Idents(object = df) <- column
phase.table <- t(table(df[[c("phase", column)]]))
write.table(phase.table, file = paste0(outdir, sampleID, "_phase.txt"), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------
# VIZ
# ------------------------------------------------------------------

umap1 <- DimPlot(df, reduction = "umap", label = FALSE, pt.size = 0.5, order = TRUE, group.by = "phase")
umap2 <- DimPlot(df, reduction = "umap", label = FALSE, pt.size = 0.5, order = TRUE, group.by = "phase", split.by = "phase")

pdf(paste0(outdir, sampleID, "_umap-phase.pdf"), width = 12, height = 12)
umap1
dev.off()

pdf(paste0(outdir, sampleID, "_umap-phase_split.pdf"), width = 12, height = 12)
umap2
dev.off()

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()