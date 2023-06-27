# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: 5_AddQCMetadata.R
# Goal: To add phase and doublet metadata to Seurat object, and filter based on doublets
# Usage: Rscript 5_AddQCMetadata.R {sampleID} {SeuratObj} {phase} {doublets} {outdir}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat("ERROR: 5 arguments expected\n")
  cat("example: Rscript 5_AddQCMetadata.R {analysisID} {SeuratObj} {phase} {doublets} {outdir}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
phase <- read.delim(args[3], sep = "\t", row.names = 1)
doublets <- read.delim(args[4], sep = "\t", row.names = 1)
outdir <- args[5]

# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "ggplot2")
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
# PRE-PROCESS DATA
# ------------------------------------------------------------------

# remove unnecessary columns produced in previous script
doublets$orig.ident <- NULL
doublets$nFeature_RNA <- NULL
# give colname to phase data
colnames(phase) <- "phase"

# ------------------------------------------------------------------
# ADD METADATA
# ------------------------------------------------------------------

df <- AddMetaData(df, metadata = doublets)
df <- AddMetaData(df, metadata = phase)

# ------------------------------------------------------------------
# FILTER BY DOUBLET CALLS
# ------------------------------------------------------------------

# doublets must be called a doublet by at least 2 of 3 scds methods, and have at least 3000 features

df <- subset(df, subset = scds.nFiltCalls < 2)

# ------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------

saveRDS(df, file = paste0(outdir, sampleID, "_preNorm_Seurat.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()