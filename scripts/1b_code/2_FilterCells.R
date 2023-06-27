# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: 2_filterCells.R
# Goal: To filter a Seurat object as a result of 3 MAD filtering
# Usage: Rscript 2_filterCells.R {sampleID} {seuratObject} {cellsToKeep} {outdir}
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("usage: Rscript myscript.R arg1 arg2 arg3 arg4\n")
  cat("example: Rscript 2_filterCells.R {sampleID} {seuratObject} {cellsToKeep} {outdir}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
cells <- args[3]
outdir <- args[4]

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
# SUBSET CELLS
# ------------------------------------------------------------------

keep <- read.delim(cells, header = FALSE)
keep <- keep$V1
df <- subset(df, cells = keep)

# ------------------------------------------------------------------
# OUTPUT FILTERED OBJECTS
# ------------------------------------------------------------------

saveRDS(df, file = paste0(outdir, sampleID, "_postQC_Seurat.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()
