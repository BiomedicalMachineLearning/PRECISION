# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 4th June 2020
# Title: 11.2_Seurat2CellRanger.R
# Goal: To convert raw counts in a Seurat object into output like CellRanger makes
# Usage: Rscript 11.2_Seurat2CellRanger.R {sampleID} {seurat.path} {outdir}
# Where: seurat.path is the path to an RDS file with a Seurat object in it
# Where: outdir is a folder which DOES NOT ALREADY contain "CellRanger"
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  cat("ERROR: 3 arguments expected\n")
  cat("example: Rscript 11.2_Seurat2CellRanger.R {sampleID} {seurat.path} {outdir}\n")
  quit()
}

sampleID <- args[1]
seurat.path <- args[2]
outdir <- args[3]

# load packages
bioc_packages <- c("DropletUtils")
r_packages <- c("Seurat")

## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
  }
}
## function to load bioconductor packages
biocondpkgTest <- function(x) {
  if (!require(x,character.only = TRUE)) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite(x)
    if(!require(x,character.only = TRUE)) stop (paste0(x, "bioconductor package not found"))
  }
}
## load packages
for (b_pkg in bioc_packages) {
  biocondpkgTest(b_pkg)
}
for (r_pkg in r_packages) {
  baseRpkgTest(r_pkg)
}

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

# ------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------

df <- readRDS(seurat.path)
df_counts <- df@assays$RNA@counts
write10xCounts(x = df_counts, path = paste0(outdir, sampleID, "_CellRanger/"), version = "2")

