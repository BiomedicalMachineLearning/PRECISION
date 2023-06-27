# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: 6_scranNormalisation.R
# Goal: To normalise libraries with scran
# Usage: Rscript 6_scranNormalisation.R {sampleID} {seuratObj} {outdir}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  cat("ERROR: 3 arguments expected\n")
  cat("example: Rscript 6_scranNormalisation.R {sampleID} {seuratObj} {outdir}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
outdir <- args[3]

# load packages
bioc_packages <- c("scran", "scater")
r_packages <- c("Seurat", "ggplot2")

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
# NORMALISATION
# ------------------------------------------------------------------

# function to normalise count data in scran/scater
# TEMPORARY NOTE: developed in `SandboxNormalisationScran.Rmd`
func_scranNorm <- function(seuratObj) {
  # USAGE: seuratObj <- func_scranNorm(seuratObj)
  # OUTPUT: a seuratObj object with (natural log) normalised counts
  # NOTE: usually, scran normalisation produces log2counts, but here we produce seurat-compatible lncounts
  
  # convert to SCE object
  mySCE <- Seurat::as.SingleCellExperiment(seuratObj)
  # calculate size factors and perform normalisation
  scranclusters <- scran::quickCluster(mySCE)
  mySCE <- scran::computeSumFactors(mySCE, clusters = scranclusters)
  # "scran sometimes calculates negative or zero size factors which will completely distort the normalized expression matrix". Let's check
  minsizefactor <- min(sizeFactors(mySCE))
  if (minsizefactor < 0) {
    warning("ALERT! scran normalisation has produced negative or zero size factors which will distort the normalised expression matrix. Proceed with care!
          \n You can try increasing the cluster and pool sizes until they are all positive
            \n See https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/normalization-confounders-and-batch-correction.html")
    # print the same message to a warning file
    sink(paste0(outdir, sampleID, "_warninglog.txt"))
    cat("ALERT! scran normalisation has produced negative or zero size factors which will distort the normalised expression matrix. Proceed with care!
          \n You can try increasing the cluster and pool sizes until they are all positive
        \n See https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/public/normalization-confounders-and-batch-correction.html")
    sink()
  }
  mySCE <- scater::logNormCounts(mySCE, log = FALSE, name = "unlog.normcounts")
  
  # natural-log transform counts and convert back to sparse matrix format
  assay(mySCE, "ln.normcounts") <- as(log(x = assay(mySCE, "unlog.normcounts") + 1), "dgCMatrix")
  # NOTE: To convert to Seurat object from now on your must run: seuratObj <- as.Seurat(mySCE, counts = "counts", data = "ln.normcounts")
  seuratObj <- as.Seurat(mySCE, counts = "counts", data = "ln.normcounts")
  return(seuratObj)
}

df <- func_scranNorm(df)

saveRDS(df, file = paste0(outdir, sampleID, "_postNorm_Seurat.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()