# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 17th June 2020
# Title: cellranger2seuratSCE.R
# Goal: To convert Cellranger output to Seurat and/or SCE object
# Usage: Rscript cellranger2seuratSCE.R {sampleID} {data_path} {seurat?} {SCE?} {outdir}
# Where: 
# sampleID = any sample name (will be used for output file)
# data_path = path/to/cellranger/output (folder containing barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
# seurat? = whether or not to output a Seurat object (yes/no)
# SCE? = whether or not to output a SCE object (yes/no)
# outdir = where to save the object/s
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat("ERROR: 5 arguments expected\n")
  cat("example:  Rscript cellranger2seuratSCE.R {analysisID} {data_path} {seurat?} {SCE?} {outdir}\n")
  quit()
}

sampleID <- args[1]
data_path <- args[2]
runseurat <- args[3]
runsce <- args[4]
outdir <- args[5]

# load packages
bioc_packages <- c("SingleCellExperiment")
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

# Check variables are entered correctly 
## checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

if (runseurat == "yes" | runseurat == "no") {
  print(paste0("running seurat? ... ", runseurat))
} else {
  stop("the argument 'runseurat' should be yes or no")
}


if (runsce == "yes" | runsce == "no") {
  print(paste0("running SCE? ... ", runsce))
} else {
  stop("the argument 'runsce' should be yes or no")
}

# ------------------------------------------------------------------
# LOAD THE DATA AND MAKE SEURAT OBJECT
# ------------------------------------------------------------------

# nb: this is done even if you don't want to save the output
df.data <- Read10X(data.dir = data_path)
print("creating seurat object. Filtering genes found in no cells, or cells containing no features")
df <- CreateSeuratObject(counts = df.data, min.cells = 1, min.features = 1)

if (runseurat == "yes") {
  print("saving seurat object")
  saveRDS(df, file = paste0(outdir, sampleID, "_rawSeurat.RDS"))
}

# ------------------------------------------------------------------
# OPTIONAL: CONVERT TO SCE OBJECT AND SAVE
# ------------------------------------------------------------------

if (runsce == "yes") {
  df.sce <- as.SingleCellExperiment(df)
  print("saving SCE object")
  saveRDS(df.sce, file = paste0(outdir, sampleID, "_rawSCE.RDS"))
}

# ------------------------------------------------------------------
# FINISHED
# ------------------------------------------------------------------

print("finished!")