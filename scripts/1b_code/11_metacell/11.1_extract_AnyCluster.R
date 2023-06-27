# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 4th March 2020
# Title: extract_AnyCluster.R
# Goal: To extract a Seurat object based on one or more clusters (or other numeric annotations) from the same column
# Usage: Rscript 3 sampleID query.path cluster.col cluster.vals outdir
# where cluster.col is the metadata column name containing whatever metadata you want to subset by (e.g. clustering)
# where cluster.vals is a comma separated list e.g. 1,2,3,5 ()
# ------------------------------------------------------------------
# USAGE EXAMPLE:
# Rscript 1_extract_AnyCluster.R nkt_cl8 seurat.RDS integrated_snn_res.0.1 1,3,5 path/to/outdir/
# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)
sampleID <- args[1]
query.path <- args[2]
cluster.col <- args[3]
cluster.vals <- args[4]
outdir <- args[5]

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  cat("ERROR: 5 argument expected\n")
  cat("example: Rscript extract_2clusters.R sampleID query.path cluster.col cluster.vals outdir\n")
  quit()
}

if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(Seurat, quietly = TRUE))

# ------------------------------------------------------------------
# LOAD THE DATA
# ------------------------------------------------------------------

query <- readRDS(file = query.path)

# ------------------------------------------------------------------
# SUBSET THE OBJECT
# ------------------------------------------------------------------

cluster.col.nbr <- which(colnames(query[[]])==cluster.col)
myclusters <- strsplit(cluster.vals, ",")[[1]]
cluster.col.data <- query[[]][,cluster.col.nbr]
cellIDs <- rownames(query[[]][which(cluster.col.data %in% myclusters),])
query.subset <- subset(query, cells = cellIDs)
query.subset[[cluster.col]] <- droplevels(query.subset[[cluster.col]])

# ------------------------------------------------------------------
# SAVE THE OUTPUT
# ------------------------------------------------------------------

saveRDS(query.subset, paste0(outdir, sampleID, "_cluster-", cluster.vals, "_subset.RDS"))
