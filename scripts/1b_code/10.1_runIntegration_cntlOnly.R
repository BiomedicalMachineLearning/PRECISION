# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 25th June 2020
# Title: 10_runIntegration.R
# Goal: To integrate the 3 control zs-green samples
# Usage: Rscript 10_runIntegration.R {outdir} {C1} {C3} {C5}
# Run on Delta2 or change the path to the Seurat installation
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------
# C1_ZsGreen_1dpo_SCI_VEH_123
# C3_ZsGreen_7dpo_SCI_VEH_123
# C5_ZsGreen_3dpo_SCI_VEH_123
# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 10_runIntegration.R {outdir} {C1} {C3} {C5}\n")
  quit()
}

outdir <- args[1]
C1 <- args[2]
C3 <- args[3]
C5 <- args[4]


# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

# ------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------

# do on Delta2
library(Seurat, lib.loc="/home/uqlgrice/R/x86_64-pc-linux-gnu-library/3.6")

# ------------------------------------------------------------------
# PREPARE SEURAT OBJECTS
# ------------------------------------------------------------------

C1.1C <- readRDS(C1)
C3.7C <- readRDS(C3)
C5.3C <- readRDS(C5)

# assign treatment/control category
C1.1C$time <- "d1"
C1.1C$state <- "cntl"
C1.1C$timestate <- "d1_cntl"
C5.3C$time <- "d3"
C5.3C$state <- "cntl"
C5.3C$timestate <- "d3_cntl"
C3.7C$time <- "d7"
C3.7C$state <- "cntl"
C3.7C$timestate <- "d7_cntl"

# I want to keep the cell barcodes consistent between the "full" dataset and the "control only" dataset. Thus, I need to append unique suffixes to the cell barcodes first
# The naming here is a bit confusing because of the 5s and the 3s. But in the original the order was C1, T1, C3, T3, C7, T7. Thus C1 = 1, C3 = 3 and C7 = 5.
C1.1C <- Seurat::RenameCells(object = C1.1C, new.names = paste0(colnames(x = C1.1C[["RNA"]]), "_1"))
C5.3C <- Seurat::RenameCells(object = C5.3C, new.names = paste0(colnames(x = C5.3C[["RNA"]]), "_3"))
C3.7C <- Seurat::RenameCells(object = C3.7C, new.names = paste0(colnames(x = C3.7C[["RNA"]]), "_5"))

# merge the 2 datasets
# NOTE that the order matters, it determines the -X suffic added to the cell barcodes
all <- merge(x = C1.1C, y = c(C5.3C, C3.7C))
Project(object = all) <- 'cntl'

# make a list of objects
all.list <- SplitObject(all, split.by = "timestate")

# ------------------------------------------------------------------
# PERFORM INTEGRATION
# ------------------------------------------------------------------

all.anchors <- FindIntegrationAnchors(object.list = all.list, dims = 1:50, anchor.features = 2000)
all.combined <- IntegrateData(anchorset = all.anchors, dims = 1:50)

# ------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------

saveRDS(all.combined, file = paste0(outdir, "cntl_integrated.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()
