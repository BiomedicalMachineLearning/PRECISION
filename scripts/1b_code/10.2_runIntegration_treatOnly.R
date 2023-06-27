# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 25th June 2020
# Title: 10_runIntegration.R
# Goal: To integrate the 3 treatment zs-green samples
# Usage: Rscript 10_runIntegration.R {outdir} {C2} {C4} {C6}
# Run on Delta2 or change the path to the Seurat installation
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------
# C2_ZsGreen_1dpo_SCI_IVIG_456
# C4_ZsGreen_7dpo_SCI_IVIG_456
# C6_ZsGreen_3dpo_SCI_IVIG_456
# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 10_runIntegration.R {outdir} {C2} {C4} {C6}\n")
  quit()
}

outdir <- args[1]
C2 <- args[2]
C4 <- args[3]
C6 <- args[4]


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

C2.1T <- readRDS(C2)
C4.7T <- readRDS(C4)
C6.3T <- readRDS(C6)

# assign treatment/control category
C2.1T$time <- "d1"
C2.1T$state <- "treat"
C2.1T$timestate <- "d1_treat"
C6.3T$time <- "d3"
C6.3T$state <- "treat"
C6.3T$timestate <- "d3_treat"
C4.7T$time <- "d7"
C4.7T$state <- "treat"
C4.7T$timestate <- "d7_treat"

# I want to keep the cell barcodes consistent between the "full" dataset and the "control only" dataset. Thus, I need to append unique suffixes to the cell barcodes first
# The naming here is a bit confusing because of the 5s and the 3s. But in the original the order was C1, T1, C3, T3, C7, T7. Thus T1 = 2, T3 = 4 and T7 = 6.
C2.1T <- Seurat::RenameCells(object = C2.1T, new.names = paste0(colnames(x = C2.1T[["RNA"]]), "_2"))
C6.3T <- Seurat::RenameCells(object = C6.3T, new.names = paste0(colnames(x = C6.3T[["RNA"]]), "_4"))
C4.7T <- Seurat::RenameCells(object = C4.7T, new.names = paste0(colnames(x = C4.7T[["RNA"]]), "_6"))

# merge the 2 datasets
# NOTE that the order matters, it determines the -X suffic added to the cell barcodes
all <- merge(x = C2.1T, y = c(C6.3T, C4.7T))
Project(object = all) <- 'treat'

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

saveRDS(all.combined, file = paste0(outdir, "treat_integrated.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, "sessionInfo.txt"))
sessionInfo()
sink()