# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 25th June 2020
# Title: 10_runIntegration.R
# Goal: To integrate the 6 zs-green samples
# Usage: Rscript 10_runIntegration.R {outdir} {C1} {C2} {C3} {C4} {C5} {C6}
# Run on Delta2 or change the path to the Seurat installation
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------
# C1_ZsGreen_1dpo_SCI_VEH_123
# C2_ZsGreen_1dpo_SCI_IVIG_456
# C3_ZsGreen_7dpo_SCI_VEH_123
# C4_ZsGreen_7dpo_SCI_IVIG_456
# C5_ZsGreen_3dpo_SCI_VEH_123
# C6_ZsGreen_3dpo_SCI_IVIG_456
# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 7) {
  cat("ERROR: 7 arguments expected\n")
  cat("example: Rscript 10_runIntegration.R {outdir} {C1} {C2} {C3} {C4} {C5} {C6}\n")
  quit()
}

outdir <- args[1]
C1 <- args[2]
C2 <- args[3]
C3 <- args[4]
C4 <- args[5]
C5 <- args[6]
C6 <- args[7]


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
C2.1T <- readRDS(C2)
C3.7C <- readRDS(C3)
C4.7T <- readRDS(C4)
C5.3C <- readRDS(C5)
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

C1.1C$time <- "d1"
C1.1C$state <- "cntl"
C1.1C$timestate <- "d1_cntl"
C5.3C$time <- "d3"
C5.3C$state <- "cntl"
C5.3C$timestate <- "d3_cntl"
C3.7C$time <- "d7"
C3.7C$state <- "cntl"
C3.7C$timestate <- "d7_cntl"

# merge the 2 datasets
# NOTE that the order matters, it determines the -X suffic added to the cell barcodes
all <- merge(x = C1.1C, y = c(C2.1T, C5.3C, C6.3T, C3.7C, C4.7T))
Project(object = all) <- 'all'

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

saveRDS(all.combined, file = paste0(outdir, "all_integrated.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()
