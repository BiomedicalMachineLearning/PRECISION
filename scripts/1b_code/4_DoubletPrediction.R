# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23 June 2020
# Title: 4_DoubletPrediction.R
# Goal: To predict doublets
# Usage: Rscript 4_DoubletPrediction.R {sampleID} {seuratObj} {outdir}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  cat("ERROR: 3 arguments expected\n")
  cat("example: Rscript 3_CellCyclePrediction.R {analysisID} {seuratObj} {outdir}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
outdir <- args[3]
# load packages
bioc_packages <- c("scds")
r_packages <- c("Seurat", "ggplot2", "dplyr")
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
# DEFINE FUNCTION
# ------------------------------------------------------------------

to.pdf <- function(expr, filename, ...) {
  pdf(filename, ...)
  on.exit(dev.off())
  print(eval.parent(substitute(expr)))
}


# ------------------------------------------------------------------
# QC 3 - DOUBLETS
# ------------------------------------------------------------------

func_predictDoublets <- function(seuratObj){
  # USAGE: mySCE <- func_predictDoublets(mySCE)
  # OUTPUT: an SCE object with doublets predicted, and removed if appropriate
  
  # we will remove cells if they are (1) predicted to be a doublet by 2-3 scds strategies and (3) if the cell also has more than n features; specify this value here
  mySCE <- as.SingleCellExperiment(seuratObj)
  nFeature_cutoff <- 3000
  
  # calculate doublet scores (and calls)
  mySCE <- cxds(mySCE, retRes = TRUE, verb=TRUE, estNdbl=TRUE)
  mySCE <- bcds(mySCE, retRes = TRUE, verb=TRUE, estNdbl=TRUE)
  mySCE = cxds_bcds_hybrid(mySCE, estNdbl=TRUE)
  
  # work out which cells we want to filter
  mytable <- data.frame(mySCE@colData)
  
  mytable$barcode <- rownames(mytable)
  mytable <- mutate(mytable, FiltStatus.cxds = ifelse(mytable$cxds_call == "TRUE" & mytable$nFeature_RNA > nFeature_cutoff, "Filter", "NoFilter"))
  mytable <- mutate(mytable, FiltStatus.bcds = ifelse(mytable$bcds_call == "TRUE" & mytable$nFeature_RNA > nFeature_cutoff, "Filter", "NoFilter"))
  mytable <- mutate(mytable, FiltStatus.hybrid = ifelse(mytable$hybrid_call == "TRUE" & mytable$nFeature_RNA > nFeature_cutoff, "Filter", "NoFilter"))
  # do 0, 1, 2 or 3 methods predict that a doublet (csds + nFeatures > 3000) should be called?
  mytable <- mutate(mytable, scds.nFiltCalls = rowSums(mytable[,c("FiltStatus.cxds", "FiltStatus.bcds", "FiltStatus.hybrid")] == "Filter"))
  
  # output a data table
  rownames(mytable) <- mytable$barcode
  mytable$barcode <- NULL
  mytable$nCount_RNA <- NULL
  mytable$ident <- NULL
  write.table(mytable, file = paste0(outdir, sampleID, "_doubletStats.txt"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
  
  # define some graph functions which will be run with `to.pdf` later
  ## cxds vs nFeature
  fig.doubletscoreVSnFeature.cxds <- function() {
    n.FiltStatus.cxds <- sum(mytable$FiltStatus.cxds == "Filter")
    plot.cxds <- ggplot(mytable, aes(x=cxds_score, y=nFeature_RNA, col = cxds_call)) + 
      geom_point() +
      scale_colour_manual(values = c("grey", "red")) +
      ggtitle(paste0(sampleID, " - scds scores (cxds). Possible doublets where nFeature > 3k = ", n.FiltStatus.cxds)) +
      geom_hline(yintercept=c(3000), linetype="dashed", color = "red") +
      labs(color = "possible doublet?")
  }
  ## bcds vs nFeature
  fig.doubletscoreVSnFeature.bcds <- function() {
    n.FiltStatus.bcds <- sum(mytable$FiltStatus.bcds == "Filter")
    plot.bcds <- ggplot(mytable, aes(x=bcds_score, y=nFeature_RNA, col = bcds_call)) + 
      geom_point() +
      scale_colour_manual(values = c("grey", "red")) +
      ggtitle(paste0(sampleID, " - scds scores (bcds). Possible doublets where nFeature > 3k = ", n.FiltStatus.bcds)) +
      geom_hline(yintercept=c(3000), linetype="dashed", color = "red") +
      labs(color = "possible doublet?")
  }
  ## hybrid vs nFeature  
  fig.doubletscoreVSnFeature.hybrid <- function() {
    n.FiltStatus.hybrid <- sum(mytable$FiltStatus.hybrid == "Filter")
    plot.hybrid <- ggplot(mytable, aes(x=hybrid_score, y=nFeature_RNA, col = hybrid_call)) + 
      geom_point() +
      scale_colour_manual(values = c("grey", "red")) +
      ggtitle(paste0(sampleID, " - scds scores (hybrid). Possible doublets where nFeature > 3k = ", n.FiltStatus.hybrid)) +
      geom_hline(yintercept=c(3000), linetype="dashed", color = "red") +
      labs(color = "possible doublet?")
  }
  
  ## #DoubletCalls vs nFeature  
  fig.doubletscoreVSnFeature.nCalls <- function() {
    n.FiltStatus.nCalls <- sum(mytable$scds.nFiltCalls > 1)
    plot.hybrid <- ggplot(mytable, aes(x=hybrid_score, y=nFeature_RNA, col = as.factor(scds.nFiltCalls))) + 
      geom_point() +
      scale_colour_manual(values = c("grey", "darkgreen", "orange", "red")) +
      ggtitle(paste0(sampleID, " - scds scores (hybrid). Possible doublets where nFeature > 3k = ", n.FiltStatus.nCalls)) +
      geom_hline(yintercept=c(3000), linetype="dashed", color = "red") +
      labs(color = "possible doublet?")
  }
  
  # Run the figure functions and save graphs as PDFs
  to.pdf(fig.doubletscoreVSnFeature.cxds(), paste0(outdir, sampleID, "_doublets_cxds.pdf"))
  to.pdf(fig.doubletscoreVSnFeature.bcds(), paste0(outdir, sampleID, "_doublets_bcds.pdf"))
  to.pdf(fig.doubletscoreVSnFeature.hybrid(), paste0(outdir, sampleID, "_doublets_hybrid.pdf"))
  to.pdf(fig.doubletscoreVSnFeature.nCalls(), paste0(outdir, sampleID, "_doublets_nCalls.pdf"))
  
  # filter the SCE object
  ## apply the new metadata to mySCE (i.e. # of methods that say filtering is required)
  nFilter <- sum(mytable$scds.nFiltCalls > 1)
  print(paste0("filtering ", nFilter, " potential doublets"))
  mySCE@colData$scds.nFiltCalls <- mytable$scds.nFiltCalls
  # NOTE: Though the following command looks wrong (extra ,) but it works (`SingleCellExperiment` manual, function `SCE-combine`, p12)
  mySCE <- subset(mySCE, , scds.nFiltCalls < 2) # remove if 2-3 of 3 tools call a doublet (doublet + 3000 sequences) 
  # return mySCE to main R environment 
  return(mySCE)
}

# ------------------------------------------------------------------
# PREDICT DOUBLETS
# ------------------------------------------------------------------

# input = seurat object
df.sce <- func_predictDoublets(df)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()