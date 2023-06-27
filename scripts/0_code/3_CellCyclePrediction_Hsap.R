# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: 3_CellCyclePrediction.R
# Goal: To predict cell types for scRNASeq data
# Usage: Rscript 3_CellCyclePrediction.R {sampleID} {seuratObj} {outdir}
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
#bioc_packages <- c("DESeq2", "pheatmap", "apeglm", "genefilter")
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
# DEFINE FUNCTION
# ------------------------------------------------------------------

to.pdf <- function(expr, filename, ...) {
  pdf(filename, ...)
  on.exit(dev.off())
  print(eval.parent(substitute(expr)))
}

func_predictCellCycle <- function(seuratObj, myspecies="mouse"){
  # USAGE: mySCE <- func_predictCellCycle(seuratObj, "mouse")
  # OUTPUT: a seuratObj object with S/G2M-phase scores and cell stage (G1, S, G2M) calls
  
  # specify the gene set used for Cell Cycle Scoring (human or mouse)
  if (identical(myspecies, "mouse")) {
    load("/scripts/0_raw_data/mouse.cc.genes.Rdata")
    geneset <- mouse.cc.genes
  } else if (identical(myspecies, "human")) {
    geneset <- cc.genes.updated.2019
  } else {
    stop("The 'species' argument must be mouse or human")
  }
  
  # make a Seurat object, normalise, run prediction
  seuratObj <- NormalizeData(seuratObj,
                             normalization.method = "LogNormalize",
                             scale.factor = 10000)
  seuratObj <- CellCycleScoring(seuratObj,
                                s.features = geneset$s.genes,
                                g2m.features = geneset$g2m.genes,
                                set.ident = TRUE)
  
  # define some graph functions which will be run with `to.pdf` later
  fig.cellcycle.bar <- function() {
    myscale <- round(max(table(seuratObj$Phase)), -3) #scale
    mybar <- barplot(table(seuratObj$Phase),
                     ylim = (c(0, myscale)),
                     main = paste0("Cell Phases in ", sampleID),
                     xlab = "cell phase",
                     ylab = "# cells", 
                     col = "white")
    text(mybar,
         table(seuratObj$Phase)+100,
         paste("n: ", table(seuratObj$Phase), sep=""), cex = 1) 
  }
  
  fig.cellcycle.pie <- function() {
    pie(table(seuratObj$Phase),
        labels = table(seuratObj$Phase),
        col = c("bisque", "cornflowerblue", "cadetblue2"),
        main = paste0("Cell phases in ", sampleID))
    legend("topright", c("G1", "G2M", "S"), cex = 0.8, fill = c("bisque", "cornflowerblue", "cadetblue2"))
  }
  
  # Run the figure functions and save graphs as PDFs
  to.pdf(fig.cellcycle.bar(), paste0(outdir, sampleID, "_CellCycle_bar.pdf"))
  to.pdf(fig.cellcycle.pie(), paste0(outdir, sampleID, "_CellCycle_pie.pdf"))
  
  # return the updated seuratObj
  return(seuratObj)
}

# ------------------------------------------------------------------
# RUN FUNCTION
# ------------------------------------------------------------------

df <- func_predictCellCycle(df, "human")

# ------------------------------------------------------------------
# SAVE THE NECESSARY OUTPUT
# ------------------------------------------------------------------

S.Score <- df$S.Score
G2M.Score <- df$G2M.Score
Phase <- df$Phase

write.table(S.Score, file = paste0(outdir, sampleID, "_Sscore.txt"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
write.table(G2M.Score, file = paste0(outdir, sampleID, "_G2MScore.txt"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)
write.table(Phase, file = paste0(outdir, sampleID, "_Phase.txt"), sep = "\t", quote = FALSE, col.names = NA, row.names = TRUE)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()