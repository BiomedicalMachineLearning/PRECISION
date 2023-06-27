# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 9th July 2020
# Title: 11.6_GeneEnrichmentPlots.R
# Goal: To plot lfp enrichment of genes of interest in Metacell network diagram
# Usage: Rscript 11.6_GeneEnrichmentPlots.R {sampleID} {MetacellDir_objects} {outdir} {genesOfInterest}
# WHERE:
# MetacellDir_objects is a Metacell directory (e.g. out_obj) containing the Rda objects
# genesOfInterest is a list of genes, separated by a comma: Tnf,Fas,Icam1
# outdir is the desired location of the output
# Note, this won't print the legends because it makes it hard to fit multiple plots on one page. If you want legends, change "show_legend = F" to T
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 11.6_GeneEnrichmentPlots.R {sampleID} {MetacellDir_objects} {genesOfInterest} {outdir}\n")
  quit()
}

sampleID <- args[1]
MetacellDir_objects <- args[2]
outdir <- args[3]
genesOfInterest <- args[4]
genes <- unique(strsplit(genesOfInterest, ",")[[1]])

# load packages
bioc_packages <- c("metacell")
r_packages <- c("ggplot2")
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
# RELOAD THE DATA
# ------------------------------------------------------------------

print(paste0("Log 1: Preparing the data for MC sample ", sampleID, " at ", date()))

#' # Reload the data
scdb_init(MetacellDir_objects, force_reinit = T)
scfigs_init(outdir)

#' Change the ID so that I can re-run graphs without over-writing the old ones
mc_play = scdb_mc("test_mc_f")
scdb_add_mc("playing", mc_play)

# ------------------------------------------------------------------
# PLOT GENES
# ------------------------------------------------------------------

# granulocytes
for (i in 1:length(genes)) {
  if (genes[i] %in% rownames(mc_play@mc_fp)) {
    mcell_mc2d_plot_gene(mc2d_id = "test_2dproj", gene = genes[i], show_mc_ids = T, show_legend = F, neto_points = F)
  } else {
    print(paste0(genes[i], " not found"))
  }}

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()