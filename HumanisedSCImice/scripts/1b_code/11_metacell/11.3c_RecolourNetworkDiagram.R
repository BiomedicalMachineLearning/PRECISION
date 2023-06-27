# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 1st July 2020
# Title: 11.3c_RecolourNetworkDiagram.R
# Goal: To take the Metacell network diagram and colour the cells based on Seurat clusters.
#       Metacell generates a network graph showing connectivity between Metacells;
#       this plot is coloured by a user-provided set of markers. However, it may be more informative
#       to colour the cells by other metadata, such as cluster ID or experimental condition.
# Usage: Rscript 11.3c_RecolourNetworkDiagram.R {sampleID} {seuratObj} {MetacellDir} {recolour_outdir}
# WHERE:
# seuratObject = the subsetted Seurat object used before Metacell analysis
# MetacellDir = the top-level Metacell directory which contains out_fig and out_obj
# recolour_outdir = where to save the results of this script
# cntl + f for "# NOTE: change to necessary cluster ID" and add in the metadata column you want to colour by
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 11.3c_RecolourNetworkDiagram.R {sampleID} {seuratObj} {MetacellDir} {recolour_outdir}\n")
  quit()
}

sampleID <- args[1]
seurat <- args[2]
outdir_metacell <- args[3]
recolour_outdir <- args[4]

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir_metacell, "/") == FALSE) {
  outdir_metacell <- paste0(outdir_metacell, "/", sep="")
}

if (endsWith(recolour_outdir, "/") == FALSE) {
  recolour_outdir <- paste0(recolour_outdir, "/", sep="")
}

# Get some other variables
outdir_metacell_figs <- paste0(outdir_metacell, "/out_fig/")
outdir_metacell_obj <- paste0(outdir_metacell, "/out_obj/")
colour.table.path <- "/scripts/0_raw_data/r_colour_conversion"

# load packages
bioc_packages <- c("metacell")
r_packages <- c("RColorBrewer", "ggplot2", "pheatmap", "randomcoloR")
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

# ------------------------------------------------------------------
# RELOAD THE DATA
# ------------------------------------------------------------------

print(paste0("Log 1: Preparing the data for MC sample ", sampleID, " at ", date()))

#' # Reload the data
scdb_init(outdir_metacell_obj, force_reinit = T)
scfigs_init(outdir_metacell_figs)

#' Change the ID so that I can re-run graphs without over-writing the old ones
mc_play = scdb_mc("test_mc_f")
scdb_add_mc("playing", mc_play)

# ------------------------------------------------------------------
# PREPARE THE ANNOTATIONS
# ------------------------------------------------------------------

print(paste0("Log 2: Preparing the annotations for MC sample ", sampleID, " at ", date()))

#' # Make an annotation table containing the data
# first get the order of cells we need to colour
mc2d = scdb_mc2d("test_2dproj")
cellnames <- names(mc2d@sc_x)
mc = scdb_mc(mc2d@mc_id)
# NOTE: change to necessary cluster ID
# get the cellID-to-cluster conversion, and filter it so it is the same length and order as the cellID-to-metacell conversion
df <- readRDS(seurat)
clusterannotations <- df[["integrated_snn_res.0.3"]]
clusterannotations <-clusterannotations[cellnames, , drop = FALSE]

# get the node and edge data
fr = mc2d@graph$mc1
to = mc2d@graph$mc2

# ------------------------------------------------------------------
# BUILD THE GRAPH
# ------------------------------------------------------------------

print(paste0("Log 3: Building the graph for MC sample ", sampleID, " at ", date()))

#' # Build the plot in ggplot2
# first need to merge all the data into a table
networkdata <- cbind(as.data.frame(mc2d@sc_x), as.data.frame(mc2d@sc_y), clusterannotations)
colnames(networkdata) <- c("mc2d_xval", "mc2d_yval", "cluster")

ggplot() + 
  geom_point(aes(x=networkdata$mc2d_xval, y=networkdata$mc2d_yval, color=networkdata$cluster), size = 1) +
  geom_point(aes(x=networkdata$mc2d_xval, y=networkdata$mc2d_yval, color=networkdata$cluster), size = 1) +
  geom_segment(aes(x = mc2d@mc_x[fr], y = mc2d@mc_y[fr], xend = mc2d@mc_x[to], yend = mc2d@mc_y[to]), col = "gray43") +
  geom_point(aes(x = mc2d@mc_x, y = mc2d@mc_y), shape = 21, fill = mc@colors, color = "black", size = 6) +
  geom_text(aes(x = mc2d@mc_x, y = mc2d@mc_y, label = 1:length(mc2d@mc_x)), size = 3.5) +
  theme_void() +
  theme(legend.position = "bottom")
ggsave(filename = paste0(recolour_outdir, sampleID, "_networkdiagram.pdf"), width = 10, height = 10)

#' # or the same plot, but coloured by cell clusters
myclusts <- unique(networkdata$cluster)
colour.tab = read.delim(colour.table.path)
rownames(colour.tab) <- colour.tab$cluster
colour.tab$cluster <- NULL
colour.tab <- subset(colour.tab, rownames(colour.tab) %in% myclusts)
colour.tab <- colour.tab[order(as.numeric(row.names(colour.tab))), , drop = F]
my.cols <- colour.tab$hex
my.cols <- as.character(colour.tab$hex)
#droplevels(my.cols)

ggplot() + 
  geom_point(aes(x=networkdata$mc2d_xval, y=networkdata$mc2d_yval, color=networkdata$cluster), size = 1) +
  scale_color_manual(values = my.cols) + #i got these colours from the plot - need t
  geom_point(aes(x=networkdata$mc2d_xval, y=networkdata$mc2d_yval, color=networkdata$cluster), size = 1) +
  geom_segment(aes(x = mc2d@mc_x[fr], y = mc2d@mc_y[fr], xend = mc2d@mc_x[to], yend = mc2d@mc_y[to]), col = "gray43") +
  geom_point(aes(x = mc2d@mc_x, y = mc2d@mc_y), shape = 21, fill = mc@colors, color = "black", size = 6) +
  geom_text(aes(x = mc2d@mc_x, y = mc2d@mc_y, label = 1:length(mc2d@mc_x)), size = 3.5) +
  theme_void() +
  theme(legend.position = "bottom")
ggsave(filename = paste0(recolour_outdir, sampleID, "_networkdiagram_coloured.pdf"), width = 10, height = 10)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(recolour_outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()