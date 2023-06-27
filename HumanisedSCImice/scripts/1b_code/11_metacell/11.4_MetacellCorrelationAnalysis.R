# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 1st July 2020
# Title: 11.4_MetacellCorrelationAnalysis.R
# Goal: To run my custom correlation analysis on R output
# Usage: Rscript 11.4_MetacellCorrelationAnalysis.R {sampleID} {seuratObj} {MetacellDir} {correlation_outdir}
# WHERE:
# seuratObject = the subsetted Seurat object used before Metacell analysis
# MetacellDir = the top-level Metacell directory which contains out_fig and out_obj
# correlation_outdir = where to save the results of the correlation analysis
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 11.4_MetacellCorrelationAnalysis.R {sampleID} {seuratObj} {MetacellDir} {correlation_outdir}\n")
  quit()
}

sampleID <- args[1]
seurat <- readRDS(args[2])
outdir_metacell <- args[3]
outdir_corr <- args[4]

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir_metacell, "/") == FALSE) {
  outdir_metacell <- paste0(outdir_metacell, "/", sep="")
}

if (endsWith(outdir_corr, "/") == FALSE) {
  outdir_corr <- paste0(outdir_corr, "/", sep="")
}

# Get some other variables
outdir_metacell_figs <- paste0(outdir_metacell, "/out_fig/")
outdir_metacell_obj <- paste0(outdir_metacell, "/out_obj/")
metacellAnnotations <- paste0(outdir_metacell, "/out_fig/metacell_metadataGroups.txt")

# load packages
bioc_packages <- c("metacell")
r_packages <- c("Seurat", "pheatmap", "ggplot2", "RColorBrewer", "randomcoloR", "dplyr")
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
# PREPARE THE DATA
# ------------------------------------------------------------------

print(paste0("Log 1: Preparing the data for MC sample ", sampleID, " at ", date()))

#' # Reload the Metacell run
scdb_init(outdir_metacell_obj, force_reinit = T)
scfigs_init(outdir_metacell_figs)

#' Change the ID so that I can re-run graphs without over-writing the old ones
mc_play = scdb_mc("test_mc_f")
scdb_add_mc("playing", mc_play)

#' # add in the metacell annotations into seurat
ann <- read.delim(metacellAnnotations, header = FALSE, row.names = 1, col.names = c("cell", "metacell", "metacell_group_colour", "metacell_group"))
head(ann)

# add into Seurat data
seurat <- AddMetaData(seurat, metadata = ann)
# filter to only keep cells with a metacell ID
cellIDs <- rownames(na.omit(seurat[["metacell"]]))
seurat <- subset(seurat, cells = cellIDs)
# visualise
DimPlot(seurat, reduction = "umap", group.by = "metacell") + 
  xlim(-10,10) +
  ylim(-10,10) +
  ggtitle("Cells coloured by metacell ID")
ggsave(paste0(outdir_corr, sampleID, "_SeuratWithMetacellData.pdf"))

# ------------------------------------------------------------------
# CORRELATION ANALYSIS
# ------------------------------------------------------------------

print(paste0("Log 2: Performing correlation analysis for MC sample ", sampleID, " at ", date()))

#' # Perform correlation analysis
Idents(object = seurat) <- "metacell"
metacelllist <- sort(as.numeric(levels(Idents(seurat))))
av.exp <- AverageExpression(seurat)$RNA
cor.exp <- as.data.frame(cor(av.exp, method = c("pearson")))
cor.exp$x <- rownames(cor.exp)
cor.df <- tidyr::gather(data = cor.exp, y, correlation, min(metacelllist):max(metacelllist))
# this step is unnecessary, you could just use the cor.exp output...
cor.mat <- reshape2::dcast(cor.df, x ~ y)
rownames(cor.mat) <- cor.mat$x
cor.mat$x <- NULL

# re-order the correlation matrix by the order shown in Metacell's confusion matrix
mc_hc = mcell_mc_hclust_confu(mc_id="playing", graph_id="test_graph")
mc_order = as.character(mc_hc$order)
cor.mat = cor.mat[mc_order, mc_order]

#' Check the min and max correlation values - these plots assume it's 0-1
## --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(paste0("min correlation value is ", min(cor.mat)))
print(paste0("max correlation value is ", max(cor.mat)))

# ------------------------------------------------------------------
# HEATMAP
# ------------------------------------------------------------------

print(paste0("Log 3: Preparing to make heatmap for MC sample ", sampleID, " at ", date()))

#' # Before making the heatmap, get some cluster annotations ready
# first get the cellID-to-metacell conversion
metacell_list <- as.data.frame(mc_play@mc)
length(unique(metacell_list[,1]))

# NOTE: change the cluster name as appropriate
# then get the cellID-to-cluster conversion, and filter it so it is the same length and order as the cellID-to-metacell conversion
clusterannotations <- seurat[["integrated_snn_res.0.3"]]
clusterannotations <-clusterannotations[rownames(metacell_list), , drop = FALSE]

# now merge the two together
metaclusterannotations <- cbind(metacell_list, clusterannotations)
rownames(metaclusterannotations) <- NULL
names(metaclusterannotations)[1] <- "metacell"

#' There are multiple clusters per metacell in some cases, so append the unique cluster sets
# NOTE: change the cluster name as appropriate
metaclusterannotations <- as.data.frame(metaclusterannotations %>% group_by(metacell) %>% summarise(integrated_snn_res.0.3 = paste(unique(sort(integrated_snn_res.0.3)), collapse = "-")))
rownames(metaclusterannotations) <- metaclusterannotations$metacell
metaclusterannotations$metacell <- NULL
metaclusterannotations <- metaclusterannotations[mc_order, , drop = FALSE]

#In this case, there's only two clusters so it doesn't matter that there's some "multiple cluster" metacells (here, cluster 1 AND cluster 6), but I'll plot it anyway
metaclusterannotations_2 <- metaclusterannotations
metaclusterannotations_2$integrated_snn_res.0.3 <- gsub(".*-.*", "multiple", metaclusterannotations_2$integrated_snn_res.0.3)
head(metaclusterannotations_2)


#' # Now we can plot the correlation analysis
# First just a basic plot
pheatmap(cor.mat, cluster_rows = FALSE, cluster_cols = FALSE, filename = paste0(outdir_corr, sampleID, "_CorrelationHeat_basic.pdf"))

# Now we'll make it pretty
palette <- colorRampPalette(brewer.pal(9,"Purples"))(100)
breaksList = seq(0, 1, by = 0.01)
pheatmap(cor.mat, cluster_rows = FALSE, cluster_cols = FALSE, color = palette, annotation_row = metaclusterannotations_2, annotation_col = metaclusterannotations_2, legend = TRUE, annotation_legend = TRUE, fontsize_row = 6, fontsize_col = 6, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA, cellwidth = 5, cellheight = 5, breaks = breaksList, width = 15, height = 15, filename = paste0(outdir_corr, sampleID, "_CorrelationHeat.pdf"))

# Now we'll make it pretty, with a different colour scheme
palette <- rev(colorRampPalette(brewer.pal(9,"PRGn"))(100))
breaksList = seq(0, 1, by = 0.01)
pheatmap(cor.mat, cluster_rows = FALSE, cluster_cols = FALSE, color = palette, annotation_row = metaclusterannotations_2, annotation_col = metaclusterannotations_2, legend = TRUE, annotation_legend = TRUE, fontsize_row = 6, fontsize_col = 6, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA, cellwidth = 5, cellheight = 5, breaks = breaksList, width = 15, height = 15, filename = paste0(outdir_corr, sampleID, "_CorrelationHeatv2.pdf"))

# Now we'll make do the same thing, but clustered by rows and columns
palette <- colorRampPalette(brewer.pal(9,"Purples"))(100)
breaksList = seq(0, 1, by = 0.01)
pheatmap(cor.mat, cluster_rows = TRUE, cluster_cols = TRUE, color = palette, annotation_row = metaclusterannotations_2, annotation_col = metaclusterannotations_2, legend = TRUE, annotation_legend = TRUE, fontsize_row = 6, fontsize_col = 6, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA, cellwidth = 5, cellheight = 5, breaks = breaksList, width = 15, height = 15, filename = paste0(outdir_corr, sampleID, "_CorrelationHeat_reclustered.pdf"))

# Now we'll make do the same thing, but clustered by rows and columns
palette <- rev(colorRampPalette(brewer.pal(9,"PRGn"))(100))
breaksList = seq(0, 1, by = 0.01)
pheatmap(cor.mat, cluster_rows = TRUE, cluster_cols = TRUE, color = palette, annotation_row = metaclusterannotations_2, annotation_col = metaclusterannotations_2, legend = TRUE, annotation_legend = TRUE, fontsize_row = 6, fontsize_col = 6, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA, cellwidth = 5, cellheight = 5, breaks = breaksList, width = 15, height = 15, filename = paste0(outdir_corr, sampleID, "_CorrelationHeat_reclusteredv2.pdf"))

# ------------------------------------------------------------------
# EXPORT RESULTS
# ------------------------------------------------------------------

write.table(av.exp, file = paste0(outdir_corr, sampleID, "_Correlation_aveExp.txt"), quote = FALSE, sep = "\t")
write.table(cor.df, file = paste0(outdir_corr, sampleID, "_Correlation_corValues.txt"), quote = FALSE, sep = "\t", row.names = FALSE)
# export order of heatmaps
## for the clustered-in-R heatmap, first save the data behind the plot to a variable
myheatmap <- pheatmap(cor.mat, cluster_rows = TRUE, cluster_cols = TRUE, color = palette, annotation_row = metaclusterannotations_2, annotation_col = metaclusterannotations_2, legend = TRUE, annotation_legend = TRUE, fontsize_row = 6, fontsize_col = 6, annotation_names_row = FALSE, annotation_names_col = FALSE, border_color = NA, cellwidth = 2, cellheight = 2, breaks = breaksList)
## get the rownames in the right order 
hmorder <- as.data.frame(myheatmap$tree_row$labels[myheatmap$tree_row$order])
write.table(hmorder, file = paste0(outdir_corr, sampleID, "_CorrelationHeat_reclustered_order.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(as.data.frame(mc_order), file = paste0(outdir_corr, sampleID, "_CorrelationHeat_MCorder.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir_corr, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()