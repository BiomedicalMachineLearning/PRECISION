# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 14th July 2020
# Title: 11.7_ConfusionAndDynamicTree.R
# Goal: To take the Metacell object, replot the confusion matrix as a heatmap, and split the resulting hc object with the DynamicCutree package
# Usage: Rscript 11.7_ConfusionAndDynamicTree.R {sampleID} {MetacellDir} {seurat} {outdir}
# WHERE:
# MetacellDir = the top-level Metacell directory which contains out_fig and out_obj
# outdir = where to save the results of this script
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 11.7_ConfusionAndDynamicTree.R {sampleID} {MetacellDir} {seurat} {outdir}\n")
  quit()
}

sampleID <- args[1]
outdir_metacell_obj <- args[2]
df <- readRDS(args[3])
outdir <- args[4]

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir_metacell_obj, "/") == FALSE) {
  outdir_metacell_obj <- paste0(outdir_metacell_obj, "/", sep="")
}

if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

# load packages
bioc_packages <- c("metacell")
r_packages <- c("RColorBrewer", "ggplot2", "pheatmap", "dynamicTreeCut", "randomcoloR", "Seurat")
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
scfigs_init(outdir)

#' Change the ID so that I can re-run graphs without over-writing the old ones
mc_play = scdb_mc("test_mc_f")
scdb_add_mc("playing", mc_play)

# ------------------------------------------------------------------
# GET THE CONFUSION MATRIX DATA
# ------------------------------------------------------------------

graph_id <- "test_graph"
mc_id <- "test_mc_f"
cgraph = scdb_cgraph(graph_id)
max_deg = nrow(cgraph@edges)
confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg, 
                               ignore_mismatch = T)
r_confu = rowSums(confu)
c_confu = colSums(confu)
norm = r_confu %*% t(c_confu)
confu_n = confu/norm
confu_nodiag = confu_n
diag(confu_nodiag) = 0
confu_n = pmin(confu_n, max(confu_nodiag))
confu_n = pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n)))
epsilon = quantile(confu_n[confu_n != 0], 0.02)
dist <- as.dist(-log10(epsilon + confu_n))
hc = hclust(dist, "average")
distmat <- -log10(epsilon + confu_n)

# from the function for the confusion matrix
# I've broken down how to do the confusion matrix and super-metacell plots here
mc_order=hc$order 
confu_n = confu_n[mc_order, mc_order]
colnames(confu_n) = (1:ncol(confu_n))[mc_order]
rownames(confu_n) = (1:ncol(confu_n))[mc_order]

# ------------------------------------------------------------------
# PLOT THE BASIC DENDROGRAM
# ------------------------------------------------------------------

png(file = paste0(outdir, sampleID, "_dendrogram.png"), width = 1000)
plot(hc)
dev.off()

# ------------------------------------------------------------------
# RUN dynamicTreeCut
# ------------------------------------------------------------------

# to deal with excessively large clusters
mc_count <- length(unique(rownames(distmat)))
if (mc_count > 100) {
  clustsize <- mc_count/20
} else {
  clustsize <- 3
}
cut <- cutreeHybrid(dendro = hc, distM = distmat, minClusterSize = clustsize, deepSplit = 0)

# plot the labels
png(file = paste0(outdir, sampleID, "_dendrogram_CutGroups.png"), width = 1000)
plot(hc, labels = cut$labels)
dev.off()

# output the groups
groups <- as.data.frame(as.vector(cut$labels[hc$order]))
metacells <- as.data.frame(hc$order)
metagroups <- cbind(groups, metacells)
colnames(metagroups) <- c("groups", "metacells")
write.table(metagroups, file = paste0(outdir, sampleID, "_metagroups.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------
# PLOT THE HEATMAP
# ------------------------------------------------------------------

# get the annotations from above
rownames(metagroups) <- metagroups$metacells
metagroups$metacells <- NULL
metagroups$groups <- as.character(metagroups$groups)

# get the colours
uniquegroups <- unique(metagroups$groups)
# get the palette - groups > 34 need special treatment
if (length(uniquegroups) > 34) {
  ann_colors <- Seurat::DiscretePalette(36, palette = "polychrome")
  ann_colors <- append(ann_colors, Seurat::DiscretePalette(32, palette = "glasbey"))
  ann_colors <- append(ann_colors, Seurat::DiscretePalette(26, palette = "alphabet2"))
  ann_colors <- unique(ann_colors)
} else {
  ann_colors <- Seurat::DiscretePalette(length(uniquegroups)+2, palette = "polychrome") #+2 because i dont like the first 2 colours
}
ann_colors <- ann_colors[3:length(ann_colors)]
names(ann_colors) <- uniquegroups
ann_colors <- list(groups = ann_colors)

con <- apply(t(confu_n),2,rev) #rotate 90 degrees CCW
# with legend, free size
pheatmap(con, color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100),
         cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_row = metagroups, annotation_col = metagroups, annotation_colors = ann_colors, border_color = NA,
         filename = paste0(outdir, sampleID, "_heatmap.png"))
# without legend
pheatmap(con, color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100),
         cluster_cols = FALSE, cluster_rows = FALSE,
         legend = FALSE, annotation_legend = FALSE,
         height = 9, width = 9,
         annotation_row = metagroups, annotation_col = metagroups, annotation_colors = ann_colors, border_color = NA,
         filename = paste0(outdir, sampleID, "_heatmap_nolegend.png"))


# ------------------------------------------------------------------
# PLOT A NETWORK DIAGRAM WITH NODES COLOURED BY GROUP
# ------------------------------------------------------------------

#' # Make an annotation table containing the data
# first get the order of cells we need to colour
mc2d = scdb_mc2d("test_2dproj")
cellnames <- names(mc2d@sc_x)
mc = scdb_mc(mc2d@mc_id)

# get the node and edge data
fr = mc2d@graph$mc1
to = mc2d@graph$mc2


#' # Build the plot in ggplot2
# first need to merge all the NODE data into a table  
sorted_groups <- metagroups[ order(as.numeric(row.names(metagroups))), , drop = FALSE]
sorted_colours <- as.vector(ann_colors$groups[order(as.numeric(names(ann_colors$groups)))])
nodedata <- cbind(as.data.frame(mc2d@mc_x), as.data.frame(mc2d@mc_y), sorted_groups)
colnames(nodedata) <- c("nodes_x", "nodes_y", "group")

ggplot() +
  geom_point(aes(x = mc2d@sc_x, y = mc2d@sc_y), size = 1) + 
  geom_segment(aes(x = mc2d@mc_x[fr], y = mc2d@mc_y[fr], xend = mc2d@mc_x[to], yend = mc2d@mc_y[to]), col = "gray43") +
  geom_point(aes(x = nodedata$nodes_x, y = nodedata$nodes_y, color = nodedata$group), size = 6) +
  scale_color_manual(values = sorted_colours) +
  geom_text(aes(x = mc2d@mc_x, y = mc2d@mc_y, label = 1:length(mc2d@mc_x)), size = 3.5) +
  theme_void() +
  theme(legend.position = "bottom")
ggsave(filename = paste0(outdir, sampleID, "_network_ColourByGroup.png"), width = 10, height = 11)

# ------------------------------------------------------------------
# PLOT A UMAP WITH THESE GROUPS MAPPED
# ------------------------------------------------------------------

# get the metadata
conversion <- as.data.frame(mc_play@mc)
conversion$cells <- rownames(conversion)
colnames(conversion) <- c("metacells", "cells")
conversion <- merge(conversion, metagroups, by.x = 1, by.y = "row.names")
rownames(conversion) <- conversion$cells
conversion$cells <- NULL
conversion$dynamicTreeCut_groups <- conversion$groups
conversion$groups <- NULL
write.table(conversion, file = paste0(outdir, sampleID, "_metagroups_conversion.txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# add to the Seurat object
df <- Seurat::AddMetaData(df, metadata = conversion)
df <- subset(df, cells = names(which(!is.na(df$dynamicTreeCut_groups))))

# Make plots
plot1 <- DimPlot(df, group.by = "dynamicTreeCut_groups", cols = sorted_colours)
ggsave(plot = plot1, paste0(outdir, sampleID, "_DimPlot.jpg"), dpi = 200, width = 8, height = 7)

plot2 <- DimPlot(df, group.by = "dynamicTreeCut_groups", split.by = "time", cols = sorted_colours)
ggsave(plot = plot2, paste0(outdir, sampleID, "_DimPlot_time.jpg"), dpi = 200, width = 10, height = 3)

plot3 <- DimPlot(df, group.by = "dynamicTreeCut_groups", split.by = "state", cols = sorted_colours)
ggsave(plot = plot3, paste0(outdir, sampleID, "_DimPlot_state.jpg"), dpi = 200, width = 7, height = 3)

plot4 <- DimPlot(df, group.by = "dynamicTreeCut_groups", split.by = "timestate", cols = sorted_colours)
ggsave(plot = plot4, paste0(outdir, sampleID, "_DimPlot_timestate.jpg"), dpi = 200, width = 20, height = 3)

# Bar plots
# time
data_time <- as.data.frame(table(plot2$data[,3:4]))
colnames(data_time) <- c("group", "condition", "count")
ggplot(data_time, aes(fill=condition, y=count, x=group)) + 
  geom_bar(position="fill", stat="identity")
ggsave(paste0(outdir, sampleID, "_barplot_time.pdf"))

# state
data_state <- as.data.frame(table(plot3$data[,3:4]))
colnames(data_state) <- c("group", "condition", "count")
ggplot(data_state, aes(fill=condition, y=count, x=group)) + 
  geom_bar(position="fill", stat="identity")
ggsave(paste0(outdir, sampleID, "_barplot_state.pdf"))

# timestate
data_timestate <- as.data.frame(table(plot4$data[,3:4]))
colnames(data_timestate) <- c("group", "condition", "count")
ggplot(data_timestate, aes(fill=condition, y=count, x=group)) + 
  geom_bar(position="fill", stat="identity")
ggsave(paste0(outdir, sampleID, "_barplot_timestate.pdf"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()