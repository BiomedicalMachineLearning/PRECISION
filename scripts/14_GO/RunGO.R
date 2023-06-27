# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 28th May 2021
# Title: RunPathGO.R
# Goal: To find the enriched GO terms associated with the top N DE gene that change along a trajectory
# Usage: Rscript RunPathGO.R {clusterID} {all/top100} {DE list} {seuratObj} {outdir}
# ------------------------------------------------------------------
# USAGE
# clusterID = name for output files, e.g. the cell type (neut_all)
# all/top100 = how many gene to test? must be "all" or "top100"
# DE list = text file of DE gene, in format: {gene | direction | p_val_adj | cluster}
# seurat object = path to seurat file
# outdir = path to save the output
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

print("initialising analysis")


args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat("ERROR: 5 arguments expected\n")
  cat("example: Usage: Rscript RunPathGO.R {clusterID} {all/top100} {DE list} {seuratObj} {outdir}\n")
  quit()
}

nToTest <- args[2]
# check it's in the right format
if (nToTest != "all" & nToTest != "top100" & nToTest != "All" & nToTest != "Top100" & nToTest != "ALL" & nToTest != "TOP100") {
  cat("ERROR: arg2 should be either all or top100")
  quit()
}
clusterID <- args[1]
DE <- read.delim(args[3])
df <- readRDS(args[4])
outdir <- args[5]

# load packages
library(presto)
library(clusterProfiler)
library(ComplexHeatmap)
library(dplyr)
library(org.Hs.eg.db)
#library(org.Hs.eg.db)
library(tidyverse)
library(ggpubr)
library(stringr)
library(aplot)
library(patchwork)

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

# ------------------------------------------------------------------
# STEP 1: PREPARE DE gene, GENE UNIVERSE, ETC.
# ------------------------------------------------------------------

print("commencing step 1: Prepare DE gene, gene universe, etc.")

# # filter the DE list
# #DE <- mutate(DE, cluster_dir = paste0(cluster, "_", direction))
# if (nToTest == "all" | nToTest == "ALL" | nToTest == "All") {
#   DE <- DE %>% drop_na(p_val_adj) %>% filter(p_val_adj < 0.05) %>% group_by(cluster)
#   print(table(DE$cluster))
# } else if (nToTest == "top100" | nToTest == "Top100" | nToTest == "TOP100") {
#   DE <- DE %>% drop_na(p_val_adj) %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>% slice_min(order_by = avg_log2FC, n = 100)
#   print(table(DE$cluster))
# } else {
#   print("nToTest does not match expected options")
#   quit()
# }
# DE <- as.data.frame(DE)

# Get Entrez codes for DE gene
prep_IDs <- function(list_to_test) {
  geneList <- pull(DE %>% dplyr::filter(cluster == list_to_test), "gene")
  myIDs = mapIds(org.Hs.eg.db, column = "ENTREZID", keytype = "SYMBOL", keys = geneList)
  myIDs = myIDs[!(is.na(myIDs))]
  return(myIDs)
}
all_entrez <- lapply(unique(DE$cluster), prep_IDs)
names(all_entrez) <- unique(DE$cluster)
# Then prepare the gene universe
mycounts <- df@assays$RNA@data
# filter non-expressed
expgene <- names(which(Matrix::rowSums(mycounts) != 0))
expcounts <- mycounts[which(rownames(mycounts) %in% expgene),]
universe <- rownames(expcounts)
universe_entrez <- mapIds(org.Hs.eg.db, keys=universe, column="ENTREZID", keytype="SYMBOL")
universe_entrez <- universe_entrez[!(is.na(universe_entrez))]
# ------------------------------------------------------------------
# STEP 2: PERFORM GO ANALYSIS
# ------------------------------------------------------------------

print("commencing step 2: performing GO analysis")

# Run the GO analysis and filter to (1) remove giant clusters, (2) reduce simplicity
ck <- compareCluster(geneCluster = all_entrez, fun = "enrichGO", OrgDb = "org.Hs.eg.db", universe = universe_entrez, readable = TRUE, ont = "BP", pvalueCutoff = 0.01, pAdjustMethod = "BH")
ck_filt <- gsfilter(ck, max = 200)
ck_filt_simple <- clusterProfiler::simplify(ck_filt)

# save output as RDS, txt file, simple dotplot
## RDS
saveRDS(ck_filt_simple, paste0(outdir, clusterID, "_top100GO_FiltSimp.RDS"))
## txt file
tosave <- ck_filt_simple@compareClusterResult
tosave$GeneRatio <- gsub("\\/", "|", tosave$GeneRatio)
tosave$BgRatio <- gsub("\\/", "|", tosave$BgRatio)
write.table(tosave, file = paste0(outdir, clusterID, "_top100GO_FiltSimp.txt"), sep = "\t", quote = FALSE, col.names = NA)
## dotplot
pdf(paste0(outdir, clusterID, "_ClusterMarker_dotplot.pdf"))
dotplot(ck_filt_simple, font.size = 8, showCategory = 7)
dev.off()

# ------------------------------------------------------------------
# STEP 3a: DEFINE FUNCTION TO MAKE PLOTS
# ------------------------------------------------------------------

print("commencing step 3a: defining function to make plots")

MakeDotplots <- function(n) {
  # where n = number of gene to pick
  analysisID <- paste0("top", n)
  
  # get the top N of each group
  topN <- ck_filt_simple@compareClusterResult %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    dplyr::slice(1:n)
  topNuniq <- unique(topN$Description)
  topN_allhits <- ck_filt_simple@compareClusterResult[which(ck_filt_simple@compareClusterResult$Description %in% topNuniq),]
  
  # prepare the data table
  table <- topN_allhits[,c(1,3,4,5,7)]
  
  # add in the various ratios to control dot size
  table <- tidyr::separate(data = table, col = GeneRatio, into = c("k_GOIinSet", "n_GOI"), sep = "\\/")
  table <- tidyr::separate(data = table, col = BgRatio, into = c("M_Set", "N_background"), sep = "\\/")
  ## make numeric
  table$k_GOIinSet <- as.numeric(table$k_GOIinSet)
  table$n_GOI <- as.numeric(table$n_GOI)
  table$M_Set <- as.numeric(table$M_Set)
  table$N_background <- as.numeric(table$N_background)  
  ## make a ratio
  table <- mutate(table, GeneRatio = k_GOIinSet / n_GOI)
  table <- mutate(table, Ratio_kM = k_GOIinSet / M_Set)
  
  # for the sake of visualisation, convert p-values to log10 and invert the numbers of the downregulated hits
  table <- mutate(table, padjdir = -log10(p.adjust + 1e-300))
  
  # simplify the "traj_0_" cluster names to just the number to aid with the ordering/legend
  #table <- mutate(table, Path = as.numeric(word(Cluster, 2, sep = "_")))
  #table$Path <- as.factor(table$Path) # convert to factor data type
  #levels(table$Path) <- sort(as.numeric(levels(table$Path)))
  
  # work out which gene to plot and order them by path ID
  markers <- table$Description %>% unique()
  table$Description <- factor(table$Description, levels = rev(topNuniq))
  
  # generate the plots
  ## size is kmratio, colour is up/downreg
  # dotplot_kmp_topN <- table %>% filter(Description %in% markers) %>% 
  #   ggplot(aes(x=Path, y = Description, color = padjdir, size = Ratio_kM)) + 
  #   geom_point() +
  #   scale_color_gradientn(colours = c("blue", "lightskyblue", "white", "orange", "red"),
  #                         values = scales::rescale(c(min(table$padjdir),
  #                                                    ((abs(min(table$padjdir)))/2)*-1, 0, max(table$padjdir)/2, max(table$padjdir)))) +
  #   scale_y_discrete(position = "right", labels = function(x) str_trunc(x, width = 80, side = "center"))  +
  # 
  #   theme_light(base_size = 10)
  
  ## size is gene ratio, colour is up/downreg - this is the best one
  dotplot_grp_topN <- table %>% filter(Description %in% markers) %>% 
    ggplot(aes(x=Cluster, y = Description, color = padjdir, size = GeneRatio)) + 
    geom_point(aes(fill=padjdir), colour="grey50",pch=21) +
    scale_fill_gradientn(colours = c("white", "gold", "orange", "red"),
                         values = scales::rescale(c(0, max(table$padjdir)/3, (max(table$padjdir)/3)*2, max(table$padjdir)))) +
    scale_y_discrete(position = "right", labels = function(x) str_trunc(x, width = 80, side = "center")) +
    theme_light(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ## size is gene ratio but no differential colouring for up/downreg
  # dotplot_ratios_topN <- table %>% filter(Description %in% markers) %>% 
  #   ggplot(aes(x=Path, y = Description, color = Ratio_kM, size = GeneRatio)) + 
  #   geom_point() +
  #   scale_color_gradient(low = "gold", high = "red") +
  #   scale_y_discrete(position = "right", labels = function(x) str_trunc(x, width = 80, side = "center")) +
  #   theme_light(base_size = 10)
  
  # save the plots
  # labels <- table %>% filter(Description %in% markers) %>% 
  # ggplot(aes(x=Path, y = 1, fill = Path)) +
  # geom_tile() + 
  # ggmap::theme_nothing() +
  # xlim2(dotplot_ratios_topN)
  
  #legend <- plot_grid(get_legend(labels + theme(legend.position="bottom")))
  
  # export the plots
  # pdf(paste0(outdir, clusterID, "_dotplot_generatios_", analysisID, ".pdf"), width = 8, height = 10)
  # #labels + dotplot_ratios + plot_layout(ncol = 1, heights = c(0.9, 20))
  # dotplot_ratios_topN
  # dev.off()
  # 
  # pdf(paste0(outdir, clusterID, "_dotplot_kmratio+p_val_adj_", analysisID, ".pdf"), width = 8, height = 10)
  # #labels + dotplot_kmp + plot_layout(ncol = 1, heights = c(0.9, 20))
  # dotplot_kmp_topN
  # dev.off()
  
  pdf(paste0(outdir, "dotplot_generatio+p_val_adj_", analysisID, ".pdf"), width = 8, height = 10)
  #labels + dotplot_grp + plot_layout(ncol = 1, heights = c(0.9, 20))
  print(dotplot_grp_topN)
  dev.off()
}

# ------------------------------------------------------------------
# STEP 3b: MAKE PLOTS
# ------------------------------------------------------------------

print("commencing step 3b: making plots")

MakeDotplots(10)
MakeDotplots(20)

# ------------------------------------------------------------------
# STEP 4: CALCULATE JACCARD DISTANCE BETWEEN GENE LISTS
# ------------------------------------------------------------------

# print("commencing step 4: calculating Jaccard distance between gene lists")
# 
# # calculate pairwise jaccard distance
# prep_gene_IDs <- function(list_to_test) {
#   geneList <- pull(DE %>% dplyr::filter(cluster == list_to_test), "gene")
#   return(geneList)
# }
# 
# all_gene_IDs <- sapply(unique(DE$cluster), prep_gene_IDs)
# 
# dist <- unlist(lapply(combn(all_gene_IDs, 2, simplify = FALSE), function(x) {
#   length(intersect(x[[1]], x[[2]]))/length(union(x[[1]], x[[2]])) }))
# dist <- as.data.frame(cbind(t(combn(names(all_gene_IDs),2)), dist))
# dist$dist <- as.numeric(dist$dist)
# dist$dist <- round(dist$dist, digits = 3)
# dist$dist <- gsub(0.000, 0, dist$dist)
# dist <- as.data.frame(pivot_wider(dist, names_from = V2, values_from = dist ))
# rownames(dist) <- dist$V1
# dist$V1 <- NULL
# colnames(dist) <- gsub("traj_", "", colnames(dist))
# colnames(dist) <- gsub("de_", "", colnames(dist))
# rownames(dist) <- gsub("traj_", "", rownames(dist))
# rownames(dist) <- gsub("de_", "", rownames(dist))  
# write.table(dist, file = paste0(outdir, clusterID, "_pathDistances.txt"), sep = "\t", quote = FALSE, col.names = NA)

# ------------------------------------------------------------------
# MAKE A HEATMAP TO SHOW GO TERM OVERLAP
# ------------------------------------------------------------------

print("commencing step 5: calculating GO overlap heatmap")

# convert the gene into a table
genehits <- tosave[,c("Cluster", "Description", "geneID")]
s <- strsplit(genehits$geneID, split = "/")
genehits <- data.frame(Cluster = rep(genehits$Cluster, sapply(s, length)),
                       Description = rep(genehits$Description, sapply(s, length)),
                       geneID = unlist(s))
# add in the DE information
genehits <- mutate(genehits, cluster_gene = paste0(Cluster, "_", geneID))
DE_mod <-  mutate(DE, cluster_gene = paste0(cluster, "_", gene))
DE_mod <- DE_mod[DE_mod$cluster_gene %in% genehits$cluster_gene,]
DE_mod <- DE_mod[,c("cluster_gene"), drop = FALSE]
genehits <- left_join(x = genehits, y = DE_mod, by = "cluster_gene")
#genehits <- mutate(genehits, direction2 = ifelse(direction == "up", 1, -1))
#genehits <- mutate(genehits, pdir = -log10(p_val_adj + 1e-300))
write.table(genehits, file = paste0(outdir, "genehits.txt"), sep = "\t", quote = FALSE, col.names = NA)

# define a function to make the heatmap
PlotGOOverlap <- function(path_to_test) {
  genehits_traj <- genehits[genehits$Cluster == path_to_test,]
  genehits_traj$on <- 1
  my_matrix <- reshape2::acast(genehits_traj, Description~geneID, value.var = "on")
  my_matrix[is.na(my_matrix)] <- 0
  # plot the heatmap
  #minval <- min(min(my_matrix), -0.01) #whatever is smaller, the min-val or -0.01
  #maxval <- max(max(my_matrix), 0.01) #whatever is smaller, the min-val or -0.01
  col_fun = circlize::colorRamp2(c(0, 1), c("white", "red"))
  ht <- ComplexHeatmap::Heatmap(my_matrix,
                                row_names_gp = gpar(fontsize = 7),
                                column_names_gp = gpar(fontsize = 7),
                                rect_gp = gpar(col = "black", lwd = 0.1),
                                border_gp = gpar(col = "black", lty = 1),
                                use_raster = TRUE, raster_quality = 1,
                                col = col_fun)
  pdf(paste0(outdir, path_to_test, "_GOgeneheatmap.pdf"))
  draw(ht, heatmap_legend_side = "bottom")
  dev.off()
}

lapply(unique(genehits$Cluster), PlotGOOverlap)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

print("finished!")

sink(paste0(outdir, clusterID, "_sessionInfo.txt"))
date()
sessionInfo()
sink()