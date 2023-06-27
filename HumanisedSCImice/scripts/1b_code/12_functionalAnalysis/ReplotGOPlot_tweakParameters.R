# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 28th May 2021
# Title: RunPathGO.R
# Goal: If you've run RunGO.R and the plots are too small/big, this code allows you to tweak them a bit
# Usage: Rscript ReplotGOPlot_tweakParameters.R {sampleID} {ck_filt_simple} {outdir} {nTerms} {width} {height}
# ------------------------------------------------------------------
# USAGE
# sampleID = name for output files, e.g. the cell type (neut_all)
# ck_filt_simple = the GO output from ClusterProfiler
# outdir = path to save the output
# nTerms = how many GO terms to plot
# width = how wide the plot should be
# height = how tall the plot should be
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

print("initialising analysis")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  cat("ERROR: 3 arguments expected\n")
  cat("example: Rscript ReplotGOPlot_tweakParameters.R {sampleID} {ck_filt_simple} {outdir}")
  quit()
}

sampleID <- args[1]
ck_filt_simple <- readRDS(args[2])
outdir <- args[3]
nTerms <- args[4]
myWidth <- args[5]
myHeight <- args[6]

# load packages
library(clusterProfiler)
library(dplyr)
library(org.Mm.eg.db)
library(ggplot2)
library(ComplexHeatmap)
library(stringr)

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}

# ------------------------------------------------------------------
# STEP 3a: DEFINE FUNCTION TO MAKE PLOTS
# ------------------------------------------------------------------

print("commencing step 3a: defining function to make plots")

MakeDotplots <- function(n, myWidth, myHeight) {
  # where n = number of gene to pick
  # where myWidth = desired width of saved ggplot
  # where myHeight = desired height of saved ggplot
  analysisID <- paste0("top", n)
  
  # get the top N of each group
  topN <- ck_filt_simple@compareClusterResult %>%
    group_by(Cluster) %>%
    arrange(p.adjust) %>%
    dplyr::slice(1:n, with_ties = FALSE)
  topN <- mutate(topN, ID_desc = paste0(ID, ": ", Description))
  
  topNuniq <- unique(topN$Description)
  topNuniq_IDdesc <- unique(topN$ID_desc)
  
  topN_allhits <- ck_filt_simple@compareClusterResult[which(ck_filt_simple@compareClusterResult$Description %in% topNuniq),]
  
  # prepare the data table
  table <- topN_allhits[,c(1,2,3,4,5,7)]
  
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
  
  table <- mutate(table, ID_desc = paste0(ID, ": ", Description))
  
  table$Description <- factor(table$Description, levels = rev(topNuniq))
  table$ID_desc <- factor(table$ID_desc, levels = rev(topNuniq_IDdesc))
  
  # generate the plots
  ## size is gene ratio, colour is up/downreg - this is the best one
  dotplot_grp_topN <- table %>% filter(Description %in% markers) %>% 
    ggplot(aes(x=Cluster, y = Description, color = padjdir, size = GeneRatio)) + 
    geom_point(aes(fill=padjdir), colour="grey50",pch=21) +
    scale_fill_gradientn(colours = c("white", "gold", "orange", "red"),
                         values = scales::rescale(c(0, max(table$padjdir)/3, (max(table$padjdir)/3)*2, max(table$padjdir)))) +
    scale_y_discrete(position = "right", labels = function(x) stringr::str_trunc(x, width = 80, side = "center")) +
    theme_light(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf(paste0(outdir, "dotplot_generatio+p_val_adj_", analysisID, "_REPLOT.pdf"), width = myWidth, height = myHeight)
  #labels + dotplot_grp + plot_layout(ncol = 1, heights = c(0.9, 20))
  print(dotplot_grp_topN)
  dev.off()
  
  ## same as above but include the GO ID in the name
  dotplot_grp_topN_withGOLabels <- table %>% filter(Description %in% markers) %>% 
    ggplot(aes(x=Cluster, y = ID_desc, color = padjdir, size = GeneRatio)) + 
    geom_point(aes(fill=padjdir), colour="grey50",pch=21) +
    scale_fill_gradientn(colours = c("white", "gold", "orange", "red"),
                         values = scales::rescale(c(0, max(table$padjdir)/3, (max(table$padjdir)/3)*2, max(table$padjdir)))) +
    scale_y_discrete(position = "right", labels = function(x) stringr::str_trunc(x, width = 80, side = "center")) +
    theme_light(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  pdf(paste0(outdir, "dotplot_generatio+p_val_adj_", analysisID, "_WithGOLabels_REPLOT.pdf"), width = (myWidth + 1), height = (myHeight + 0.5)) #extra width to allow for GO ID; idk why height varies
  #labels + dotplot_grp + plot_layout(ncol = 1, heights = c(0.9, 20))
  print(dotplot_grp_topN_withGOLabels)
  dev.off()
}

# ------------------------------------------------------------------
# STEP 3b: MAKE PLOTS
# ------------------------------------------------------------------

print("commencing step 3b: making plots")

# MakeDotplots(nGenes, Width, Height)
MakeDotplots(3, 5.5, 3.5)
MakeDotplots(10, 6, 7.5)
MakeDotplots(20, 8, 10)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

print("finished!")

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
date()
sessionInfo()
sink()