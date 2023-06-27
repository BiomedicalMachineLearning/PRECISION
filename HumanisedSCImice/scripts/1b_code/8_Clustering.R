# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th July 2020
# Title: 8_Clustering.R
# Goal: To cluster individual samples
# ------------------------------------------------------------------
# USAGE
# Usage: Rscript 8_Clustering.R {sampleID} {seuratObj} {outdir_figs} {outdir_obj} {res.override} {ndims}
# Where:
# ndims = # dimensions used for clustering (should be the same as used for UMAP) (use 50 unless you changed it in script 7b)
# res.override = no or a number - use if you know which resolution you want to use
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 6) {
  cat("ERROR: 6 arguments expected\n")
  cat("example: Rscript 8_Clustering.R {sampleID} {seuratObj} {outdir_figs} {outdir_obj} {res.override} {ndims}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
outdir <- args[3]
outdir_obj <- args[4]
res.override <- args[5]
ndims <- as.numeric(args[6])

# check that the res.override parameter is formatted correctly
if (is.numeric(res.override) | res.override == "no") {
  print(paste0("clustering resolution override? ... ", res.override))
} else {
  stop("the argument 'res.override' should be no or some value")
}

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}
if (endsWith(outdir_obj, "/") == FALSE) {
  outdir_obj <- paste0(outdir_obj, "/", sep="")
}

# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "ggplot2", "clustree")
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

# ------------------------------------------------------------------
# RUN CLUSTERING WITH USER-SPECIFIED PARAMETERS
# ------------------------------------------------------------------

if (is.numeric(res.override)) {
  df <- FindClusters(df, resolution = as.numeric(res.override))
  DimPlot(df, reduction = "umap") + ggtitle(paste0("res = ", res.override))
  ggsave(paste0(outdir, sampleID, "_UserChosenClust_res", res.override, ".pdf"))
}

# ------------------------------------------------------------------
# RUN CLUSTERING AND SELECT MOST STABLE PARAMETER
# ------------------------------------------------------------------

if (res.override == "no") {
  # perform clustering across a range of resolution values
  df <- FindNeighbors(df, reduction = "pca", dims = 1:ndims)
  df.2 <- FindClusters(df, resolution = 0)
  df.2 <- FindClusters(df.2, resolution = 0.1)
  df.2 <- FindClusters(df.2, resolution = 0.2)
  df.2 <- FindClusters(df.2, resolution = 0.3)
  df.2 <- FindClusters(df.2, resolution = 0.4)
  df.2 <- FindClusters(df.2, resolution = 0.5)
  df.2 <- FindClusters(df.2, resolution = 0.6)
  df.2 <- FindClusters(df.2, resolution = 0.7)
  df.2 <- FindClusters(df.2, resolution = 0.8)
  df.2 <- FindClusters(df.2, resolution = 0.9)
  df.2 <- FindClusters(df.2, resolution = 1)
  df.2 <- FindClusters(df.2, resolution = 1.2)
  df.2 <- FindClusters(df.2, resolution = 1.4)
  df.2 <- FindClusters(df.2, resolution = 1.6)
  
  # plot clustering results for the tested resolution values
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0") + ggtitle("candidate res = 0")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.1") + ggtitle("candidate res = 0.1")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.1.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.2") + ggtitle("candidate res = 0.2")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.2.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.3") + ggtitle("candidate res = 0.3")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.3.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.4") + ggtitle("candidate res = 0.4")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.4.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.5") + ggtitle("candidate res = 0.5")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.5.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.6") + ggtitle("candidate res = 0.6")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.6.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.7") + ggtitle("candidate res = 0.7")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.7.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.8") + ggtitle("candidate res = 0.8")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.8.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.0.9") + ggtitle("candidate res = 0.9")
  ggsave(paste0(outdir, sampleID, "_possClust_res0.9.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.1") + ggtitle("rcandidate res = 1")
  ggsave(paste0(outdir, sampleID, "_possClust_res1.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.1.2") + ggtitle("candidate res = 1.2")
  ggsave(paste0(outdir, sampleID, "_possClust_res1.2.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.1.4") + ggtitle("candidate res = 1.4")
  ggsave(paste0(outdir, sampleID, "_possClust_res1.4.pdf"))
  DimPlot(df.2, reduction = "umap", group.by = "RNA_snn_res.1.6") + ggtitle("candidate res = 1.6")
  ggsave(paste0(outdir, sampleID, "_possClust_res1.6.pdf"))
  
  # assess clustering results with clustree
  clust <- clustree(df.2, prefix = "RNA_snn_res.", node_colour = "sc3_stability")
  clust
  ggsave(plot = clust, file = paste0(outdir, sampleID, "_clustree.pdf"), height = 12, width = 12)
  # extract the stability values for the different resolutions
  stability <- clust$data[,c("RNA_snn_res.", "sc3_stability")]
  write.table(stability, file = paste0(outdir, "clustree_stability.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # select the most stable clustering value
  stability <- stability[stability$RNA_snn_res. %in% names(which(table(stability$RNA_snn_res.) > 1)), ]
  stability.ave <- aggregate(as.numeric(stability$sc3_stability), list(stability$RNA_snn_res.), mean)
  rownames(stability.ave) <- stability.ave$Group.1
  stability.ave$Group.1 <- NULL
  stability.ave.no0 <- stability.ave[2:nrow(stability.ave), , drop = FALSE]
  bestres <- as.numeric(rownames(stability.ave.no0)[which.max(stability.ave.no0$x)])
  stability.ave
  
  sink(paste0(outdir, sampleID, "_bestRes.txt"))
  cat(paste0("the chosen most stable resolution parameter is ", bestres))
  sink()
  
  # re-run the clustering for df with the most stable value
  rm(df.2)
  df <- FindClusters(df, resolution = bestres)
  DimPlot(df, reduction = "umap") + ggtitle(paste0("res = ", bestres))
  ggsave(paste0(outdir, sampleID, "_chosenClust_res", bestres, ".pdf"))
}

# ------------------------------------------------------------------
# SAVE THE R OBJECT WITH CLUSTERING
# ------------------------------------------------------------------

saveRDS(df, file = paste0(outdir_obj, sampleID, "clustered_Seurat.RDS"))
# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()