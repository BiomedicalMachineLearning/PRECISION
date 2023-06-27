# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: 7_seuratProcessing.R
# Goal: Perform Seurat scaling, PCA, TSNE and UMAP
# Usage: Rscript 7_seuratProcessing.R {sampleID} {seuratObj} {outdir_figs} {outdir_obj}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 7_seuratProcessing.R {sampleID} {seuratObj} {outdir_figs} {outdir_obj}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
outdir <- args[3] #for figures
outdir_obj <- args[4] #for objects

# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "ggplot2", "dplyr", "tibble")
## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
  }
}
# # function to load bioconductor packages
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

if (endsWith(outdir_obj, "/") == FALSE) {
  outdir_obj <- paste0(outdir_obj, "/", sep="")
}

# ------------------------------------------------------------------
# DEFINE FUNCTIONS
# ------------------------------------------------------------------

# -------------------
# function: PDFs
# -------------------

# function to help write cleaner code for generating images
# from https://nicercode.github.io/blog/2013-07-09-figure-functions/
# you must put your figure-generating code into a function
# USAGE: to.pdf(some.fig.making.function(), "figs/trend.pdf", width=6, height=4)
to.pdf <- function(expr, filename, ...) {
  pdf(filename, ...)
  on.exit(dev.off())
  print(eval.parent(substitute(expr)))
}

# -------------------------------------
# function: VARIABLE GENES + SCALE
# -------------------------------------

# function to find variable features and scale data using Seurat
func_ScaleData <- function(seuratObj) {
  # USAGE: seuratObj <- func_ScaleData(seuratObj)
  # OUTPUT: a Seurat object with scaled normalised counts
  
  # find variable features
  seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(seuratObj), 10)
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(seuratObj)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  ggsave(plot = plot1, filename = paste0(outdir, sampleID, "_variableFeatures.pdf"))
  ggsave(plot = plot2, filename = paste0(outdir, sampleID, "_variableFeatures_labelled.pdf"))
  
  # perform scaling
  seuratObj <- ScaleData(seuratObj)
  return(seuratObj)
}

# ------------------
# function: PCA
# ------------------

# function to perform PCA (linear dimensionality reduction)
func_runPCA <- function(seuratObj, runJackstraw = "TRUE") {
  # USAGE: seuratObj <- func_runPCA(seuratObj, runJackstraw = "TRUE" or "FALSE")
  # OUTPUT: a Seurat object with PCA run
  
  # Run PCA
  seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj), npcs = 50)
  
  # List top 10 genes associated with each PC
  sink(paste0(outdir, sampleID, "_PCgenes.txt"))
  print(seuratObj[["pca"]], dims = 1:50, nfeatures = 10)
  sink()
  
  # Loadings of the top 5 PCs
  VizDimLoadings(seuratObj, dims = 1:2, reduction = "pca", balanced = TRUE, nfeatures = 40)  
  ggsave(paste0(outdir, sampleID, "_top20UpDowngenes_PCs.pdf"))
  
  # Basic PCA plot
  DimPlot(seuratObj, reduction = "pca")
  ggsave(paste0(outdir, sampleID, "_PCA_raw.pdf"))
  
  # PCA heatmap
  pdf(paste0(outdir, sampleID, "_PCA_heat.pdf"), width = 12, height = 12)
  DimHeatmap(seuratObj, dims = 1:10, cells = 500, balanced = TRUE)
  dev.off()
  
  # calculate variance explained by each PC
  total_variance <- seuratObj@reductions$pca@misc$total.variance
  eigValues <- (seuratObj[["pca"]]@stdev)^2
  varExplained <- eigValues / total_variance
  varExplained.cum <- cumsum(varExplained)
  
  ### how many PCs before 20 % of the variance is explained?
  var.20pc <- sum(varExplained.cum <= 0.2)
  ### how much variance do 50 PCs explain?
  varpc.50PCA <- 100*(varExplained.cum[50])
  sink(paste0(outdir, sampleID, "_explanatoryPCs.txt"))
  cat(paste0("The first 50 PCs explain ", round(varpc.50PCA), "% of the variance. 20% of the variance is explained by the first ", var.20pc, " PCs"))
  sink()
  
  # define some graph functions which will be run with `to.pdf` later
  ## scree plot
  fig.scree <- function() {
    varExplained %>% tibble::enframe(name = "PC", value = "varExplained" ) %>%
      ggplot(aes(x = PC, y = varExplained)) + 
      theme_bw() +
      geom_bar(stat = "identity") +
      theme_classic() +
      ggtitle(paste0(sampleID, ": scree plot")) +
      ylab("explained variance")
  }
  ## cumulative variance
  fig.cumulativeVar <- function() {
    ggplot(as.data.frame(varExplained.cum), aes(y = varExplained.cum, x = seq(1, length(varExplained.cum)))) +
      geom_point(size = 1) +
      theme_bw() +
      ggtitle("cumulative variance explained by increasing PCs") +
      xlab("PCs") +
      ylab("cumulative explained variance") +
      geom_hline(yintercept = c(0.2), linetype = "dashed", color = "blue") +
      geom_vline(xintercept = c(20), linetype = "dashed", color = "blue")
  }
  
  # Make an elbow plot with elbow point annotated (adapted from Seurat's ElbowPlot() but to show all tested PCs)
  fig.elbow <- function() {
    ElbowPlot(seuratObj, ndims = 50, reduction = "pca") +
      theme_bw() +
      ggtitle(paste0(sampleID, ": elbow plot of standard deviations of principal components"))
  }
  
  # Perform JackStraw analysis
  if (runJackstraw == "TRUE") {
    seuratObj <- JackStraw(seuratObj, num.replicate = 100, dims = 50)
    seuratObj <- ScoreJackStraw(seuratObj, dims = 1:50) # because `RunPCA` calculates 50x PCs by defalt (you can change this)
    fig.jackstraw <- function() {
      JackStrawPlot(seuratObj, dims = 1:50) +
        ggtitle(paste0("PCA JackStraw - last significant PC is ", chosen.jack)) + 
        NoLegend()
    }
    
    fig.jackstraw2 <- function() {
      JackStrawPlot(seuratObj, dims = 1:50) +
        ggtitle(paste0("PCA JackStraw - last significant PC is ", chosen.jack))
    }
    # the PC p-vals are in seuratObj@reductions$pca@jackstraw$overall.p.values
    # get the PC number of the last PC before one is not significant
    jscores <- as.data.frame(seuratObj@reductions$pca@jackstraw$overall.p.values > 0.05)
    chosen.jack <- as.numeric(rownames(jscores[jscores$Score == "TRUE", ][1,])) - 1
    sink(paste0(outdir, sampleID, "_lastSigPC.txt"))
    cat(paste0("The PC number of the last PC before one is not significant is ", chosen.jack))
    sink()
    to.pdf(fig.jackstraw(), paste0(outdir, sampleID, "_PCA_jackstraw.pdf"))
    to.pdf(fig.jackstraw2(), paste0(outdir, sampleID, "_PCA_jackstraw_legend.pdf"))
    
  } else {
    if (runJackstraw == "FALSE") {
      print("skipping Jackstraw analysis")
    } else {
      stop("runJackstraw must be TRUE or FALSE")
    }
  }
  
  # Run the figure functions and save graphs as PDFs
  to.pdf(fig.scree(), paste0(outdir, sampleID, "_scree.pdf"))
  to.pdf(fig.cumulativeVar(), paste0(outdir, sampleID, "_cumulativeVariance.pdf"))
  to.pdf(fig.elbow(), paste0(outdir, sampleID, "_PCA_elbow.pdf"))
  # to.pdf(fig.jackstraw(), paste0(outdir, sampleID, "_PCA_jackstraw.pdf")) #run above in if/else bit
  
  # for now, just return the Seurat object
  return(seuratObj)
}

# --------------------------
# function: UMAP + tSNE
# --------------------------

func_runNonLinearDR <- function(seuratObj, runTSNE = "TRUE", dims = 20, neigh = 30, min.dist = 0.3) {
  # USAGE: seuratObj <- func_runNonLinearDR(seuratObj, runTSNE = "TRUE" or "FALSE")
  # OUTPUT: a Seurat object with tSNE and UMAP coordinates
  
  # Ensure that the parameters are numeric
  dims <- as.numeric(dims)
  neigh <- as.numeric(neigh)
  min.dist <- as.numeric(min.dist)
  # Run UMAP
  seuratObj <- Seurat::RunUMAP(seuratObj, dims = 1:dims, n.neighbors = neigh, min.dist = min.dist)
  fig.umap.raw <- function() {
    # alternative to DimPlot(seuratObj, reduction = "umap")
    Embeddings(seuratObj, reduction = "umap") %>%
      as.data.frame() %>%
      ggplot(aes(x = UMAP_1, y = UMAP_2)) +
      geom_point(size = 0.3) +
      theme_bw(base_size = 14) +
      ggtitle(paste0(sampleID, ": UMAP"))
  }
  to.pdf(fig.umap.raw(), paste0(outdir, sampleID, "_UMAP_raw.pdf"))
  
  # Run tSNE
  if (runTSNE == "TRUE") {
    seuratObj <- Seurat::RunTSNE(seuratObj, dims = 1:dims)
    
    fig.tSNE.raw <- function() {
      # alternative to DimPlot(seuratObj, reduction = "tSNE")
      Embeddings(seuratObj, reduction = "tsne") %>%
        as.data.frame() %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(size = 0.3) +
        theme_bw(base_size = 14) +
        ggtitle(paste0(sampleID, ": tSNE"))
    }
    to.pdf(fig.tSNE.raw(), paste0(outdir, sampleID, "_tSNE_raw.pdf"))
  } else {
    if (runTSNE == "FALSE") {
      print("skipping tSNE plot")
    } else {
      stop("runTSNE must be TRUE or FALSE")
    }
  }
  return(seuratObj)
}

# ------------------------------------------------------------------
# RUN FUNCTIONS
# ------------------------------------------------------------------

df <- func_ScaleData(df)
df <- func_runPCA(df, runJackstraw = "FALSE") 
df <- func_runNonLinearDR(df, runTSNE = "FALSE", dims = 50, neigh = 30, min.dist = 0.1)

# ------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------

saveRDS(df, file = paste0(outdir_obj, sampleID, "_processed_Seurat.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()