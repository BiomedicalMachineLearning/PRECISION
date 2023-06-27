# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: 7b_optional_rerunUMAP.R
# Goal: Run if, based on the results of script #7, the user decides to change:
#       * The number of PCA dims for UMAP
#       * The UMAP parameters (min dist, n. neighbours)
# Usage: Rscript 7b_optional_rerunUMAP.R {sampleID} {seuratObj} {outdir_figs} {outdir_obj} {dims} {neighbours} {min.dist} {spread}
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 8) {
  cat("ERROR: 8 arguments expected\n")
  cat("example: R Rscript 7b_optional_rerunUMAP.R {sampleID} {seuratObj} {outdir_figs} {outdir_obj} {dims} {neighbours} {min.dist} {spread}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
outdir <- args[3]
outdir_obj <- args[4]
dims <- as.numeric(args[5])
neighbours <- as.numeric(args[6])
dist <- as.numeric(args[7])
spread <- as.numeric(args[7])

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
# RE-RUN UMAP
# ------------------------------------------------------------------

df <- Seurat::RunUMAP(df, dims = 1:dims, n.neighbors = neighbours, min.dist = dist, spread = spread)

# alternative to DimPlot(seuratObj, reduction = "umap")
plot <- Embeddings(df, reduction = "umap") %>%
    as.data.frame() %>%
    ggplot(aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = 0.3) +
    theme_bw(base_size = 14) +
    ggtitle(paste0(sampleID, ": UMAP"))
ggsave(plot = plot, paste0(outdir, sampleID, "_reUMAP_raw.pdf"))

# ------------------------------------------------------------------
# SAVE OBJECT
# ------------------------------------------------------------------

saveRDS(df, file = paste0(outdir_obj, sampleID, "_reprocessed_Seurat.RDS"))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()