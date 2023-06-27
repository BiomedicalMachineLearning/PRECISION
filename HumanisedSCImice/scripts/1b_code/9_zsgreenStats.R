# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: 9_zsgreenStats.R
# Goal: To visualise the zs-green expression in our samples
# Usage: Rscript 9_zsgreenStats.R {sampleID} {seuratObj} {outdir}
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  cat("ERROR: 3 arguments expected\n")
  cat("example: Rscript 9_zsgreenStats.R {sampleID} {seuratObj} {outdir}\n")
  quit()
}

sampleID <- args[1]
df <- readRDS(args[2])
outdir <- args[3]

# load packages
#bioc_packages <- c()
r_packages <- c("Seurat", "ggplot2")
## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    instdf.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
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
# ZS-GREEN STATS
# ------------------------------------------------------------------

DefaultAssay(df) <- "RNA"
#Idents(object = df) <- "integrated_snn_res.0.1"
my.features <- c("zsgreen", "Ptprc")

pdf(paste0(outdir, sampleID, "_zsgreen+cd45_ridge.pdf"))
RidgePlot(df, features = my.features)
dev.off()

pdf(paste0(outdir, sampleID, "_zsgreen+cd45_violin.pdf"))
VlnPlot(df, features = my.features)
dev.off()

plot1 <- FeaturePlot(df, features = c("zsgreen"), max.cutoff = "q1", pt.size = 0.5) +
  ggtitle("zsgreen +/- (on = 1st percentile upwards)") +
  theme(legend.position = "none")
plot1$data <- plot1$data[order(plot1$data$zsgreen),]
ggsave(paste0(outdir, sampleID, "_zsgreen_feature_onoff.pdf"), plot = plot1)


plot2 <- FeaturePlot(df, features = c("Ptprc"), max.cutoff = "q1", pt.size = 0.5) +
  ggtitle("Ptprc +/- (on = 1st percentile upwards)") +
  theme(legend.position = "none")
plot2$data <- plot2$data[order(plot2$data$Ptprc),]
ggsave(paste0(outdir, sampleID, "_cd45_feature_onoff.pdf"), plot = plot2)

my.dotplot <- DotPlot(df, features = my.features) + RotatedAxis()
pdf(paste0(outdir, sampleID, "_zsgreen+cd45_dotplot.pdf"))
my.dotplot
dev.off()
write.table(my.dotplot$data, file = paste0(outdir, sampleID, "_zsgreen+cd45_percentageCluster.txt"), sep = "\t", quote = FALSE)

my.dotplot2 <- DotPlot(df, features = my.features, group.by = "orig.ident") + RotatedAxis()
write.table(my.dotplot2$data, file = paste0(outdir, sampleID, "_zsgreen+cd45_percentageOverall.txt"), sep = "\t", quote = FALSE)

pdf(paste0(outdir, sampleID, "_zsgreen+cd45_blend_feature.pdf"), width = 12)
FeaturePlot(object = df, features = my.features, cols = c("grey", "red", "blue"), blend = TRUE, sort.cell = TRUE)
dev.off()

pdf(paste0(outdir, sampleID, "_zsgreen+cd45_blend_feature_onoff.pdf"), width = 12, height = 3)
FeaturePlot(object = df, features = my.features, cols= c("grey", "red", "blue"), blend = TRUE, sort.cell = TRUE, max.cutoff = "q1")
dev.off()

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()