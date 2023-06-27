# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 22nd June 2020
# Title: 1_QC_3MAD.R
# Goal: To perform 
# Output: QC plots/tables, and a list of filtered cells
# Usage: Rscript 1_QC_3MAD.R {sampleID} {path/to/SCE} {path/to/outdir}
# ------------------------------------------------------------------
# TO DO:
# 1. Make Seurat comparison optional
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

# load arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 3) {
  cat("ERROR: 3 arguments expected\n")
  cat("example: Rscript QC_3MAD_scater.R {sampleID} {path/to/SCE} {path/to/outdir}")
  quit()
}

sampleID <- args[1]
sce.path <- args[2]
outdir <- args[3]

# load packages
bioc_packages <- c("scater", "edgeR")
r_packages <- c("Seurat", "ggplot2", "VennDiagram")
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

print("If you are not working with mouse data, change the mitochondrial/ribosomal gene codes")

# ------------------------------------------------------------------
# SETUP
# ------------------------------------------------------------------

print("loading data")

# load data
df <- readRDS(sce.path)
df.seurat <- as.Seurat(df)

# ------------------------------------------------------------------
# QC CALCULATION
# ------------------------------------------------------------------

print("QC metrics")

is.mito <- rownames(df)[grep("^MT-",rownames(df))]
is.ribo <- rownames(df)[grep("^RPS|^RPL",rownames(df))]

qc <- perCellQCMetrics(df, subsets=list(Mito=is.mito, Ribo=is.ribo))

# ------------------------------------------------------------------
# OUTLIER DETECTION
# ------------------------------------------------------------------

print("outlier detection")

# log-transformed library size
qc.lib2 <- isOutlier(qc$sum, log=TRUE, type="lower")
## get some stats
nfilt_libsize <- sum(qc.lib2)
thresh_libsize <- attr(qc.lib2, "thresholds")
thresh_libsize <- thresh_libsize[1]  #lower bound
dim(thresh_libsize) <- NULL

# log-transformed library size
# here we just do lower-end library sizes; large cells are filtered later in the droplet filtration step
qc.nexprs2 <- isOutlier(qc$detected, log=TRUE, type="lower")
## get some stats
nfilt_nexprs <- sum(qc.nexprs2)
thresh_nexprs <- attr(qc.nexprs2, "thresholds")
thresh_nexprs <- thresh_nexprs[1] #lower bound
dim(thresh_nexprs) <- NULL

# log-transformed library size
qc.mito2 <- isOutlier(qc$subsets_Mito_percent, type="higher")
## get some stats
nfilt_mito <- sum(qc.mito2)
thresh_mito <- attr(qc.mito2, "thresholds")
thresh_mito <- thresh_mito[2] #upper bound
dim(thresh_mito) <- NULL

# log-transformed library size
qc.ribo2 <- isOutlier(qc$subsets_Ribo_percent, type="higher")
## get some stats
nfilt_ribo <- sum(qc.ribo2)
thresh_ribo <- attr(qc.ribo2, "thresholds")
thresh_ribo <- thresh_ribo[2] #upper bound
dim(thresh_ribo) <- NULL

#A cell that is an outlier for any of these metrics is considered to be of low quality and discarded.
discard2 <- qc.lib2 | qc.nexprs2 | qc.ribo2 | qc.mito2
nfilt_total <- sum(discard2)

# summarise discarded cells
reasons <- quickPerCellQC(qc, percent_subsets=c("subsets_Mito_percent", "subsets_Ribo_percent"))
reasons.tab <- as.data.frame(colSums(as.matrix(reasons)))
colnames(reasons.tab) <- "nDiscard"

# ------------------------------------------------------------------
# DATA VIZ
# ------------------------------------------------------------------

print("data viz")

# add metadata to df
colData(df) <- cbind(colData(df), qc)
df$discard <- reasons$discard
df$low_lib_size <- reasons$low_lib_size
df$low_n_features <- reasons$low_n_features
df$high_subsets_Mito_percent <- reasons$high_subsets_Mito_percent
df$high_subsets_Ribo_percent <- reasons$high_subsets_Ribo_percent

# library size stats
## jitter
plotColData(df, x="orig.ident", y="sum", colour_by="discard",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Total count (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_libsize_allfilt_jitter.pdf"))

plotColData(df, x="orig.ident", y="sum", colour_by="low_lib_size",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Total count (filtered due to libsize)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_libsize_libfilt_jitter.pdf"))

## swarm
plotColData(df, x="orig.ident", y="sum", colour_by="discard",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) + 
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Total count (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_libsize_allfilt_swarm.pdf"))

plotColData(df, x="orig.ident", y="sum", colour_by="low_lib_size",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") + 
  ggtitle("Total count (filtered due to libsize)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_libsize_libfilt_swarm.pdf"))

# nGenes stats
## jitter
plotColData(df, x="orig.ident", y="detected", colour_by="discard",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Detected features (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ngenes_allfilt_jitter.pdf"))

plotColData(df, x="orig.ident", y="detected", colour_by="low_n_features",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Detected features (filtered due to low # features)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ngenes_genefilt_jitter.pdf"))

## swarm
plotColData(df, x="orig.ident", y="detected", colour_by="discard",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Detected features (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ngenes_allfilt_swarm.pdf"))

plotColData(df, x="orig.ident", y="detected", colour_by="low_n_features",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Detected features (filtered due to low # features)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ngenes_genefilt_swarm.pdf"))

# mitochondrial stats
## jitter
plotColData(df, x="orig.ident", y="subsets_Mito_percent", colour_by="discard",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Mito percent (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_mito_allfilt_jitter.pdf"))

plotColData(df, x="orig.ident", y="subsets_Mito_percent", colour_by="high_subsets_Mito_percent",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Mito percent (filtered due to high % mitochondria)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_mito_mitofilt_jitter.pdf"))

## swarm
plotColData(df, x="orig.ident", y="subsets_Mito_percent", colour_by="discard",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) + 
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Detected features (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_mito_allfilt_swarm.pdf"))

plotColData(df, x="orig.ident", y="subsets_Mito_percent", colour_by="high_subsets_Mito_percent",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) + 
  scale_y_log10() + ggtitle("Detected features (filtered due to low # features)") + scale_fill_discrete(name="Discard?")
ggsave(file = paste0(outdir, sampleID, "_filtcells_mito_mitofilt_swarm.pdf"))

# ribosomal stats
## jitter
plotColData(df, x="orig.ident", y="subsets_Ribo_percent", colour_by="discard",
            point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Ribosomal percent (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ribo_allfilt_jitter.pdf"))

plotColData(df, x="orig.ident", y="subsets_Ribo_percent", colour_by="high_subsets_Ribo_percent", point_alpha = 0.5, jitter_type = "jitter", show_median = FALSE, show_violin = FALSE) + 
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Ribosomal percent (filtered due to high % ribosomes)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ribo_ribofilt_jitter.pdf"))

## swarm
plotColData(df, x="orig.ident", y="subsets_Ribo_percent", colour_by="discard",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Ribosomal percent (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ribo_allfilt_swarm.pdf"))

plotColData(df, x="orig.ident", y="subsets_Ribo_percent", colour_by="high_subsets_Ribo_percent",
            point_alpha = 0.2, jitter_type = "swarm", show_median = TRUE, show_violin = TRUE) +
  scale_y_log10() +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Ribosomal percent (filtered due to high % ribosomes)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_ribo_ribofilt_swarm.pdf"))

# mitochondrial vs library size stats
plotColData(df, x="sum", y="subsets_Mito_percent", colour_by="discard") +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Mitochondrial percent vs. library size (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_mitoVSlibsize_allfilt.pdf"))

plotColData(df, x="sum", y="subsets_Mito_percent", colour_by="high_subsets_Mito_percent") +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Mitochondrial percent vs. library size (filtered due to high % mitochondria)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_mitoVSlibsize_mitofilt.pdf"))

# ribosomal vs library size stats
plotColData(df, x="sum", y="subsets_Ribo_percent", colour_by="discard") +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Ribosomal percent vs. library size (all filtered)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_riboVSlibsize_allfilt.pdf"))

plotColData(df, x="sum", y="subsets_Ribo_percent", colour_by="high_subsets_Ribo_percent") +
  scale_fill_discrete(name="Discard?") +
  ggtitle("Ribosomal percent vs. library size (filtered due to high % ribosomes)")
ggsave(file = paste0(outdir, sampleID, "_filtcells_riboVSlibsize_ribofilt.pdf"))

# ------------------------------------------------------------------
# STATISTICAL COMPARISON OF DISCARDED VS KEPT CELLS
# ------------------------------------------------------------------

print("statistical comparison")

# get lists of kept vs discarded cells
lost <- calculateAverage(counts(df)[,!reasons$discard])
kept <- calculateAverage(counts(df)[,reasons$discard])

# work out logFC between kept vs discarded cells
logged <- cpm(cbind(lost, kept), log=TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

# Plot the logFC, and highlight mitochondrial genes
pdf(paste0(outdir, sampleID, "_keptVSlost_mt.pdf"))
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16, main = "Kept vs discarded cells (with mitochondrial genes)")
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16)
dev.off()

# Plot the logFC, and highlight ribosomes genes
pdf(paste0(outdir, sampleID, "_keptVSlost_rb.pdf"))
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16, main = "Kept vs discarded cells (with ribosomal genes)")
points(abundance[is.ribo], logFC[is.ribo], col="dodgerblue", pch=16)
dev.off()

#Let’s look at the most up- and down-regulated genes - a positive number is “up” in lost, and a negative number is “down” in lost.
fc <- as.data.frame(logFC)
fc$logFC <- as.numeric(fc$logFC)
fc <- fc[order(fc$logFC), , drop = FALSE]
top10_up <- rownames(tail(fc, n = 10L))
top10_down <- rownames(head(fc, n = 10L))

# Plot the logFC, and highlight genes that are upregulated in the discarded cells
pdf(paste0(outdir, sampleID, "_keptVSlost_top10upInDiscard.pdf"))
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[top10_up], logFC[top10_up], col="orange", pch=16)
dev.off()

# Plot the logFC, and highlight genes that are upregulated in the discarded cells
pdf(paste0(outdir, sampleID, "_keptVSlost_top10downInDiscard.pdf"))
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
points(abundance[top10_down], logFC[top10_down], col="blue", pch=16)
dev.off()

# ------------------------------------------------------------------
# FILTER OUTLIERS
# ------------------------------------------------------------------

print("filter outliers")

cells.keep <- rownames(colData(df)[!reasons$discard,])
cells.discard <- rownames(colData(df)[reasons$discard,])
df.filt <- df[,!reasons$discard]

# ------------------------------------------------------------------
# COMPARE 3 MAD FILTERING TO PREVIOUS PIPELINE SEURAT FILTERING
# ------------------------------------------------------------------

print("compare to seurat")

df.seurat[["percent.mt"]] <- PercentageFeatureSet(df.seurat, pattern = "^MT-")
df.seurat[["percent.rb"]] <- PercentageFeatureSet(df.seurat, pattern = "^RPS|^RPL")

# mitochondrial plot
VlnPlot(df.seurat, features = c("percent.mt")) +
  ggtitle(label = paste0(sampleID, " - percentage mitochondrial genes per cell"),
          subtitle = "blue line = Seurat upper threshold, red line = 3MAD upper threshold") +
  ylab("% mitochondrial genes") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = thresh_mito, linetype = "dashed", color = "red") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file = paste0(outdir, sampleID, "_seurat_mitopercent.pdf"))

# ribosomal plot
VlnPlot(df.seurat, features = c("percent.rb")) +
  ggtitle(label = paste0(sampleID, " - percentage ribosomal genes per cell"),
          subtitle = "blue line = Seurat upper threshold, red line = 3MAD upper threshold") +
  ylab("% ribosomal genes") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = thresh_ribo, linetype = "dashed", color = "red") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file = paste0(outdir, sampleID, "_seurat_ribopercent.pdf"))

# mitochondria % vs nFeatures
basicplot <- FeatureScatter(df.seurat, feature1 = "nFeature_RNA", feature2 = "percent.mt", pt.size = 0.5)
pearson <- basicplot$labels$title
basicplot +
  ggtitle(label = paste0(sampleID, " - % mitochondrial genes vs number of genes per cell (pearson = ", pearson, ")"),
          subtitle = "blue lines = Seurat upper thresholds, red line = 3MAD upper threshold,\ngreen line = 3MAD lower threshold") +
  ylab("% mitochondrial genes") +
  xlab("number of genes") +
  geom_hline(yintercept = c(20), linetype="dashed", color = "blue") +
  geom_hline(yintercept = thresh_mito, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(4000), linetype="dashed", color = "blue") +
  geom_vline(xintercept = thresh_nexprs, linetype="dashed", color = "darkgreen") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme(legend.position = "none")
ggsave(file = paste0(outdir, sampleID, "_seurat_mitoVSnGenes.pdf"))

# ribosome % vs nFeatures
basicplot <- FeatureScatter(df.seurat, feature1 = "nFeature_RNA", feature2 = "percent.rb", pt.size = 0.5)
pearson <- basicplot$labels$title
basicplot +
  ggtitle(label = paste0(sampleID, " - % ribosomal genes vs number of genes per cell (pearson = ", pearson, ")"),
          subtitle = "blue lines = Seurat upper thresholds, red line = 3MAD upper threshold,\ngreen line = 3MAD lower threshold") +
  ylab("% ribosomal genes") +
  xlab("number of genes") +
  geom_hline(yintercept = 50, linetype="dashed", color = "blue") +
  geom_hline(yintercept = thresh_ribo, linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(4000), linetype="dashed", color = "blue") +
  geom_vline(xintercept = thresh_nexprs, linetype="dashed", color = "darkgreen") +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  theme_bw() +
  theme(legend.position = "none")
ggsave(file = paste0(outdir, sampleID, "_seurat_riboVSnGenes.pdf"))

# get lists of kept and discarded cells from this method
df.seurat_keep <- subset(df.seurat, subset = percent.mt < 20 & percent.rb < 50)
df.seurat_discard <- subset(df.seurat, subset = percent.mt < 20 & percent.rb < 50, invert = TRUE)
seurat_cells.keep <- rownames(df.seurat_keep[[]])
seurat_cells.discard <- rownames(df.seurat_discard[[]])

# compare the cells discarded by previous vs 3MAD filtering
venn.diagram(x = list(seurat_cells.discard, cells.discard),
             category.names = c("seurat_discard", "scater_discard"),
             filename = paste0(outdir, sampleID, "_discard.png"),
             output = TRUE)

# ------------------------------------------------------------------
# SAVE OUTPUT
# ------------------------------------------------------------------

print("save output")

# save a list of the selected thresholds
threshtab <- as.data.frame(DataFrame(parameter=c("library size (lower bound)", "# expressed genes (lower bound)", "% mitochondrial genes (upper bound)", "% ribosomal genes (upper bound)"), threshold_values = c(
  thresh_libsize, thresh_nexprs, thresh_mito, thresh_ribo))
)
write.table(threshtab, file = paste0(outdir, sampleID, "_filterthresholds.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# save a list of the top 10 up/downregulated genes in discarded vs kept cells
top10genes <- as.data.frame(DataFrame(Up=top10_up, Down=top10_down))
write.table(top10genes, file = paste0(outdir, sampleID, "_keptVSdiscardCells_top10updown.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# save a table of counts for # cells discarded and for what reason
write.table(reasons.tab, file = paste0(outdir, sampleID, "_cellFiltrationReasons.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# save the QC metadata
write.table(colData(df), file = paste0(outdir, sampleID, "_cellQC.txt"), sep = "\t", quote = FALSE, row.names = TRUE)

# save lists of kept and discarded cells from 3MAD and previous-pipeline filtering
write.table(as.data.frame(cells.keep), file = paste0(outdir, sampleID, "_cellQC_3MADkeptcells.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(cells.discard), file = paste0(outdir, sampleID, "_cellQC_3MADdiscardedcells.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(seurat_cells.keep), file = paste0(outdir, sampleID, "_cellQC_seuratkeptcells.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(as.data.frame(seurat_cells.discard), file = paste0(outdir, sampleID, "_cellQC_seuratdiscardedcells.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()

print("done")
