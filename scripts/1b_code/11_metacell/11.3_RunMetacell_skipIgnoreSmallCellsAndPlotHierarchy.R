# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: 11.3_RunMetacell.R
# Goal: To run metacell for a given sample
# Usage: Rscript 11.3_RunMetacell.R {sampleID} {CellRangerDir} {outdir_fig} {outdir_obj}
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("ERROR: 4 arguments expected\n")
  cat("example: Rscript 11.3_RunMetacell.R {sampleID} {CellRangerDir} {outdir_fig} {outdir_obj}\n")
  quit()
}

sampleID <- args[1]
CellRangerDir <- args[2]
outdir_fig <- args[3]
outdir_obj <- args[4]

# load packages
bioc_packages <- c("metacell")
#r_packages <- c()
## function to load R packages
# baseRpkgTest <- function(x) {
#   if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
#     install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
#     if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
#   }
# }
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
#for (r_pkg in r_packages) {
#  baseRpkgTest(r_pkg)
#}

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir_fig, "/") == FALSE) {
  outdir_fig <- paste0(outdir_fig, "/", sep="")
}

if (endsWith(outdir_obj, "/") == FALSE) {
  outdir_obj <- paste0(outdir_obj, "/", sep="")
}

# ------------------------------------------------------------------
# SET UP OUTPUT DIRECTORIES
# ------------------------------------------------------------------

print(paste0("Log 1: Starting to set up output directories for MC sample ", sampleID, " at ", date()))

# setup the object library
if(!dir.exists(outdir_obj)) dir.create(outdir_obj)
scdb_init(outdir_obj, force_reinit=T)

# setup the figure library
if(!dir.exists(outdir_fig)) dir.create(outdir_fig)
scfigs_init(outdir_fig)

# ------------------------------------------------------------------
# SET UP COUNT MATRICES
# ------------------------------------------------------------------

print(paste0("Log 2: Starting to set up count matrices for MC sample ", sampleID, " at ", date()))

mcell_import_scmat_10x("test", base_dir=CellRangerDir)
mat = scdb_mat("test")
#print(dim(mat@mat))

# ------------------------------------------------------------------
# EXPLORING AND FILTERING THE UMI MATRIX
# ------------------------------------------------------------------

print(paste0("Log 3: Exploring and filtering UMI matrix for MC sample ", sampleID, " at ", date()))


#' To get a basic understanding of the new data, we will plot the distribution of UMI count per cell (the plot is thresholded at 800 UMIs by default).
#' MetaCell print the figures to the outdir_fig directory (e.g. for this step, it makes a figure called `test.total_umi_distr.png`).
#' to read this plot, each point on the axis is 2^x e.g. 2^9 = 512 UMIs

mcell_plot_umis_per_cell("test")

#' Or we can plot this ourselves without a log2 axis
hist(Matrix::colSums(mat@mat), breaks = 200)
hist(Matrix::colSums(mat@mat), breaks = 200, xlim = c(0,5000))

#' We want to clean some known issues from the matrix before starting to work with it.
#' We generate a list of mitochondrial genes that typically mark cells as being stressed or dying,
#' as well as immunoglobulin genes that may represent strong clonal signatures in plasma cells,
#' rather than cellular identity. We will ask the package to ignore these genes:
mat = scdb_mat("test")
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
# NOTE: Mouse only!
ig_genes = c(grep("^Igj", nms, v=T), 
             grep("^Igh",nms,v=T),
             grep("^Igk", nms, v=T), 
             grep("^Igl", nms, v=T))

bad_genes = unique(c(grep("^mt-", nms, v=T), grep("^Mtmr", nms, v=T), grep("^Mtnd", nms, v=T),"Neat1","Tmsb4X", "Tmsb10", ig_genes))
mcell_mat_ignore_genes(new_mat_id="test", mat_id="test", bad_genes, reverse=F) 

#' Ignored genes are kept in the matrix for reference, but all downstream analysis will disregard them.
#' This means that the number of UMIs from these genes cannot be used to distinguish between cells.
#' In the current example we will also eliminate cells with less than 800 UMIs (threshold can be set based on examination of the UMI count distribution):
#' NB: nUMIs = nCounts
#mcell_mat_ignore_small_cells("test", "test", 800) # EDIT: PossDC and B failed at this step - skip

# ------------------------------------------------------------------
# SELECTING FEATURE GENES
# ------------------------------------------------------------------

print(paste0("Log 4: Selecting feature genes for MC sample ", sampleID, " at ", date()))

#' We move on to computing statistics on the distributions of each gene in the data,
#' which are going to be our main tool for selecting feature genes for MetaCell analysis.
#' This will makes a new object (type **gstat**) under the name `test`
#' (NB: a file called `gstat.test.Rda` is saved to the `testdb` folder), by analysing the `test` count matrix.  
mcell_add_gene_stat(gstat_id="test", mat_id="test", force=T)

#' If we wanted to, we could now explore intesting genes and their distributions.
#' But here we will move directly to select a gene set for downstream analysis.
#' We create a new object of type gset (gene set), to which all genes whose scaled variance (variance divided by mean)
#' exceeds a given threshold are added. (NB: a file called `gset.test_feats.Rda` is made)
#' In `mcell_gset_filter_varmean` we create a new gene set of all genes for which scaled variance is 0.08 and higher.
#' In `mcell_gset_filter_cov` we restrict the gene set to only genes with 100+ UMIs across the whole dataset, and to genes with 3+ cells for 2+ UMIs
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.08, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=100, T_top3=2)


#' To check that our parameters were OK, we can plot all our genes on a mean vs variance plot.
#' This makes a plot called `test.varmin.png`. Selected genes are coloured red and unselected genes are plotted black.
#' We want to make sure that there's not too many black spots or they're in a weird distribution. For comparison, see the plot in the Metacell vignette.
mcell_plot_gstats(gstat_id="test", gset_id="test_feats")

#' Strategies for studying the selected genes will be discussed in another vignette ("Supervised filtering of feature genes").

# ------------------------------------------------------------------
# BUILDING THE BALANCED CELL GRAPH
# ------------------------------------------------------------------

print(paste0("Log 5: Building the balanced cell graph for MC sample ", sampleID, " at ", date()))

#' Assuming we are happy with the selected genes, we will move forward to create a similarity graph (cgraph),
#' using a construction called balanced K-nn graph. When we run this command, we get a new `cgraph` object called `test_graph`.
#' (NB: a file called `cgraph.test_graph.Rda` is saved to the `testdb` folder)
#' 
#' The `k = 100` parameter is important, it affects the size distributions of the derived metacells.
#' **NOTE** In Supplementary Table 1 of the Metacell paper, they recommend you do not change this parameter.
#' If you have more than 20-30k cells, this may get slow.
mcell_add_cgraph_from_mat_bknn(mat_id="test", 
                               gset_id = "test_feats", 
                               graph_id="test_graph",
                               K=100,
                               dsamp=T)

# ------------------------------------------------------------------
# RESAMPLING AND GENERATING THE CO-CLUSTERING GRAPH
# ------------------------------------------------------------------

print(paste0("Log 6: Resampling and generating the co-clustering graph for MC sample ", sampleID, " at ", date()))

#' Now we use the cgraph to sample 500 metacell partitions, each covering 75% of the cells and organising them into dense subgraphs.
#' The metacell size distribution of the resampled partitions will mostly be determined by the `K` parameter used above.
#' This resampling makes a new `coclust` object in the database called `test_coc500` and stores the **number of times each pair of cells was part of the same metacell**.
#' (NB: a file called `coclust.test_coc500.Rda` is saved to the `testdb` folder) 
#' This process can be slow if the graphs are very big, but you can reduce the `n_resamp` to run fewer resampling steps.
mcell_coclust_from_graph_resamp(
  coc_id="test_coc500", 
  graph_id="test_graph",
  min_mc_size=20, 
  p_resamp=0.75, n_resamp=500)


#' The co-clustering statistics we just ran are used to generate a new similarity graph,
#' based on which accurate calling of the final set of metacells is done.
#' We will create an object `test_mc` based on analysis of the co-clustering graphs.
#' The `K` parameter determines the number of neighbours we wish to minimally associate with each cell (minimum cluster size).
#' (NB: a file called `mc.test_mc.Rda` is saved to the `testdb` folder)
#' Prior to partitioning the co-cluster graph is filtered to eliminate highly unbalanced edges, with smaller alpha resulting in harsher filtering.
mcell_mc_from_coclust_balanced(
  coc_id="test_coc500", 
  mat_id= "test",
  mc_id= "test_mc", 
  K=30, min_mc_size=30, alpha=2)

# ------------------------------------------------------------------
# REMOVING OUTLIER CELLS
# ------------------------------------------------------------------

print(paste0("Log 7: Removing outlier cells for MC sample ", sampleID, " at ", date()))

#' We now have a preliminary metacell object.'
#' It is a good practice to make sure all metacells within it are homogeneous.
#' This is done by the outlier scan procedure, which splits metacells whose underlying similarity structure supports the existence of multiple sub-clusters,
#' and removes outlier cells that strongly deviate from their metacell’s expression profile.
#' In `mcell_plot_outlier_heatmap` we make a heatmap summarising the detected outlier behaviours
#' (this is only possible for datasets of modest size). This heatmap is called `test_mc.outlier.png`
#' In `mcell_mc_split_filt` we remove cells with outlier expression from each metacell
mcell_plot_outlier_heatmap(mc_id="test_mc", mat_id = "test", T_lfc=3)

mcell_mc_split_filt(new_mc_id="test_mc_f", 
                    mc_id="test_mc", 
                    mat_id="test",
                    T_lfc=3, plot_mats=F)

# ------------------------------------------------------------------
# SELECTING MARKERS AND COLOURING METACELLS
# ------------------------------------------------------------------

print(paste0("Log 8: Selecting markers and colouring metacells for MC sample ", sampleID, " at ", date()))

#' The filtered metacell object `test_mc_f` is saved to the `testdb` dir as `mc.test_mc_f.Rda`.
#' It can now be visualized. In order to do this effectively, we usually go through one or two iterations of selecting informative marker genes.
#' The package can select markers for you automatically - by simply looking for genes that are strongly enriched in any of the metacells:
mcell_gset_from_mc_markers(gset_id="test_markers", mc_id="test_mc_f")

#' It is however very useful to analyze metacell models in depth, and select genes with known or hypothesized biological significance. We will assign each of these genes with a color, fold change threshold and priority, using a table that look like this:
# NOTE: Assess these!
marks_colors = read.table(system.file("extdata", "pbmc_mc_colorize.txt", package="metacell"), sep="\t", h=T, stringsAsFactors=F)
# But the gene IDs need to be modified to be lowercase except the first letter.
marks_colors$gene <- stringr::str_to_sentence(marks_colors$gene)
marks_colors$color <- gsub("cyan2", "slateblue1", marks_colors$color, ignore.case=T)
marks_colors$color <- gsub("darkgreen", "olivedrab2", marks_colors$color, ignore.case=T)
# DNTT and FLT3 are both steel blue, is this on purpose? Change to at least a slightly different blue
marks_colors$color <- ifelse(grepl("Dntt", marks_colors$gene), gsub("steelblue", "royalblue3", marks_colors$color), marks_colors$color)
marks_colors

#' Applying this table to color metacells is done using the command `mc_colorize` as shown below. Note that there are more sophisticated ways to color/annotate metacells, and that the model’s understanding and annotation can often greatly benefit from looking at gene distributions and reading literature about possible functions and regulatory mechanisms.
mc_colorize("test_mc_f", marker_colors=marks_colors)

#' We are now equipped with some basic coloring of metacells, which can also be accessed directly:
mc = scdb_mc("test_mc_f")
table(mc@colors)

# ------------------------------------------------------------------
# CREATING A HEATMAP OF GENES AND METACELLS
# ------------------------------------------------------------------

print(paste0("Log 9: Creating a heatmap of genes and metacells for MC sample ", sampleID, " at ", date()))

#' We can use the colors to produce a labeled heat map, showing selected genes and their distributions over metacells,
#' with the colored annotation shown at the bottom. Note that the values plotted are color coded log2(fold enrichment) 
#' value of the metacell over the median of all other metacells. 
mcell_mc_plot_marks(mc_id="test_mc_f", gset_id="test_markers", mat_id="test")


#' It can be useful to explore these values directly - e.g.:
lfp = log2(mc@mc_fp)
#tail(sort(lfp["Cd8a",]))
somegene <- head(rownames(lfp), n = 1L)
somegene
tail(sort(lfp[somegene,]))

# ------------------------------------------------------------------
# PROJECTING METACELLS AND CELLS IN 2D
# ------------------------------------------------------------------

print(paste0("Log 10: Projecting metacells and cells in 2D for MC sample ", sampleID, " at ", date()))

#' Heat maps are useful but sometimes hard to interpret, and so we may want to visualize the similarity structure among metacells
#' (or among cells within metacells). To this end we construct a 2D projection of the metacells, and use it to plot the metacells 
#' and key similarities between them (shown as connecting edges), as well as the cells. This plot will use the same metacell
#' coloring we established before (and in case we improve the coloring based on additional analysis, the plot can be regenerated):
mcell_mc2d_force_knn(mc2d_id="test_2dproj",mc_id="test_mc_f", graph_id="test_graph")
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="test_2dproj")

#' Note that we changed the metacell parameters “mcell_mc2d_height/width” to get a reasonably-sized figure.
#' There are many additional parameters that can be tuned in MetaCell, and more of those meant for routine tuning
#' will be discussed in other vignettes. We obtain the following figure:

# ------------------------------------------------------------------
# VISUALISING THE MC CONFUSION MATRIX
# ------------------------------------------------------------------

print(paste0("Log 11: Visualising the MC confusion matrix for MC sample ", sampleID, " at ", date()))

#' While 2D projections are popular and intuitive (albeit sometimes misleading) ways to visualize scRNA-seq results,
#' we can also summarize the similarity structure among metacells using a “confusion matrix” which encodes the 
#' pairwise similarities between all metacells. This matrix may capture hierarchical structures or other complex organizations among metacells.
#' 
#' We first create a hierarchical clustering of metacells, based on the number of similarity relations between their cells:
mc_hc = mcell_mc_hclust_confu(mc_id="test_mc_f", 
                              graph_id="test_graph")

mc_sup = mcell_mc_hierarchy(mc_id="test_mc_f",
                            mc_hc=mc_hc, T_gap=0.04)
#mcell_mc_plot_hierarchy(mc_id="test_mc_f",  # Failed for B cell step
#                        graph_id="test_graph", 
#                        mc_order=mc_hc$order, 
#                        sup_mc = mc_sup, 
#                        width=2800, heigh=2000, min_nmc=2)

# ------------------------------------------------------------------
# SAVE METACELL PARTITIONS
# ------------------------------------------------------------------

print(paste0("Log 12: Saving metacell partitions for MC sample ", sampleID, " at ", date()))

metacell.list <- as.data.frame(mc@mc)
metacell.colourchart <- mc@color_key
metacell.groupings <- as.data.frame(mc@colors)

write.table(metacell.list, file = paste0(outdir_fig, "metacell.list.txt"),
            sep = "\t", quote = FALSE)
write.table(metacell.colourchart, file = paste0(outdir_fig, "metacell.colourchart.txt"),
            sep = "\t", quote = FALSE)
write.table(metacell.groupings, file = paste0(outdir_fig, "metacell.groupings.txt"),
            sep = "\t", quote = FALSE)

print(paste0("Log 13: Wrapping up metacell run for MC sample ", sampleID, " at ", date()))

# ------------------------------------------------------------------
# SESSION INFO
# ------------------------------------------------------------------

sink(paste0(outdir_fig, sampleID, "_sessionInfo.txt"))
sessionInfo()
sink()

print(paste0("Log 14: Finished for MC sample ", sampleID, " at ", date()))
