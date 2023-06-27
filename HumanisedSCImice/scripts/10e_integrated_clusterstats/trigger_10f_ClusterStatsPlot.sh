#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: trigger_10f_ClusterStatsPlot.sh
# Goal: To plot cluster vs annotation plots
# Usage: nohup ./trigger_10f_ClusterStatsPlot.sh > nohup_trigger_10f_ClusterStatsPlot.out 2>&1&
# Usage: Rscript 10f_ClusterStatsPlot.R {sampleID} {comparisonID} {cluster_col} {annotation_col} {seuratObj} {outdir}
# ------------------------------------------------------------------

cd /scripts/10e_integrated_clusterstats/
myscript=/scripts/1b_code/10f_ClusterStatsPlot.R

for i in sample injury tissue drug injurydrug
do
echo "starting to analyse sample "$i""
Rscript "$myscript" all clust_vs_"$i" integrated_snn_res.0.5 "$i" /scripts/10c_integrated_clustering/out_obj/allintegratedClustered_Seurat.RDS /scripts/10e_integrated_clusterstats/out_data
echo "finished analysing sample "$i""
done

