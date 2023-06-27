#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 27th August 2020
# Title: trigger_10j_ClusterMarkers.sh
# Goal: To identify all/top10 positive markers associated with each cluster of a given ident
# Usage: nohup ./trigger_10j_ClusterMarkers.sh > nohup_trigger_10j_ClusterMarkers.out 2>&1&
# Usage:Rscript 10f_ClusterStatsPlot.R {sampleID} {cluster_col} {seuratObj} {outdir}
# ------------------------------------------------------------------

cd /scripts/10f_integrated_clusterMarkers/
myscript=/scripts/0_code/10j_ClusterMarkers_newSeurat.R

echo "starting to analyse sample"
Rscript "$myscript" all integrated_snn_res.0.5 /scripts/10c_integrated_clustering/out_obj/allintegratedClustered_Seurat.RDS /scripts/10f_integrated_clusterMarkers/out_data/

echo "finished analysing sample"

