#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 27th August 2020
# Title: trigger_10j_ClusterMarkers.sh
# Goal: To identify all/top10 positive markers associated with each cluster of a given ident
# Usage: nohup ./trigger_10j_ClusterMarkers.sh > nohup_trigger_10j_ClusterMarkers.out 2>&1&
# Usage:Rscript 10f_ClusterStatsPlot.R {sampleID} {cluster_col} {seuratObj} {outdir}
# ------------------------------------------------------------------

cd /scripts/17b_Rerun_integrated_clusterMarkers_Vehicle/
myscript=/scripts/0_code/10j_ClusterMarkers_newSeurat.R

echo "starting to analyse sample"
Rscript "$myscript" all clusters_wholedataset /scripts/17_MetacellVehicleOnly/vehicle/11.1a_ReplotUMAP/out_obj/ReplotUMAP_Seurat.RDS /scripts/17b_Rerun_integrated_clusterMarkers_Vehicle/out_data/

echo "finished analysing sample"

