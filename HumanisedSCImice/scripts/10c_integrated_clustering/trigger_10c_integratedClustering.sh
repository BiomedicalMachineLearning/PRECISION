#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: 10c_integrationClustering.R.sh
# Goal: To run clustering for individual Seurat objects
# Usage: nohup ./10c_integrationClustering.R > nohup_10c_integrationClustering.R.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/10c_integrated_clustering/

myscript=/scripts/1b_code/10c_integrationClustering.R

echo "starting to analyse integrated sample"
Rscript "$myscript" all /scripts/10b_integrated_seuratProcessing/out_obj/all_processed_Seurat.RDS /scripts/10c_integrated_clustering/out_data/ /scripts/10c_integrated_clustering/out_obj/ no 20

echo "finished analysing integrated sample"


