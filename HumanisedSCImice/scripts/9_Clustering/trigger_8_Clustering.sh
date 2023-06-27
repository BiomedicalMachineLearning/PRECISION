#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: trigger_8_Clustering.sh
# Goal: To run clustering for individual Seurat objects
# Usage: nohup ./trigger_8_Clustering.sh > nohup_trigger_8_Clustering.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/9_Clustering/
myscript=/scripts/1b_code/8_Clustering.R

for i in B1 B2 B3 B4 C2 C3
do
echo "starting to analyse sample "$i""
Rscript "$myscript" "$i"_filtered /scripts/8_SeuratProcessing/out_obj/"$i"_filtered_processed_Seurat.RDS /scripts/9_Clustering/out_data/ /scripts/9_Clustering/out_obj/ no 20
echo "finished analysing sample "$i""
done

