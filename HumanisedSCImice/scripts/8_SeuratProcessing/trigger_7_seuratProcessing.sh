#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: trigger_7_seuratProcessing.sh
# Goal: To convert Cellranger output to Seurat objects
# Usage: nohup ./trigger_7_seuratProcessing.sh > nohup_trigger_7_seuratProcessing.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/8_SeuratProcessing/

myscript=/scripts/0_code/7_SeuratProcessing_modParams.R

for i in B1 B2 B3 B4 C2 C3
do
echo "starting to analyse sample "$i""
Rscript "$myscript" "$i"_filtered /scripts/7_scranNormalisation/out_obj/"$i"_filtered_postNorm_Seurat.RDS /scripts/8_SeuratProcessing/out_data/ /scripts/8_SeuratProcessing/out_obj/
echo "finished analysing sample "$i""
done

