#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: trigger_10b_scr7_seuratProcessing.sh
# Goal: To convert Cellranger output to Seurat objects
# Usage: nohup ./trigger_10b_scr7_seuratProcessing.sh > nohup_trigger_10b_scr7_seuratProcessing.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/10b_integrated_seuratProcessing/
myscript=/scripts/0_code/7_SeuratProcessing_modParams.R

echo "starting to analyse integrated sample"
Rscript "$myscript" all /scripts/10_integration/out_obj/all_integrated.RDS /scripts/10b_integrated_seuratProcessing/out_data/ /scripts/10b_integrated_seuratProcessing/out_obj/
echo "finished analysing integrated sample"

