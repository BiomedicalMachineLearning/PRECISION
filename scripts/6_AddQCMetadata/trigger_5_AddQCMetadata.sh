#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: trigger_5_AddQCMetadata.sh
# Goal: To add doublet and phase data to seurat object, and filter based on doublet results
# Usage: nohup ./trigger_5_AddQCMetadata.sh > nohup_trigger_5_AddQCMetadata.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/6_AddQCMetadata/

myscript=/scripts/1b_code/5_AddQCMetadata.R

# ------------------------------------------------------------------
# Sample B1
# ------------------------------------------------------------------

echo "B1"
Rscript "$myscript" B1_filtered /scripts/3_FilterCells/out_obj/B1_filtered_postQC_Seurat.RDS /scripts/4_CellCyclePrediction/out_data/B1_filtered_Phase.txt /scripts/5_doublets/out_data/B1_filtered_doubletStats.txt /scripts/6_AddQCMetadata/out_obj/

# ------------------------------------------------------------------
# Sample B2
# ------------------------------------------------------------------

echo "B2"
Rscript "$myscript" B2_filtered /scripts/3_FilterCells/out_obj/B2_filtered_postQC_Seurat.RDS /scripts/4_CellCyclePrediction/out_data/B2_filtered_Phase.txt /scripts/5_doublets/out_data/B2_filtered_doubletStats.txt /scripts/6_AddQCMetadata/out_obj/


# ------------------------------------------------------------------
# Sample B3
# ------------------------------------------------------------------

echo "B3"
Rscript "$myscript" B3_filtered /scripts/3_FilterCells/out_obj/B3_filtered_postQC_Seurat.RDS /scripts/4_CellCyclePrediction/out_data/B3_filtered_Phase.txt /scripts/5_doublets/out_data/B3_filtered_doubletStats.txt /scripts/6_AddQCMetadata/out_obj/


# ------------------------------------------------------------------
# Sample B4
# ------------------------------------------------------------------

echo "B4"
Rscript "$myscript" B4_filtered /scripts/3_FilterCells/out_obj/B4_filtered_postQC_Seurat.RDS /scripts/4_CellCyclePrediction/out_data/B4_filtered_Phase.txt /scripts/5_doublets/out_data/B4_filtered_doubletStats.txt /scripts/6_AddQCMetadata/out_obj/


# ------------------------------------------------------------------
# Sample C2
# ------------------------------------------------------------------
echo "C2"
Rscript "$myscript" C2_filtered /scripts/3_FilterCells/out_obj/C2_filtered_postQC_Seurat.RDS /scripts/4_CellCyclePrediction/out_data/C2_filtered_Phase.txt /scripts/5_doublets/out_data/C2_filtered_doubletStats.txt /scripts/6_AddQCMetadata/out_obj/


# ------------------------------------------------------------------
# Sample C3
# ------------------------------------------------------------------

echo "C3"
Rscript "$myscript" C3_filtered /scripts/3_FilterCells/out_obj/C3_filtered_postQC_Seurat.RDS /scripts/4_CellCyclePrediction/out_data/C3_filtered_Phase.txt /scripts/5_doublets/out_data/C3_filtered_doubletStats.txt /scripts/6_AddQCMetadata/out_obj/



