#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: trigger_4_DoubletPrediction.sh
# Goal: To predict doublets
# Usage: nohup ./trigger_4_DoubletPrediction.sh > nohup_trigger_4_DoubletPrediction.out 2>&1&
# ------------------------------------------------------------------

myscript=/scripts/1b_code/4_DoubletPrediction.R

# ------------------------------------------------------------------
# Sample B1
# ------------------------------------------------------------------
echo "B1"
Rscript "$myscript" B1_filtered /scripts/3_FilterCells/out_obj/B1_filtered_postQC_Seurat.RDS /scripts/5_doublets/out_data

# ------------------------------------------------------------------
# Sample B2
# ------------------------------------------------------------------

echo "B2"
Rscript "$myscript" B2_filtered /scripts/3_FilterCells/out_obj/B2_filtered_postQC_Seurat.RDS /scripts/5_doublets/out_data


# ------------------------------------------------------------------
# Sample B3
# ------------------------------------------------------------------

echo "B3"
Rscript "$myscript" B3_filtered /scripts/3_FilterCells/out_obj/B3_filtered_postQC_Seurat.RDS /scripts/5_doublets/out_data

# ------------------------------------------------------------------
# Sample B4
# ------------------------------------------------------------------

echo "B4"
Rscript "$myscript" B4_filtered /scripts/3_FilterCells/out_obj/B4_filtered_postQC_Seurat.RDS /scripts/5_doublets/out_data

# ------------------------------------------------------------------
# Sample C2
# ------------------------------------------------------------------

echo "C2"
Rscript "$myscript" C2_filtered /scripts/3_FilterCells/out_obj/C2_filtered_postQC_Seurat.RDS /scripts/5_doublets/out_data

# ------------------------------------------------------------------
# Sample C3
# ------------------------------------------------------------------

echo "C3"
Rscript "$myscript" C3_filtered /scripts/3_FilterCells/out_obj/C3_filtered_postQC_Seurat.RDS /scripts/5_doublets/out_data


