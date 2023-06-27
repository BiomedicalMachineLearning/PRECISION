#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: trigger_2_FilterCells.sh
# Goal: To convert Cellranger output to Seurat objects
# Usage: nohup ./trigger_2_FilterCells.sh > nohup_trigger_2_FilterCells.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/3_FilterCells/

myscript=/scripts/1b_code/2_FilterCells.R

# ------------------------------------------------------------------
# Sample B1
# ------------------------------------------------------------------

echo "B1"
# filtered
Rscript "$myscript" B1_filtered /scripts/1_ChangeHumanNames/outdir/B1.RDS /scripts/2_QC_3MAD/out_data/B1_filtered_cellQC_3MADkeptcells.txt /scripts/3_FilterCells/out_obj/

# ------------------------------------------------------------------
# Sample B2
# ------------------------------------------------------------------
echo "B2"
# filtered
Rscript "$myscript" B2_filtered /scripts/1_ChangeHumanNames/outdir/B2.RDS /scripts/2_QC_3MAD/out_data/B2_filtered_cellQC_3MADkeptcells.txt /scripts/3_FilterCells/out_obj/

# ------------------------------------------------------------------
# Sample B3
# ------------------------------------------------------------------
echo "B3"
# filtered
Rscript "$myscript" B3_filtered /scripts/1_ChangeHumanNames/outdir/B3.RDS /scripts/2_QC_3MAD/out_data/B3_filtered_cellQC_3MADkeptcells.txt /scripts/3_FilterCells/out_obj/

# ------------------------------------------------------------------
# Sample B4
# ------------------------------------------------------------------
echo "B4"
# filtered
Rscript "$myscript" B4_filtered /scripts/1_ChangeHumanNames/outdir/B4.RDS /scripts/2_QC_3MAD/out_data/B4_filtered_cellQC_3MADkeptcells.txt /scripts/3_FilterCells/out_obj/

# ------------------------------------------------------------------
# Sample C3
# ------------------------------------------------------------------
echo "C3"
# filtered
Rscript "$myscript" C3_filtered /scripts/1_ChangeHumanNames/outdir/C3.RDS /scripts/2_QC_3MAD/out_data/C3_filtered_cellQC_3MADkeptcells.txt /scripts/3_FilterCells/out_obj/

# ------------------------------------------------------------------
# Sample C2
# ------------------------------------------------------------------
echo "C2"
# filtered
Rscript "$myscript" C2_filtered /scripts/1_ChangeHumanNames/outdir/C2.RDS /scripts/2_QC_3MAD/out_data/C2_filtered_cellQC_3MADkeptcells.txt /scripts/3_FilterCells/out_obj/


