#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 23rd June 2020
# Title: trigger_1_QC_3MAD.sh
# Goal: To convert Cellranger output to Seurat objects
# Usage: nohup ./trigger_1_QC_3MAD.sh > nohup_trigger_1_QC_3MAD.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/2_QC_3MAD/

myscript=/scripts/0_code/1_QC_3MAD_Hsap.R

# ------------------------------------------------------------------
# Sample B1
# ------------------------------------------------------------------

echo "B1"

# filtered
Rscript "$myscript" B1_filtered /scripts/1b_ConvertToSCE/outdir/B1_sce.RDS /scripts/2_QC_3MAD/out_data

# ------------------------------------------------------------------
# Sample B2
# ------------------------------------------------------------------
echo "B2"
Rscript "$myscript" B2_filtered /scripts/1b_ConvertToSCE/outdir/B2_sce.RDS /scripts/2_QC_3MAD/out_data

# ------------------------------------------------------------------
# Sample B3
# ------------------------------------------------------------------
echo "B3"
Rscript "$myscript" B3_filtered /scripts/1b_ConvertToSCE/outdir/B3_sce.RDS /scripts/2_QC_3MAD/out_data

# ------------------------------------------------------------------
# Sample B4
# ------------------------------------------------------------------
echo "B4"
Rscript "$myscript" B4_filtered /scripts/1b_ConvertToSCE/outdir/B4_sce.RDS /scripts/2_QC_3MAD/out_data

# ------------------------------------------------------------------
# Sample C2
# ------------------------------------------------------------------
echo "C2"
Rscript "$myscript" C2_filtered /scripts/1b_ConvertToSCE/outdir/C2_sce.RDS /scripts/2_QC_3MAD/out_data

# ------------------------------------------------------------------
# Sample C3
# ------------------------------------------------------------------
echo "C3"
Rscript "$myscript" C3_filtered /scripts/1b_ConvertToSCE/outdir/C3_sce.RDS /scripts/2_QC_3MAD/out_data

