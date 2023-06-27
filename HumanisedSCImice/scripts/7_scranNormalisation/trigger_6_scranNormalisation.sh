#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: trigger_6_scranNormalisation.sh
# Goal: To convert Cellranger output to Seurat objects
# Usage: nohup ./trigger_6_scranNormalisation.sh > nohup_trigger_6_scranNormalisation.out 2>&1&
# ------------------------------------------------------------------

myscript=/scripts/0_code/6_scranNormalisation_RenameAssay.R
cd /scripts/7_scranNormalisation/out_obj

for i in B1 B2 B3 B4 C2 C3
do
Rscript "$myscript" "$i"_filtered /scripts/6_AddQCMetadata/out_obj/"$i"_filtered_preNorm_Seurat.RDS /scripts/7_scranNormalisation/out_obj
done

