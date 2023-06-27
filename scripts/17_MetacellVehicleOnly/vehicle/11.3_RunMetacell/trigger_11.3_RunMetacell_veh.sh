#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: trigger_11.3_RunMetacell.sh
# Goal: To convert the count matrix to a fake CellRanger-style folder hierarchy
# Usage: nohup ./trigger_11.3_RunMetacell.sh > nohup_trigger_11.3_RunMetacell.out 2>&1&
# Example:
# Rscript 11.3_RunMetacell.R {sampleID} {CellRangerDir} {outdir_fig} {outdir_obj}
# ------------------------------------------------------------------

cd /scripts/17_Metacell/vehicle/11.3_RunMetacell/

# Run with default settings
## vehicle

echo "starting vehicle default at $(date)"
Rscript /scripts/1b_code/11_metacell/11.3_RunMetacell.R vehicle /scripts/17_Metacell/vehicle/11.2_Seurat2CellRanger/out_obj/vehicle_CellRanger/ /scripts/17_Metacell/vehicle/11.3_RunMetacell/out_fig/ /scripts/17_Metacell/vehicle/11.3_RunMetacell/out_obj/
echo "ending vehicle default at $(date)"
