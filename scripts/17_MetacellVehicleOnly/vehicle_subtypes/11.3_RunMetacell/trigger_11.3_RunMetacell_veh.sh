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

## T cells (all)

echo "starting T cells at $(date)"
Rscript /scripts/1b_code/11_metacell/11.3_RunMetacell.R Tcells /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/Tcells_CellRanger/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Tcells/out_fig/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Tcells/out_obj/
echo "ending T cells default at $(date)"

## T cells (CD4)

echo "starting CD4 T cells at $(date)"
Rscript /scripts/1b_code/11_metacell/11.3_RunMetacell.R Tcells_CD4 /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/CD4Tcells_CellRanger/  /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Tcells_CD4/out_fig/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Tcells_CD4/out_obj/
echo "ending T cells default at $(date)"

## T cells (CD8)

echo "starting CD8 T cells at $(date)"
Rscript /scripts/1b_code/11_metacell/11.3_RunMetacell.R Tcells_CD8 /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/CD8Tcells_CellRanger/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Tcells_CD8/out_fig/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Tcells_CD8/out_obj/
echo "ending T cells default at $(date)"

## B cells (all)

echo "starting B cells at $(date)"
Rscript /scripts/1b_code/11_metacell/11.3_RunMetacell.R Bcells /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/Bcells_CellRanger/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Bcells/out_fig/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Bcells/out_obj/
echo "ending B cells default at $(date)"

## Granulocytes

echo "starting Granulocytes at $(date)"
Rscript /scripts/1b_code/11_metacell/11.3_RunMetacell.R granulocytes /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/Granulocytes_CellRanger/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Granulocytes/out_fig/ /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.3_RunMetacell/Granulocytes/out_obj/
echo "ending Granulocytes default at $(date)"