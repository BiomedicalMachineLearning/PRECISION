#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 30th June 2020
# Title: trigger_11.2_Seurat2CellRanger.sh
# Goal: To convert the count matrix to a fake CellRanger-style folder hierarchy
# Usage: nohup ./trigger_11.2_Seurat2CellRanger.sh > nohup_trigger_11.2_Seurat2CellRanger.out 2>&1&
# Example:
# Rscript 11.2_Seurat2CellRanger.R {sampleID} {seurat.path} {output.dir}

# ------------------------------------------------------------------

cd /scripts/17_Metacell/vehicle/11.2_Seurat2CellRanger/

# Vehicle
Rscript /scripts/1b_code/11_metacell/11.2_Seurat2CellRanger.R vehicle /scripts/17_Metacell/vehicle/11.1_ExtractClusters/out/vehicle_subset.RDS /scripts/17_Metacell/vehicle/11.2_Seurat2CellRanger/out_obj/
