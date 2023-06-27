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

cd /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/

# T cells all
Rscript /scripts/1b_code/11_metacell/11.2_Seurat2CellRanger.R Tcells /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.1_ExtractClusters/out/Tcells_subset.RDS /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/

# T cells (CD4)
Rscript /scripts/1b_code/11_metacell/11.2_Seurat2CellRanger.R CD4Tcells /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.1_ExtractClusters/out/CD4_Tcells_subset.RDS /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/

# T cells (CD8)
Rscript /scripts/1b_code/11_metacell/11.2_Seurat2CellRanger.R CD8Tcells /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.1_ExtractClusters/out/CD8_Tcells_subset.RDS /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/

# B cells
Rscript /scripts/1b_code/11_metacell/11.2_Seurat2CellRanger.R Bcells /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.1_ExtractClusters/out/Bcells_subset.RDS /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/

# Granulocytes
Rscript /scripts/1b_code/11_metacell/11.2_Seurat2CellRanger.R Granulocytes /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.1_ExtractClusters/out/granulocytes_subset.RDS /scripts/17_MetacellVehicleOnly/vehicle_subtypes/11.2_Seurat2CellRanger/out_obj/