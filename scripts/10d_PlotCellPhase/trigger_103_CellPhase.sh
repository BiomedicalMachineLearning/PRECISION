#!/bin/bash
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 24th June 2020
# Title: trigger_10e_PlotCellPhase.sh
# Goal: To plot cell phase data for integrated object
# Usage: nohup ./trigger_10e_PlotCellPhase.sh > nohup_trigger_10e_PlotCellPhase.out 2>&1&
# ------------------------------------------------------------------

cd /scripts/2_output/10e_integrated_phase/
myscript=/scripts/1b_code/10e_PlotCellPhase.R

echo "starting to analyse sample"
Rscript "$myscript" all integrated_snn_res.0.5 /scripts/10c_integrated_clustering/out_obj/allintegratedClustered_Seurat.RDS /scripts/10d_PlotCellPhase/out_data/
echo "finished analysing sample"

