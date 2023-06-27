#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R

# analysis 4.2 vs 4.3 (ivig only)
Rscript "$myscript" cf_4.2_4.3_ivigOnly all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.2v4.3_IvigOnly_ComparisonVenn_Markers_RemoveRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.2vs4.3_ivigOnly

echo "done 4.2 vs 4.3 (ivig only)"

# analysis 4.2 vs 4.3 (veh only)
Rscript "$myscript" cf_4.2_4.3_vehOnly all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.2v4.3_VehOnly_ComparisonVenn_Markers_RemoveRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.2vs4.3_vehOnly

echo "done 4.2 vs 4.3 (veh only)"
