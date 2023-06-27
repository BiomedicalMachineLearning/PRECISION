#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R


# analysis 4.2 vs 4.3 (ivig only)
Rscript "$myscript" cf_4.2_4.3_ivig all /scripts/13b_ProcessForGO/Samplewise/outdir/4.2v4.3_ivigonly_ComparisonVenn_Markers_RemoveRbMt_ivigOnly.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Samplewise_AllGenes_RemoveRbMt/outdir_4.2vs4.3_ivigOnly

echo "done 4.2 vs 4.3 (ivig only)"

# analysis 4.2 vs 4.3 (veh only)
Rscript "$myscript" cf_4.2_4.3_veh all /scripts/13b_ProcessForGO/Samplewise/outdir/4.2v4.3_vehonly_ComparisonVenn_Markers_RemoveRbMt_vehOnly.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Samplewise_AllGenes_RemoveRbMt/outdir_4.2vs4.3_vehOnly

echo "done 4.2 vs 4.3 (veh only)"