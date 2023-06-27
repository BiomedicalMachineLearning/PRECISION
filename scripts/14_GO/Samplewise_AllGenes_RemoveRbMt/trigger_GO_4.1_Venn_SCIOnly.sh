#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R


# analysis 4.1 (SCI only)
Rscript "$myscript" cf_4.1_SCI all /scripts/13b_ProcessForGO/Samplewise/outdir/4.1_SCIvsshamnaive_blood_ComparisonVenn_Markers_RemoveRbMt_SCIonly.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Samplewise_AllGenes_RemoveRbMt/outdir_4.1_VennOverlaps_SCIOnly

echo "done 4.1 (SCI only)"
