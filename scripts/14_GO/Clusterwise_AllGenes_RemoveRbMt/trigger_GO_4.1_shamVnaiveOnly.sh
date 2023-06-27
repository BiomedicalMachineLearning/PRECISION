#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R

# analysis 4.1 (sham v naive only)
Rscript "$myscript" cf_4.1_shamVnaive all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.1_shamVnaive_blood_Markers_KeepSeparate_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.1_shamVnaiveOnly

echo "done 4.1 (sham v naive only)"
