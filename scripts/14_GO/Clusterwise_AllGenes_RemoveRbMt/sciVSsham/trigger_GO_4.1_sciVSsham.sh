#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R

# analysis 4.2 vs 4.3 (ivig only)
Rscript "$myscript" cf_4.1_sciVSsham all /scripts/13b_ProcessForGO/Clusterwise/sciVSsham/outdir/Clusterwise_SCIvsSham.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/sciVSsham/outdir

echo "done traumagenes"
