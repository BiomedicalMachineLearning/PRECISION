#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R


# analysis 4.1 (SCI only)
Rscript "$myscript" cf_traumaGenes all /scripts/13b_ProcessForGO/Samplewise/TraumaInducedGenes/outdir/TraumaGenes.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Samplewise_AllGenes_RemoveRbMt/TraumaGenes/outdir

echo "done trauma genes (samplewise)"
