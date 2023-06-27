#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R


# analysis 4.6 (SCI only)
Rscript "$myscript" BothSCIvsEitherCntl all /scripts/13b_ProcessForGO/Samplewise/4.6/outdir/BothSCIvsEitherCntl_MarkersForGO.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Samplewise_AllGenes_RemoveRbMt/4.6/outdir_BothSCIvsEitherCntl


Rscript "$myscript" EitherSCIvsBothCntl all /scripts/13b_ProcessForGO/Samplewise/4.6/outdir/EitherSCIvsBothCntl_MarkersForGO.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Samplewise_AllGenes_RemoveRbMt/4.6/outdir_EitherSCIvsBothCntl

echo "done 4.6"
