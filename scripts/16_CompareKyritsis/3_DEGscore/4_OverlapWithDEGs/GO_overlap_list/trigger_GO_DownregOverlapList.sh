#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R


# analysis 4.1 (SCI only)
Rscript "$myscript" downreg all /scripts/16_CompareKyritsis/3_DEGscore/4_OverlapWithDEGs/outdir/sharedGenes.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/16_CompareKyritsis/3_DEGscore/4_OverlapWithDEGs/GO_overlap_list/outdir

echo "done 4.1 (SCI only)"
