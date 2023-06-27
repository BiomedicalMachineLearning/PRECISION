#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO_AllNAOneClust.R

# analysis 4.5
Rscript "$myscript" clust_4.5 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.5_vehVSivig_cord_Markers_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.5

echo "done 4.5"
