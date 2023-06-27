#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R

# analysis 4.5
Rscript "$myscript" clust_4.5 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.5_vehVSivig_cord_Markers_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.5

echo "done 4.5"

# analysis 4.4
Rscript "$myscript" clust_4.4 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.4_vehVSivig_blood_Markers_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.4

echo "done 4.4"

# analysis 4.3
Rscript "$myscript" clust_4.3 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.3_bloodVScord_ivig_Markers_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.3

echo "done 4.3"

# analysis 4.2
Rscript "$myscript" clust_4.2 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.2_bloodVScord_veh_Markers_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.2

echo "done 4.2"

# analysis 4.1 (each analysis kept separate)
Rscript "$myscript" clust_4.1_separate all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.1_SCIvsshamnaive_blood_Markers_KeepSeparate_dropRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.1_3ComparisonsSeparate

echo "done 4.1 (keep separate)"

# analysis 4.1 (make venn of overlaps)
Rscript "$myscript" clust_4.1_overlaps all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.1_SCIvsshamnaive_blood_ComparisonVenn_Markers_RemoveRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.1_VennOverlaps

echo "done 4.1 (venn overlaps)"

# analysis 4.2 vs 4.3
Rscript "$myscript" cf_4.2_4.3 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.2v4.3_ComparisonVenn_Markers_RemoveRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.2vs4.3

echo "done 4.2 vs 4.3"

# analysis 4.4 vs 4.5
Rscript "$myscript" cf_4.4_4.5 all /scripts/13b_ProcessForGO/Clusterwise/outdir/4.4v4.5_ComparisonVenn_Markers_RemoveRbMt.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/outdir_4.4vs4.5

echo "done 4.4 vs 4.5"