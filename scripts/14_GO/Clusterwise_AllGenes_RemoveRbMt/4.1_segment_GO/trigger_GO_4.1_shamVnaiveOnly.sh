#! /bin/bash

# trigger all 

#Rscript RunPathGO.R {sampleID} {all/top100} {DE list} {seuratObj} {outdir}

myscript=/scripts/14_GO/RunGO.R

# analysis 1: SCI vs sham exclusive
Rscript "$myscript" cf_4.1_SCIvSham_exclusive all /scripts/13b_ProcessForGO/Clusterwise/4.1_splitIntoVennSegments/PutClustersTogether_outdir/SegmentsForGO_SCIvsham.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/4.1_segment_GO/outdir_SCIvShamExclusive

echo "done 4.1 (SCI vs sham exclusive)"


# analysis 2: SCI vs naive exclusive
Rscript "$myscript" cf_4.1_SCIvNaive_exclusive all /scripts/13b_ProcessForGO/Clusterwise/4.1_splitIntoVennSegments/PutClustersTogether_outdir/SegmentsForGO_SCIvnaive.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/4.1_segment_GO/outdir_SCIvNaiveExclusive

echo "done 4.1 (SCI vs naive exclusive)"

# analysis 2: SCI vs sham and SCI vs naive overlap
Rscript "$myscript" cf_4.1_SCIvSham+SCIvNaive_overlap all /scripts/13b_ProcessForGO/Clusterwise/4.1_splitIntoVennSegments/PutClustersTogether_outdir/SegmentsForGO_SCIvnaive_SCIvsham.txt /scripts/11b_RunClusteringRes0.4/out_obj/allintegratedClustered_Seurat2.RDS /scripts/14_GO/Clusterwise_AllGenes_RemoveRbMt/4.1_segment_GO/outdir_SCIvSham+SCIvNaiveOverlap

echo "done 4.1 (SCI vs sham and SCI vs naive overlap)"