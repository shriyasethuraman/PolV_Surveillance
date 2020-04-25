#!/bin/sh 

cd ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/ddm1_keith/

intersectBed -wa -a ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d10_q01.bed -b /ssd/shriya_data/remap_Masa_RIP/HMM_200bp_150slide/output3_geneNORM_BothContainConvert/postHMM_filter/PVS_filtered1_loci.bed | uniq > PVS_ddm1Col_DMR_overlap.bed
intersectBed -wa -a ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d10_q01.bed -b /ssd/shriya_data/remap_Masa_RIP/HMM_200bp_150slide/output3_geneNORM_BothContainConvert/postHMM_filter/RdDM_nonOtherRegions_postFilter.bed | uniq > RdDM_ddm1Col_DMR_overlap.bed 
intersectBed -wa -a ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d10_q01.bed -b /ssd/shriya_data/remap_Masa_RIP/HMM_200bp_150slide/output3_geneNORM_BothContainConvert/postHMM_filter/other_regions_postFilter.bed | uniq > other_ddm1Col_DMR_overlap.bed

awk '{OFS="\t"}{print "chr"$1,$2,$3}' PVS_ddm1Col_DMR_overlap.bed > PVS_ddm1Col_DMR_overlap_methylKit.bed
awk '{OFS="\t"}{print "chr"$1,$2,$3}' RdDM_ddm1Col_DMR_overlap.bed > RdDM_ddm1Col_DMR_overlap_methylKit.bed
awk '{OFS="\t"}{print "chr"$1,$2,$3}' other_ddm1Col_DMR_overlap.bed > other_ddm1Col_DMR_overlap_methylKit.bed


