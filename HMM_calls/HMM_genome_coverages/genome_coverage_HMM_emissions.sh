#!/bin/sh

cd ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/postHMM_filter/genome_coverage_data/

multiIntersectBed -i ../FILTER1_RdDM_1read.bed ../FILTER1_other.bed ../FILTER1_PVS_antisense_1read.bed > RdDM_other_PVS_intersect.bed

awk '{if((($6==1)&&($4==1))||(($6==1)&&($8==1))) print $0}' RdDM_other_PVS_intersect.bed > RdDM_nonOtherRegions_postFilter.bed
awk '{if(($4==1)&&($8==1)) print $0}' RdDM_other_PVS_intersect.bed > PVS_onlyRegions_postFilter.bed

awk '{if(($5==2)||($5=="2,3")||($5=="1,2")) print $0}' RdDM_other_PVS_intersect.bed > other_regions_postFilter.bed

sortBed -i RdDM_nonOtherRegions_postFilter.bed | genomeCoverageBed -i stdin -g ~/cifs-lab/Shriya/Sequences/chrom_sizes.txt | awk '($1=="genome") && ($2==0)' 
#genome	0	94898830	119146348	0.79649
##20.351%
sortBed -i PVS_onlyRegions_postFilter.bed | genomeCoverageBed -i stdin -g ~/cifs-lab/Shriya/Sequences/chrom_sizes.txt | awk '($1=="genome") && ($2==0)' 
#genome	0	92830348	119146348	0.779129
##22.087%
sortBed -i other_regions_postFilter.bed | genomeCoverageBed -i stdin -g ~/cifs-lab/Shriya/Sequences/chrom_sizes.txt | awk '($1=="genome") && ($2==0)' 
#genome	0	50563518	119146348	0.424382
##57.562%
