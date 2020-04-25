#!/bin/sh


# cd ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling

bedtools makewindows -g ~/cifs-lab/Shriya/Sequences/chrom_sizes.txt -w 200 -s 150 > 200bp_150bp_sliding_windows.bed ## Creating bins 200bp wide and 150bp sliding


awk '{OFS="\t"}{print $0,"0","0","+"}' 200bp_150bp_sliding_windows.bed | intersectBed -c -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | intersectBed -c -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/nrpe1_3reps_sort_dedup_filtChr.bed | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,($7+0.1)/($8+0.1)}' > 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio.bed ## Calculating the number of combined Col-0 and nrpe1 reads in each bin irrespective of strands. Also, calculating the ration of col0/nrpe1 in each bin.

coverageBed -c -s -a 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio.bed -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | intersectBed -c -S -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | awk '{OFS="\t"}{if($10>$11) print $0,($10+0.1)/($11+0.1); else print $0,($11+0.1)/($10+0.1)}' > 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio.bed ## Calculating the number of reads in Col-0 on the +strand and -strand. Also, calculating the positive ratio of the strand bias of reads in each bin.

intersectBed -wao -a 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio.bed -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_col0_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8,9,12 -c 26 -o mean | awk '{OFS="\t"}{if($11==".") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,"NA"; else print $0}' > temp1 ## Calculating the average Col-0 CHHme in each bin. If there is no CHHme in a bin, reporting NA. METHYL-DATA: our lab.. 2reps.. methylation needs to be present in both.. C considered methylated if it has atleast 5reads at the site..
intersectBed -wao -a temp1 -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_nrpe1_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8,9,10,11 -c 25 -o mean | awk '{OFS="\t"}{if($12==".") print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,"NA"; else print $0}' > 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_CHHcolNrpe_BothContainConvert.bed ## Calculating the average nrpe1 CHHme in each bin. If there is no CHHme in a bin, reporting NA.

intersectBed -c -a 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio.bed -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646727_col0_adapter.bed | intersectBed -c -a stdin -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646746_nrpd1_adapter.bed > 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_smRNACol0Nrpd.bed ## Calculating the number of Col-0 and nrpd1 smRNA overlapping each bin. DATASETS: Jacobsen.
awk '{OFS="\t"}{print $13,$14}' 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_smRNACol0Nrpd.bed > bins_smRNA

paste 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_CHHcolNrpe_BothContainConvert.bed bins_smRNA > 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_CHHcolNrpe_BothContainConvert_smRNACol0Nrpd.bed 

intersectBed -wa -a ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed -b ~/cifs-lab/Shriya/Sequences/genes_list.bed | uniq | wc -l ## Finding the number of reads in Col-0 that overlap genes
##2684910
intersectBed -wa -a ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/nrpe1_3reps_sort_dedup_filtChr.bed -b ~/cifs-lab/Shriya/Sequences/genes_list.bed | uniq | wc -l ## Finding the number of reads in nrpe1 that overlap genes
##2052106

#gene-normalization ratio: 2052106/2684910=0.77
#gene-normalization ratio: 2023850/2416589=0.84

awk '{OFS="\t"}{print $9*0.84,$10,$11-$12}' 200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_CHHcolNrpe_BothContainConvert.bed > features_geneNORMcolNrpe_antiSense_methDiff_BothContainConvert.bed

