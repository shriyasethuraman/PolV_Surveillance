#!/bin/sh

cd ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/
paste ../200bp_150bp_sliding_windows.bed ../hmm_1_output_full_3feature_4st_op3_geneNORM_BothContainConvert.txt > HMM_locus_plus_minus_emissions_3feature_4st_op3_geneNORM_BothContainConvert.bed


awk '{OFS="\t"} $4==0' HMM_locus_plus_minus_emissions_3feature_4st_op3_geneNORM_BothContainConvert.bed > output_em0.bed
awk '{OFS="\t"} $4==1' HMM_locus_plus_minus_emissions_3feature_4st_op3_geneNORM_BothContainConvert.bed > output_em1.bed
awk '{OFS="\t"} $4==2' HMM_locus_plus_minus_emissions_3feature_4st_op3_geneNORM_BothContainConvert.bed > output_em2.bed
awk '{OFS="\t"} $4==3' HMM_locus_plus_minus_emissions_3feature_4st_op3_geneNORM_BothContainConvert.bed > output_em3.bed

mergeBed -d 3 -i output_em0.bed | awk '{OFS="\t"}{print $0,"0","0","+"}' > output_em0_2.bed
mergeBed -d 3 -i output_em1.bed | awk '{OFS="\t"}{print $0,"0","0","+"}' > output_em1_2.bed
mergeBed -d 3 -i output_em2.bed | awk '{OFS="\t"}{print $0,"0","0","+"}' > output_em2_2.bed
mergeBed -d 3 -i output_em3.bed | awk '{OFS="\t"}{print $0,"0","0","+"}' > output_em3_2.bed


intersectBed -c -a output_em0_2.bed -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | intersectBed -c -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/nrpe1_3reps_sort_dedup_filtChr.bed > em0_col0_nrpe1_RIP.bed
intersectBed -c -a output_em1_2.bed -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | intersectBed -c -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/nrpe1_3reps_sort_dedup_filtChr.bed > em1_col0_nrpe1_RIP.bed
intersectBed -c -a output_em2_2.bed -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | intersectBed -c -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/nrpe1_3reps_sort_dedup_filtChr.bed > em2_col0_nrpe1_RIP.bed
intersectBed -c -a output_em3_2.bed -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/col0_3reps_sort_dedup_filtChr.bed | intersectBed -c -a stdin -b ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/mapped_combined_reads/nrpe1_3reps_sort_dedup_filtChr.bed > em3_col0_nrpe1_RIP.bed


~/Downloads/bedtools2/bin/intersectBed -wao -a em0_col0_nrpe1_RIP.bed -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_col0_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 22 -o mean | awk '{OFS="\t"}{if($9==".") print $1,$2,$3,$4,$5,$6,$7,$8,"NA"; else print $0}' > temp1
~/Downloads/bedtools2/bin/intersectBed -wao -a temp1 -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_nrpe1_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8,9 -c 23 -o mean | awk '{OFS="\t"}{if($10==".") print $1,$2,$3,$4,$5,$6,$7,$8,$9,"NA"; else print $0}' > em0_col0_nrpe1_RIP_col0-nrpe1_CHH.bed
~/Downloads/bedtools2/bin/intersectBed -wao -a em1_col0_nrpe1_RIP.bed -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_col0_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 22 -o mean | awk '{OFS="\t"}{if($9==".") print $1,$2,$3,$4,$5,$6,$7,$8,"NA"; else print $0}' > temp1
~/Downloads/bedtools2/bin/intersectBed -wao -a temp1 -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_nrpe1_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8,9 -c 23 -o mean | awk '{OFS="\t"}{if($10==".") print $1,$2,$3,$4,$5,$6,$7,$8,$9,"NA"; else print $0}' > em1_col0_nrpe1_RIP_col0-nrpe1_CHH.bed
~/Downloads/bedtools2/bin/intersectBed -wao -a em2_col0_nrpe1_RIP.bed -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_col0_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 22 -o mean | awk '{OFS="\t"}{if($9==".") print $1,$2,$3,$4,$5,$6,$7,$8,"NA"; else print $0}' > temp1
~/Downloads/bedtools2/bin/intersectBed -wao -a temp1 -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_nrpe1_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8,9 -c 23 -o mean | awk '{OFS="\t"}{if($10==".") print $1,$2,$3,$4,$5,$6,$7,$8,$9,"NA"; else print $0}' > em2_col0_nrpe1_RIP_col0-nrpe1_CHH.bed
intersectBed -wao -a em3_col0_nrpe1_RIP.bed -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_col0_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8 -c 22 -o mean | awk '{OFS="\t"}{if($9==".") print $1,$2,$3,$4,$5,$6,$7,$8,"NA"; else print $0}' > temp1
intersectBed -wao -a temp1 -b ~/cifs-lab/Shriya/Samples/methylation_data/CHH_nrpe1_ALL_2rep_2_BothContainConvert.bed | ~/Downloads/bedtools2/bin/groupBy -i stdin -g 1,2,3,4,5,6,7,8,9 -c 23 -o mean | awk '{OFS="\t"}{if($10==".") print $1,$2,$3,$4,$5,$6,$7,$8,$9,"NA"; else print $0}' > em3_col0_nrpe1_RIP_col0-nrpe1_CHH.bed

~/Downloads/bedtools2/bin/intersectBed -c -a em0_col0_nrpe1_RIP.bed -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646727_col0_adapter.bed | intersectBed -c -a stdin -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646746_nrpd1_adapter.bed > em0_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed
~/Downloads/bedtools2/bin/intersectBed -c -a em1_col0_nrpe1_RIP.bed -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646727_col0_adapter.bed | intersectBed -c -a stdin -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646746_nrpd1_adapter.bed > em1_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed
~/Downloads/bedtools2/bin/intersectBed -c -a em2_col0_nrpe1_RIP.bed -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646727_col0_adapter.bed | intersectBed -c -a stdin -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646746_nrpd1_adapter.bed > em2_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed
~/Downloads/bedtools2/bin/intersectBed -c -a em3_col0_nrpe1_RIP.bed -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646727_col0_adapter.bed | intersectBed -c -a stdin -b ~/cifs-lab/Masayuki/RIP_Pol_V_Terminator/sRNA_Jacobsen/SRR5646746_nrpd1_adapter.bed > em3_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed

## for heatmaps
paste ../200bp_150bp_sliding_windows_UNstranded_col0_nrpe_ratio_plus_minus_SASRatio_CHHcolNrpe_BothContainConvert_smRNACol0Nrpd.bed ../hmm_1_output_full_3feature_4st_op3_geneNORM_BothContainConvert.txt | awk '{OFS="\t"}{print $15,$7,$8,$7+$8,$9,$10,$13,$14,$11,$12,$11-$12}' | sort -k1 -n > outputStates_inputs_list.txt

paste ../hmm_1_output_full_3feature_4st_op3_geneNORM_BothContainConvert.txt ../log_likelihood_enission_values.txt | sort -k1 -n > outputStates_logLikelihood_list.txt

