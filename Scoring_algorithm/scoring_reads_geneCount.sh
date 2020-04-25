#!/bin/sh

#cd /ssd/shriya_data/remap_Masa_RIP/scoring_reads/

intersectBed -s -c -a /ssd/shriya_data/remap_Masa_RIP/self_trim2/col0_3reps_sort_dedup_filtChr.bed -b ~/cifs-lab/Shriya/Sequences/genes_list.bed > col0_minusTer_sort_dedup_geneCount.bed 
## mapped dedup sorted bed file with 7columns:1-chr#,2-start,3-stop,4-read name,5-random,6-strand,7-gene overlap binary

awk '$1!="Pt" && $1!="Mt"' col0_minusTer_sort_dedup_geneCount.bed | sort -k1 -k2 -n > sorted_col0_minusTer_sort_dedup_geneCount.bed ## producing sorted input file


awk '$1==1' sorted_col0_minusTer_sort_dedup_geneCount.bed > chr1_col0_input.bed
awk '$1==2' sorted_col0_minusTer_sort_dedup_geneCount.bed > chr2_col0_input.bed
awk '$1==3' sorted_col0_minusTer_sort_dedup_geneCount.bed > chr3_col0_input.bed
awk '$1==4' sorted_col0_minusTer_sort_dedup_geneCount.bed > chr4_col0_input.bed
awk '$1==5' sorted_col0_minusTer_sort_dedup_geneCount.bed > chr5_col0_input.bed


coverageBed -s -d -a ~/cifs-lab/Shriya/Sequences/chrom_lengths.bed -b /ssd/shriya_data/remap_Masa_RIP/self_trim2/col0_3reps_sort_dedup_filtChr.bed | awk '{OFS="\t"}{print $1,$2+$7-1,$2+$7,$4,$5,$6,$8}' > sorted_input_count_temp.bed
intersectBed -s -c -a sorted_input_count_temp.bed -b /ssd/shriya_data/remap_Masa_RIP/self_trim2/nrpe1_3reps_sort_dedup_filtChr.bed | awk '{OFS="\t"}{if($7>0 || $8>0) print $1,$2,$6,$7,$8}' > sorted_input_chr_loci_counts_nonzero.bed
## file with the genome-wide coverage of each dataset (Col0,Nrpe1): 1-chr#, 2-position, 3-strand, 4-Col0 count, 5-nrpe1 count


### Running the scoring algorithm for Col-0-wildtype input
python identify_nonPol2_reads_opt4.py chr1_col0_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_chr1.bed
python identify_nonPol2_reads_opt4.py chr2_col0_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_chr2.bed
python identify_nonPol2_reads_opt4.py chr3_col0_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_chr3.bed
python identify_nonPol2_reads_opt4.py chr4_col0_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_chr4.bed
python identify_nonPol2_reads_opt4.py chr5_col0_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_chr5.bed

## Adding the scores of each read to each read line of the input
paste chr1_col0_input.bed scoring_4_chr1.bed > chr1_reads_scores.bed
paste chr2_col0_input.bed scoring_4_chr2.bed > chr2_reads_scores.bed
paste chr3_col0_input.bed scoring_4_chr3.bed > chr3_reads_scores.bed
paste chr4_col0_input.bed scoring_4_chr4.bed > chr4_reads_scores.bed
paste chr5_col0_input.bed scoring_4_chr5.bed > chr5_reads_scores.bed

cat chr1_reads_scores.bed chr2_reads_scores.bed chr3_reads_scores.bed chr4_reads_scores.bed chr5_reads_scores.bed > chr#_read_scores.bed


#### Redoing the same for the nrpe1 mutant as we did for Col-0 wild type
intersectBed -s -c -a /ssd/shriya_data/remap_Masa_RIP/self_trim2/nrpe1_3reps_sort_dedup_filtChr.bed -b ~/cifs-lab/Shriya/Sequences/genes_list.bed > nrpe1_minusTer_sort_dedup_geneCount.bed

awk '$1!="Pt" && $1!="Mt"' nrpe1_minusTer_sort_dedup_geneCount.bed | sort -k1 -k2 -n > sorted_nrpe1_minusTer_sort_dedup_geneCount.bed 


awk '$1==1' sorted_nrpe1_minusTer_sort_dedup_geneCount.bed > chr1_nrpe1_input.bed
awk '$1==2' sorted_nrpe1_minusTer_sort_dedup_geneCount.bed > chr2_nrpe1_input.bed
awk '$1==3' sorted_nrpe1_minusTer_sort_dedup_geneCount.bed > chr3_nrpe1_input.bed
awk '$1==4' sorted_nrpe1_minusTer_sort_dedup_geneCount.bed > chr4_nrpe1_input.bed
awk '$1==5' sorted_nrpe1_minusTer_sort_dedup_geneCount.bed > chr5_nrpe1_input.bed

python identify_nonPol2_reads_opt4.py chr1_nrpe1_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_nrpe1_chr1.bed
python identify_nonPol2_reads_opt4.py chr2_nrpe1_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_nrpe1_chr2.bed
python identify_nonPol2_reads_opt4.py chr3_nrpe1_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_nrpe1_chr3.bed
python identify_nonPol2_reads_opt4.py chr4_nrpe1_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_nrpe1_chr4.bed
python identify_nonPol2_reads_opt4.py chr5_nrpe1_input.bed sorted_input_chr_loci_counts_nonzero.bed scoring_4_nrpe1_chr5.bed

paste chr1_nrpe1_input.bed scoring_4_nrpe1_chr1.bed > chr1_nrpe1_reads_scores.bed
paste chr2_nrpe1_input.bed scoring_4_nrpe1_chr2.bed > chr2_nrpe1_reads_scores.bed
paste chr3_nrpe1_input.bed scoring_4_nrpe1_chr3.bed > chr3_nrpe1_reads_scores.bed
paste chr4_nrpe1_input.bed scoring_4_nrpe1_chr4.bed > chr4_nrpe1_reads_scores.bed
paste chr5_nrpe1_input.bed scoring_4_nrpe1_chr5.bed > chr5_nrpe1_reads_scores.bed

cat chr1_nrpe1_reads_scores.bed chr2_nrpe1_reads_scores.bed chr3_nrpe1_reads_scores.bed chr4_nrpe1_reads_scores.bed chr5_nrpe1_reads_scores.bed > chr#_nrpe1_read_scores.bed


mkdir select_subset
cd select_subset/
awk '$8>-2' ../chr#_read_scores.bed | intersectBed -s -v -a stdin -b ~/cifs-lab/Shriya/Sequences/genes_list.bed > scores_greater_than_minus2.bed

intersectBed -s -c -a ~/cifs-lab/Shriya/Samples/Masa_RIP_terminator/window_read_test/strandSpecific_windows.bed -b scores_greater_than_minus2.bed > window_wise_col0_scored_count.bed

awk '$8>-2' ../chr#_nrpe1_read_scores.bed | intersectBed -s -v -a stdin -b ~/cifs-lab/Shriya/Sequences/genes_list.bed > scores_nrpe1_greater_than_minus2.bed

intersectBed -s -c -a window_wise_col0_scored_count.bed -b scores_nrpe1_greater_than_minus2.bed > window_wise_col0_nrpe1_scored_count.bed
intersectBed -v -a window_wise_col0_nrpe1_scored_count.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_regions_Col0GreaterThan0.bed | intersectBed -s -v -a stdin -b ~/cifs-lab/Shriya/Sequences/genes_list.bed > windows_nongene_nonPVT_counts.bed
intersectBed -v -a window_wise_col0_nrpe1_scored_count.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_regions_Col0GreaterThan0.bed | intersectBed -v -a stdin -b ~/cifs-lab/Shriya/Sequences/genes_list.bed > windows_nongene_nonPVT_counts_nonstranded.bed



#awk '($3-$2)>=20' scores_greater_than_minus2.bed | intersectBed -wa -a stdin -b ~/cifs-lab/Shriya/Samples/Masa_RIP_terminator/methyltrans_PVT/all_merged_PVT_merge50.bed | uniq | awk '{print $3-$2}' > col0_minLt_20_overlap_PVT.bed
#awk '($3-$2)>=20' scores_nrpe1_greater_than_minus2.bed | intersectBed -wa -a stdin -b ~/cifs-lab/Shriya/Samples/Masa_RIP_terminator/methyltrans_PVT/all_merged_PVT_merge50.bed | uniq | awk '{print $3-$2}' > nrpe1_minLt_20_overlap_PVT.bed


intersectBed -wa -a scores_greater_than_minus2.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_regions_Col0GreaterThan0.bed | uniq | awk '{print $3-$2}' > col0_overlap_PVT.bed
intersectBed -wa -a scores_nrpe1_greater_than_minus2.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_regions_Col0GreaterThan0.bed | uniq | awk '{print $3-$2}' > nrpe1_overlap_PVT.bed

