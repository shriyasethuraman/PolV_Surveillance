#!/bin/sh

#cd /ssd/shriya_data/remap_Masa_RIP/self_trim/
#cd ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/trimmed_reads/self_trim/

## UMI from each read extracted from the read2 of the sequenced fastq files using umi_tools_extract command of UMI_tools
## Append the UMI info extracted to the read names of read1. This can then be used to remove the duplicates

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118708_Col0_R2.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_R2.log --stdout=118708_Col0_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118708_Col0_R1.fastq 118708_Col0_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118708_Col0_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118709_nrpe1_R2.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_R2.log --stdout=118709_nrpe1_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118709_nrpe1_R1.fastq 118709_nrpe1_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118709_nrpe1_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/118708_Col0_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_R2_2ndSeq.log --stdout=118708_Col0_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/118708_Col0_R1_2ndSeq.fastq 118708_Col0_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118708_Col0_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/118709_nrpe1_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_R2_2ndSeq.log --stdout=118709_nrpe1_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/118709_nrpe1_R1_2ndSeq.fastq 118709_nrpe1_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118709_nrpe1_R12_2ndSeq_umiheader.fastq


umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118710_spt5l_R2.fastq --bc-pattern=NNNNNNNN --log=umi_spt5l_minusTer_R2.log --stdout=118710_spt5l_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118710_spt5l_R1.fastq 118710_spt5l_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118710_spt5l_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/118710_spt5l_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_spt5l_minusTer_R2_2ndSeq.log --stdout=118710_spt5l_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/118710_spt5l_R1_2ndSeq.fastq 118710_spt5l_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118710_spt5l_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118711_cmt3_R2.fastq --bc-pattern=NNNNNNNN --log=umi_cmt3_minusTer_R2.log --stdout=118711_cmt3_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118711_cmt3_R1.fastq 118711_cmt3_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118711_cmt3_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118712_cmt2_R2.fastq --bc-pattern=NNNNNNNN --log=umi_cmt2_minusTer_R2.log --stdout=118712_cmt2_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/118712_cmt2_R1.fastq 118712_cmt2_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 118712_cmt2_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119205_col02_R2.fastq --bc-pattern=NNNNNNNN --log=umi_col02_minusTer_R2.log --stdout=119205_col02_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119205_col02_R1.fastq 119205_col02_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 119205_col02_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119206_ago4_R2.fastq --bc-pattern=NNNNNNNN --log=umi_ago4_minusTer_R2.log --stdout=119206_ago4_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119206_ago4_R1.fastq 119206_ago4_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 119206_ago4_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119207_drm2_R2.fastq --bc-pattern=NNNNNNNN --log=umi_drm2_minusTer_R2.log --stdout=119207_drm2_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119207_drm2_R1.fastq 119207_drm2_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 119207_drm2_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119208_met13_R2.fastq --bc-pattern=NNNNNNNN --log=umi_met13_minusTer_R2.log --stdout=119208_met13_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2553/119208_met13_R1.fastq 119208_met13_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 119208_met13_R12_umiheader.fastq


#umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/123313_Col0_rep3_R2.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_rep3_R2.log --stdout=123313_Col0_rep3_R2_umiextract.fastq
#paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/123313_Col0_rep3_R1.fastq 123313_Col0_rep3_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123313_Col0_rep3_R12_umiheader.fastq

#umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/123314_nrpe1_rep3_R2.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_rep3_R2.log --stdout=123314_nrpe1_rep3_R2_umiextract.fastq
#paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/123314_nrpe1_rep3_R1.fastq 123314_nrpe1_rep3_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123314_nrpe1_rep3_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/123315_spt5l_rep3_R2.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_rep3_R2.log --stdout=123315_spt5l_rep3_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_rep1_2652/123315_spt5l_rep3_R1.fastq 123315_spt5l_rep3_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123315_spt5l_rep3_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123313_Col0_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_rep3_R2_2ndSeq.log --stdout=123313_Col0_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123313_Col0_rep3_R1_2ndSeq.fastq 123313_Col0_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123313_Col0_rep3_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123314_nrpe1_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_rep3_R2_2ndSeq.log --stdout=123314_nrpe1_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123314_nrpe1_rep3_R1_2ndSeq.fastq 123314_nrpe1_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123314_nrpe1_rep3_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123315_spt5l_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_spt5l_minusTer_rep3_R2_2ndSeq.log --stdout=123315_spt5l_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123315_spt5l_rep3_R1_2ndSeq.fastq 123315_spt5l_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123315_spt5l_rep3_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123325_ago4_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_ago4_minusTer_rep3_R2_2ndSeq.log --stdout=123325_ago4_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123325_ago4_rep3_R1_2ndSeq.fastq 123325_ago4_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123325_ago4_rep3_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123326_drm2_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_drm2_minusTer_rep3_R2_2ndSeq.log --stdout=123326_drm2_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123326_drm2_rep3_R1_2ndSeq.fastq 123326_drm2_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123326_drm2_rep3_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123328_cmt3_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_cmt3_minusTer_rep3_R2_2ndSeq.log --stdout=123328_cmt3_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123328_cmt3_rep3_R1_2ndSeq.fastq 123328_cmt3_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123328_cmt3_rep3_R12_2ndSeq_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123329_cmt2_rep3_R2_2ndSeq.fastq --bc-pattern=NNNNNNNN --log=umi_cmt2_minusTer_rep3_R2_2ndSeq.log --stdout=123329_cmt2_rep3_R2_2ndSeq_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123329_cmt2_rep3_R1_2ndSeq.fastq 123329_cmt2_rep3_R2_2ndSeq_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123329_cmt2_rep3_R12_2ndSeq_umiheader.fastq



umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123330_Col0_rep5_R2.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_rep5_R2.log --stdout=123330_Col0_rep5_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123330_Col0_rep5_R1.fastq 123330_Col0_rep5_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123330_Col0_rep5_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123331_nrpe1_rep5_R2.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_rep5_R2.log --stdout=123331_nrpe1_rep5_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123331_nrpe1_rep5_R1.fastq 123331_nrpe1_rep5_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123331_nrpe1_rep5_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123332_drd1_R2.fastq --bc-pattern=NNNNNNNN --log=umi_drd1_minusTer_R2.log --stdout=123332_drd1_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123332_drd1_R1.fastq 123332_drd1_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123332_drd1_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123333_dms3_R2.fastq --bc-pattern=NNNNNNNN --log=umi_dms3_minusTer_R2.log --stdout=123333_dms3_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123333_dms3_R1.fastq 123333_dms3_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123333_dms3_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123334_suvh2_R2.fastq --bc-pattern=NNNNNNNN --log=umi_suvh2_minusTer_R2.log --stdout=123334_suvh2_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123334_suvh2_R1.fastq 123334_suvh2_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123334_suvh2_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123335_suvh9_R2.fastq --bc-pattern=NNNNNNNN --log=umi_suvh9_minusTer_R2.log --stdout=123335_suvh9_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123335_suvh9_R1.fastq 123335_suvh9_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123335_suvh9_R12_umiheader.fastq

umi_tools extract --stdin=/ssd/shriya_data/remap_Masa_RIP/reads_2666/123336_suvh29_R2.fastq --bc-pattern=NNNNNNNN --log=umi_suvh29_minusTer_R2.log --stdout=123336_suvh29_R2_umiextract.fastq
paste /ssd/shriya_data/remap_Masa_RIP/reads_2666/123336_suvh29_R1.fastq 123336_suvh29_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 123336_suvh29_R12_umiheader.fastq


umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127881/127881_CGTCTGCG_S19_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_rep6_R2.log --stdout=127881_Col0_rep6_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127881/127881_CGTCTGCG_S19_R1_001.fastq 127881_Col0_rep6_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127881_Col0_rep6_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127882/127882_TACTCATA_S20_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_rep6_R2.log --stdout=127882_nrpe1_rep6_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127882/127882_TACTCATA_S20_R1_001.fastq 127882_nrpe1_rep6_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127882_nrpe1_rep6_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127883/127883_ACGCACCT_S21_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_dms5_rep1_minusTer_R2.log --stdout=127883_dms5_rep1_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127883/127883_ACGCACCT_S21_R1_001.fastq 127883_dms5_rep1_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127883_dms5_rep1_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127884/127884_GTATGTTC_S22_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_drd3_rep1_minusTer_R2.log --stdout=127884_drd3_rep1_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127884/127884_GTATGTTC_S22_R1_001.fastq 127884_drd3_rep1_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127884_drd3_rep1_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127885/127885_CGCTATGT_S23_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_col0_minusTer_rep7_R2.log --stdout=127885_Col0_rep7_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127885/127885_CGCTATGT_S23_R1_001.fastq 127885_Col0_rep7_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127885_Col0_rep7_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127886/127886_CCAACAGA_S24_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_nrpe1_minusTer_rep7_R2.log --stdout=127886_nrpe1_rep7_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127886/127886_CCAACAGA_S24_R1_001.fastq 127886_nrpe1_rep7_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127886_nrpe1_rep7_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127887/127887_TTCAGGTC_S25_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_dms5_rep2_minusTer_R2.log --stdout=127887_dms5_rep2_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127887/127887_TTCAGGTC_S25_R1_001.fastq 127887_dms5_rep2_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127887_dms5_rep2_R12_umiheader.fastq

umi_tools extract --stdin=/home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127888/127888_TCTGTTGG_S26_R2_001.fastq --bc-pattern=NNNNNNNN --log=umi_drd3_rep2_minusTer_R2.log --stdout=127888_drd3_rep2_R2_umiextract.fastq
paste /home/shriyas/cifs-lab/sequencing_dataset_archive/RIP_Pol_V_Terminator/Run_2866/wierzbicki/Sample_127888/127888_TCTGTTGG_S26_R1_001.fastq 127888_drd3_rep2_R2_umiextract.fastq | awk '{OFS="\t"}{if($1 ~ /@/) print $3" "$4; else print $1}' > 127888_drd3_rep2_R12_umiheader.fastq



