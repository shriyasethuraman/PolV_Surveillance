#!/bin/sh

#cd /ssd/shriya_data/remap_Masa_RIP/self_trim2/
# cd ~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/trimmed_reads/self_trim2/

## The 3' adapters and polyA sequences are removed using cutadapt.
## Removed additional polyA sequences and a minimum length cut-off of 20bps.
## Trimmed reads are then mapped to bowtie2 allowing 1 mismatch 
## Converting sam -> bam -> remove PCR duplicates using UMI information of each read -> bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118708_Col0_R12_umiheader_c1.fq ../self_trim/118708_Col0_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118708_Col0_R12_umiheader_c2.fq 118708_Col0_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118708_Col0_R12_umiheader_c2.fq -S col0_R12_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118709_nrpe1_R12_umiheader_c1.fq ../self_trim/118709_nrpe1_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118709_nrpe1_R12_umiheader_c2.fq 118709_nrpe1_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118709_nrpe1_R12_umiheader_c2.fq -S nrpe1_R12_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118708_Col0_R12_2ndSeq_umiheader_c1.fq ../self_trim/118708_Col0_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118708_Col0_R12_2ndSeq_umiheader_c2.fq 118708_Col0_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118708_Col0_R12_2ndSeq_umiheader_c2.fq -S col0_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118709_nrpe1_R12_2ndSeq_umiheader_c1.fq ../self_trim/118709_nrpe1_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118709_nrpe1_R12_2ndSeq_umiheader_c2.fq 118709_nrpe1_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118709_nrpe1_R12_2ndSeq_umiheader_c2.fq -S nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS col0_R12_uminame_c2_MasaStyleMap.sam > col0_R12_uminame_c2_MasaStyleMap.bam ## Converting mapped sam files to bam format
samtools sort col0_R12_uminame_c2_MasaStyleMap.bam col0_R12_uminame_c2_MasaStyleMap_sort ## sort the bam files
samtools index col0_R12_uminame_c2_MasaStyleMap_sort.bam ## index the bam files
umi_tools dedup -I col0_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_R12_uminame_c2_MasaStyleMap_sort_dedup.bam ## removed the PCR duplicates utilizing UMI information in the read name
bamToBed -i col0_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_R12_uminame_c2_MasaStyleMap_sort_dedup.bed ## conversing the bam files to bed files

samtools view -bS nrpe1_R12_uminame_c2_MasaStyleMap.sam > nrpe1_R12_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_R12_uminame_c2_MasaStyleMap.bam nrpe1_R12_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

samtools view -bS col0_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > col0_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort col0_R12_2ndSeq_uminame_c2_MasaStyleMap.bam col0_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index col0_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col0_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col0_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

samtools view -bS nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap.bam nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed


## rep1 other mutants## 

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118710_spt5l_R12_umiheader_c1.fq ../self_trim/118710_spt5l_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118710_spt5l_R12_umiheader_c2.fq 118710_spt5l_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118710_spt5l_R12_umiheader_c2.fq -S spt5l_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS spt5l_R12_uminame_c2_MasaStyleMap.sam > spt5l_R12_uminame_c2_MasaStyleMap.bam
samtools sort spt5l_R12_uminame_c2_MasaStyleMap.bam spt5l_R12_uminame_c2_MasaStyleMap_sort
samtools index spt5l_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I spt5l_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S spt5l_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i spt5l_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > spt5l_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118710_spt5l_R12_2ndSeq_umiheader_c1.fq ../self_trim/118710_spt5l_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118710_spt5l_R12_2ndSeq_umiheader_c2.fq 118710_spt5l_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118710_spt5l_R12_2ndSeq_umiheader_c2.fq -S spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap.bam spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > spt5l_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118711_cmt3_R12_umiheader_c1.fq ../self_trim/118711_cmt3_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118711_cmt3_R12_umiheader_c2.fq 118711_cmt3_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118711_cmt3_R12_umiheader_c2.fq -S cmt3_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS cmt3_R12_uminame_c2_MasaStyleMap.sam > cmt3_R12_uminame_c2_MasaStyleMap.bam
samtools sort cmt3_R12_uminame_c2_MasaStyleMap.bam cmt3_R12_uminame_c2_MasaStyleMap_sort
samtools index cmt3_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I cmt3_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S cmt3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i cmt3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > cmt3_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 118712_cmt2_R12_umiheader_c1.fq ../self_trim/118712_cmt2_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 118712_cmt2_R12_umiheader_c2.fq 118712_cmt2_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 118712_cmt2_R12_umiheader_c2.fq -S cmt2_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS cmt2_R12_uminame_c2_MasaStyleMap.sam > cmt2_R12_uminame_c2_MasaStyleMap.bam
samtools sort cmt2_R12_uminame_c2_MasaStyleMap.bam cmt2_R12_uminame_c2_MasaStyleMap_sort
samtools index cmt2_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I cmt2_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S cmt2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i cmt2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > cmt2_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 119205_col02_R12_umiheader_c1.fq ../self_trim/119205_col02_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 119205_col02_R12_umiheader_c2.fq 119205_col02_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 119205_col02_R12_umiheader_c2.fq -S col02_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS col02_R12_uminame_c2_MasaStyleMap.sam > col02_R12_uminame_c2_MasaStyleMap.bam
samtools sort col02_R12_uminame_c2_MasaStyleMap.bam col02_R12_uminame_c2_MasaStyleMap_sort
samtools index col02_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col02_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col02_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col02_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > col02_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 119206_ago4_R12_umiheader_c1.fq ../self_trim/119206_ago4_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 119206_ago4_R12_umiheader_c2.fq 119206_ago4_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 119206_ago4_R12_umiheader_c2.fq -S ago4_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS ago4_R12_uminame_c2_MasaStyleMap.sam > ago4_R12_uminame_c2_MasaStyleMap.bam
samtools sort ago4_R12_uminame_c2_MasaStyleMap.bam ago4_R12_uminame_c2_MasaStyleMap_sort
samtools index ago4_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I ago4_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S ago4_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i ago4_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > ago4_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 119207_drm2_R12_umiheader_c1.fq ../self_trim/119207_drm2_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 119207_drm2_R12_umiheader_c2.fq 119207_drm2_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 119207_drm2_R12_umiheader_c2.fq -S drm2_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS drm2_R12_uminame_c2_MasaStyleMap.sam > drm2_R12_uminame_c2_MasaStyleMap.bam
samtools sort drm2_R12_uminame_c2_MasaStyleMap.bam drm2_R12_uminame_c2_MasaStyleMap_sort
samtools index drm2_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I drm2_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S drm2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i drm2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > drm2_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 119208_met13_R12_umiheader_c1.fq ../self_trim/119208_met13_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 119208_met13_R12_umiheader_c2.fq 119208_met13_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 119208_met13_R12_umiheader_c2.fq -S met13_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS met13_R12_uminame_c2_MasaStyleMap.sam > met13_R12_uminame_c2_MasaStyleMap.bam
samtools sort met13_R12_uminame_c2_MasaStyleMap.bam met13_R12_uminame_c2_MasaStyleMap_sort
samtools index met13_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I met13_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S met13_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i met13_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > met13_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


#### rep2 

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123313_Col0_rep3_R12_umiheader_c1.fq ../self_trim/123313_Col0_rep3_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123313_Col0_rep3_R12_umiheader_c2.fq 123313_Col0_rep3_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123313_Col0_rep3_R12_umiheader_c2.fq -S col0_rep3_R12_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123314_nrpe1_rep3_R12_umiheader_c1.fq ../self_trim/123314_nrpe1_rep3_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123314_nrpe1_rep3_R12_umiheader_c2.fq 123314_nrpe1_rep3_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123314_nrpe1_rep3_R12_umiheader_c2.fq -S nrpe1_rep3_R12_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123313_Col0_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123313_Col0_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123313_Col0_rep3_R12_2ndSeq_umiheader_c2.fq 123313_Col0_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123313_Col0_rep3_R12_2ndSeq_umiheader_c2.fq -S col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123314_nrpe1_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123314_nrpe1_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123314_nrpe1_rep3_R12_2ndSeq_umiheader_c2.fq 123314_nrpe1_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123314_nrpe1_rep3_R12_2ndSeq_umiheader_c2.fq -S nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam


samtools view -bS col0_rep3_R12_uminame_c2_MasaStyleMap.sam > col0_rep3_R12_uminame_c2_MasaStyleMap.bam
samtools sort col0_rep3_R12_uminame_c2_MasaStyleMap.bam col0_rep3_R12_uminame_c2_MasaStyleMap_sort
samtools index col0_rep3_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col0_rep3_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col0_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

samtools view -bS nrpe1_rep3_R12_uminame_c2_MasaStyleMap.sam > nrpe1_rep3_R12_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_rep3_R12_uminame_c2_MasaStyleMap.bam nrpe1_rep3_R12_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_rep3_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_rep3_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

samtools view -bS col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

samtools view -bS nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123315_spt5l_rep3_R12_umiheader_c1.fq ../self_trim/123315_spt5l_rep3_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123315_spt5l_rep3_R12_umiheader_c2.fq 123315_spt5l_rep3_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123315_spt5l_rep3_R12_umiheader_c2.fq -S spt5l_rep3_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS spt5l_rep3_R12_uminame_c2_MasaStyleMap.sam > spt5l_rep3_R12_uminame_c2_MasaStyleMap.bam
samtools sort spt5l_rep3_R12_uminame_c2_MasaStyleMap.bam spt5l_rep3_R12_uminame_c2_MasaStyleMap_sort
samtools index spt5l_rep3_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I spt5l_rep3_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S spt5l_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i spt5l_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > spt5l_rep3_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123315_spt5l_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123315_spt5l_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123315_spt5l_rep3_R12_2ndSeq_umiheader_c2.fq 123315_spt5l_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123315_spt5l_rep3_R12_2ndSeq_umiheader_c2.fq -S spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > spt5l_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123325_ago4_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123325_ago4_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123325_ago4_rep3_R12_2ndSeq_umiheader_c2.fq 123325_ago4_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123325_ago4_rep3_R12_2ndSeq_umiheader_c2.fq -S ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > ago4_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123326_drm2_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123326_drm2_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123326_drm2_rep3_R12_2ndSeq_umiheader_c2.fq 123326_drm2_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123326_drm2_rep3_R12_2ndSeq_umiheader_c2.fq -S drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > drm2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123328_cmt3_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123328_cmt3_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123328_cmt3_rep3_R12_2ndSeq_umiheader_c2.fq 123328_cmt3_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123328_cmt3_rep3_R12_2ndSeq_umiheader_c2.fq -S cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > cmt3_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123329_cmt2_rep3_R12_2ndSeq_umiheader_c1.fq ../self_trim/123329_cmt2_rep3_R12_2ndSeq_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123329_cmt2_rep3_R12_2ndSeq_umiheader_c2.fq 123329_cmt2_rep3_R12_2ndSeq_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123329_cmt2_rep3_R12_2ndSeq_umiheader_c2.fq -S cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam

samtools view -bS cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.sam > cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam
samtools sort cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap.bam cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort
samtools index cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bam > cmt2_rep3_R12_2ndSeq_uminame_c2_MasaStyleMap_sort_dedup.bed



##### rep3 #############

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123330_Col0_rep5_R12_umiheader_c1.fq ../self_trim/123330_Col0_rep5_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123330_Col0_rep5_R12_umiheader_c2.fq 123330_Col0_rep5_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123330_Col0_rep5_R12_umiheader_c2.fq -S col0_rep5_R12_uminame_c2_MasaStyleMap.sam

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123331_nrpe1_rep5_R12_umiheader_c1.fq ../self_trim/123331_nrpe1_rep5_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123331_nrpe1_rep5_R12_umiheader_c2.fq 123331_nrpe1_rep5_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123331_nrpe1_rep5_R12_umiheader_c2.fq -S nrpe1_rep5_R12_uminame_c2_MasaStyleMap.sam


samtools view -bS col0_rep5_R12_uminame_c2_MasaStyleMap.sam > col0_rep5_R12_uminame_c2_MasaStyleMap.bam
samtools sort col0_rep5_R12_uminame_c2_MasaStyleMap.bam col0_rep5_R12_uminame_c2_MasaStyleMap_sort
samtools index col0_rep5_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col0_rep5_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_rep5_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col0_rep5_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_rep5_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

samtools view -bS nrpe1_rep5_R12_uminame_c2_MasaStyleMap.sam > nrpe1_rep5_R12_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_rep5_R12_uminame_c2_MasaStyleMap.bam nrpe1_rep5_R12_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_rep5_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_rep5_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_rep5_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_rep5_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_rep5_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123332_drd1_R12_umiheader_c1.fq ../self_trim/123332_drd1_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123332_drd1_R12_umiheader_c2.fq 123332_drd1_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123332_drd1_R12_umiheader_c2.fq -S drd1_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS drd1_R12_uminame_c2_MasaStyleMap.sam > drd1_R12_uminame_c2_MasaStyleMap.bam
samtools sort drd1_R12_uminame_c2_MasaStyleMap.bam drd1_R12_uminame_c2_MasaStyleMap_sort
samtools index drd1_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I drd1_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S drd1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i drd1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > drd1_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123333_dms3_R12_umiheader_c1.fq ../self_trim/123333_dms3_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123333_dms3_R12_umiheader_c2.fq 123333_dms3_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123333_dms3_R12_umiheader_c2.fq -S dms3_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS dms3_R12_uminame_c2_MasaStyleMap.sam > dms3_R12_uminame_c2_MasaStyleMap.bam
samtools sort dms3_R12_uminame_c2_MasaStyleMap.bam dms3_R12_uminame_c2_MasaStyleMap_sort
samtools index dms3_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I dms3_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S dms3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i dms3_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > dms3_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123334_suvh2_R12_umiheader_c1.fq ../self_trim/123334_suvh2_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123334_suvh2_R12_umiheader_c2.fq 123334_suvh2_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123334_suvh2_R12_umiheader_c2.fq -S suvh2_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS suvh2_R12_uminame_c2_MasaStyleMap.sam > suvh2_R12_uminame_c2_MasaStyleMap.bam
samtools sort suvh2_R12_uminame_c2_MasaStyleMap.bam suvh2_R12_uminame_c2_MasaStyleMap_sort
samtools index suvh2_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I suvh2_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S suvh2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i suvh2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > suvh2_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123335_suvh9_R12_umiheader_c1.fq ../self_trim/123335_suvh9_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123335_suvh9_R12_umiheader_c2.fq 123335_suvh9_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123335_suvh9_R12_umiheader_c2.fq -S suvh9_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS suvh9_R12_uminame_c2_MasaStyleMap.sam > suvh9_R12_uminame_c2_MasaStyleMap.bam
samtools sort suvh9_R12_uminame_c2_MasaStyleMap.bam suvh9_R12_uminame_c2_MasaStyleMap_sort
samtools index suvh9_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I suvh9_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S suvh9_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i suvh9_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > suvh9_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 123336_suvh29_R12_umiheader_c1.fq ../self_trim/123336_suvh29_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 123336_suvh29_R12_umiheader_c2.fq 123336_suvh29_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 123336_suvh29_R12_umiheader_c2.fq -S suvh29_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS suvh29_R12_uminame_c2_MasaStyleMap.sam > suvh29_R12_uminame_c2_MasaStyleMap.bam
samtools sort suvh29_R12_uminame_c2_MasaStyleMap.bam suvh29_R12_uminame_c2_MasaStyleMap_sort
samtools index suvh29_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I suvh29_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S suvh29_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i suvh29_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > suvh29_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


### point mutants #####

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127881_Col0_rep6_R12_umiheader_c1.fq ../self_trim/127881_Col0_rep6_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127881_Col0_rep6_R12_umiheader_c2.fq 127881_Col0_rep6_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127881_Col0_rep6_R12_umiheader_c2.fq -S col0_rep6_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS col0_rep6_R12_uminame_c2_MasaStyleMap.sam > col0_rep6_R12_uminame_c2_MasaStyleMap.bam
samtools sort col0_rep6_R12_uminame_c2_MasaStyleMap.bam col0_rep6_R12_uminame_c2_MasaStyleMap_sort
samtools index col0_rep6_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col0_rep6_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_rep6_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col0_rep6_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_rep6_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127882_nrpe1_rep6_R12_umiheader_c1.fq ../self_trim/127882_nrpe1_rep6_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127882_nrpe1_rep6_R12_umiheader_c2.fq 127882_nrpe1_rep6_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127882_nrpe1_rep6_R12_umiheader_c2.fq -S nrpe1_rep6_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS nrpe1_rep6_R12_uminame_c2_MasaStyleMap.sam > nrpe1_rep6_R12_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_rep6_R12_uminame_c2_MasaStyleMap.bam nrpe1_rep6_R12_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_rep6_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_rep6_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_rep6_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_rep6_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_rep6_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127885_Col0_rep7_R12_umiheader_c1.fq ../self_trim/127885_Col0_rep7_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127885_Col0_rep7_R12_umiheader_c2.fq 127885_Col0_rep7_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127885_Col0_rep7_R12_umiheader_c2.fq -S col0_rep7_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS col0_rep7_R12_uminame_c2_MasaStyleMap.sam > col0_rep7_R12_uminame_c2_MasaStyleMap.bam
samtools sort col0_rep7_R12_uminame_c2_MasaStyleMap.bam col0_rep7_R12_uminame_c2_MasaStyleMap_sort
samtools index col0_rep7_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I col0_rep7_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S col0_rep7_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i col0_rep7_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > col0_rep7_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127886_nrpe1_rep7_R12_umiheader_c1.fq ../self_trim/127886_nrpe1_rep7_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127886_nrpe1_rep7_R12_umiheader_c2.fq 127886_nrpe1_rep7_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127886_nrpe1_rep7_R12_umiheader_c2.fq -S nrpe1_rep7_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS nrpe1_rep7_R12_uminame_c2_MasaStyleMap.sam > nrpe1_rep7_R12_uminame_c2_MasaStyleMap.bam
samtools sort nrpe1_rep7_R12_uminame_c2_MasaStyleMap.bam nrpe1_rep7_R12_uminame_c2_MasaStyleMap_sort
samtools index nrpe1_rep7_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I nrpe1_rep7_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S nrpe1_rep7_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i nrpe1_rep7_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > nrpe1_rep7_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127883_dms5_rep1_R12_umiheader_c1.fq ../self_trim/127883_dms5_rep1_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127883_dms5_rep1_R12_umiheader_c2.fq 127883_dms5_rep1_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127883_dms5_rep1_R12_umiheader_c2.fq -S dms5_rep1_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS dms5_rep1_R12_uminame_c2_MasaStyleMap.sam > dms5_rep1_R12_uminame_c2_MasaStyleMap.bam
samtools sort dms5_rep1_R12_uminame_c2_MasaStyleMap.bam dms5_rep1_R12_uminame_c2_MasaStyleMap_sort
samtools index dms5_rep1_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I dms5_rep1_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S dms5_rep1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i dms5_rep1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > dms5_rep1_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127884_drd3_rep1_R12_umiheader_c1.fq ../self_trim/127884_drd3_rep1_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127884_drd3_rep1_R12_umiheader_c2.fq 127884_drd3_rep1_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127884_drd3_rep1_R12_umiheader_c2.fq -S drd3_rep1_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS drd3_rep1_R12_uminame_c2_MasaStyleMap.sam > drd3_rep1_R12_uminame_c2_MasaStyleMap.bam
samtools sort drd3_rep1_R12_uminame_c2_MasaStyleMap.bam drd3_rep1_R12_uminame_c2_MasaStyleMap_sort
samtools index drd3_rep1_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I drd3_rep1_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S drd3_rep1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i drd3_rep1_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > drd3_rep1_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127887_dms5_rep2_R12_umiheader_c1.fq ../self_trim/127887_dms5_rep2_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127887_dms5_rep2_R12_umiheader_c2.fq 127887_dms5_rep2_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127887_dms5_rep2_R12_umiheader_c2.fq -S dms5_rep2_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS dms5_rep2_R12_uminame_c2_MasaStyleMap.sam > dms5_rep2_R12_uminame_c2_MasaStyleMap.bam
samtools sort dms5_rep2_R12_uminame_c2_MasaStyleMap.bam dms5_rep2_R12_uminame_c2_MasaStyleMap_sort
samtools index dms5_rep2_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I dms5_rep2_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S dms5_rep2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i dms5_rep2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > dms5_rep2_R12_uminame_c2_MasaStyleMap_sort_dedup.bed

cutadapt -e 0.05 -a AAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNNNCTGTCTCTTATACACATCTCCGAGCCCACGAG -o 127888_drd3_rep2_R12_umiheader_c1.fq ../self_trim/127888_drd3_rep2_R12_umiheader.fastq
cutadapt -m 20 -e 0.05 -a AAAAAAAAAAAAAAAAAAAA -o 127888_drd3_rep2_R12_umiheader_c2.fq 127888_drd3_rep2_R12_umiheader_c1.fq
bowtie2 -q -N 1 -p 8 -x ~/cifs-lab/Shriya/Reference/Whole_genome/bowtie2_index/bowtie2_indx -U 127888_drd3_rep2_R12_umiheader_c2.fq -S drd3_rep2_R12_uminame_c2_MasaStyleMap.sam

samtools view -bS drd3_rep2_R12_uminame_c2_MasaStyleMap.sam > drd3_rep2_R12_uminame_c2_MasaStyleMap.bam
samtools sort drd3_rep2_R12_uminame_c2_MasaStyleMap.bam drd3_rep2_R12_uminame_c2_MasaStyleMap_sort
samtools index drd3_rep2_R12_uminame_c2_MasaStyleMap_sort.bam
umi_tools dedup -I drd3_rep2_R12_uminame_c2_MasaStyleMap_sort.bam --output-stats=deduplicated -S drd3_rep2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam
bamToBed -i drd3_rep2_R12_uminame_c2_MasaStyleMap_sort_dedup.bam > drd3_rep2_R12_uminame_c2_MasaStyleMap_sort_dedup.bed


