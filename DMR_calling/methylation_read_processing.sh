#!/bin/sh

#cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/
#fastqc *fastq

#cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/
#fastqc *fastq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/
trim_galore --length 100 --fastqc -o trimmed_reads/ SRR3316941.fastq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/
trim_galore --length 100 --fastqc -o trimmed_reads/ SRR3316949.fastq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/
trim_galore --length 100 --fastqc -o trimmed_reads/ SRR3316952.fastq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/
trim_galore --length 100 --fastqc -o trimmed_reads/ SRR3316959.fastq



cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/
bismark -q --bowtie2 --non_directional -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/se_mapped_reads/ --un --multicore 2 -p 3 -N 1 /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ trimmed_reads/*.fq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/
bismark -q --bowtie2 --non_directional -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/se_mapped_reads/--un --multicore 2 -p 3 -N 1 /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ trimmed_reads/*.fq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/
bismark -q --bowtie2 --non_directional -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/se_mapped_reads/ --un --multicore 2 -p 3 -N 1 /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ trimmed_reads/*.fq

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/
bismark -q --bowtie2 --non_directional -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/se_mapped_reads/ --un --multicore 2 -p 3 -N 1 /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ trimmed_reads/*.fq



cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/
bismark_methylation_extractor -s -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/meth_extr/ --report --multicore 3 --bedGraph --counts --CX_context --cytosine_report --CX --genome_folder /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/se_mapped_reads/*.bam 

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/
bismark_methylation_extractor -s -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/meth_extr/ --report --multicore 3 --bedGraph --counts --CX_context --cytosine_report --CX --genome_folder /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/se_mapped_reads/*.bam 

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/
bismark_methylation_extractor -s -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/meth_extr/ --report --multicore 3 --bedGraph --counts --CX_context --cytosine_report --CX --genome_folder /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/se_mapped_reads/*.bam 

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/
bismark_methylation_extractor -s -o /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/meth_extr/ --report --multicore 3 --bedGraph --counts --CX_context --cytosine_report --CX --genome_folder /home/shriyas/cifs-lab/Shriya/Reference/Whole_genome/ /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/se_mapped_reads/*.bam 

### processing

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/meth_extr/

awk '{OFS="\t"} $6=="CHH" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316941_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > Col0_CHHmeth.bed
awk '{OFS="\t"} $6=="CHG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316941_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > Col0_CHGmeth.bed
awk '{OFS="\t"} $6=="CG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316941_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > Col0_CGmeth.bed
awk '{OFS="\t"}{if(($4+$5)>4) print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316941_trimmed_bismark_bt2.CX_report.txt > Col0_allCs.bed

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/meth_extr/

awk '{OFS="\t"} $6=="CHH" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316949_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > nrpe1_CHHmeth.bed
awk '{OFS="\t"} $6=="CHG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316949_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > nrpe1_CHGmeth.bed
awk '{OFS="\t"} $6=="CG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316949_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > nrpe1_CGmeth.bed
awk '{OFS="\t"}{if(($4+$5)>4) print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316949_trimmed_bismark_bt2.CX_report.txt > nrpe1_allCs.bed

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/meth_extr/

awk '{OFS="\t"} $6=="CHH" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316952_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > ddm1_CHHmeth.bed
awk '{OFS="\t"} $6=="CHG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316952_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > ddm1_CHGmeth.bed
awk '{OFS="\t"} $6=="CG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316952_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > ddm1_CGmeth.bed
awk '{OFS="\t"}{if(($4+$5)>4) print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316952_trimmed_bismark_bt2.CX_report.txt > ddm1_allCs.bed

cd /home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/meth_extr/

awk '{OFS="\t"} $6=="CHH" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316959_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > ddm1nrpe1_CHHmeth.bed
awk '{OFS="\t"} $6=="CHG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316959_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > ddm1nrpe1_CHGmeth.bed
awk '{OFS="\t"} $6=="CG" && $4+$5>4 {print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316959_trimmed_bismark_bt2.CX_report.txt | awk '$1==1 || $1==2 || $1==3 || $1==4 || $1==5' > ddm1nrpe1_CGmeth.bed
awk '{OFS="\t"}{if(($4+$5)>4) print $1,$2-1,$2,$4,$4+$5,$4/($4+$5)}' SRR3316959_trimmed_bismark_bt2.CX_report.txt > ddm1nrpe1_allCs.bed


