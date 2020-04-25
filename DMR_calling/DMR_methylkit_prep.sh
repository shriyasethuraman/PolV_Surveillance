#!/bin/sh

cd ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/

awk '{OFS="\t"}{if(($3=="+") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/meth_extr/SRR3316941_trimmed_bismark_bt2.CX_report.txt > Col0_labMap_methylKit.txt 
awk '{OFS="\t"}{if(($3=="+") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/meth_extr/SRR3316949_trimmed_bismark_bt2.CX_report.txt > nrpe1_labMap_methylKit.txt 
awk '{OFS="\t"}{if(($3=="+") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/meth_extr/SRR3316952_trimmed_bismark_bt2.CX_report.txt > ddm1_labMap_methylKit.txt 
awk '{OFS="\t"}{if(($3=="+") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CHH") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/meth_extr/SRR3316959_trimmed_bismark_bt2.CX_report.txt > ddm1nrpe1_labMap_methylKit.txt 


sed 's/"//g' mydiff_reg_ddm1_CHH_col0_d10_q01.txt | sed 's/chr//g' > mydiff_reg_ddm1_CHH_col0_d10_q01.bed

intersectBed -wa -a mydiff_reg_ddm1_CHH_col0_d10_q01.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/PVS_loci_ColGreaterThan0.bed | uniq | wc -l
3302
intersectBed -wa -a mydiff_reg_ddm1_CHH_col0_d10_q01.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/prospectivePVS_100bp_80slidebins_nostr_RIPrep1_rep2_rep3_mergedRIP_scoredRIP_plusTer_labCNme_ALLjacobsen_methylation_h3-h3k9m2-hafiz_histone_sRNA_UMI.bed | uniq | wc -l
3353
intersectBed -wa -a mydiff_reg_ddm1_CHH_col0_d10_q01.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_loci_ColGreaterThan0.bed | uniq | wc -l
23877

wc -l mydiff_reg_ddm1_CHH_col0_d10_q01.bed 
24354 mydiff_reg_ddm1_CHH_col0_d10_q01.bed

intersectBed -wa -a mydiff_reg_ddm1_CHH_col0_d10_q01.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/PVS_loci_ColGreaterThan0.bed | uniq | intersectBed -v -a stdin -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_loci_ColGreaterThan0.bed | wc -l
398
intersectBed -wa -b mydiff_reg_ddm1_CHH_col0_d10_q01.bed -a ~/cifs-lab/Shriya/PolV_surveillance_regions/PVS_loci_ColGreaterThan0.bed | uniq | intersectBed -v -a stdin -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_loci_ColGreaterThan0.bed | wc -l
132


intersectBed -wa -f 1 -r -a mydiff_reg_ddm1_CHH_col0_d10_q01.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/PVS_loci_ColGreaterThan0.bed | uniq | intersectBed -v -a stdin -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_loci_ColGreaterThan0.bed | wc -l
276

intersectBed -wa -f 1 -r -a mydiff_reg_ddm1_CHH_col0_d10_q01.bed -b ~/cifs-lab/Shriya/PolV_surveillance_regions/PVS_loci_ColGreaterThan0.bed | uniq | intersectBed -v -a stdin -b ~/cifs-lab/Shriya/PolV_surveillance_regions/RdDM_loci_ColGreaterThan0.bed | awk '{OFS="\t"}{print "chr"$1,$2,$3}' > just_PVS_DMR_methylkit.bed



awk '{OFS="\t"}{if(($3=="+") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/Col0/meth_extr/SRR3316941_trimmed_bismark_bt2.CX_report.txt > Col0_labMap_CG_methylKit.txt 
awk '{OFS="\t"}{if(($3=="+") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/nrpe1/meth_extr/SRR3316949_trimmed_bismark_bt2.CX_report.txt > nrpe1_labMap_CG_methylKit.txt 
awk '{OFS="\t"}{if(($3=="+") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1/meth_extr/SRR3316952_trimmed_bismark_bt2.CX_report.txt > ddm1_labMap_CG_methylKit.txt 
awk '{OFS="\t"}{if(($3=="+") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"F",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5); if(($3=="-") && ($6=="CG") && (($4+$5)>4)) print "chr"$1"."$2,"chr"$1,$2,"R",$4+$5,($4*100)/($4+$5),($5*100)/($4+$5)}' ~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/ddm1_nrpe1/meth_extr/SRR3316959_trimmed_bismark_bt2.CX_report.txt > ddm1nrpe1_labMap_CG_methylKit.txt 

