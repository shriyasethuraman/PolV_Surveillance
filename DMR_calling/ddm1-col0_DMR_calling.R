source('~/cifs-lab/RIP_manuscript/Revised_Figures/Fig7/Fig_7A/hp-profile-heatmap_modGB.r')
source('~/cifs-lab/Shriya/Codes/window_generator_mod.R')
library('methylKit')
library('gplots')
setwd("/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/")
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/ddm1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/col0_labMap_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("ddm1","col0"),
               assembly="tair10",
               treatment=c(0,1),
               context="CHH")
tiles=tileMethylCounts(myobj,win.size=200,step.size=150,cov.bases = 10)

meth_reg=unite(tiles, destrand=FALSE)
myDiff_reg=calculateDiffMeth(meth_reg)
myDiff10p.hypo_reg=getMethylDiff(myDiff_reg,difference=10,qvalue=0.01,type="hypo")
Qn = myDiff_reg$qvalue<0.01 & myDiff_reg$meth.diff < -10 #alternate way of doing the getMethylDiff function but with a usable index. difference in methylation is more than 10, but it is stated in the negative compared to the control i.e. Col-0.

perc.meth_reg=percMethylation(meth_reg)
perc.meth_reg_qn = perc.meth_reg[Qn,]

hmcol1 = colorRampPalette(c("red","white","blue"))(299)
hmcol2 = colorRampPalette(c("red", "yellow", "green"))(n = 299)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)
col_breaks = c(seq(quantile(perc.meth_reg_qn, 0.1),quantile(perc.meth_reg_qn,0.9), length = 300))

b <- perc.meth_reg_qn[seq(1, nrow(perc.meth_reg_qn), 2),] #### plot only every 20th DMR
pdf("ddm1-col0_CHH_10percDiff_heatmap_qn_every2nd_DMR.pdf", width=8, height=5)
heatmap.2(b,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(10,5), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",
          main = "CHH-DMR (3572): ddm1-col0")
dev.off()


b <- perc.meth_reg_qn[seq(1, nrow(perc.meth_reg_qn), 1),] #### plot only every 20th DMR
pdf("ddm1-col0_CHH_10percDiff_heatmap_qn_every_DMR.pdf", width=8, height=5)
heatmap.2(b,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(10,5), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",
          main = "CHH-DMR (3572): ddm1-col0")
dev.off()

write.table(myDiff10p.hypo_reg, file="mydiff_reg_ddm1-col0_CHH_col0_d10_q01.txt", row.names = FALSE,col.names = F, sep="\t")


myDiff10p.hypo_reg=getMethylDiff(myDiff_reg,difference=20,qvalue=0.01,type="hypo")
Qn = myDiff_reg$qvalue<0.01 & myDiff_reg$meth.diff < -20 #alternate way of doing the getMethylDiff function but with a usable index. difference in methylation is more than 10, but it is stated in the negative compared to the control i.e. Col-0.

perc.meth_reg=percMethylation(meth_reg)
perc.meth_reg_qn = perc.meth_reg[Qn,]

hmcol1 = colorRampPalette(c("red","white","blue"))(299)
hmcol2 = colorRampPalette(c("red", "yellow", "green"))(n = 299)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)
col_breaks = c(seq(quantile(perc.meth_reg_qn, 0.1),quantile(perc.meth_reg_qn,0.9), length = 300))



b <- perc.meth_reg_qn[seq(1, nrow(perc.meth_reg_qn), 1),] #### plot only every 20th DMR
pdf("ddm1-col0_CHH_20percDiff_heatmap_qn_every_DMR.pdf", width=8, height=5)
heatmap.2(b,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(10,5), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",
          main = "CHH-DMR (379): ddm1-col0")
dev.off()

write.table(myDiff10p.hypo_reg, file="mydiff_reg_ddm1-col0_CHH_col0_d20_q01.txt", row.names = FALSE,col.names = F, sep="\t")





bed.obj_PVT = readBed(file = "~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d10_q01_methylkit.bed")

file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/Col0_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/nrpe1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CHH")

dmr.drm2_PVT = regionCounts(myobj,bed.obj_PVT)

meth_PVT=unite(dmr.drm2_PVT, destrand=FALSE)
perc.meth_PVT=percMethylation(meth_PVT)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)

# color break settings for heatmaps
# the length value must be equivalent to the color palette n value
# in this case being 299, the lengths are split 100 each
# the seq (start,stop) values determine the exact breaking point of the different colours
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("ddm1-col0_heatmap_10percDiff_allGenotypes.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "ddm1-col0 CHH-DMRs (#3530)")
dev.off()

pdf("boxplots at ddm1-col0 CHH DMRs percMeth.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at ddm1-col0 CHH DMRs overlapping PVSloci (#3530)", outline=F, ylab="Average CHHme")
dev.off()

perc.meth_PVT_2 <- data.frame(perc.meth_PVT)
temp1 <- perc.meth_PVT_2[which(perc.meth_PVT_2$ddm1>(perc.meth_PVT_2$ddm1_nrpe1*2)),]
temp2 <- temp1[which(temp1$col0<10),]
onlyPVS_noCol <- temp2[which(temp2$nrpe1<10),]

write.table(onlyPVS_noCol, file="mydiff_reg_ddm1-col0_CHH_ddm1MoreThan2DN_col0_and_nrpe1_lessThan10_d10_q01.txt", row.names = FALSE,col.names = F, sep="\t")





bed.obj_PVT = readBed(file = "~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d10_q01_methylkit.bed")

file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/Col0_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/nrpe1_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_CG_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CG")

dmr.drm2_PVT = regionCounts(myobj,bed.obj_PVT)

meth_PVT=unite(dmr.drm2_PVT, destrand=FALSE)
perc.meth_PVT=percMethylation(meth_PVT)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)

# color break settings for heatmaps
# the length value must be equivalent to the color palette n value
# in this case being 299, the lengths are split 100 each
# the seq (start,stop) values determine the exact breaking point of the different colours
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("ddm1-col0_CHH_DMR_CGme_heatmap_10percDiff_allGenotypes.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "ddm1-col0 CHH-DMRs CGme% (#3530)")
dev.off()

pdf("boxplots of CGme at ddm1-col0 CHH DMRs percMeth.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "CG DNA Methylation (%) at ddm1-col0 CHH DMRs overlapping PVSloci (#3530)", outline=F, ylab="Average CHHme")
dev.off()



#### diff20%

bed.obj_PVT = readBed(file = "~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d20_q01_methylkit.bed")

file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/Col0_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/nrpe1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CHH")

dmr.drm2_PVT = regionCounts(myobj,bed.obj_PVT)

meth_PVT=unite(dmr.drm2_PVT, destrand=FALSE)
perc.meth_PVT=percMethylation(meth_PVT)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)

# color break settings for heatmaps
# the length value must be equivalent to the color palette n value
# in this case being 299, the lengths are split 100 each
# the seq (start,stop) values determine the exact breaking point of the different colours
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("ddm1-col0_heatmap_20percDiff_allGenotypes.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "ddm1-col0 CHH-DMRs (#371)")
dev.off()

pdf("boxplots at ddm1-col0 CHH DMRs 20percDiff percMeth.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at ddm1-col0 CHH DMRs 20percDiff overlapping PVSloci (#371)", outline=F, ylab="Average CHHme")
dev.off()

perc.meth_PVT_2 <- data.frame(perc.meth_PVT)
onlyPVS_noCol <- perc.meth_PVT_2[which(perc.meth_PVT_2$ddm1>(perc.meth_PVT_2$ddm1_nrpe1*2)),]




bed.obj_PVT = readBed(file = "~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d20_q01_methylkit.bed")

file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/Col0_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/nrpe1_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_CG_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CG")

dmr.drm2_PVT = regionCounts(myobj,bed.obj_PVT)

meth_PVT=unite(dmr.drm2_PVT, destrand=FALSE)
perc.meth_PVT=percMethylation(meth_PVT)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)

# color break settings for heatmaps
# the length value must be equivalent to the color palette n value
# in this case being 299, the lengths are split 100 each
# the seq (start,stop) values determine the exact breaking point of the different colours
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("ddm1-col0_CHH_DMR_CGme_heatmap_20percDiff_allGenotypes.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "ddm1-col0 CHH-DMRs CGme% (#371)")
dev.off()

pdf("boxplots of CGme at ddm1-col0 CHH DMRs 20percDiff percMeth.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "CG DNA Methylation (%) at ddm1-col0 CHH DMRs 20percDiff overlapping PVSloci (#371)", outline=F, ylab="Average CHHme")
dev.off()




bed.obj_PVT = readBed(file = "~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/mydiff_reg_ddm1-col0_CHH_col0_d10_q01_overlapPVS_methylkit.bed")

file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/Col0_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/nrpe1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1_labMap_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CHH")

dmr.drm2_PVT = regionCounts(myobj,bed.obj_PVT)

meth_PVT=unite(dmr.drm2_PVT, destrand=FALSE)
perc.meth_PVT=percMethylation(meth_PVT)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)

# color break settings for heatmaps
# the length value must be equivalent to the color palette n value
# in this case being 299, the lengths are split 100 each
# the seq (start,stop) values determine the exact breaking point of the different colours
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("ddm1-col0_CHH_DMR_CHme_heatmap_10percDiff_allGenotypes_overlapPVS.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CHH DNA Methylation (%)",                 
          main = "ddm1-col0 CHH-DMRs overlapping PVS CHHme% (#560)")
dev.off()

pdf("boxplots of CHHme at ddm1-col0 CHH DMRs 10percDiff percMeth overlapping PVS.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "CHH DNA Methylation (%) at ddm1-col0 CHH DMRs 10percDiff overlapping PVSloci (#560)", outline=F, ylab="Average CHHme")
dev.off()




file.list_ddm1=list("../../methylation/DMR_calls_genomeWide/Col0_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/nrpe1_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1_labMap_CG_methylKit.txt","../../methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_CG_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CG")

dmr.drm2_PVT = regionCounts(myobj,bed.obj_PVT)

meth_PVT=unite(dmr.drm2_PVT, destrand=FALSE)
perc.meth_PVT=percMethylation(meth_PVT)
hmcol3 <- colorRampPalette(c("#fee8c8", "#fdbb84", "#e34a33"))(n=299)

# color break settings for heatmaps
# the length value must be equivalent to the color palette n value
# in this case being 299, the lengths are split 100 each
# the seq (start,stop) values determine the exact breaking point of the different colours
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("ddm1-col0_CHH_DMR_CGme_heatmap_10percDiff_allGenotypes_overlapPVS.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "ddm1-col0 CHH-DMRs overlapping PVS CGme% (#436)")
dev.off()

pdf("boxplots of CGme at ddm1-col0 CHH DMRs percMeth overlapping PVS.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "CG DNA Methylation (%) at ddm1-col0 CHH DMRs overlapping PVSloci (#436)", outline=F, ylab="Average CHHme")
dev.off()

