source('~/cifs-lab/RIP_manuscript/Revised_Figures/Fig7/Fig_7A/hp-profile-heatmap_modGB.r')
source('~/cifs-lab/Shriya/Codes/window_generator_mod.R')
library('methylKit')
library('gplots')
setwd("~/cifs-lab/Shriya/Samples/Keith_Panda_datasets/new_PVS_2020/col0_ddm1_CHH_DMR/each_state_DMRoverlap/")
library(devtools)
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
library(genomation)

bed.obj_PVT = readBed(file = "PVS_ddm1Col_DMR_overlap_methylKit.bed")

file.list_ddm1=list("/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/Col0_labMap_methylKit.txt","/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/nrpe1_labMap_methylKit.txt","/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/ddm1_labMap_methylKit.txt","/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CHH")
#tiles=tileMethylCounts(myobj,win.size=100,step.size=80,cov.bases = 10)

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
pdf("heatmap_PVS_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "PVS overlap ddm1-col0 CHH-DMR (#285)")
dev.off()

pdf("boxplots_PVS_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at PVS overlap ddm1-col0 CHH-DMR (#285)", outline=F, ylab="Average CHHme")
dev.off()

pdf("boxplots_PVS_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at PVS overlap ddm1-col0 CHH-DMR (#285)", outline=F, ylab="Average CHHme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()


hmcol3 <- colorRampPalette(c("#fbdfdf", "#ee6363", "#8e3b3b"))(n=299)
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

pdf("colored_heatmap_PVS_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "PVS overlap ddm1-col0 CHH-DMR (#285)")
dev.off()

pdf("boxplots_colored_PVS_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at PVS overlap ddm1-col0 CHH-DMR (#285)", outline=F, ylab="Average CHHme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()

wilcox.test(perc.meth_PVT[,3],perc.meth_PVT[,4])
wilcox.test(perc.meth_PVT[,3],perc.meth_PVT[,1])
wilcox.test(perc.meth_PVT[,3],perc.meth_PVT[,2])


bed.obj_PVT = readBed(file = "RdDM_ddm1Col_DMR_overlap_methylKit.bed")

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
pdf("heatmap_RdDM_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "RdDM overlap ddm1-col0 CHH-DMR (#3418)")
dev.off()

pdf("boxplots_RdDM_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at RdDM overlap ddm1-col0 CHH-DMR (#3418)", outline=F, ylab="Average CHHme")
dev.off()

pdf("boxplots_RdDM_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at RdDM overlap ddm1-col0 CHH-DMR (#3418)", outline=F, ylab="Average CHHme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()



hmcol3 <- colorRampPalette(c("#d8f8f5", "#40e0d0", "#195953"))(n=299)

col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

pdf("colored_heatmap_RdDM_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "RdDM overlap ddm1-col0 CHH-DMR (#3418)")
dev.off()

pdf("boxplots_colored_RdDM_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at RdDM overlap ddm1-col0 CHH-DMR (#3418)", outline=F, ylab="Average CHHme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()



bed.obj_PVT = readBed(file = "other_ddm1Col_DMR_overlap_methylKit.bed")

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
pdf("heatmap_other_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "other overlap ddm1-col0 CHH-DMR (#57)")
dev.off()

pdf("boxplots_other_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at other overlap ddm1-col0 CHH-DMR (#57)", outline=F, ylab="Average CHHme")
dev.off()

pdf("boxplots_other_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at other overlap ddm1-col0 CHH-DMR (#57)", outline=F, ylab="Average CHHme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()



hmcol3 <- colorRampPalette(c("#f8eade", "#eecbad", "#473c33"))(n=299)

col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

pdf("colored_heatmap_other_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "DNA Methylation (%)",                 
          main = "other overlap ddm1-col0 CHH-DMR (#57)")
dev.off()

pdf("boxplots_colored_other_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CHH context at other overlap ddm1-col0 CHH-DMR (#57)", outline=F, ylab="Average CHHme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()




#### CGme

bed.obj_PVT = readBed(file = "PVS_ddm1Col_DMR_overlap_methylKit.bed")

file.list_ddm1=list("/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/Col0_labMap_CG_methylKit.txt","/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/nrpe1_labMap_CG_methylKit.txt","/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/ddm1_labMap_CG_methylKit.txt","/home/shriyas/cifs-lab/Shriya/Samples/Keith_Panda_datasets/methylation/DMR_calls_genomeWide/ddm1nrpe1_labMap_CG_methylKit.txt")
myobj=methRead(file.list_ddm1,
               sample.id=list("col0","nrpe1","ddm1","ddm1_nrpe1"),
               assembly="tair10",
               treatment=c(0,1,2,3),
               context="CG")
#tiles=tileMethylCounts(myobj,win.size=100,step.size=80,cov.bases = 10)

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
pdf("CG_heatmap_PVS_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "CG DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "PVS overlap ddm1-col0 CHH-DMR CGme (#285)")
dev.off()

pdf("CG_boxplots_PVS_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at PVS overlap ddm1-col0 CHH-DMR (#285)", outline=F, ylab="Average CGme")
dev.off()

pdf("CG_boxplots_PVS_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at PVS overlap ddm1-col0 CHH-DMR (#285)", outline=F, ylab="Average CGme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()



hmcol3 <- colorRampPalette(c("#fbdfdf", "#ee6363", "#8e3b3b"))(n=299)
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

pdf("CG_colored_heatmap_PVS_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "PVS overlap ddm1-col0 CHH-DMR CGme (#285)")
dev.off()

pdf("CG_boxplots_colored_PVS_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at PVS overlap ddm1-col0 CHH-DMR (#285)", outline=F, ylab="Average CGme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()

wilcox.test(perc.meth_PVT[,3],perc.meth_PVT[,4]) #p-value = 2.215e-11
wilcox.test(perc.meth_PVT[,1],perc.meth_PVT[,3]) #p-value = 5.553e-13

bed.obj_PVT = readBed(file = "RdDM_ddm1Col_DMR_overlap_methylKit.bed")

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
pdf("CG_heatmap_RdDM_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "CG DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "RdDM overlap ddm1-col0 CHH-DMR CGme (#3418)")
dev.off()

pdf("CG_boxplots_RdDM_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at RdDM overlap ddm1-col0 CHH-DMR (#3418)", outline=F, ylab="Average CGme")
dev.off()

pdf("CG_boxplots_RdDM_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at RdDM overlap ddm1-col0 CHH-DMR (#3418)", outline=F, ylab="Average CGme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()


hmcol3 <- colorRampPalette(c("#d8f8f5", "#40e0d0", "#195953"))(n=299)
col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

#generate heatmap with heatmap function (perc.meth_100)
pdf("CG_colored_heatmap_RdDM_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "CG DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "RdDM overlap ddm1-col0 CHH-DMR CGme (#3418)")
dev.off()

pdf("CG_boxplots_colored_RdDM_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at RdDM overlap ddm1-col0 CHH-DMR (#3418)", outline=F, ylab="Average CGme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()


bed.obj_PVT = readBed(file = "other_ddm1Col_DMR_overlap_methylKit.bed")

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
pdf("CG_heatmap_other_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "CG DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "other overlap ddm1-col0 CHH-DMR CGme (#57)")
dev.off()

pdf("CG_boxplots_other_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at other overlap ddm1-col0 CHH-DMR (#57)", outline=F, ylab="Average CGme")
dev.off()

pdf("CG_boxplots_other_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at other overlap ddm1-col0 CHH-DMR (#57)", outline=F, ylab="Average CGme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()


hmcol3 <- colorRampPalette(c("#f8eade", "#eecbad", "#473c33"))(n=299)

col_breaks = c(seq(quantile(perc.meth_PVT, 0.1),quantile(perc.meth_PVT,0.9), length = 300))

pdf("CG_colored_heatmap_other_ddm1Col_DMR_overlap_loci.pdf", width=8, height=5)
heatmap.3(perc.meth_PVT,
          col = hmcol3, breaks = col_breaks,
          cexCol = 2, trace = "none", density.info= "none",key.xlab = "CG DNA Methylation (%)",
          margins = c(11,3), symm = F,Colv = F,dendrogram = "none",
          labRow = F, xlab = "genotypes", ylab =  "CG DNA Methylation (%)",                 
          main = "other overlap ddm1-col0 CHH-DMR CGme (#57)")
dev.off()

pdf("CG_boxplots_colored_other_ddm1Col_DMR_overlap_loci_2.pdf", width=8, height=5)
boxplot(perc.meth_PVT, main = "DNA Methylation (%) CG context at other overlap ddm1-col0 CHH-DMR (#57)", outline=F, ylab="Average CGme", col=c(hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,1])))==min(abs(col_breaks-(mean(perc.meth_PVT[,1])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,2])))==min(abs(col_breaks-(mean(perc.meth_PVT[,2])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,3])))==min(abs(col_breaks-(mean(perc.meth_PVT[,3])))))],hmcol3[which(abs(col_breaks-(mean(perc.meth_PVT[,4])))==min(abs(col_breaks-(mean(perc.meth_PVT[,4])))))]))
dev.off()
