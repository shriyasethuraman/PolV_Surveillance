setwd("~/cifs-lab/PolV_surveillance_manuscript/Figures/Fig2/2C/")

table1 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em0_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed")
table2 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em1_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed")
table3 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em2_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed")
table4 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em3_col0_nrpe1_RIP_col0_nrpd1_smRNA.bed")

# png("col0 and nrpe1 RIP RPM and CHHme at every emission.png", width=1500, height=1200)
# boxplot(table1$V9/5.95,table1$V10/3.63,table2$V9/5.95,table2$V10/3.63,table3$V9/5.95,table3$V10/3.63,table4$V9/5.95,table4$V10/3.63, outline=F, names=c("col0-em0","nrpd1-em0","col0-em1","nrpd1-em1","col0-em2","nrpd1-em2","col0-em3","nrpd1-em3"),las=2,lwd=3, main="col0 and nrpd1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="RPMs")
# abline(h=0,col="grey")
# 
# pdf("col0 over nrpd1 siRNA RPM at every emission.pdf", width=10, height=8)
# boxplot(((table1$V9+0.1)/5.95)/((table1$V10+0.1)/3.63),((table2$V9+0.1)/5.95)/((table2$V10+0.1)/3.63),((table3$V9+0.1)/5.95)/((table3$V10+0.1)/3.63),((table4$V9+0.1)/5.95)/((table4$V10+0.1)/3.63), outline=F, names=c("em0","em1","em2","em3"),las=2,lwd=3, main="col0/nrpd1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="col0/nrpd1 siRNA RPMs", col="turquoise3")
# abline(h=1,col="grey")
# dev.off()
# 
# pdf("col0 over nrpd1 minNorm siRNA RPM at every emission.pdf", width=10, height=8)
# boxplot(((table1$V9+min(table1[(which(table1$V10>0)),10], na.rm=T))/5.95)/((table1$V10+min(table1[(which(table1$V10>0)),10], na.rm=T))/3.63),((table2$V9+min(table2[(which(table2$V10>0)),10], na.rm=T))/5.95)/((table2$V10+min(table2[(which(table2$V10>0)),10], na.rm=T))/3.63),((table3$V9+min(table3[(which(table3$V10>0)),10], na.rm=T))/5.95)/((table3$V10+min(table3[(which(table3$V10>0)),10], na.rm=T))/3.63),((table4$V9+min(table4[(which(table4$V10>0)),10], na.rm=T))/5.95)/((table4$V10+min(table4[(which(table4$V10>0)),10], na.rm=T))/3.63), outline=F, names=c("em0","em1","em2","em3"),las=2,lwd=3, main="col0/nrpd1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="col0/nrpd1 siRNA RPMs", col="turquoise3")
# abline(h=1,col="grey")
# dev.off()

pdf("log2 col0 over nrpd1 siRNA RPM at every emission.pdf", width=10, height=8)
boxplot(log2(((table1$V9+0.1)/5.95)/((table1$V10+0.1)/3.63)),log2(((table2$V9+0.1)/5.95)/((table2$V10+0.1)/3.63)),log2(((table3$V9+0.1)/5.95)/((table3$V10+0.1)/3.63)),log2(((table4$V9+0.1)/5.95)/((table4$V10+0.1)/3.63)), outline=F, names=c("em0","em1","em2","em3"),las=2,lwd=3, main="col0/nrpd1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="log2(col0/nrpd1) siRNA RPMs", col=c("turquoise3","indianred2","peachpuff2","peachpuff2"))
abline(h=0,col="grey")
dev.off()

# pdf("log2 col0 over nrpd1 minNorm siRNA RPM at every emission.pdf", width=10, height=8)
# boxplot(log2(((table1$V9+min(table1[(which(table1$V10>0)),10], na.rm=T))/5.95)/((table1$V10+min(table1[(which(table1$V10>0)),10], na.rm=T))/3.63)),log2(((table2$V9+min(table2[(which(table2$V10>0)),10], na.rm=T))/5.95)/((table2$V10+min(table2[(which(table2$V10>0)),10], na.rm=T))/3.63)),log2(((table3$V9+min(table3[(which(table3$V10>0)),10], na.rm=T))/5.95)/((table3$V10+min(table2[(which(table2$V10>0)),10], na.rm=T))/3.63)),log2(((table4$V9+min(table4[(which(table4$V10>0)),10], na.rm=T))/5.95)/((table4$V10+min(table2[(which(table2$V10>0)),10], na.rm=T))/3.63)), outline=F, names=c("em0","em1","em2","em3"),las=2,lwd=3, main="col0/nrpd1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="log2(col0/nrpd1) siRNA RPMs", col=c("turquoise3","peachpuff2","peachpuff2","indianred2"))
# abline(h=0,col="grey")
# dev.off()

wilcox.test(log2(((table1$V9+0.1)/5.95)/((table1$V10+0.1)/3.63)),log2(((table2$V9+0.1)/5.95)/((table2$V10+0.1)/3.63)))
wilcox.test(log2(((table1$V9+0.1)/5.95)/((table1$V10+0.1)/3.63)),log2(((table3$V9+0.1)/5.95)/((table3$V10+0.1)/3.63)))
wilcox.test(log2(((table1$V9+0.1)/5.95)/((table1$V10+0.1)/3.63)),log2(((table4$V9+0.1)/5.95)/((table4$V10+0.1)/3.63)))

wilcox.test(log2(((table2$V9+0.1)/5.95)/((table2$V10+0.1)/3.63)),log2(((table3$V9+0.1)/5.95)/((table3$V10+0.1)/3.63)))
wilcox.test(log2(((table2$V9+0.1)/5.95)/((table2$V10+0.1)/3.63)),log2(((table4$V9+0.1)/5.95)/((table4$V10+0.1)/3.63)))
