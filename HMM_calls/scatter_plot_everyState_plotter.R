setwd("~/cifs-lab/Shriya/Samples/manuscript_prep/PVS/Figures/Fig1/1E/")

setwd("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/")

table1 <- read.table("em0_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table2 <- read.table("em1_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table3 <- read.table("em2_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table4 <- read.table("em3_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")


library("ggplot2")
png("plot_emissions_RIP_scatter.png", width=800, height=600, units="px")
ggplot(data=final_table, aes(x = log2(V7/2.685), y = log2(V8/2.052), color=factor(V11))) + xlab("log2(Col0)") + ylab("log2(nrpe1)") + geom_point(shape=1) + geom_smooth(method="lm", se=FALSE) + #+ geom_abline(lm((log2(V7)) ~ (log2(V8)), data = final_table))#+ geom_abline(intercept = log2(18.0987), slope = log2(2.37))
  geom_abline(intercept = 0, slope = 1,col="grey")
dev.off()


# png("col0 and nrpe1 RIP read counts at every emission.png", width=1500, height=1200)
# boxplot(table1$V7,table1$V8,table2$V7,table2$V8,table3$V7,table3$V8,table4$V7,table4$V8,table5$V7,table5$V8,table6$V7,table6$V8,table7$V7,table7$V8,table8$V7,table8$V8, outline=F, names=c("col0-em0","nrpe1-em0","col0-em1","nrpe1-em1","col0-em2","nrpe1-em2","col0-em3","nrpe1-em3","col0-em4","nrpe1-em4","col0-em5","nrpe1-em5","col0-em6","nrpe1-em6","col0-em7","nrpe1-em7"),las=2,lwd=3, main="col0 and nrpe1 RIP read counts at every emission of 60bp-40bp slide output3", ylab="read counts")
# abline(h=1,col="grey")
# dev.off()
# 
# png("col0 over nrpe1 RIP read counts at every emission.png", width=1500, height=1200)
# boxplot((table1$V7+0.1)/(table1$V8+0.1),(table2$V7+0.1)/(table2$V8+0.1),(table3$V7+0.1)/(table3$V8+0.1),(table4$V7+0.1)/(table4$V8+0.1), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="Col0/nrpe1 read count ratio", main="col0/nrpe1 RIP read counts at every emission of 40bp output1")
# dev.off()

png("col0 and nrpe1 RIP RPM and CHHme at every emission.png", width=1500, height=1200)
par(mfrow=c(2,1))
boxplot(table1$V7/10.11,table1$V8/5.27,table2$V7/10.11,table2$V8/5.27,table3$V7/10.11,table3$V8/5.27,table4$V7/10.11,table4$V8/5.27, outline=F, names=c("col0-em0","nrpe1-em0","col0-em1","nrpe1-em1","col0-em2","nrpe1-em2","col0-em3","nrpe1-em3"),las=2,lwd=3, main="col0 and nrpe1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="RPMs")
abline(h=0,col="grey")
boxplot(table1$V9,table1$V10,table2$V9,table2$V10,table3$V9,table3$V10,table4$V9,table4$V10, outline=F, names=c("col0-em0","nrpe1-em0","col0-em1","nrpe1-em1","col0-em2","nrpe1-em2","col0-em3","nrpe1-em3"),las=2,lwd=3, main="col0 and nrpe1 RIP read counts at every emission of 200bp-150bp slide output3 BothContainConvert", ylab="CHH methylation")
dev.off()


png("col0 over nrpe1 RIP read counts and CHHme at every emission.png", width=1500, height=1200)
par(mfrow=c(2,1))
boxplot(((table1$V7+0.1)*5.27)/((table1$V8+0.1)*10.11),((table2$V7+0.1)*5.27)/((table2$V8+0.1)*10.11),((table3$V7+0.1)*5.27)/((table3$V8+0.1)*10.11),((table4$V7+0.1)*5.27)/((table4$V8+0.1)*10.11), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="Col0/nrpe1 read count ratio", main="col0/nrpe1 RIP read counts at every emission of 200bp-150bp slide output3 (normalized HMM input ratio)")
abline(h=1,col="grey")
boxplot((table1$V9-table1$V10),(table2$V9-table2$V10),(table3$V9-table3$V10),(table4$V9-table4$V10), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="Col0-nrpe1 CHHme", main="(col0-nrpe1) CHH methylation at every emission of 200bp-150bp slide output3 BothContainConvert")
dev.off()

png("log2 col0 over nrpe1 RIP read counts and CHHme at every emission.png", width=1500, height=1200)
par(mfrow=c(2,1))
boxplot(log2(((table1$V7+0.1)*5.27)/((table1$V8+0.1)*10.11)),log2(((table2$V7+0.1)*5.27)/((table2$V8+0.1)*10.11)),log2(((table3$V7+0.1)*5.27)/((table3$V8+0.1)*10.11)),log2(((table4$V7+0.1)*5.27)/((table4$V8+0.1)*10.11)), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="log2(Col0/nrpe1) read count ratio", main="log2(col0/nrpe1) RIP read counts at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain")
abline(h=0,col="grey")
boxplot((table1$V9-table1$V10),(table2$V9-table2$V10),(table3$V9-table3$V10),(table4$V9-table4$V10), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="Col0-nrpe1 CHHme", main="(col0-nrpe1) CHH methylation at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain")
dev.off()

png("log2 col0 over nrpe1 RIP read counts at every emission.png", width=1500, height=1200)
pdf("~/cifs-lab/Shriya/Samples/manuscript_prep/PVS/Figures/Fig1/1E/log2 col0 over nrpe1 RIP read counts at every emission.pdf", width=10, height=8)
boxplot(log2(((table1$V7+0.1)*5.27)/((table1$V8+0.1)*10.11)),log2(((table2$V7+0.1)*5.27)/((table2$V8+0.1)*10.11)),log2(((table3$V7+0.1)*5.27)/((table3$V8+0.1)*10.11)),log2(((table4$V7+0.1)*5.27)/((table4$V8+0.1)*10.11)), lwd=4,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="log2(Col0/nrpe1) read count ratio", main="log2(col0/nrpe1) RIP read counts at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain",cex.lab=2,cex.axis=2)
abline(h=0,col="grey")
dev.off()

table1$V11 <- "em0"
table2$V11 <- "em1"
table3$V11 <- "em2"
table4$V11 <- "em3"
final_table <- rbind(table1,table2,table3,table4)

library("ggplot2")
png("plot_emissions_RIP_scatter.png", width=800, height=600, units="px")
ggplot(data=final_table, aes(x = log2(V7/2.685), y = log2(V8/2.052), color=factor(V11))) + xlab("log2(Col0)") + ylab("log2(nrpe1)") + geom_point(shape=1) + geom_smooth(method="lm", se=FALSE) + #+ geom_abline(lm((log2(V7)) ~ (log2(V8)), data = final_table))#+ geom_abline(intercept = log2(18.0987), slope = log2(2.37))
  geom_abline(intercept = 0, slope = 1,col="grey")
dev.off()



#### ONLY PVS and RdDM
RdDM <- read.table("em0_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
PVS <- read.table("em3_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")

png("col0 and nrpe1 RIP read counts (RPM) and CHHme at PVS and RdDMs.png", width=1500, height=1200)
par(mfrow=c(2,1))
boxplot(PVS$V7/10.11,PVS$V8/5.27,RdDM$V7/10.11,RdDM$V8/5.27, outline=F, names=c("PVS-Col0","PVS-nrpe1","RdDM-Col0","RdDM-nrpe1"),las=2,lwd=3, main="col0 and nrpe1 RIP read counts (RPM) at PVS and RdDM loci for output3-200bp-150bp slide HMM calls", ylab="read counts (RPM)")
abline(h=0,col="grey")
boxplot(PVS$V9,PVS$V10,RdDM$V9,RdDM$V10, outline=F, names=c("PVS-Col0","PVS-nrpe1","RdDM-Col0","RdDM-nrpe1"),las=2,lwd=3, main="col0 and nrpe1 CHHme at PVS and RdDM loci for output3-200bp-150bp slide HMM calls", ylab="CHH methylation")
dev.off()

png("col0 over nrpe1 RIP read counts (RPM) and CHHme at PVS and RdDMs.png", width=1500, height=1200)
par(mfrow=c(2,1))
boxplot(((PVS$V7+0.1)*5.27)/((PVS$V8+0.1)*10.11),((RdDM$V7+0.1)*5.27)/((RdDM$V8+0.1)*10.11), lwd=3,outline=F, names=c("PVS","RdDMs"),las=2, ylab="Col0/nrpe1 read count (RPM) ratio", main="col0/nrpe1 RIP read counts (RPM) at PVS and RdDM loci for output3-200bp-150bp HMM calls")
abline(h=1,col="grey")
boxplot((PVS$V9-PVS$V10),(RdDM$V9-RdDM$V10), lwd=3,outline=F, names=c("PVS","RdDMs"),las=2, ylab="Col0-nrpe1 CHHme", main="(col0-nrpe1) CHH methylation at PVS and RdDM loci for output3-200bp-150bp HMM calls")
dev.off()


###### PolV vs non-PolV

PolV <- rbind(table1[,-1],table4[,-1])
nonPolV <- rbind(table2[,-1],table3[,-1])

png("log2 col0 over nrpe1 RIP read counts at PolV_nonPV.png", width=1500, height=1200)
pdf("log2 col0 over nrpe1 RIP read counts at PolV_nonPV.pdf", width=10, height=8)
boxplot(log2(((PolV$V7+0.1)*5.27)/((PolV$V8+0.1)*10.11)),log2(((nonPolV$V7+0.1)*5.27)/((nonPolV$V8+0.1)*10.11)), lwd=4,outline=F, names=c("PolV regions","nonPolV regions"), ylab="log2(Col0/nrpe1) RPM ratio", main="log2(col0/nrpe1) RIP RPM at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain",cex.lab=2,cex.axis=2,col = c("#91bca6","peachpuff2"))
abline(h=0,col="grey")
dev.off()
