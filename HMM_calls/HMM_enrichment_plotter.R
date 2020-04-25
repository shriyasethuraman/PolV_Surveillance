setwd("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/")

table1 <- read.table("em0_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table2 <- read.table("em1_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table3 <- read.table("em2_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table4 <- read.table("em3_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")

table1 <- read.table("em0_col0_nrpe1_RIP.bed")
table2 <- read.table("em1_col0_nrpe1_RIP.bed")
table3 <- read.table("em2_col0_nrpe1_RIP.bed")
table4 <- read.table("em3_col0_nrpe1_RIP.bed")


# png("log2 col0 over nrpe1 RIP read counts and CHHme at every emission.png", width=1500, height=1200)
# par(mfrow=c(2,1))
# boxplot(log2(((table1$V7+0.1)*3.93)/((table1$V8+0.1)*7.93)),log2(((table2$V7+0.1)*3.93)/((table2$V8+0.1)*7.93)),log2(((table3$V7+0.1)*3.93)/((table3$V8+0.1)*7.93)),log2(((table4$V7+0.1)*3.93)/((table4$V8+0.1)*7.93)), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="log2(Col0/nrpe1) read count ratio", main="log2(col0/nrpe1) RIP read counts at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain")
# abline(h=0,col="grey")
# boxplot((table1$V9-table1$V10),(table2$V9-table2$V10),(table3$V9-table3$V10),(table4$V9-table4$V10), lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="Col0-nrpe1 CHHme", main="(col0-nrpe1) CHH methylation at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain")
# dev.off()

png("log2 col0 over nrpe1 RIP read counts at every emission.png", width=1500, height=1200)
pdf("log2 col0 over nrpe1 RIP read counts at every emission.pdf", width=10, height=8)
#pdf("~/cifs-lab/Shriya/Samples/manuscript_prep/PVS/Figures/Fig1/1E/log2 col0 over nrpe1 RIP read counts at every emission.pdf", width=10, height=8)
boxplot(log2(((table1$V7+0.1)*3.93)/((table1$V8+0.1)*7.93)),log2(((table2$V7+0.1)*3.93)/((table2$V8+0.1)*7.93)),log2(((table3$V7+0.1)*3.93)/((table3$V8+0.1)*7.93)),log2(((table4$V7+0.1)*3.93)/((table4$V8+0.1)*7.93)), lwd=4,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="log2(Col0/nrpe1) read count ratio", main="log2(col0/nrpe1) RIP read counts at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain",cex.lab=2,cex.axis=2,col=c("turquoise3","indianred2","peachpuff2","peachpuff2","indianred2"))
abline(h=0,col="grey")
dev.off()

###### PolV vs non-PolV

PolV <- rbind(table1,table2)
nonPolV <- rbind(table3,table4)

png("log2 col0 over nrpe1 RIP read counts at PolV_nonPV.png", width=1500, height=1200)
pdf("log2 col0 over nrpe1 RIP read counts at PolV_nonPV.pdf", width=10, height=8)
boxplot(log2(((PolV$V7+0.1)*3.93)/((PolV$V8+0.1)*7.93)),log2(((nonPolV$V7+0.1)*3.93)/((nonPolV$V8+0.1)*7.93)), lwd=4,outline=F, names=c("PolV regions","nonPolV regions"), ylab="log2(Col0/nrpe1) RPM ratio", main="log2(col0/nrpe1) RIP RPM at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain",cex.lab=2,cex.axis=2,col = c("#91bca6","peachpuff2"))
abline(h=0,col="grey")
dev.off()

noNA1 <- PolV[-which(is.na(PolV$V9)),]
PolV_2 <- noNA1[-which(is.na(noNA1$V10)),]
noNA1 <- nonPolV[-which(is.na(nonPolV$V9)),]
nonPolV_2 <- noNA1[-which(is.na(noNA1$V10)),]



boxplot(((PolV_2$V9-PolV_2$V10)*100),((nonPolV_2$V9-nonPolV_2$V10)*100), lwd=3,outline=F, names=c("PolV regions","nonPolV regions"), ylab="Col0-nrpe1 CHHme", main="(col0-nrpe1) CHH methylation at every emission of 200bp-150bp slide output3_geneNorm HMM bothContain")
