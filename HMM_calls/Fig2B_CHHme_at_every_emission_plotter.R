setwd("~/cifs-lab/PolV_surveillance_manuscript/Figures/Fig2/2B/")

table1 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em0_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table2 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em1_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table3 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em2_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")
table4 <- read.table("~/cifs-lab/PolV_surveillance_manuscript/Bioinformatics/HMM_calling/our_fav_HMM/em3_col0_nrpe1_RIP_col0-nrpe1_CHH.bed")

png("col0 minus nrpe1 CHHme at every emission.png", width=1500, height=1200)
pdf("col0 minus nrpe1 CHHme at every emission.pdf", width=10, height=8)
boxplot((table1$V9-table1$V10)*100,(table2$V9-table2$V10)*100,(table3$V9-table3$V10)*100,(table4$V9-table4$V10)*100, lwd=3,outline=F, names=c("em0","em1","em2","em3"),las=2, ylab="Col0-nrpe1 CHHme%", col=c("turquoise3","peachpuff2","peachpuff2","indianred2"), main="(col0-nrpe1) CHH methylation at every emission of 200bp-150bp slide output3 BothContainConvert")
abline(h=0,col="grey")
dev.off()

wilcox.test((table1$V9-table1$V10)*100,(table2$V9-table2$V10)*100)
wilcox.test((table1$V9-table1$V10)*100,(table3$V9-table3$V10)*100)
wilcox.test((table1$V9-table1$V10)*100,(table4$V9-table4$V10)*100)

wilcox.test((table2$V9-table2$V10)*100,(table3$V9-table3$V10)*100) #p-value = 0.4164
wilcox.test((table2$V9-table2$V10)*100,(table4$V9-table4$V10)*100) #p-value = 0.3474
wilcox.test((table3$V9-table3$V10)*100,(table4$V9-table4$V10)*100) #p-value = 0.2275
