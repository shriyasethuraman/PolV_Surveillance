setwd("~/cifs-lab/PolV_surveillance_manuscript/Figures/Fig2/S2A/")

percentages_table <- c(20.351,22.087,57.562)

pdf("pieChart_percent_postFilter_coverage_HMM_states.pdf")
pie(percentages_table,labels = c("RdDM","PVS","others"), col = c("turquoise3","indianred2","peachpuff2"), main="postfilter coverage")
dev.off()

