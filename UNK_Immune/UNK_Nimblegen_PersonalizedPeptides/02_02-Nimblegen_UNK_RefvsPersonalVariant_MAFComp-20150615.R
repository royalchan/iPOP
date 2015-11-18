############################################################
# TO-DO LIST.
############################################################
DONE. # 1. This attempt is to compare the behavior of peptides containing personal variants with their reference counterparts.





#############################################################################################
#############################################################################################
#############################################################################################
### QN U-TEST SUMMARY (AFTER SHAPIRO-WILK TEST):                                          ###
###                                                                                       ###
### HRV1 vs HEALTHY: P < 0.005;  63 Hits (No Significant Hits Post FDR).                  ###
### HRV2 vs HEALTHY: P < 0.005; 151 Hits (No Significant Hits Post FDR).                  ###
### HRV3 vs HEALTHY: P < 0.005; 132 Hits (No Significant Hits Post FDR).                  ###
### RSV  vs HEALTHY: P < 0.005; 366 Hits (No Significant Hits Post FDR).                  ###
### ADV1 vs HEALTHY: P < 0.005; 380 Hits (No Significant Hits Post FDR).                  ###
### ADV2 vs HEALTHY: P < 0.005; 217 Hits (FDR_P < 0.05; 3 Hits).                          ###
### HRVALL vs HEALTHY: FDR_P < 0.05; 10 Hits.                                             ###
### ADVALL vs HEALTHY: FDR_P < 0.05; 207 Hits.                                            ###
### INFALL vs HEALTHY: P < 0.005; 295 Hits (No Significant Hits Post FDR).                ###
### AUTOCORRELATED (Cluster 1, Pre-SVA): 3290 Hits.                                       ###
### AUTOCORRELATED (Cluster 2, Pre-SVA): 2502 Hits.                                       ###
#############################################################################################
#############################################################################################
#############################################################################################





############################################################
# A. Library Used.
############################################################
library("gplots")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library("ggplot2")
############
# LOAD APPENDIX FUNCTIONS.
############
# Package A2R:
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
# colored dendrogram op = par(bg = "#EFEFEF") A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("#FF6B6B","#4ECDC4", "#556270"))
# another colored dendrogram
op = par(bg = "gray15") cols = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = cols)







############################################################
# B. U TEST HITS UPDATED WITH UNIFIED CUT OFF (P < 0.005 OR FDR-P < 0.05).
############################################################
############################################################
# B.1. Raw Data Loading (All Quantiles).
############################################################
# Load data.

dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/UNK_NimbleGen_OnArray_RareVar_MAF.txt",  header =TRUE, sep = "\t")

# Column Names.
# CHR_COORD_REF_ALT	ESP_EUR_MAF	1000G_EUR_MAF	PROBE_ID	P_Value	FDR_P	Difference_MeanRef-MeanPsn	MeanRef	MeanPsn	D_.123_Ref	D_0_Ref	D_4_Ref	D_21_Ref	D_116_Ref	D_185_Ref	D_255_Ref	D_289_Ref	D_290_Ref	D_292_Ref	D_294_Ref	D_297_Ref	D_301_Ref	D_307_Ref	D_311_Ref	D_322_Ref	D_369_Ref	D_380_Ref	D_400_Ref	D_476_Ref	D_532_Ref	D_546_Ref	D_602_Ref	D_615_Ref	D_618_Ref	D_620_Ref	D_625_Ref	D_630_Ref	D_647_Ref	D_679_Ref	D_680_Ref	D_683_Ref	D_688_Ref	D_694_Ref	D_700_Ref	D_711_Ref	D_735_Ref	D_796_Ref	D_840_Ref	D_912_Ref	D_944_Ref	D_945_Ref	D_948_Ref	D_959_Ref	D_966_Ref	D_984_Ref	D_1029_Ref	D_1030_Ref	D_1032_Ref	D_1038_Ref	D_1045_Ref	D_.123_Psn	D_0_Psn	D_4_Psn	D_21_Psn	D_116_Psn	D_185_Psn	D_255_Psn	D_289_Psn	D_290_Psn	D_292_Psn	D_294_Psn	D_297_Psn	D_301_Psn	D_307_Psn	D_311_Psn	D_322_Psn	D_369_Psn	D_380_Psn	D_400_Psn	D_476_Psn	D_532_Psn	D_546_Psn	D_602_Psn	D_615_Psn	D_618_Psn	D_620_Psn	D_625_Psn	D_630_Psn	D_647_Psn	D_679_Psn	D_680_Psn	D_683_Psn	D_688_Psn	D_694_Psn	D_700_Psn	D_711_Psn	D_735_Psn	D_796_Psn	D_840_Psn	D_912_Psn	D_944_Psn	D_945_Psn	D_948_Psn	D_959_Psn	D_966_Psn	D_984_Psn	D_1029_Psn	D_1030_Psn	D_1032_Psn	D_1038_Psn	D_1045_Psn	Chr	Coordinate	QUERY_KEY	rsID	NimbleGen_Record	UNK_Nonsynonymous_Var


row.names(dat) <- dat$PROBE_ID
datrownames <- row.names(dat)
n <- length(dat[,1])


############################################################
# B.2. DIVIDE DATA BY ESP MAF.
############################################################
# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, Max].
dat_data_a <- as.matrix(dat[(1:n)[dat[,2] < 0.01],])
dat_data_b <- as.matrix(dat[(1:n)[dat[,2] >= 0.01 & dat[,2] < 0.05],])
dat_data_c <- as.matrix(dat[(1:n)[dat[,2] >= 0.05],])

# Mean Distribution.
par(mfrow = c(1,3))
datdataameanforbp <- data.frame(datmean = c(as.numeric(as.character(dat_data_a[,8])), as.numeric(as.character(dat_data_a[,9]))), category = c(rep("Reference Peptide", length(dat_data_a[,8])), rep("Personal Peptide", length(dat_data_a[,9]))))
boxplot(datdataameanforbp[,1]~ datdataameanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, \nMAF < 0.01", col = c("gold", "turquoise"), ylim = c(6,12))
# Add Lines to Box Plot.
for (i in 1:length(dat_data_a[,8])) {
	refpsn <- c(as.numeric(as.character(dat_data_a[i,9])), as.numeric(as.character(dat_data_a[i,8])))
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}
datdatabmeanforbp <- data.frame(datmean = c(as.numeric(as.character(dat_data_b[,8])), as.numeric(as.character(dat_data_b[,9]))), category = c(rep("Reference Peptide", length(dat_data_b[,8])), rep("Personal Peptide", length(dat_data_b[,9]))))
boxplot(datdatabmeanforbp[,1]~ datdatabmeanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, \n0.01 <= MAF < 0.05", col = c("gold", "turquoise"), ylim = c(6,12))
# Add Lines to Box Plot.
for (i in 1:length(dat_data_b[,8])) {
	refpsn <- c(as.numeric(as.character(dat_data_b[i,9])), as.numeric(as.character(dat_data_b[i,8])))
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}
datdatacmeanforbp <- data.frame(datmean = c(as.numeric(as.character(dat_data_c[,8])), as.numeric(as.character(dat_data_c[,9]))), category = c(rep("Reference Peptide", length(dat_data_c[,8])), rep("Personal Peptide", length(dat_data_c[,9]))))
boxplot(datdatacmeanforbp[,1]~ datdatacmeanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, \nMAF >= 0.05", col = c("gold", "turquoise"), ylim = c(6,12))
# Add Lines to Box Plot.
for (i in 1:length(dat_data_c[,8])) {
	refpsn <- c(as.numeric(as.character(dat_data_c[i,9])), as.numeric(as.character(dat_data_c[i,8])))
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,5]))), -log10(as.numeric(as.character(dat_data_b[,5]))), -log10(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF > 0.05", length(dat_data_c[,5]))))
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = "Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = 0.03727)", col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 6.5792, df = 2, p-value = 0.03727
############

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,6]))), -log10(as.numeric(as.character(dat_data_b[,6]))), -log10(as.numeric(as.character(dat_data_c[,6])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,6])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,6])), rep("MAF > 0.05", length(dat_data_c[,6]))))
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = "Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = 0.03997)", col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 6.4392, df = 2, p-value = 0.03997
############


############################################################
# B.3. DIVIDE DATA BY 1000G MAF.
############################################################
# Divide Peptides by 1000G MAF [0, 0.01), [0.01, 0.05), [0.05, Max].
dat_data_a <- as.matrix(dat[(1:n)[dat[,3] < 0.01],])
dat_data_b <- as.matrix(dat[(1:n)[dat[,3] >= 0.01 & dat[,3] < 0.05],])
dat_data_c <- as.matrix(dat[(1:n)[dat[,3] >= 0.05],])

# Mean Distribution.
par(mfrow = c(1,3))
datdataameanforbp <- data.frame(datmean = c(as.numeric(as.character(dat_data_a[,8])), as.numeric(as.character(dat_data_a[,9]))), category = c(rep("Reference Peptide", length(dat_data_a[,8])), rep("Personal Peptide", length(dat_data_a[,9]))))
boxplot(datdataameanforbp[,1]~ datdataameanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, \nMAF < 0.01", col = c("gold", "turquoise"), ylim = c(6,12))
# Add Lines to Box Plot.
for (i in 1:length(dat_data_a[,8])) {
	refpsn <- c(as.numeric(as.character(dat_data_a[i,9])), as.numeric(as.character(dat_data_a[i,8])))
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}
datdatabmeanforbp <- data.frame(datmean = c(as.numeric(as.character(dat_data_b[,8])), as.numeric(as.character(dat_data_b[,9]))), category = c(rep("Reference Peptide", length(dat_data_b[,8])), rep("Personal Peptide", length(dat_data_b[,9]))))
boxplot(datdatabmeanforbp[,1]~ datdatabmeanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, \n0.01 <= MAF < 0.05", col = c("gold", "turquoise"), ylim = c(6,12))
# Add Lines to Box Plot.
for (i in 1:length(dat_data_b[,8])) {
	refpsn <- c(as.numeric(as.character(dat_data_b[i,9])), as.numeric(as.character(dat_data_b[i,8])))
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}
datdatacmeanforbp <- data.frame(datmean = c(as.numeric(as.character(dat_data_c[,8])), as.numeric(as.character(dat_data_c[,9]))), category = c(rep("Reference Peptide", length(dat_data_c[,8])), rep("Personal Peptide", length(dat_data_c[,9]))))
boxplot(datdatacmeanforbp[,1]~ datdatacmeanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, \nMAF >= 0.05", col = c("gold", "turquoise"), ylim = c(6,12))
# Add Lines to Box Plot.
for (i in 1:length(dat_data_c[,8])) {
	refpsn <- c(as.numeric(as.character(dat_data_c[i,9])), as.numeric(as.character(dat_data_c[i,8])))
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,5]))), -log10(as.numeric(as.character(dat_data_b[,5]))), -log10(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF > 0.05", length(dat_data_c[,5]))))
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = "Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = 0.0005274)", col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 15.0951, df = 2, p-value = 0.0005274
############

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,6]))), -log10(as.numeric(as.character(dat_data_b[,6]))), -log10(as.numeric(as.character(dat_data_c[,6])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,6])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,6])), rep("MAF > 0.05", length(dat_data_c[,6]))))
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = "Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = 0.0005692)", col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 14.9425, df = 2, p-value = 0.0005692
############




