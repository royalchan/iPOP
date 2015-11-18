############################################################
# TO-DO LIST.
############################################################
# This is an extended version of the 03-Nimblegen_UNK_RefvsPersonalVariant_MAFComp_MoreCat-20151003 file.
# SVA background correction was not performed due to between-analyte comparisons.
DONE. # 1. This attempt is to plot mean difference as well as different categorization schemes of the UNK reactivities.





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
library("vioplot")


############################################################
# B. U TEST P-VALUE DISTRIBUTION BY MAF IN DIFFERENT GROUPS.
############################################################
############################################################
# B.1. RAW DATA LOADING.
############################################################
# Load data.

dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/UNK_NimbleGen_OnArray_RareVar_MAF.txt",  header =TRUE, sep = "\t")

# Column Names.
# CHR_COORD_REF_ALT	ESP_EUR_MAF	1000G_EUR_MAF	PROBE_ID	P_Value	FDR_P	Difference_MeanRef-MeanPsn	MeanRef	MeanPsn	D_.123_Ref	D_0_Ref	D_4_Ref	D_21_Ref	D_116_Ref	D_185_Ref	D_255_Ref	D_289_Ref	D_290_Ref	D_292_Ref	D_294_Ref	D_297_Ref	D_301_Ref	D_307_Ref	D_311_Ref	D_322_Ref	D_369_Ref	D_380_Ref	D_400_Ref	D_476_Ref	D_532_Ref	D_546_Ref	D_602_Ref	D_615_Ref	D_618_Ref	D_620_Ref	D_625_Ref	D_630_Ref	D_647_Ref	D_679_Ref	D_680_Ref	D_683_Ref	D_688_Ref	D_694_Ref	D_700_Ref	D_711_Ref	D_735_Ref	D_796_Ref	D_840_Ref	D_912_Ref	D_944_Ref	D_945_Ref	D_948_Ref	D_959_Ref	D_966_Ref	D_984_Ref	D_1029_Ref	D_1030_Ref	D_1032_Ref	D_1038_Ref	D_1045_Ref	D_.123_Psn	D_0_Psn	D_4_Psn	D_21_Psn	D_116_Psn	D_185_Psn	D_255_Psn	D_289_Psn	D_290_Psn	D_292_Psn	D_294_Psn	D_297_Psn	D_301_Psn	D_307_Psn	D_311_Psn	D_322_Psn	D_369_Psn	D_380_Psn	D_400_Psn	D_476_Psn	D_532_Psn	D_546_Psn	D_602_Psn	D_615_Psn	D_618_Psn	D_620_Psn	D_625_Psn	D_630_Psn	D_647_Psn	D_679_Psn	D_680_Psn	D_683_Psn	D_688_Psn	D_694_Psn	D_700_Psn	D_711_Psn	D_735_Psn	D_796_Psn	D_840_Psn	D_912_Psn	D_944_Psn	D_945_Psn	D_948_Psn	D_959_Psn	D_966_Psn	D_984_Psn	D_1029_Psn	D_1030_Psn	D_1032_Psn	D_1038_Psn	D_1045_Psn	Chr	Coordinate	QUERY_KEY	rsID	NimbleGen_Record	UNK_Nonsynonymous_Var

datrownames <- row.names(dat) <- dat$PROBE_ID
n <- length(dat[,1])
n
[1] 3386
dim(dat)
[1] 3386  117


############################################################
# B.2. P-VALUE DISTRIBUTION BY MAF.
############################################################
# Columns 1 - 9.
# CHR_COORD_REF_ALT	ESP_EUR_MAF	1000G_EUR_MAF	PROBE_ID	P_Value	FDR_P	Difference_MeanRef-MeanPsn	MeanRef	MeanPsn

# UNK.
dat_a <- dat[,c(1:9,10:60,61:111)]


############################################################
# B.2.1. DIVIDE DATA BY ESP MAF.
############################################################
############################################################
# B.2.1.1. DIVIDE DATA BY MAF (3 CATEGORIES).
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames

# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,1] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.01 & dat_data[,1] < 0.05],])
dat_data_c <-  as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.05],])

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("MAF >= 0.05", length(dat_data_c[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 6.5792, df = 2, p-value = 0.03727
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,3]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("MAF >= 0.05", length(dat_data_c[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 6.4392, df = 2, p-value = 0.03997
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,4]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# Mean Difference Distribution.
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF >= 0.05", length(dat_data_c[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 7.9571, df = 2, p-value = 0.01871
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "gray", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_c[,5]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "gray", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_c[,5])), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5]))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF >= 0.05", length(dat_data_c[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "gray", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_c[,5])))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")


############################################################
# B.2.1.2. DIVIDE DATA BY MAF (5 CATEGORIES).
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames

# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, 0.10), [0.10, 0.50), [0.50, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,1] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.01 & dat_data[,1] < 0.05],])
dat_data_c <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.05 & dat_data[,1] < 0.10],])
dat_data_d <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.10 & dat_data[,1] < 0.50],])
dat_data_e <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.50],])

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))), -log10(as.numeric(as.character(dat_data_d[,3]))), -log10(as.numeric(as.character(dat_data_e[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,3])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,3])), rep("MAF >= 0.50", length(dat_data_e[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 12.141, df = 4, p-value = 0.01633
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))), -log10(as.numeric(as.character(dat_data_d[,3]))), -log10(as.numeric(as.character(dat_data_e[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,3]))), at = 3, col = "turquoise", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_d[,3]))), at = 4, col = "red", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_e[,3]))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# P_Value Distribution (MAF < 0.10).
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 11.798, df = 2, p-value = 0.002743
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,3]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))), -log10(as.numeric(as.character(dat_data_d[,4]))), -log10(as.numeric(as.character(dat_data_e[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,4])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,4])), rep("MAF >= 0.50", length(dat_data_e[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 11.908, df = 4, p-value = 0.01805
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))), -log10(as.numeric(as.character(dat_data_d[,4]))), -log10(as.numeric(as.character(dat_data_e[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,4]))), at = 3, col = "turquoise", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_d[,4]))), at = 4, col = "red", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_e[,4]))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of FDR-Ajusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(FDR_P_Value)")

# FDR_P_Value Distribution (MAF < 0.10).
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 11.63, df = 2, p-value = 0.002982
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,4]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(FDR_P_Value)")

# Mean Difference Distribution.
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))), abs(as.numeric(as.character(dat_data_d[,5]))), abs(as.numeric(as.character(dat_data_e[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,5])), rep("MAF >= 0.50", length(dat_data_e[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 13.677, df = 4, p-value = 0.008401
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))), abs(as.numeric(as.character(dat_data_d[,5]))), abs(as.numeric(as.character(dat_data_e[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "gray", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_c[,5]))), at = 3, col = "turquoise", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_d[,5]))), at = 4, col = "red", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_e[,5]))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])), -as.numeric(as.character(dat_data_d[,5])), -as.numeric(as.character(dat_data_e[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "gray", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_c[,5])), at = 3, col = "turquoise", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_d[,5])), at = 4, col = "red", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_e[,5])), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])), -as.numeric(as.character(dat_data_d[,5])), -as.numeric(as.character(dat_data_e[,5]))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,5])), rep("MAF >= 0.50", length(dat_data_e[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))), abs(as.numeric(as.character(dat_data_d[,5]))), abs(as.numeric(as.character(dat_data_e[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "gray", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_c[,5])))), at = 3, col = "turquoise", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_d[,5])))), at = 4, col = "red", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_e[,5])))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")

# Mean Difference Distribution (MAF < 0.10).
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 13.022, df = 2, p-value = 0.001487
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "gray", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_c[,5]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "gray", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_c[,5])), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5]))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "gray", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_c[,5])))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")


############################################################
# B.2.1.3. DIVIDE DATA BY MAF (100 CATEGORIES).
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
k <- length(dat_data[,1])

# Divide Peptides by ESP MAF Percentage.
# P_Value / FDR_P_Value Distributions.
pv <- c()
fdrp <- c()
category <- c()
for (i in 1:100) {
	dat_data_a <- as.matrix(dat_data[(1:k)[dat_data[,1] >= (i-1)/100 & dat_data[,1] < i/100],])
	pv <- c(pv, -log10(as.numeric(as.character(dat_data_a[,3]))))
	fdrp <- c(fdrp, -log10(as.numeric(as.character(dat_data_a[,4]))))
	category <- c(category, rep(i, length(dat_data_a[,1])))
}
# P_Value.
dat_data_pv <- data.frame(LOG10_P_VALUE = pv, CATEGORY = category)
# FDR_P_Value.
dat_data_fdrp <- data.frame(LOG10_FDR_P = fdrp, CATEGORY = category)

# Mean Difference Distribution.
mean_diff <- c()
category <- c()
for (i in 1:100) {
	dat_data_a <- as.matrix(dat_data[(1:k)[dat_data[,1] >= (i-1)/100 & dat_data[,1] < i/100],])
	mean_diff <- c(mean_diff, abs(as.numeric(as.character(dat_data_a[,5]))))
	category <- c(category, rep(i, length(dat_data_a[,1])))
}
dat_data_mean_diff <- data.frame(ABS_MEAN_DIFF = mean_diff, CATEGORY = category)

# Box Plots.
# P_Value.
boxplot(dat_data_pv[,1]~dat_data_pv[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF"))
stripchart(as.numeric(dat_data_pv[,1])~dat_data_pv[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 
# FDR_P_Value.
boxplot(dat_data_fdrp[,1]~dat_data_fdrp[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adusted P_Value for Different MAF"))
stripchart(as.numeric(dat_data_fdrp[,1])~dat_data_fdrp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 
# Mean Difference.
# Plot Absolute Difference.
boxplot(dat_data_mean_diff[,1]~dat_data_mean_diff[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "ABSOLUTE MEAN DIFFERENCE (|REF-PSN|)", main = paste("Comparison of Mean Differences for Different MAF"))
stripchart(as.numeric(dat_data_mean_diff[,1])~dat_data_mean_diff[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
boxplot(log10(dat_data_mean_diff[,1])~dat_data_mean_diff[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|REF-PSN|))", main = paste("Comparison of Mean Differences for Different MAF"))
stripchart(log10(as.numeric(dat_data_mean_diff[,1]))~dat_data_mean_diff[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 

# Violin Plots.
# P_Value.
p <- ggplot(dat_data_pv, aes(x=as.factor(CATEGORY), y=LOG10_P_VALUE))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "-Log10(P_Value)")
# FDR_P_Value.
p <- ggplot(dat_data_fdrp, aes(x=as.factor(CATEGORY), y=LOG10_FDR_P))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "-Log10(FDR_P_Value)")
# Mean Difference.
# Plot Absolute Difference.
p <- ggplot(dat_data_mean_diff, aes(x=as.factor(CATEGORY), y=ABS_MEAN_DIFF))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "ABSOLUTE MEAN DIFFERENCE (|REF-PSN|)")
# Plot log10 Absolute Difference.
p <- ggplot(dat_data_mean_diff, aes(x=as.factor(CATEGORY), y=log10(ABS_MEAN_DIFF)))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "LOG10(ABSOLUTE MEAN DIFFERENCE (|REF-PSN|))")


############################################################
# B.2.1.4. COMPARE RARE VS COMMON MAP.
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames

# Divide Peptides by ESP MAF [0, 0.05), (0.05, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,1] < 0.05],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,1] > 0.3],])

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3])))), category = c(rep("MAF < 0.05", length(dat_data_a[,3])), rep("MAF > 0.30", length(dat_data_b[,3]))))
test <- wilcox.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Wilcoxon rank sum test with continuity correction

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
W = 475730, p-value = 0.2094
alternative hypothesis: true location shift is not equal to 0
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4])))), category = c(rep("MAF < 0.05", length(dat_data_a[,4])), rep("MAF > 0.30", length(dat_data_b[,4]))))
test <- wilcox.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Wilcoxon rank sum test with continuity correction

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
W = 475390, p-value = 0.2186
alternative hypothesis: true location shift is not equal to 0
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(FDR_P_Value)")

# Mean Difference Distribution.
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5])))), category = c(rep("MAF < 0.05", length(dat_data_a[,5])), rep("MAF > 0.30", length(dat_data_b[,5]))))
test <- wilcox.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Wilcoxon rank sum test with continuity correction

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
W = 477500, p-value = 0.1657
alternative hypothesis: true location shift is not equal to 0
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a[,5])), rep("MAF > 0.30", length(dat_data_b[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2. DIVIDE DATA BY 1000G MAF.
############################################################
############################################################
# B.2.2.1. DIVIDE DATA BY MAF (3 CATEGORIES).
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames

# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.01 & dat_data[,2] < 0.05],])
dat_data_c <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.05],])

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("MAF >= 0.05", length(dat_data_c[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 15.095, df = 2, p-value = 0.0005274
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,3]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("MAF >= 0.05", length(dat_data_c[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 14.942, df = 2, p-value = 0.0005692
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,4]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# Mean Difference Distribution.
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF >= 0.05", length(dat_data_c[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 13.749, df = 2, p-value = 0.001034
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "gray", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_c[,5]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "gray", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_c[,5])), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5]))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF >= 0.05", length(dat_data_c[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "gray", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_c[,5])))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","MAF >= 0.05"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")


############################################################
# B.2.2.2. DIVIDE DATA BY MAF (5 CATEGORIES).
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames

# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, 0.10), [0.10, 0.50), [0.50, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.01 & dat_data[,2] < 0.05],])
dat_data_c <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.05 & dat_data[,2] < 0.10],])
dat_data_d <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.10 & dat_data[,2] < 0.50],])
dat_data_e <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.50],])

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))), -log10(as.numeric(as.character(dat_data_d[,3]))), -log10(as.numeric(as.character(dat_data_e[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,3])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,3])), rep("MAF >= 0.50", length(dat_data_e[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 17.677, df = 4, p-value = 0.001427
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))), -log10(as.numeric(as.character(dat_data_d[,3]))), -log10(as.numeric(as.character(dat_data_e[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,3]))), at = 3, col = "turquoise", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_d[,3]))), at = 4, col = "red", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_e[,3]))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# P_Value Distribution (MAF < 0.10).
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 12.785, df = 2, p-value = 0.001674
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))), -log10(as.numeric(as.character(dat_data_c[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,3]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))), -log10(as.numeric(as.character(dat_data_d[,4]))), -log10(as.numeric(as.character(dat_data_e[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,4])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,4])), rep("MAF >= 0.50", length(dat_data_e[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 17.503, df = 4, p-value = 0.001543
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))), -log10(as.numeric(as.character(dat_data_d[,4]))), -log10(as.numeric(as.character(dat_data_e[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,4]))), at = 3, col = "turquoise", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_d[,4]))), at = 4, col = "red", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_e[,4]))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(FDR_P_Value)")

# FDR_P_Value Distribution (MAF < 0.10).
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 12.726, df = 2, p-value = 0.001724
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))), -log10(as.numeric(as.character(dat_data_c[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "gray", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_c[,4]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(FDR_P_Value)")

# Mean Difference Distribution.
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))), abs(as.numeric(as.character(dat_data_d[,5]))), abs(as.numeric(as.character(dat_data_e[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,5])), rep("MAF >= 0.50", length(dat_data_e[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 16.679, df = 4, p-value = 0.002231
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))), abs(as.numeric(as.character(dat_data_d[,5]))), abs(as.numeric(as.character(dat_data_e[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "gray", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_c[,5]))), at = 3, col = "turquoise", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_d[,5]))), at = 4, col = "red", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_e[,5]))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])), -as.numeric(as.character(dat_data_d[,5])), -as.numeric(as.character(dat_data_e[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "gray", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_c[,5])), at = 3, col = "turquoise", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_d[,5])), at = 4, col = "red", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_e[,5])), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])), -as.numeric(as.character(dat_data_d[,5])), -as.numeric(as.character(dat_data_e[,5]))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,5])), rep("MAF >= 0.50", length(dat_data_e[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,5.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))), abs(as.numeric(as.character(dat_data_d[,5]))), abs(as.numeric(as.character(dat_data_e[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "gray", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_c[,5])))), at = 3, col = "turquoise", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_d[,5])))), at = 4, col = "red", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_e[,5])))), at = 5, col = "darkblue", add=TRUE)
axis(side=1,at=1:5,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10","0.10 <= MAF < 0.50","MAF >= 0.50"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")

# Mean Difference Distribution (MAF < 0.10).
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 11.691, df = 2, p-value = 0.002893
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "gray", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_c[,5]))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "gray", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_c[,5])), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])), -as.numeric(as.character(dat_data_c[,5]))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,3.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))), abs(as.numeric(as.character(dat_data_c[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "gray", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_c[,5])))), at = 3, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("0 <= MAF < 0.01","0.01 <= MAF < 0.05","0.05 <= MAF < 0.10"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")


############################################################
# B.2.2.3. DIVIDE DATA BY MAF (100 CATEGORIES).
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
k <- length(dat_data[,1])

# Divide Peptides by TG MAF Percentage.
# P_Value / FDR_P_Value Distributions.
pv <- c()
fdrp <- c()
category <- c()
for (i in 1:100) {
	dat_data_a <- as.matrix(dat_data[(1:k)[dat_data[,2] >= (i-1)/100 & dat_data[,2] < i/100],])
	pv <- c(pv, -log10(as.numeric(as.character(dat_data_a[,3]))))
	fdrp <- c(fdrp, -log10(as.numeric(as.character(dat_data_a[,4]))))
	category <- c(category, rep(i, length(dat_data_a[,1])))
}
# P_Value.
dat_data_pv <- data.frame(LOG10_P_VALUE = pv, CATEGORY = category)
# FDR_P_Value.
dat_data_fdrp <- data.frame(LOG10_FDR_P = fdrp, CATEGORY = category)

# Mean Difference Distribution.
mean_diff <- c()
category <- c()
for (i in 1:100) {
	dat_data_a <- as.matrix(dat_data[(1:k)[dat_data[,2] >= (i-1)/100 & dat_data[,2] < i/100],])
	mean_diff <- c(mean_diff, abs(as.numeric(as.character(dat_data_a[,5]))))
	category <- c(category, rep(i, length(dat_data_a[,1])))
}
dat_data_mean_diff <- data.frame(ABS_MEAN_DIFF = mean_diff, CATEGORY = category)

# Box Plots.
# P_Value.
boxplot(dat_data_pv[,1]~dat_data_pv[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF"))
stripchart(as.numeric(dat_data_pv[,1])~dat_data_pv[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 
# FDR_P_Value.
boxplot(dat_data_fdrp[,1]~dat_data_fdrp[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adusted P_Value for Different MAF"))
stripchart(as.numeric(dat_data_fdrp[,1])~dat_data_fdrp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 
# Mean Difference.
# Plot Absolute Difference.
boxplot(dat_data_mean_diff[,1]~dat_data_mean_diff[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "ABSOLUTE MEAN DIFFERENCE (|REF-PSN|)", main = paste("Comparison of Mean Differences for Different MAF"))
stripchart(as.numeric(dat_data_mean_diff[,1])~dat_data_mean_diff[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
boxplot(log10(dat_data_mean_diff[,1])~dat_data_mean_diff[,2], col = rainbow(150), xlab = "MAF (%)", ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|REF-PSN|))", main = paste("Comparison of Mean Differences for Different MAF"))
stripchart(log10(as.numeric(dat_data_mean_diff[,1]))~dat_data_mean_diff[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c(rev(rainbow(150))[51:80], rev(rainbow(150))[11:50], rev(rainbow(150))[121:150]), bg = "bisque", add = TRUE) 

# Violin Plots.
# P_Value.
p <- ggplot(dat_data_pv, aes(x=as.factor(CATEGORY), y=LOG10_P_VALUE))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "-Log10(P_Value)")
# FDR_P_Value.
p <- ggplot(dat_data_fdrp, aes(x=as.factor(CATEGORY), y=LOG10_FDR_P))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "-Log10(FDR_P_Value)")
# Mean Difference.
# Plot Absolute Difference.
p <- ggplot(dat_data_mean_diff, aes(x=as.factor(CATEGORY), y=ABS_MEAN_DIFF))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "ABSOLUTE MEAN DIFFERENCE (|REF-PSN|)")
# Plot log10 Absolute Difference.
p <- ggplot(dat_data_mean_diff, aes(x=as.factor(CATEGORY), y=log10(ABS_MEAN_DIFF)))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:100], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "MAF (%)", y = "LOG10(ABSOLUTE MEAN DIFFERENCE (|REF-PSN|))")


############################################################
# B.2.2.4. COMPARE RARE VS COMMON MAP.
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	X1000G_EUR_MAF	P_Value	FDR_P	Difference_MeanRef.MeanPsn	MeanRef	MeanPsn
dat_data <- as.matrix(dat_a[,c(2:3, 5:9)])
dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames

# Divide Peptides by ESP MAF [0, 0.05), (0.05, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.3],])

# P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3])))), category = c(rep("MAF < 0.05", length(dat_data_a[,3])), rep("MAF > 0.30", length(dat_data_b[,3]))))
test <- wilcox.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Wilcoxon rank sum test with continuity correction

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
W = 391290, p-value = 0.008288
alternative hypothesis: true location shift is not equal to 0
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(P_Value)", main = paste("Comparison of P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,3]))), -log10(as.numeric(as.character(dat_data_b[,3]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,3]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,3]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(P_Value)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 

# FDR_P_Value Distribution.
datdatapforbp <- data.frame(datmean = c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4])))), category = c(rep("MAF < 0.05", length(dat_data_a[,4])), rep("MAF > 0.30", length(dat_data_b[,4]))))
test <- wilcox.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Wilcoxon rank sum test with continuity correction

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
W = 391040, p-value = 0.008804
alternative hypothesis: true location shift is not equal to 0
############
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "-Log10(FDR_P_Value)", main = paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-log10(as.numeric(as.character(dat_data_a[,4]))), -log10(as.numeric(as.character(dat_data_b[,4]))))), axes=FALSE,ann=FALSE)
vioplot(-log10(as.numeric(as.character(dat_data_a[,4]))), at = 1, col = "gold", add=TRUE)
vioplot(-log10(as.numeric(as.character(dat_data_b[,4]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of FDR-Adjusted P_Value for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="-Log10(FDR_P_Value)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 

# Mean Difference Distribution.
datdatapforbp <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5])))), category = c(rep("MAF < 0.05", length(dat_data_a[,5])), rep("MAF > 0.30", length(dat_data_b[,5]))))
test <- wilcox.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Wilcoxon rank sum test with continuity correction

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
W = 391770, p-value = 0.007359
alternative hypothesis: true location shift is not equal to 0
############
# Plot Absolute Difference.
# Box Plot.
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)", main = paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Absolute Diff.
# All MAF > 0.30 Plotted.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b[,5]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 Plotted.
dim(dat_data_a)
[1] 364   7
set.seed(1)
dat_data_b_row_num <- sample(c(1:length(dat_data_b[,1])), length(dat_data_a[,1]))
dat_data_b_sub <- dat_data_b[dat_data_b_row_num,]
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b_sub[,5]))))), axes=FALSE,ann=FALSE)
vioplot(abs(as.numeric(as.character(dat_data_a[,5]))), at = 1, col = "gold", add=TRUE)
vioplot(abs(as.numeric(as.character(dat_data_b_sub[,5]))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="ABSOLUTE MEAN DIFFERENCE (|PSN-REF|)")
datdata_sub <- data.frame(datmean = c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b_sub[,5])))), category = c(rep("MAF < 0.05", length(dat_data_a[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdata_sub[,1])~ datdata_sub[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot -- Original Diff.
# All MAF > 0.30 Plotted.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a[,5])), rep("MAF > 0.30", length(dat_data_b[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 Plotted.
dim(dat_data_a)
[1] 364   7
set.seed(1)
dat_data_b_row_num <- sample(c(1:length(dat_data_b[,1])), length(dat_data_a[,1]))
dat_data_b_sub <- dat_data_b[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# All MAF > 0.30 && FDR_P < 0.05 Plotted.
# quantile(dat_data[,4], probs = seq(0,1,0.01))
summary(dat_data[,4])
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000000 0.0000005 0.0048450 0.1805000 0.2620000 1.0000000 
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.05],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.05],])
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_fdrp[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_fdrp[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_fdrp[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_fdrp[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_fdrp[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 && FDR_P < 0.05 Plotted -- Iteration 1.
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.05],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.05],])
dim(dat_data_a_fdrp)
[1] 243   7
set.seed(1)
dat_data_b_row_num <- sample(c(1:length(dat_data_b_fdrp[,1])), length(dat_data_a_fdrp[,1]))
dat_data_b_sub <- dat_data_b_fdrp[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\nFDR-Adjusted P-Value < 0.05 (Iteration 1)\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 && FDR_P < 0.05 Plotted -- Iteration 2.
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.05],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.05],])
dim(dat_data_a_fdrp)
[1] 243   7
set.seed(2)
dat_data_b_row_num <- sample(c(1:length(dat_data_b_fdrp[,1])), length(dat_data_a_fdrp[,1]))
dat_data_b_sub <- dat_data_b_fdrp[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\nFDR-Adjusted P-Value < 0.05 (Iteration 2)\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 && FDR_P < 0.05 Plotted -- Iteration 3.
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.05],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.05],])
dim(dat_data_a_fdrp)
[1] 243   7
set.seed(5)
dat_data_b_row_num <- sample(c(1:length(dat_data_b_fdrp[,1])), length(dat_data_a_fdrp[,1]))
dat_data_b_sub <- dat_data_b_fdrp[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\nFDR-Adjusted P-Value < 0.05 (Iteration 3)\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 && FDR_P < 0.0000005 Plotted -- Iteration 1.
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.0000005],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.0000005],])
dim(dat_data_a_fdrp)
[1] 243   7
set.seed(1)
dat_data_b_row_num <- sample(c(1:length(dat_data_b_fdrp[,1])), length(dat_data_a_fdrp[,1]))
dat_data_b_sub <- dat_data_b_fdrp[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\nFDR-Adjusted P-Value < 0.0000005 (1st Qt., Iteration 1)\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 && FDR_P < 0.0000005 Plotted -- Iteration 2.
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.0000005],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.0000005],])
dim(dat_data_a_fdrp)
[1] 243   7
set.seed(2)
dat_data_b_row_num <- sample(c(1:length(dat_data_b_fdrp[,1])), length(dat_data_a_fdrp[,1]))
dat_data_b_sub <- dat_data_b_fdrp[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\nFDR-Adjusted P-Value < 0.0000005 (1st Qt., Iteration 2)\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Randomly Selected (Equal Length) MAF > 0.30 && FDR_P < 0.0000005 Plotted -- Iteration 3.
dat_data_a_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.05 & dat_data[,4] < 0.0000005],])
dat_data_b_fdrp <- as.matrix(dat_data[(1:n)[dat_data[,2] > 0.30 & dat_data[,4] < 0.0000005],])
dim(dat_data_a_fdrp)
[1] 243   7
set.seed(5)
dat_data_b_row_num <- sample(c(1:length(dat_data_b_fdrp[,1])), length(dat_data_a_fdrp[,1]))
dat_data_b_sub <- dat_data_b_fdrp[dat_data_b_row_num,]
test <- wilcox.test(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])), alternative = "two.sided", )
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5])))), axes=FALSE,ann=FALSE)
vioplot(-as.numeric(as.character(dat_data_a_fdrp[,5])), at = 1, col = "gold", add=TRUE)
vioplot(-as.numeric(as.character(dat_data_b_sub[,5])), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:2,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\nFDR-Adjusted P-Value < 0.0000005 (1st Qt., Iteration 3)\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="MEAN DIFFERENCE (PSN-REF)")
datdatapforvp <- data.frame(datmean = c(-as.numeric(as.character(dat_data_a_fdrp[,5])), -as.numeric(as.character(dat_data_b_sub[,5]))), category = c(rep("MAF < 0.05", length(dat_data_a_fdrp[,5])), rep("MAF > 0.30", length(dat_data_b_sub[,5]))))
stripchart(as.numeric(datdatapforvp[,1])~ datdatapforvp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Plot log10 Absolute Difference.
# Box Plot.
boxplot(log10(datdatapforbp[,1])~ datdatapforbp[,2], ylab = "LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))", main = paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "turquoise"))
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 
# Violin Plot.
par(las=1,bty="l")
plot(0:1,0:1,type="n",xlim=c(0.5,2.5),ylim=range(log10(c(abs(as.numeric(as.character(dat_data_a[,5]))), abs(as.numeric(as.character(dat_data_b[,5])))))), axes=FALSE,ann=FALSE)
vioplot(log10(abs(as.numeric(as.character(dat_data_a[,5])))), at = 1, col = "gold", add=TRUE)
vioplot(log10(abs(as.numeric(as.character(dat_data_b[,5])))), at = 2, col = "turquoise", add=TRUE)
axis(side=1,at=1:3,labels=c("MAF < 0.05","MAF > 0.30"))
axis(side=2)
title(paste("Comparison of Mean Differences for Different MAF\n(Mann-Whitney U Test P = ", sprintf("%1.2E", test$p.value), ")"), ylab="LOG10(ABSOLUTE MEAN DIFFERENCE (|PSN-REF|))")
stripchart(log10(as.numeric(datdatapforbp[,1]))~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "gold"), bg = "bisque", add = TRUE) 




