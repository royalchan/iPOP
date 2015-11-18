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

dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_All_UNK_Log2_QN_lm_id_refpsn_concat.txt",  header =TRUE, sep = "\t")

# SAMPLE Names.
# PROBE_ID_Ref	UNK-0_Ref	UNK-1_Ref	UNK-2_Ref	UNK-3_Ref	UNK-5_Ref	UNK-6_Ref	UNK-8_Ref	UNK-9_Ref	UNK-10_Ref	UNK-11_Ref	UNK-12_Ref	UNK-13_Ref	UNK-14_Ref	UNK-15_Ref	UNK-16_Ref	UNK-17_Ref	UNK-19_Ref	UNK-20_Ref	UNK-21_Ref	UNK-23_Ref	UNK-24_Ref	UNK-25_Ref	UNK-26_Ref	UNK-27_Ref	UNK-29_Ref	UNK-30_Ref	UNK-31_Ref	UNK-32_Ref	UNK-33_Ref	UNK-34_Ref	UNK-35_Ref	UNK-36_Ref	UNK-37_Ref	UNK-38_Ref	UNK-39_Ref	UNK-40_Ref	UNK-41_Ref	UNK-42_Ref	UNK-43_Ref	UNK-45_Ref	UNK-46_Ref	UNK-47_Ref	UNK-48_Ref	UNK-49_Ref	UNK-50_Ref	UNK-51_Ref	UNK-52_Ref	UNK-53_Ref	UNK-54_Ref	UNK-55_Ref	UNK-56_Ref	PROBE_ID_Psn	UNK-0_Psn	UNK-1_Psn	UNK-2_Psn	UNK-3_Psn	UNK-5_Psn	UNK-6_Psn	UNK-8_Psn	UNK-9_Psn	UNK-10_Psn	UNK-11_Psn	UNK-12_Psn	UNK-13_Psn	UNK-14_Psn	UNK-15_Psn	UNK-16_Psn	UNK-17_Psn	UNK-19_Psn	UNK-20_Psn	UNK-21_Psn	UNK-23_Psn	UNK-24_Psn	UNK-25_Psn	UNK-26_Psn	UNK-27_Psn	UNK-29_Psn	UNK-30_Psn	UNK-31_Psn	UNK-32_Psn	UNK-33_Psn	UNK-34_Psn	UNK-35_Psn	UNK-36_Psn	UNK-37_Psn	UNK-38_Psn	UNK-39_Psn	UNK-40_Psn	UNK-41_Psn	UNK-42_Psn	UNK-43_Psn	UNK-45_Psn	UNK-46_Psn	UNK-47_Psn	UNK-48_Psn	UNK-49_Psn	UNK-50_Psn	UNK-51_Psn	UNK-52_Psn	UNK-53_Psn	UNK-54_Psn	UNK-55_Psn	UNK-56_Psn
# PROBE_ID_Ref	D_-123_Ref	D_0_Ref	D_4_Ref	D_21_Ref	D_116_Ref	D_185_Ref	D_255_Ref	D_289_Ref	D_290_Ref	D_292_Ref	D_294_Ref	D_297_Ref	D_301_Ref	D_307_Ref	D_311_Ref	D_322_Ref	D_369_Ref	D_380_Ref	D_400_Ref	D_476_Ref	D_532_Ref	D_546_Ref	D_602_Ref	D_615_Ref	D_618_Ref	D_620_Ref	D_625_Ref	D_630_Ref	D_647_Ref	D_679_Ref	D_680_Ref	D_683_Ref	D_688_Ref	D_694_Ref	D_700_Ref	D_711_Ref	D_735_Ref	D_796_Ref	D_840_Ref	D_912_Ref	D_944_Ref	D_945_Ref	D_948_Ref	D_959_Ref	D_966_Ref	D_984_Ref	D_1029_Ref	D_1030_Ref	D_1032_Ref	D_1038_Ref	D_1045_Ref	PROBE_ID_Psn	D_-123_Psn	D_0_Psn	D_4_Psn	D_21_Psn	D_116_Psn	D_185_Psn	D_255_Psn	D_289_Psn	D_290_Psn	D_292_Psn	D_294_Psn	D_297_Psn	D_301_Psn	D_307_Psn	D_311_Psn	D_322_Psn	D_369_Psn	D_380_Psn	D_400_Psn	D_476_Psn	D_532_Psn	D_546_Psn	D_602_Psn	D_615_Psn	D_618_Psn	D_620_Psn	D_625_Psn	D_630_Psn	D_647_Psn	D_679_Psn	D_680_Psn	D_683_Psn	D_688_Psn	D_694_Psn	D_700_Psn	D_711_Psn	D_735_Psn	D_796_Psn	D_840_Psn	D_912_Psn	D_944_Psn	D_945_Psn	D_948_Psn	D_959_Psn	D_966_Psn	D_984_Psn	D_1029_Psn	D_1030_Psn	D_1032_Psn	D_1038_Psn	D_1045_Psn

row.names(dat) <- dat$PROBE_ID_Psn
datrownames <- row.names(dat)
# Reference peptides.
dat_data1 <- as.matrix(dat[,2:52])
row.names(dat_data1) <- datrownames
nrow1 <- length(dat_data1[,1])
ncol1 <- length(dat_data1[1,])
datcolnames1 <- colnames(dat)[2:52]
# Peptides with personal variants.
dat_data2 <- as.matrix(dat[,54:length(dat[1,])])
row.names(dat_data2) <- datrownames
nrow2 <- length(dat_data2[,1])
ncol2 <- length(dat_data2[1,])
datcolnames2 <- colnames(dat)[54:length(dat[1,])]

# Histogram.
hist(dat_data1, nclass = 200, col = rgb(0,0,1,1/2), main = "Histogram Comparison of Reactivity to Reference vs Personal Peptides", ylim = c(0, 7000))
hist(dat_data2, nclass = 200, col = rgb(0,1,0,1/4), add = TRUE)
legendtext <- c("Reference Peptides", "Personal Peptides")
legend("topright", legendtext, fill = c(rgb(0,0,1,1/2), rgb(0,1,0,1/4)), cex =1.2)

# Density Plot.
allsignals <- c(as.vector(dat_data1), as.vector(dat_data2))
dendat <- data.frame(INTENSITY = allsignals, CATEGORIES = c(rep("Reference Peptides", length(as.vector(dat_data1))), rep("Personal Peptides", length(as.vector(dat_data2)))))
ggplot(dendat, aes(x = INTENSITY, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle("Comparison of Reactivity to Reference vs Personal Peptides")


############################################################
# B.2. Reactivity Comparison between Reference Peptides and Personal Peptides.
############################################################
# Look for Significant Differences.
# Mean Differences.
dat_data1_mean <- apply(dat_data1, 1, mean)
dat_data2_mean <- apply(dat_data2, 1, mean)
test <- wilcox.test(dat_data1_mean, dat_data2_mean, paired = TRUE, alternative = "two.sided")
######OUTPUT_START######
test
	Wilcoxon signed rank test with continuity correction

data:  dat_data1_mean and dat_data2_mean
V = 2940838, p-value = 0.1949
alternative hypothesis: true location shift is not equal to 0
######OUTPUT_END######
# Box Plot.
dat12meanforbp <- data.frame(datmean = c(dat_data1_mean, dat_data2_mean), category = c(rep("Reference Peptide", length(dat_data1_mean)), rep("Personal Peptide", length(dat_data2_mean))))
boxplot(dat12meanforbp[,1]~dat12meanforbp[,2], ylab = "Log2(Reactivity)", main = "Comparison of Mean Reactivity to Reference vs Personal Peptides, p = 0.19", col = c("gold", "turquoise"))
# Add Lines to Box Plot.
for (i in 1:length(dat_data1_mean)) {
	refpsn <- c(dat_data2_mean[i], dat_data1_mean[i])
	lines(refpsn, type = "o", pch = ".", col = i, lwd = 0.1, add = TRUE)
}

# Per Peptide Differences.
n <- length(dat_data1[,1])
pv <- c()
diff12 <- c()
mean1 <- c()
mean2 <- c()
for (i in 1:n)
   {
    	t1 <- c(dat_data1[i,])
    	t2 <- c(dat_data2[i,])
    	mean1[i] <- mean(t1)
    	mean2[i] <- mean(t2)
    	diff12[i] <- mean1[i] - mean2[i]
		test <- wilcox.test(t1, t2, paired = TRUE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   
summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000000 0.0000001 0.0024270 0.1578000 0.1966000 1.0000000 
# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000000 0.0000005 0.0048450 0.1805000 0.2620000 1.0000000 
# Export pv and new_p for all peptides.
writepeptide <- cbind(pv, new_p, diff12, mean1, mean2, dat_data1, dat_data2)
write(c("PROBE_ID","P_Value","FDR_P","Difference_MeanRef-MeanPsn","MeanRef","MeanPsn",colnames(writepeptide)[6:length(colnames(writepeptide))]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_Reactivity_All_ts.txt", sep = "\t", ncolumn = 108)
write.table(writepeptide, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_Reactivity_All_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export peptides with one of the means >= 8.13612 (Shapiro Wilk Cutoff).
writepeptide <- cbind(diff12, mean1, mean2, dat_data1, dat_data2)
writepeptidesw <- writepeptide[(1:n)[writepeptide[,2] >= 8.13612 | writepeptide[,3] >= 8.13612],]
write(c("PROBE_ID","Difference_MeanRef-MeanPsn","MeanRef","MeanPsn",colnames(writepeptidesw)[4:length(colnames(writepeptidesw))]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_ReactivityGrEqSW_All_ts.txt", sep = "\t", ncolumn = 106)
write.table(writepeptidesw, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_ReactivityGrEqSW_All_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 2156
length((1:n)[pv < 0.005])
[1] 1773
length((1:n)[new_p < 0.10])
[1] 2235
length((1:n)[new_p < 0.05])
[1] 2047
# Export hits with FDR-adjusted p-value < 0.05.
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], diff12[(1:n)[new_p < 0.05]], mean1[(1:n)[new_p < 0.05]], mean2[(1:n)[new_p < 0.05]], dat_data1[(1:n)[new_p < 0.05],], dat_data2[(1:n)[new_p < 0.05],])
write(c("PROBE_ID","FDR_p-value","Difference_MeanRef-MeanPsn","MeanRef","MeanPsn",colnames(writecontentall)[5:length(colnames(writecontentall))]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_Reactivity_FDRp0_05_ts.txt", sep = "\t", ncolumn = 107)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_Reactivity_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,5:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,5:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(5,7.5,length=167), seq(7.5,8.5,length=167), seq(8.5,15,length=167)), keysize = 0.5, cexRow = 0.5)
total_heatmap <- heatmap.2(writecontentall[,5:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "row", symkey = TRUE, keysize = 0.5, cexRow = 0.5)

# Export hits with FDR-adjusted p-value < 0.05 AND one of the means >= 8.13612 (Shapiro Wilk Cutoff).
m <- length(writecontentall[,1])
writecontentall9 <- writecontentall[(1:m)[writecontentall[,3] >= 8.13612 | writecontentall[,4] >= 8.13612],]
write(c("PROBE_ID","FDR_p-value","Difference_MeanRef-MeanPsn","MeanRef","MeanPsn",colnames(writecontentall)[5:length(colnames(writecontentall))]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_ReactivityGrEqSW_FDRp0_05_ts.txt", sep = "\t", ncolumn = 107)
write.table(writecontentall9, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150403-Pass19_PersonalizedPeptides/Nimblegen_UNK_RefvsPsnPeptides_ReactivityGrEqSW_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall9[,5:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall9[,5:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(5,7.5,length=167), seq(7.5,8.5,length=167), seq(8.5,15,length=167)), keysize = 0.5, cexRow = 0.5)
total_heatmap <- heatmap.2(writecontentall9[,5:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "row", symkey = TRUE, keysize = 0.5, cexRow = 0.5)




