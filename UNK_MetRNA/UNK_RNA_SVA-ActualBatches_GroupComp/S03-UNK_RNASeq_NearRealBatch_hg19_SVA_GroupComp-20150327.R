############################################################
# TO-DO LIST.
############################################################
DONE. 1. Analyze SVA Adjusted UNK RNA-Seq data (group comparison + autocorrelated).





##################################################################################
##################################################################################
##################################################################################
### QN U-TEST SUMMARY:                                                         ###
###                                                                            ###
### HRV1 vs HEALTHY: FDR_P < 0.10; 171 Hits.                                   ###
### HRV2 vs HEALTHY: P < 0.005; 114 Hits.                                      ###
### HRV3 vs HEALTHY: FDR_P < 0.05; 153 Hits.                                   ###
### RSV  vs HEALTHY: FDR_P < 0.05; 2388 Hits.                                  ###
### ADV1 vs HEALTHY: FDR_P < 0.05; 698 Hits.                                   ###
### ADV2 vs HEALTHY: FDR_P < 0.05; 536 Hits.                                   ###
### HRVALL vs HEALTHY: FDR_P < 0.05; 606 Hits.                                 ###
### ADVALL vs HEALTHY: FDR_P < 0.05; 1056 Hits.                                ###
### INFALL vs HEALTHY: FDR_P < 0.05; 1223 Hits.                                ###
### LSP (Cluster 1): FDR_P < 0.05; 123 Hits.                                   ###
### LSP (Cluster 2): FDR_P < 0.05; 165 Hits.                                   ###
### AUTOCORRELATED (Cluster 1): 1.645 sigma (One-Tailed); 1422 Hits.           ###
### AUTOCORRELATED (Cluster 2): 1.645 sigma (One-Tailed); 633 Hits.            ###
##################################################################################
##################################################################################
##################################################################################





############################################################
# A. LIBRARIES USED
############################################################
library("gplots")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library("MASS")
library("geneplotter")
library("scatterplot3d")
library("affy")
library("ggplot2")
# Package A2R:
source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
# colored dendrogram op = par(bg = "#EFEFEF") A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("#FF6B6B","#4ECDC4", "#556270"))
# another colored dendrogram
op = par(bg = "gray15") cols = hsv(c(0.2, 0.57, 0.95), 1, 1, 0.8) A2Rplot(hc, k = 3, boxes = FALSE, col.up = "gray50", col.down = cols)





############################################################
# B. LOAD DATA AND INITIAL PROCESSING
############################################################
# Load stitched data with linear fitting.

# Entrez Gene Symbol used.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_NearRealBatch_ComBat_Rescaled_no_ac.txt",  header =TRUE, sep = "\t")


# Colume Names.
# Column 1: GeneID;
# Columns 2-58: UNK (1-60);
# Columns 59-60: UNK-07 & 08 RiboZero;
# GeneID	U01RNA	U02RNA	U03RNA	U05RNA	U06RNA	U07RNA	U08RNA	U09RNA	U10RNA	U11RNA	U12RNA	U13RNA	U14RNA	U15RNA	U16RNA	U17RNA	U18RNA	U19RNA	U20RNA	U21RNA	U23RNA	U24RNA	U25RNA	U26RNA	U27RNA	U28RNA	U29RNA	U30RNA	U31RNA	U32RNA	U33RNA	U34RNA	U35RNA	U36RNA	U37RNA	U38RNA	U39RNA	U40RNA	U41RNA	U42RNA	U43RNA	U45RNA	U46RNA	U47RNA	U48RNA	U49RNA	U50RNA	U51RNA	U52RNA	U53RNA	U54RNA	U55RNA	U56RNA	U58RNA	U58_5RNA	U59RNA	U60RNA	U07RNARZ	U08RNARZ

# GeneID	D_0	D_4	D_21	D_116	D_185	D_186	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_329	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045	D_1051	D_1060	D_1109	D_1124	D_186RZ	D_255RZ

dat_data <- as.matrix(dat[,c(2:58)])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
row.names(dat_data) <- dat$GeneID
# Alternative Column Names (by Day):
colID <- c("D_0", "D_4", "D_21", "D_116", "D_185", "D_186", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_329", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
colnames(dat_data) <- colID


############################################################
# C. U TEST
############################################################
############################################################
# AS THE DATA WERE ALREADY RESCALED SO THIS WAS NOT DONE.
# FOR EASY ADAPTATION THE VARIABLE NAME WAS UNCHANGED.
############################################################
# Input this time is already quantile normalized and rescaled.
total_matrix_unk_qn <- dat_data

############################################################
# SAME BATCH:
# 	D: 0 - 602;
# 	D: 615 - 1124;

# GROUP DETERMINATION:
# 	HRV1 (Days 0-21; 3 Columns): Columns 1-3;
# 	HRV2 (Days 615-630; 6 Columns): Columns 25-30;
# 	HRV3 (Days 1029-1060; 7 Columns): Columns 49-55;
# 	RSV (Days 289-311; 8 Columns): Columns 8-15;
# 	ADV1 (Days 679-700; 6 Columns): Columns 32-37;
# 	ADV2 (Days 944-984; 6 Columns): Columns 43-48;
# 	HLT (Days 116-255, 322-602, 647, 711-912, 1109-1124; 21 Columns): Columns 4-7, 16-24, 31, 38-42, 56-57.
############################################################
total_matrix_unk_qn_hrv1 <- total_matrix_unk_qn[,1:3]
total_matrix_unk_qn_hrv2 <- total_matrix_unk_qn[,25:30]
total_matrix_unk_qn_hrv3 <- total_matrix_unk_qn[,49:55]
total_matrix_unk_qn_rsv <- total_matrix_unk_qn[,8:15]
total_matrix_unk_qn_adv1 <- total_matrix_unk_qn[,32:37]
total_matrix_unk_qn_adv2 <- total_matrix_unk_qn[,43:48]
total_matrix_unk_qn_hlt <- total_matrix_unk_qn[,c(4:7,16:24,31,38:42,56:57)]

############################################################
# Mann-Whitney U Test:
# HRV1 vs HEALTHY.
# HRV2 vs HEALTHY.
# HRV3 vs HEALTHY.
# RSV vs HEALTHY.
# ADV1 vs HEALTHY.
# ADV2 vs HEALTHY.
############################################################
n <- length(total_matrix_unk_qn[,1])

############################################################
# B.2.1. HRV1 vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_hrv1[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0009881 0.0971500 0.3103000 0.3831000 0.6196000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
0.0009881423 0.0009881423 0.0019762846 0.0039525692 0.0069169960 0.0069169960 0.0108695652 0.0158102767 0.0158102767 0.0227272727 0.0227272727 0.0232100355 0.0306324111 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
0.0306324111 0.0405138340 0.0405138340 0.0446389171 0.0523715415 0.0523715415 0.0662055336 0.0662055336 0.0662055336 0.0807896235 0.0820158103 0.0820158103 0.0971527508 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
0.1007905138 0.1007905138 0.1160641222 0.1215415020 0.1215415020 0.1377547183 0.1452569170 0.1452569170 0.1624441848 0.1719367589 0.1719367589 0.1903338738 0.2015810277 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
0.2015810277 0.2033609929 0.2341897233 0.2341897233 0.2341897233 0.2563837860 0.2707509881 0.2707509881 0.2707509881 0.3102766798 0.3102766798 0.3102766798 0.3153686920 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
0.3537549407 0.3537549407 0.3537549407 0.3592922939 0.4011857708 0.4011857708 0.4011857708 0.4068734016 0.4515810277 0.4515810277 0.4515810277 0.4580254346 0.5049407115 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
0.5049407115 0.5049407115 0.5049407115 0.5612648221 0.5612648221 0.5612648221 0.5612648221 0.6003920370 0.6195652174 0.6195652174 0.6195652174 0.6624517383 0.6798418972 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
0.6798418972 0.6798418972 0.7269207361 0.7420948617 0.7420948617 0.7420948617 0.7598803137 0.8053359684 0.8053359684 0.8053359684 0.8613699736 0.8695652174 0.8695652174 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
0.8695652174 0.9193187654 0.9347826087 0.9347826087 0.9347826087 0.9347826087 1.0000000000 1.0000000000 1.0000000000 1.0000000000 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05801 0.36470 0.60970 0.59190 0.82010 1.00000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1183
length((1:n)[pv < 0.005])
[1] 262
length((1:n)[new_p < 0.10])
[1] 171
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with FDR-Adjusted p-value < 0.10.
hrv1vshlt <- cbind(total_matrix_unk_qn_hrv1[(1:n)[new_p < 0.10],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.10],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.10]],total_matrix_unk_qn_hrv1[(1:n)[new_p < 0.10],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.10],])
write(c("GeneID","FDR_p-value",colnames(hrv1vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10.txt", sep = "\t", ncolumn = 26)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with FDR-Adjusted p-value < 0.1 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.10]], total_matrix_unk_qn[(1:n)[new_p < 0.10],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_hrv1[(1:n)[new_p < 0.10],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.10],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,9), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV1vsHLT_FDRp0_10_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,10), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# B.2.2. HRV2 vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_hrv2[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000473 0.1752000 0.4412000 0.4533000 0.7119000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
0.0000472957 0.0033647512 0.0067295024 0.0124590385 0.0150940847 0.0217154826 0.0258031823 0.0305124827 0.0358636533 0.0419580420 0.0488159184 0.0565589000 0.0651869869 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
0.0748285531 0.0852059278 0.0854971116 0.0973075234 0.1085507022 0.1102665450 0.1245093071 0.1245093071 0.1400087835 0.1447710426 0.1569001047 0.1569001047 0.1751630012 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
0.1751630012 0.1949123340 0.1949123340 0.2161075639 0.2161075639 0.2388703084 0.2388703084 0.2552832004 0.2631194892 0.2631194892 0.2889699672 0.2889699672 0.3069242046 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
0.3163406642 0.3163406642 0.3453059018 0.3453059018 0.3453059018 0.3757575758 0.3757575758 0.3816003962 0.4077835208 0.4077835208 0.4307423919 0.4412283369 0.4412283369 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
0.4412283369 0.4761595892 0.4761595892 0.4761595892 0.5124218776 0.5124218776 0.5124218776 0.5402295390 0.5500422283 0.5500422283 0.5594436085 0.5888449715 0.5888449715 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
0.5888449715 0.6200314710 0.6288638897 0.6288638897 0.6407059714 0.6698760177 0.6698760177 0.6698760177 0.6698760177 0.7045797682 0.7119016249 0.7119016249 0.7547245026 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
0.7547245026 0.7547245026 0.7547245026 0.7983243809 0.7983243809 0.7983243809 0.8382321357 0.8424715381 0.8424715381 0.8610452806 0.8871659741 0.8871659741 0.8871659741 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
0.9302779551 0.9321374278 0.9321374278 0.9321374278 0.9773723861 0.9773723861 0.9773723861 1.0000000000 1.0000000000 1.0000000000 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1476  0.6708  0.8478  0.7848  0.9313  1.0000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 729
length((1:n)[pv < 0.005])
[1] 114
length((1:n)[new_p < 0.10])
[1] 0
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
hrv2vshlt <- cbind(total_matrix_unk_qn_hrv2[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_hrv2[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("GeneID","p-value",colnames(hrv2vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005.txt", sep = "\t", ncolumn = 29)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],1:57])
write(c("GeneID","p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_hrv2[(1:n)[pv < 0.005],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV2vsHLT_p0_005_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.7)


############################################################
# B.2.3. HRV3 vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_hrv3[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000017 0.1009000 0.3478000 0.4051000 0.6781000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
1.689132e-06 3.057329e-04 9.526705e-04 1.988109e-03 3.831528e-03 4.918205e-03 6.998074e-03 1.012297e-02 1.207561e-02 1.692679e-02 1.989967e-02 2.328300e-02 2.712577e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
3.146515e-02 3.635181e-02 4.182629e-02 4.794602e-02 5.474984e-02 6.229688e-02 6.328163e-02 7.063106e-02 7.980980e-02 8.950431e-02 8.987534e-02 1.008885e-01 1.008885e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
1.128847e-01 1.259214e-01 1.259214e-01 1.400358e-01 1.518530e-01 1.552752e-01 1.716699e-01 1.716699e-01 1.892656e-01 1.892656e-01 2.080808e-01 2.080808e-01 2.281511e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
2.281511e-01 2.494899e-01 2.494899e-01 2.721175e-01 2.721175e-01 2.960356e-01 2.960356e-01 3.133592e-01 3.212577e-01 3.212577e-01 3.477687e-01 3.477687e-01 3.755667e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
3.755667e-01 3.812798e-01 4.046282e-01 4.046282e-01 4.349313e-01 4.349313e-01 4.569024e-01 4.664386e-01 4.664386e-01 4.903074e-01 4.991200e-01 5.050759e-01 5.329178e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
5.329178e-01 5.417126e-01 5.677849e-01 5.677849e-01 6.036570e-01 6.036570e-01 6.329580e-01 6.404665e-01 6.404665e-01 6.781376e-01 6.781376e-01 6.906501e-01 7.165957e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
7.165957e-01 7.165957e-01 7.557481e-01 7.557481e-01 7.704037e-01 7.955086e-01 7.955086e-01 8.318875e-01 8.357826e-01 8.357826e-01 8.735198e-01 8.764704e-01 8.764704e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.154823e-01 9.174707e-01 9.174707e-01 9.365597e-01 9.586838e-01 9.586838e-01 9.788334e-01 1.000000e+00 1.000000e+00 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.003074 0.402000 0.691900 0.635600 0.893400 1.000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1225
length((1:n)[pv < 0.005])
[1] 367
length((1:n)[new_p < 0.10])
[1] 367
length((1:n)[new_p < 0.05])
[1] 153
# Export hits with FDR-Adjusted p-value < 0.05.
hrv3vshlt <- cbind(total_matrix_unk_qn_hrv3[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_hrv3[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(hrv3vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 30)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_hrv3[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRV3vsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.7)


############################################################
# B.2.4. RSV vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_rsv[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000005 0.0058920 0.0927000 0.2487000 0.4286000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
4.659675e-07 4.659675e-07 1.863870e-06 8.853382e-06 2.096854e-05 4.473288e-05 6.337158e-05 1.202196e-04 1.616907e-04 2.779704e-04 3.635526e-04 4.776167e-04 6.118153e-04 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
7.758359e-04 8.291690e-04 9.824418e-04 1.222699e-03 1.519520e-03 1.874587e-03 2.108813e-03 2.473556e-03 2.805590e-03 3.404358e-03 4.107037e-03 4.932266e-03 5.891693e-03 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
7.008151e-03 8.295619e-03 8.407416e-03 9.780657e-03 1.148004e-02 1.342359e-02 1.563135e-02 1.794024e-02 1.813778e-02 2.096434e-02 2.414970e-02 2.771854e-02 3.171281e-02 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.374084e-02 3.615861e-02 4.110206e-02 4.657019e-02 5.261192e-02 5.925569e-02 5.925569e-02 6.655134e-02 7.452731e-02 8.319596e-02 8.762976e-02 9.270237e-02 1.029793e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
1.127318e-01 1.140912e-01 1.260880e-01 1.389902e-01 1.431856e-01 1.528462e-01 1.676709e-01 1.835092e-01 1.876284e-01 2.003697e-01 2.182909e-01 2.372744e-01 2.372744e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
2.573534e-01 2.785190e-01 2.785190e-01 3.007965e-01 3.241689e-01 3.411887e-01 3.538127e-01 3.742157e-01 4.008690e-01 4.285713e-01 4.285713e-01 4.573196e-01 4.870632e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
4.870632e-01 5.177882e-01 5.494330e-01 5.819752e-01 6.153408e-01 6.153408e-01 6.494981e-01 6.494981e-01 6.843655e-01 6.961896e-01 7.199016e-01 7.560155e-01 7.560155e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
7.926601e-01 8.297366e-01 8.297366e-01 8.671930e-01 9.049256e-01 9.416386e-01 9.428768e-01 9.809408e-01 9.809408e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000365 0.0232200 0.1843000 0.3161000 0.5699000 1.0000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 3114
length((1:n)[pv < 0.005])
[1] 1765
length((1:n)[new_p < 0.10])
[1] 3015
length((1:n)[new_p < 0.05])
[1] 2388
# Export hits with FDR-Adjusted p-value < 0.05.
rsvvshlt <- cbind(total_matrix_unk_qn_rsv[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_rsv[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(rsvvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 31)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_rsv[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_RSVvsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# B.2.5. ADV1 vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_adv1[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000068 0.0419600 0.2389000 0.3332000 0.5888000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
6.756528e-06 4.729570e-05 2.026959e-04 4.324178e-04 8.445661e-04 1.141853e-03 1.533732e-03 2.416548e-03 3.268065e-03 4.283639e-03 5.384953e-03 6.729502e-03 8.317287e-03 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
1.022263e-02 1.245904e-02 1.313007e-02 1.509408e-02 1.814804e-02 2.171548e-02 2.291352e-02 2.580318e-02 3.051248e-02 3.091126e-02 3.586365e-02 4.195804e-02 4.195804e-02 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
4.881592e-02 5.500709e-02 5.655890e-02 6.518699e-02 6.518699e-02 7.482855e-02 8.013512e-02 8.549711e-02 9.643110e-02 9.730752e-02 1.102665e-01 1.102665e-01 1.245093e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
1.245093e-01 1.400088e-01 1.475805e-01 1.569001e-01 1.704507e-01 1.751630e-01 1.893761e-01 1.949123e-01 1.949123e-01 2.161076e-01 2.161076e-01 2.388703e-01 2.388703e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
2.631195e-01 2.631195e-01 2.889700e-01 2.889700e-01 3.073603e-01 3.163407e-01 3.213871e-01 3.453059e-01 3.453059e-01 3.757576e-01 3.757576e-01 4.054579e-01 4.077835e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
4.141449e-01 4.412283e-01 4.482074e-01 4.761596e-01 4.761596e-01 5.124219e-01 5.124219e-01 5.500422e-01 5.500422e-01 5.534807e-01 5.888450e-01 5.888450e-01 6.288639e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
6.288639e-01 6.407570e-01 6.698760e-01 6.698760e-01 7.119016e-01 7.119016e-01 7.119016e-01 7.547245e-01 7.547245e-01 7.983244e-01 7.983244e-01 8.424715e-01 8.424715e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
8.590022e-01 8.871660e-01 8.871660e-01 9.321374e-01 9.321374e-01 9.321374e-01 9.773724e-01 9.773724e-01 1.000000e+00 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.001587 0.167200 0.466200 0.474500 0.769300 1.000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1946
length((1:n)[pv < 0.005])
[1] 698
length((1:n)[new_p < 0.10])
[1] 1192
length((1:n)[new_p < 0.05])
[1] 698
# Export hits with FDR-Adjusted p-value < 0.05.
adv1vshlt <- cbind(total_matrix_unk_qn_adv1[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_adv1[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(adv1vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 29)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_adv1[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV1vsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# B.2.6. ADV2 vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_adv2[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000068 0.0473300 0.2161000 0.3250000 0.5500000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
6.756528e-06 8.107834e-05 2.972873e-04 6.080876e-04 1.141853e-03 1.533732e-03 2.621533e-03 3.364751e-03 4.283639e-03 5.384953e-03 6.729502e-03 8.317287e-03 1.022263e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
1.245904e-02 1.317407e-02 1.509408e-02 1.814804e-02 2.171548e-02 2.171548e-02 2.580318e-02 3.051248e-02 3.051248e-02 3.586365e-02 3.586365e-02 4.195804e-02 4.733408e-02 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
4.881592e-02 5.069437e-02 5.655890e-02 5.799604e-02 6.518699e-02 7.482855e-02 7.482855e-02 8.520307e-02 8.549711e-02 9.643110e-02 9.730752e-02 1.086970e-01 1.102665e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
1.245093e-01 1.245093e-01 1.400088e-01 1.400088e-01 1.569001e-01 1.569001e-01 1.751630e-01 1.751630e-01 1.949123e-01 1.949123e-01 2.161076e-01 2.161076e-01 2.388703e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
2.388703e-01 2.553558e-01 2.631195e-01 2.805402e-01 2.889700e-01 2.889700e-01 3.163407e-01 3.163407e-01 3.453059e-01 3.453059e-01 3.757576e-01 3.757576e-01 3.816004e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
4.077835e-01 4.077835e-01 4.412283e-01 4.412283e-01 4.761596e-01 4.761596e-01 5.124219e-01 5.124219e-01 5.334226e-01 5.500422e-01 5.500422e-01 5.888450e-01 5.888450e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
6.200315e-01 6.288639e-01 6.407570e-01 6.698760e-01 7.045703e-01 7.119016e-01 7.119016e-01 7.547245e-01 7.547245e-01 7.983244e-01 7.983244e-01 8.382321e-01 8.424715e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
8.424715e-01 8.871660e-01 8.871660e-01 9.302780e-01 9.321374e-01 9.321374e-01 9.773724e-01 9.773724e-01 1.000000e+00 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.004099 0.182000 0.427900 0.462600 0.729400 1.000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1953
length((1:n)[pv < 0.005])
[1] 608
length((1:n)[new_p < 0.10])
[1] 1112
length((1:n)[new_p < 0.05])
[1] 536
# Export hits with FDR-Adjusted p-value < 0.05.
adv2vshlt <- cbind(total_matrix_unk_qn_adv2[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_adv2[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(adv2vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 29)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_adv2[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADV2vsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,9), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# B.3. Grouped Comparison (Group All HRVs together, All ADVs together, and all infections together).
############################################################
# HRV All: 16 Columns.
total_matrix_unk_qn_hrvall <- total_matrix_unk_qn[,c(1:3,25:30,49:55)]
# ADV All: 12 Columns.
total_matrix_unk_qn_advall <- total_matrix_unk_qn[,c(32:37,43:48)]
# Infection All: 36 Columns.
total_matrix_unk_qn_infectionall <- total_matrix_unk_qn[,c(1:3,8:15,25:30,32:37,43:48,49:55)]
# Healthy All: 21 Columns.
total_matrix_unk_qn_hlt <- total_matrix_unk_qn[,c(4:7,16:24,31,38:42,56:57)]


############################################################
# B.3.1. HRVALL vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_hrvall[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.04795 0.24110 0.33600 0.57500 1.00000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
3.028944e-08 9.225946e-05 2.652145e-04 5.018407e-04 9.291382e-04 1.468648e-03 2.283876e-03 3.022857e-03 4.091476e-03 5.057251e-03 6.359767e-03 7.633252e-03 9.250516e-03 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
1.119828e-02 1.348981e-02 1.477913e-02 1.767650e-02 1.929797e-02 2.233505e-02 2.493936e-02 2.814692e-02 3.190866e-02 3.456513e-02 3.740276e-02 4.365745e-02 4.795300e-02 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
5.155478e-02 5.732237e-02 6.311280e-02 6.812072e-02 7.536946e-02 8.052391e-02 8.678840e-02 9.474069e-02 1.014359e-01 1.081384e-01 1.151771e-01 1.225614e-01 1.303004e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
1.384029e-01 1.468778e-01 1.557331e-01 1.649770e-01 1.746168e-01 1.846599e-01 1.951127e-01 1.978612e-01 2.087265e-01 2.172719e-01 2.289889e-01 2.411369e-01 2.537198e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
2.667407e-01 2.763922e-01 2.829009e-01 2.941057e-01 3.084525e-01 3.232425e-01 3.384754e-01 3.541495e-01 3.657652e-01 3.702628e-01 3.868120e-01 4.037933e-01 4.125760e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
4.212017e-01 4.390317e-01 4.572765e-01 4.759287e-01 4.949798e-01 5.097906e-01 5.196828e-01 5.342413e-01 5.544306e-01 5.648302e-01 5.749765e-01 5.958667e-01 6.170874e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
6.170874e-01 6.386246e-01 6.604630e-01 6.825872e-01 7.049805e-01 7.244063e-01 7.359310e-01 7.505054e-01 7.736011e-01 7.944012e-01 8.062472e-01 8.203642e-01 8.439925e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
8.660851e-01 8.677586e-01 8.916416e-01 9.156210e-01 9.267168e-01 9.396755e-01 9.637838e-01 9.879243e-01 1.000000e+00 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0002205 0.1908000 0.4779000 0.4866000 0.7652000 1.0000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1840
length((1:n)[pv < 0.005])
[1] 641
length((1:n)[new_p < 0.10])
[1] 1106
length((1:n)[new_p < 0.05])
[1] 606
# Export hits with FDR-Adjusted p-value < 0.05.
hrvallvshlt <- cbind(total_matrix_unk_qn_hrvall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_hrvall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(hrvallvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 39)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_hrvall[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_HRVALLvsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# B.3.2. ADVALL vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_advall[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.03302 0.22760 0.32560 0.56720 1.00000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
5.636703e-09 1.462161e-05 5.886409e-05 1.637011e-04 2.879228e-04 4.901903e-04 8.084695e-04 9.524394e-04 1.374159e-03 2.047104e-03 2.729202e-03 3.303006e-03 4.124934e-03 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
5.369907e-03 6.674483e-03 7.876224e-03 9.996023e-03 1.124863e-02 1.263341e-02 1.577614e-02 1.769062e-02 2.001435e-02 2.241996e-02 2.700763e-02 2.988835e-02 3.301920e-02 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
3.954261e-02 4.407537e-02 4.837037e-02 5.299933e-02 5.797914e-02 6.332839e-02 6.906425e-02 7.520570e-02 8.176988e-02 8.857986e-02 9.581109e-02 1.004323e-01 1.076175e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
1.203376e-01 1.249031e-01 1.393077e-01 1.443835e-01 1.516362e-01 1.627875e-01 1.745255e-01 1.868640e-01 1.910308e-01 1.998135e-01 2.133857e-01 2.275888e-01 2.310455e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
2.426533e-01 2.579205e-01 2.740610e-01 2.908555e-01 2.908555e-01 3.083075e-01 3.212851e-01 3.264158e-01 3.451803e-01 3.591604e-01 3.645963e-01 3.846601e-01 4.053633e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
4.053633e-01 4.266983e-01 4.486530e-01 4.486530e-01 4.712161e-01 4.903904e-01 4.943712e-01 5.181035e-01 5.423927e-01 5.643110e-01 5.672202e-01 5.925622e-01 6.183962e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
6.265748e-01 6.446951e-01 6.714331e-01 6.714331e-01 6.985796e-01 7.221491e-01 7.261061e-01 7.539791e-01 7.789195e-01 7.821672e-01 8.106347e-01 8.106347e-01 8.393481e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
8.515541e-01 8.682695e-01 8.973637e-01 8.973637e-01 9.265915e-01 9.552297e-01 9.559166e-01 9.852988e-01 9.883860e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000319 0.1319000 0.4476000 0.4573000 0.7559000 1.0000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 2065
length((1:n)[pv < 0.005])
[1] 932
length((1:n)[new_p < 0.10])
[1] 1599
length((1:n)[new_p < 0.05])
[1] 1056
# Export hits with FDR-Adjusted p-value < 0.05.
advallvshlt <- cbind(total_matrix_unk_qn_advall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_advall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(advallvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 35)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_advall[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_ADVALLvsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)


############################################################
# B.3.3. INFALL vs HEALTHY.
############################################################
pv <- c()
for (i in 1:n)
   {
    	t1 <- c(total_matrix_unk_qn_infectionall[i,])
    	t2 <- c(total_matrix_unk_qn_hlt[i,])
#    	test <- wilcoxsign_test(t1 ~ t2, zero.method = "Pratt", ties.method = NULL, distribution = exact(), alternative="two.sided")
#    	pv[i] <- pvalue(test)
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "two.sided")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00000 0.02608 0.17910 0.29850 0.52730 1.00000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
4.043491e-09 1.115055e-05 6.182974e-05 1.567342e-04 2.827917e-04 4.387141e-04 6.180118e-04 8.618569e-04 1.267935e-03 1.627541e-03 2.076637e-03 2.483441e-03 3.136523e-03 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
3.904445e-03 4.792845e-03 6.109786e-03 7.360856e-03 8.814228e-03 1.026278e-02 1.191492e-02 1.379365e-02 1.592373e-02 1.833169e-02 2.010530e-02 2.292334e-02 2.607607e-02 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
2.874254e-02 3.216449e-02 3.559403e-02 3.869366e-02 4.356154e-02 4.743885e-02 5.136611e-02 5.775445e-02 6.236293e-02 6.742343e-02 7.247510e-02 7.800401e-02 8.386503e-02 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
9.007079e-02 9.663391e-02 1.035669e-01 1.108820e-01 1.185915e-01 1.267073e-01 1.352409e-01 1.431710e-01 1.500677e-01 1.584759e-01 1.697205e-01 1.791050e-01 1.901260e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
2.016295e-01 2.136232e-01 2.261138e-01 2.365456e-01 2.502423e-01 2.641336e-01 2.785137e-01 2.886165e-01 3.039265e-01 3.197569e-01 3.278673e-01 3.444775e-01 3.542246e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
3.703628e-01 3.851073e-01 3.974013e-01 4.160593e-01 4.319801e-01 4.449781e-01 4.616224e-01 4.749899e-01 4.906472e-01 5.080895e-01 5.273431e-01 5.459558e-01 5.600861e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
5.737938e-01 5.937672e-01 6.138542e-01 6.313951e-01 6.518312e-01 6.671037e-01 6.913304e-01 7.096819e-01 7.243203e-01 7.366632e-01 7.615472e-01 7.785255e-01 7.993246e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
8.232689e-01 8.426383e-01 8.632827e-01 8.881715e-01 9.021129e-01 9.275017e-01 9.411747e-01 9.672991e-01 9.934001e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000294 0.1041000 0.3578000 0.4121000 0.7006000 1.0000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 2317
length((1:n)[pv < 0.005])
[1] 1034
length((1:n)[new_p < 0.10])
[1] 1793
length((1:n)[new_p < 0.05])
[1] 1223
# Export hits with FDR-Adjusted p-value < 0.05.
infectionallvshlt <- cbind(total_matrix_unk_qn_infectionall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_infectionall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("GeneID","FDR_p-value",colnames(infectionallvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05.txt", sep = "\t", ncolumn = 59)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],1:57])
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts.txt", sep = "\t", ncolumn = 59)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Output Up/Dn Hits.
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59)
write(c("GeneID","FDR_p-value",colnames(total_matrix_unk_qn)[1:57]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59)
for (i in 1:length(writecontentall[,1])) {
	if (mean(total_matrix_unk_qn_infectionall[(1:n)[new_p < 0.05],][i,]) > mean(total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],][i,])) {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts_UP.txt", sep = "\t", ncolumn = 59, append = TRUE)
	} else {
		write(c(row.names(writecontentall)[i], writecontentall[i,]), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts_DN.txt", sep = "\t", ncolumn = 59, append = TRUE)
	}
}
# Heatmap of Up/Dn Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
# Up.
up <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts_UP.txt",  header =TRUE, sep = "\t")
up_cluster <- data.matrix(up[,3:length(up[1,])])
row.names(up_cluster) <- up[,1]
hc_ts <- hclust(dist(up_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(up_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Down.
dn <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-PASS08-SVA-ActualBatches/UNK_RNASeq_hg19_QNSVA_INFALLvsHLT_FDRp0_05_ts_DN.txt",  header =TRUE, sep = "\t")
dn_cluster <- data.matrix(dn[,3:length(up[1,])])
row.names(dn_cluster) <- dn[,1]
hc_ts <- hclust(dist(dn_cluster, method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dn_cluster, Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)




