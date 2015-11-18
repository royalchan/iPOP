############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Analyze Filtered Data post Shapiro-Wilk Test (One-Tailed) + ACF Hits Removed + SVA, and Replot hit heatmaps with quantile normalized, rescaled data with both my_palette color scheme and 4 color (human, hrv, rsv, adv) color scheme.





#############################################################################################
#############################################################################################
#############################################################################################
### QN U-TEST SUMMARY (AFTER SHAPIRO-WILK TEST):                                          ###
###                                                                                       ###
### HRV1 vs HEALTHY: P < 0.005;  49 Hits (No Significant Hits Post FDR).                  ###
### HRV2 vs HEALTHY: P < 0.005; 125 Hits (No Significant Hits Post FDR).                  ###
### HRV3 vs HEALTHY: P < 0.005; 104 Hits (No Significant Hits Post FDR).                  ###
### RSV  vs HEALTHY: P < 0.005; 278 Hits (No Significant Hits Post FDR).                  ###
### ADV1 vs HEALTHY: P < 0.005; 316 Hits (No Significant Hits Post FDR).                  ###
### ADV2 vs HEALTHY: P < 0.005; 173 Hits (No Significant Hits Post FDR).                  ###
### HRVALL vs HEALTHY: FDR_P < 0.05; 7 Hits (P < 0.005; 302 Hits).                        ###
### ADVALL vs HEALTHY: FDR_P < 0.05; 147 Hits (P < 0.005; 538 Hits).                      ###
### INFALL vs HEALTHY: P < 0.005; 244 Hits (No Significant Hits Post FDR).                ###
### AUTOCORRELATED (Cluster 1, Pre-SVA): 3934 Hits.                                       ###
### AUTOCORRELATED (Cluster 2, Pre-SVA): 3776 Hits.                                       ###
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

dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_Log2_QN_id_NoGrAC_ComBat_Rescaled.txt",  header =TRUE, sep = "\t")



# attach(dat)
# SAMPLE Numbers.

# PROBE_ID	UNK-0	UNK-1	UNK-2	UNK-3	UNK-5	UNK-6	UNK-8	UNK-9	UNK-10	UNK-11	UNK-12	UNK-13	UNK-14	UNK-15	UNK-16	UNK-17	UNK-19	UNK-20	UNK-21	UNK-23	UNK-24	UNK-25	UNK-26	UNK-27	UNK-28	UNK-29	UNK-30	UNK-31	UNK-32	UNK-33	UNK-34	UNK-35	UNK-36	UNK-37	UNK-38	UNK-39	UNK-40	UNK-41	UNK-42	UNK-43	UNK-45	UNK-46	UNK-47	UNK-48	UNK-49	UNK-50	UNK-51	UNK-52	UNK-53	UNK-54	UNK-55	UNK-56	UNKM-1	UNKM-2	UNK6-1	UNK7-1	UNK8-1	UNK9-1	UNK3-1	UNK4-1	UNK5-1
# PROBE_ID	D_-123	D_0	D_4	D_21	D_116	D_185	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045	UNKM-1	UNKM-2	UNK6-1	UNK7-1	UNK8-1	UNK9-1	UNK3-1	UNK4-1	UNK5-1


row.names(dat) <- dat$PROBE_ID
dat_data <- as.matrix(dat[,2:length(dat[1,])])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)
datcolnames <- colnames(dat)[2:length(dat[1,])]


############################################################
# B.2. U TEST (All Quantiles, UNK, No D_616)
############################################################
# Input this time is already quantile normalized and rescaled.
total_matrix_unk_qn <- dat_data

############################################################
# SAME BATCH:
# 	D:-123 - 297;
# 	D: 301 - 615;
# 	D: 616 - 700;
# 	D: 711 - 1029;
# 	D: 1030 - 1045; UNKM-1; UNKM-2; UNK6-1; UNK7-1; UNK8-1; UNK9-1; UNK3-1; UNK4-1;
# 	UNK5-1.

# GROUP DETERMINATION:
# 	HRV1 (Days 0-21; 3 Columns): Columns 2-4;
# 	HRV2 (Days 615-630; 5 Columns): Columns 24-28;
# 	HRV3 (Days 1029-1045; 5 Columns): Columns 47-51;
# 	RSV (Days 289-311; 8 Columns): Columns 8-15;
# 	ADV1 (Days 679-700; 6 Columns): Columns 30-35;
# 	ADV2 (Days 944-984; 6 Columns): Columns 41-46;
# 	HLT (Days -123, 116-255, 322-602, 647, 711-912; 18 Columns): Columns 1, 5-7, 16-23, 29, 36-40.
############################################################
total_matrix_unk_qn_hrv1 <- total_matrix_unk_qn[,2:4]
total_matrix_unk_qn_hrv2 <- total_matrix_unk_qn[,24:28]
total_matrix_unk_qn_hrv3 <- total_matrix_unk_qn[,47:51]
total_matrix_unk_qn_rsv <- total_matrix_unk_qn[,8:15]
total_matrix_unk_qn_adv1 <- total_matrix_unk_qn[,30:35]
total_matrix_unk_qn_adv2 <- total_matrix_unk_qn[,41:46]
total_matrix_unk_qn_hlt <- total_matrix_unk_qn[,c(1,5:7,16:23,29,36:40)]

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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0007519 0.2567000 0.5188000 0.5069000 0.7647000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
0.0007518797 0.0120300752 0.0233082707 0.0308270677 0.0398496241 0.0503759398 0.0624060150 0.0706126404 0.0766917293 0.0924812030 0.0924812030 0.1105263158 0.1105263158 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
0.1308270677 0.1308270677 0.1533834586 0.1533834586 0.1781954887 0.1781954887 0.1781954887 0.2060150376 0.2060150376 0.2353383459 0.2353383459 0.2353383459 0.2567212530 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
0.2669172932 0.2669172932 0.2877077083 0.3000000000 0.3000000000 0.3000000000 0.3345864662 0.3345864662 0.3345864662 0.3345864662 0.3699248120 0.3699248120 0.3699248120 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
0.4067669173 0.4067669173 0.4067669173 0.4067669173 0.4436090226 0.4436090226 0.4436090226 0.4812030075 0.4812030075 0.4812030075 0.4812030075 0.5187969925 0.5187969925 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
0.5187969925 0.5187969925 0.5563909774 0.5563909774 0.5563909774 0.5932330827 0.5932330827 0.5932330827 0.5932330827 0.6300751880 0.6300751880 0.6300751880 0.6654135338 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
0.6654135338 0.6654135338 0.6745145926 0.7000000000 0.7000000000 0.7000000000 0.7330827068 0.7330827068 0.7330827068 0.7646616541 0.7646616541 0.7646616541 0.7939849624 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
0.7939849624 0.7939849624 0.8218045113 0.8218045113 0.8426392647 0.8466165414 0.8466165414 0.8691729323 0.8691729323 0.8894736842 0.8894736842 0.9075187970 0.9075187970 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
0.9233082707 0.9375939850 0.9496240602 0.9496240602 0.9601503759 0.9691729323 0.9766917293 0.9879699248 0.9947368421 1.0000000000 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.7102  0.9542  0.9781  0.9673  1.0000  1.0000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 697
length((1:n)[pv < 0.005])
[1] 49
length((1:n)[new_p < 0.10])
[1] 0
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
hrv1vshlt <- cbind(total_matrix_unk_qn_hrv1[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_hrv1[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(hrv1vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV1vsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 23)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV1vsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV1vsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV1vsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV1vsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV1vsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000297 0.2230000 0.5144000 0.5025000 0.7770000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
2.971857e-05 5.854557e-03 1.212517e-02 1.866326e-02 2.769770e-02 3.334423e-02 4.701458e-02 5.441673e-02 5.545484e-02 6.481619e-02 7.521769e-02 8.677821e-02 9.946804e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
9.946804e-02 1.134060e-01 1.285328e-01 1.285328e-01 1.449672e-01 1.449672e-01 1.625903e-01 1.625903e-01 1.815210e-01 1.854876e-01 2.016405e-01 2.168672e-01 2.230081e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
2.230081e-01 2.454753e-01 2.454753e-01 2.691016e-01 2.691016e-01 2.936789e-01 2.936789e-01 3.139797e-01 3.192368e-01 3.192368e-01 3.455675e-01 3.455675e-01 3.726708e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.726708e-01 3.726708e-01 4.003091e-01 4.003091e-01 4.284823e-01 4.284823e-01 4.284823e-01 4.569229e-01 4.569229e-01 4.856311e-01 4.856311e-01 5.143689e-01 5.143689e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.143689e-01 5.430771e-01 5.430771e-01 5.715177e-01 5.715177e-01 5.715177e-01 5.996909e-01 5.996909e-01 6.273292e-01 6.273292e-01 6.313991e-01 6.544325e-01 6.544325e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
6.807632e-01 6.807632e-01 7.063211e-01 7.063211e-01 7.063211e-01 7.308984e-01 7.308984e-01 7.545247e-01 7.545247e-01 7.769919e-01 7.769919e-01 7.983595e-01 7.983595e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
8.184790e-01 8.184790e-01 8.374097e-01 8.374097e-01 8.550328e-01 8.714672e-01 8.714672e-01 8.865940e-01 9.005320e-01 9.132218e-01 9.132218e-01 9.247823e-01 9.351838e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.445452e-01 9.528366e-01 9.602068e-01 9.666558e-01 9.723023e-01 9.792696e-01 9.848733e-01 9.924515e-01 9.959355e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1796  0.8575  0.9859  0.9116  1.0000  1.0000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1054
length((1:n)[pv < 0.005])
[1] 125
length((1:n)[new_p < 0.10])
[1] 0
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
hrv2vshlt <- cbind(total_matrix_unk_qn_hrv2[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_hrv2[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(hrv2vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV2vsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 25)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV2vsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV2vsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV2vsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,17), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,17), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV2vsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,17), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,17), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,17), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV2vsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,17), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000297 0.2455000 0.5144000 0.5006000 0.7545000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
2.971857e-05 9.191800e-03 1.644697e-02 2.769770e-02 3.334423e-02 4.672575e-02 5.545484e-02 6.481619e-02 7.521769e-02 8.677821e-02 8.677821e-02 9.946804e-02 1.134060e-01 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
1.239250e-01 1.285328e-01 1.449672e-01 1.625903e-01 1.625903e-01 1.815210e-01 1.815210e-01 2.016405e-01 2.016405e-01 2.168672e-01 2.230081e-01 2.230081e-01 2.454753e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
2.454753e-01 2.691016e-01 2.691016e-01 2.936789e-01 2.936789e-01 2.936789e-01 3.192368e-01 3.192368e-01 3.455675e-01 3.455675e-01 3.455675e-01 3.726708e-01 3.726708e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
4.003091e-01 4.003091e-01 4.003091e-01 4.284823e-01 4.284823e-01 4.284823e-01 4.569229e-01 4.569229e-01 4.851321e-01 4.856311e-01 4.856311e-01 5.143689e-01 5.143689e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.143689e-01 5.430771e-01 5.430771e-01 5.430771e-01 5.715177e-01 5.715177e-01 5.715177e-01 5.996909e-01 5.996909e-01 6.273292e-01 6.273292e-01 6.273292e-01 6.544325e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
6.544325e-01 6.544325e-01 6.807632e-01 6.807632e-01 7.063211e-01 7.063211e-01 7.063211e-01 7.308984e-01 7.308984e-01 7.545247e-01 7.545247e-01 7.769919e-01 7.769919e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
7.769919e-01 7.983595e-01 7.983595e-01 8.184790e-01 8.184790e-01 8.374097e-01 8.550328e-01 8.550328e-01 8.714672e-01 8.714672e-01 8.865940e-01 9.005320e-01 9.101975e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.132218e-01 9.247823e-01 9.351838e-01 9.445452e-01 9.528366e-01 9.666558e-01 9.723023e-01 9.827885e-01 9.924515e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.2566  0.9323  0.9817  0.9526  0.9954  1.0000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 885
length((1:n)[pv < 0.005])
[1] 104
length((1:n)[new_p < 0.10])
[1] 0
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
hrv3vshlt <- cbind(total_matrix_unk_qn_hrv3[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_hrv3[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(hrv3vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV3vsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 25)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV3vsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV3vsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV3vsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV3vsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRV3vsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000077 0.1835000 0.4675000 0.4924000 0.8014000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
7.681106e-06 2.079019e-03 5.551519e-03 9.366469e-03 1.512538e-02 2.036837e-02 2.349394e-02 3.088061e-02 3.519995e-02 4.524812e-02 5.103583e-02 5.735930e-02 6.426141e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
7.175177e-02 7.987006e-02 8.862140e-02 9.804420e-02 1.081353e-01 1.189285e-01 1.304156e-01 1.426234e-01 1.555347e-01 1.555347e-01 1.691712e-01 1.835055e-01 1.835055e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
1.985521e-01 2.142747e-01 2.142747e-01 2.306783e-01 2.306783e-01 2.477176e-01 2.653912e-01 2.653912e-01 2.836415e-01 2.987955e-01 3.024596e-01 3.217814e-01 3.217814e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.415890e-01 3.486546e-01 3.618095e-01 3.824205e-01 3.824205e-01 4.033413e-01 4.033413e-01 4.245443e-01 4.459452e-01 4.459452e-01 4.675112e-01 4.675112e-01 4.891552e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.108448e-01 5.108448e-01 5.324888e-01 5.332192e-01 5.540548e-01 5.754557e-01 5.754557e-01 5.966587e-01 6.175795e-01 6.175795e-01 6.381905e-01 6.410149e-01 6.584110e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
6.782186e-01 6.791094e-01 6.975404e-01 7.163585e-01 7.346088e-01 7.346088e-01 7.522824e-01 7.693217e-01 7.734113e-01 7.857253e-01 8.014479e-01 8.164945e-01 8.308288e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
8.345713e-01 8.444653e-01 8.573766e-01 8.695844e-01 8.810715e-01 8.918647e-01 9.019558e-01 9.113786e-01 9.201299e-01 9.295535e-01 9.426407e-01 9.489642e-01 9.600160e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.648001e-01 9.730156e-01 9.796316e-01 9.848746e-01 9.906335e-01 9.944485e-01 9.968635e-01 9.986584e-01 9.996102e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1161  0.7284  0.9337  0.8364  1.0000  1.0000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1469
length((1:n)[pv < 0.005])
[1] 278
length((1:n)[new_p < 0.10])
[1] 0
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
rsvvshlt <- cbind(total_matrix_unk_qn_rsv[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_rsv[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(rsvvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_RSVvsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 25)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_RSVvsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_RSVvsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_RSVvsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_RSVvsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_RSVvsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000074 0.1716000 0.5000000 0.5038000 0.8284000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
7.429641e-06 2.214033e-03 4.710393e-03 8.185855e-03 1.121876e-02 1.650123e-02 1.978513e-02 2.792802e-02 3.288359e-02 3.846325e-02 4.476359e-02 5.177717e-02 5.962287e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
6.826354e-02 7.781806e-02 8.825671e-02 8.825671e-02 9.967607e-02 1.120316e-01 1.216213e-01 1.254346e-01 1.397961e-01 1.430088e-01 1.552201e-01 1.716173e-01 1.716173e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
1.890621e-01 2.074430e-01 2.074430e-01 2.268418e-01 2.471024e-01 2.471024e-01 2.682918e-01 2.682918e-01 2.902612e-01 2.902612e-01 3.130405e-01 3.199261e-01 3.364587e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.605605e-01 3.605605e-01 3.851303e-01 3.851303e-01 4.102054e-01 4.102054e-01 4.355850e-01 4.355850e-01 4.612693e-01 4.612693e-01 4.870501e-01 5.000000e-01 5.129499e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.387307e-01 5.387307e-01 5.644150e-01 5.644150e-01 5.897946e-01 5.897946e-01 6.148697e-01 6.179363e-01 6.394395e-01 6.635413e-01 6.635413e-01 6.869595e-01 6.869595e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
7.097388e-01 7.151822e-01 7.317082e-01 7.528976e-01 7.528976e-01 7.731582e-01 7.731582e-01 7.925570e-01 8.109379e-01 8.109379e-01 8.283827e-01 8.447799e-01 8.447799e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
8.602039e-01 8.745654e-01 8.879684e-01 8.915745e-01 9.003239e-01 9.117433e-01 9.221819e-01 9.317365e-01 9.403771e-01 9.482228e-01 9.552364e-01 9.615367e-01 9.671164e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.764035e-01 9.802149e-01 9.834988e-01 9.887812e-01 9.926001e-01 9.940786e-01 9.965432e-01 9.983135e-01 9.995245e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05614 0.68360 0.99430 0.82940 1.00000 1.00000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1585
length((1:n)[pv < 0.005])
[1] 316
length((1:n)[new_p < 0.10])
[1] 34
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
adv1vshlt <- cbind(total_matrix_unk_qn_adv1[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_adv1[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(adv1vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV1vsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 25)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV1vsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV1vsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV1vsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV1vsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV1vsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000074 0.2074000 0.4871000 0.4945000 0.7732000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
7.429641e-06 4.710393e-03 9.484187e-03 1.650123e-02 2.359654e-02 2.792802e-02 3.318937e-02 4.476359e-02 5.177717e-02 5.962287e-02 6.826354e-02 7.781806e-02 8.581832e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
8.825671e-02 9.967607e-02 1.120316e-01 1.254346e-01 1.254346e-01 1.397961e-01 1.506714e-01 1.552201e-01 1.716173e-01 1.716173e-01 1.890621e-01 1.890621e-01 2.074430e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
2.074430e-01 2.268418e-01 2.471024e-01 2.471024e-01 2.682918e-01 2.682918e-01 2.682918e-01 2.902612e-01 2.902612e-01 3.130405e-01 3.130405e-01 3.364587e-01 3.364587e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.605605e-01 3.605605e-01 3.851303e-01 3.851303e-01 4.102054e-01 4.102054e-01 4.355850e-01 4.355850e-01 4.612693e-01 4.612693e-01 4.867015e-01 4.870501e-01 4.870501e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.129499e-01 5.129499e-01 5.387307e-01 5.387307e-01 5.644150e-01 5.644150e-01 5.897946e-01 5.897946e-01 6.148697e-01 6.148697e-01 6.305860e-01 6.394395e-01 6.430959e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
6.635413e-01 6.676479e-01 6.869595e-01 7.097388e-01 7.097388e-01 7.317082e-01 7.317082e-01 7.528976e-01 7.528976e-01 7.731582e-01 7.731582e-01 7.925570e-01 8.109379e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
8.109379e-01 8.283827e-01 8.334498e-01 8.447799e-01 8.602039e-01 8.745654e-01 8.745654e-01 8.879684e-01 9.003239e-01 9.117433e-01 9.221819e-01 9.317365e-01 9.403771e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.482228e-01 9.552364e-01 9.671164e-01 9.720720e-01 9.802149e-01 9.863517e-01 9.926001e-01 9.963000e-01 9.990713e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.05614 0.79750 0.95250 0.87920 1.00000 1.00000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1167
length((1:n)[pv < 0.005])
[1] 173
length((1:n)[new_p < 0.10])
[1] 2
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
adv2vshlt <- cbind(total_matrix_unk_qn_adv2[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_adv2[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(adv2vshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV2vsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 25)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV2vsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV2vsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV2vsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV2vsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADV2vsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()


############################################################
# B.3. Grouped Comparison (Group All HRVs together, All ADVs together, and all infections together).
############################################################
# HRV All: 13 Columns.
total_matrix_unk_qn_hrvall <- total_matrix_unk_qn[,c(2:4,24:28,47:51)]
# ADV All: 12 Columns.
total_matrix_unk_qn_advall <- total_matrix_unk_qn[,c(30:35,41:46)]
# Infection All: 33 Columns.
total_matrix_unk_qn_infectionall <- total_matrix_unk_qn[,c(2:4,24:28,47:51,30:35,41:46,8:15)]
# Healthy All: 18 Columns.
total_matrix_unk_qn_hlt <- total_matrix_unk_qn[,c(1,5:7,16:23,29,36:40)]


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000018 0.1761000 0.4921000 0.4936000 0.8025000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
1.808458e-06 1.575133e-03 5.189045e-03 8.773658e-03 1.235169e-02 1.707582e-02 2.098667e-02 2.603084e-02 3.132269e-02 3.740681e-02 4.452477e-02 5.281952e-02 6.225195e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
6.742027e-02 7.870982e-02 8.484888e-02 9.132793e-02 1.053363e-01 1.128781e-01 1.207850e-01 1.290616e-01 1.377105e-01 1.467346e-01 1.561341e-01 1.659098e-01 1.760597e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
1.865823e-01 1.974731e-01 2.087281e-01 2.203405e-01 2.323037e-01 2.446083e-01 2.572453e-01 2.702026e-01 2.834689e-01 2.970297e-01 3.108713e-01 3.249769e-01 3.375051e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.393305e-01 3.539133e-01 3.687076e-01 3.836924e-01 3.988484e-01 4.141533e-01 4.295861e-01 4.295861e-01 4.451236e-01 4.607438e-01 4.764224e-01 4.921370e-01 4.969334e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.078630e-01 5.235776e-01 5.392562e-01 5.548764e-01 5.704139e-01 5.811647e-01 5.858467e-01 6.011516e-01 6.163076e-01 6.312924e-01 6.460867e-01 6.606695e-01 6.701751e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
6.750231e-01 6.891287e-01 7.029703e-01 7.165311e-01 7.297974e-01 7.427547e-01 7.553917e-01 7.676963e-01 7.796595e-01 7.912719e-01 8.025269e-01 8.161558e-01 8.340902e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
8.438659e-01 8.532654e-01 8.622895e-01 8.709384e-01 8.812114e-01 8.946637e-01 9.018449e-01 9.151511e-01 9.212902e-01 9.325797e-01 9.377480e-01 9.471805e-01 9.539252e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.609305e-01 9.672375e-01 9.726802e-01 9.790133e-01 9.846415e-01 9.889564e-01 9.930886e-01 9.963399e-01 9.984249e-01 9.999987e-01 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.02218 0.69670 0.96500 0.83140 1.00000 1.00000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1600
length((1:n)[pv < 0.005])
[1] 302
length((1:n)[new_p < 0.10])
[1] 33
length((1:n)[new_p < 0.05])
[1] 7
# Export hits with FDR-Adjusted p-value < 0.05.
hrvallvshlt <- cbind(total_matrix_unk_qn_hrvall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_hrvall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("PROBE_ID","FDR_p-value",colnames(hrvallvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRVALLvsHLT_FDRp0_05Gr.txt", sep = "\t", ncolumn = 28)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRVALLvsHLT_FDRp0_05Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with FDR-Adjusted p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],])
write(c("PROBE_ID","FDR_p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRVALLvsHLT_FDRp0_05Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRVALLvsHLT_FDRp0_05Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRVALLvsHLT_FDRp0_05Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_HRVALLvsHLT_FDRp0_05Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000002 0.1426000 0.4751000 0.4935000 0.8476000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
2.196704e-07 5.727501e-04 1.839971e-03 3.850706e-03 5.798859e-03 9.656051e-03 1.228052e-02 1.547500e-02 1.932764e-02 2.393272e-02 2.938926e-02 3.579980e-02 3.939577e-02 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
4.743287e-02 5.189981e-02 5.668341e-02 6.724646e-02 7.304899e-02 7.921455e-02 9.267122e-02 9.998013e-02 1.076878e-01 1.157990e-01 1.243210e-01 1.332564e-01 1.426102e-01 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
1.523822e-01 1.703946e-01 1.842134e-01 1.956526e-01 2.074999e-01 2.197459e-01 2.323845e-01 2.454031e-01 2.587923e-01 2.725366e-01 2.866236e-01 2.982236e-01 3.070464e-01 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
3.157542e-01 3.307610e-01 3.460372e-01 3.615586e-01 3.773051e-01 3.932503e-01 4.093723e-01 4.256426e-01 4.420378e-01 4.578526e-01 4.585284e-01 4.750897e-01 4.916914e-01 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
5.083086e-01 5.249103e-01 5.414716e-01 5.579622e-01 5.743574e-01 5.906277e-01 6.067497e-01 6.226949e-01 6.384414e-01 6.539628e-01 6.692390e-01 6.842458e-01 6.989656e-01 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
7.133764e-01 7.274634e-01 7.412077e-01 7.545969e-01 7.676155e-01 7.802541e-01 7.925001e-01 8.043474e-01 8.186879e-01 8.374249e-01 8.476178e-01 8.573898e-01 8.667436e-01 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
8.756790e-01 8.903681e-01 9.000199e-01 9.073288e-01 9.155873e-01 9.269510e-01 9.327535e-01 9.433166e-01 9.481002e-01 9.567313e-01 9.642002e-01 9.706107e-01 9.768654e-01 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
9.826865e-01 9.861986e-01 9.903439e-01 9.933878e-01 9.955745e-01 9.972344e-01 9.986594e-01 9.995221e-01 9.998734e-01 1.000000e+00 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
0.002621 0.567000 0.947800 0.775700 1.000000 1.000000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 2060
length((1:n)[pv < 0.005])
[1] 538
length((1:n)[new_p < 0.10])
[1] 377
length((1:n)[new_p < 0.05])
[1] 147
# Export hits with FDR-Adjusted p-value < 0.05.
advallvshlt <- cbind(total_matrix_unk_qn_advall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
writecontent <- cbind(new_p[(1:n)[new_p < 0.05]],total_matrix_unk_qn_advall[(1:n)[new_p < 0.05],],total_matrix_unk_qn_hlt[(1:n)[new_p < 0.05],])
write(c("PROBE_ID","FDR_p-value",colnames(advallvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADVALLvsHLT_FDRp0_05Gr.txt", sep = "\t", ncolumn = 28)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADVALLvsHLT_FDRp0_05Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with FDR-Adjusted p-value < 0.05 (whole time course).
writecontentall <- cbind(new_p[(1:n)[new_p < 0.05]], total_matrix_unk_qn[(1:n)[new_p < 0.05],])
write(c("PROBE_ID","FDR_p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADVALLvsHLT_FDRp0_05Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADVALLvsHLT_FDRp0_05Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADVALLvsHLT_FDRp0_05Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_ADVALLvsHLT_FDRp0_05Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.6)
dev.off()


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
		test <- wilcox.test(t1, t2, paired = FALSE, alternative = "greater")
		pv[i] <- test$p.value
   }
pv[is.na(pv)] <- 1   

summary(pv)
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0000134 0.1902000 0.4805000 0.4913000 0.7915000 1.0000000 

quantile(pv, probs = seq(0,1,0.01))
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
0.0000134396 0.0025932150 0.0065548778 0.0103079335 0.0157574371 0.0202612662 0.0267004833 0.0324845107 0.0388014671 0.0460614225 0.0521763828 0.0599640229 0.0688888752 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
0.0743366658 0.0831033515 0.0926088227 0.1018020071 0.1101597343 0.1177959625 0.1299252786 0.1384665859 0.1473747951 0.1566511800 0.1670641020 0.1814515050 0.1901963592 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
0.1974256761 0.2085232722 0.2199729723 0.2317678953 0.2438999254 0.2563597151 0.2691366947 0.2822190904 0.2888708167 0.3023866875 0.3161734496 0.3287041094 0.3373267828 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
0.3503551719 0.3589964942 0.3736992879 0.3837437712 0.3960891693 0.4112097028 0.4257309893 0.4341282053 0.4490267731 0.4572658570 0.4727751142 0.4805463984 0.4961078279 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
0.5116750649 0.5196542552 0.5349860254 0.5472419818 0.5581798653 0.5735393574 0.5887902972 0.5973092350 0.6114155160 0.6263007121 0.6336761488 0.6482803694 0.6626732172 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
0.6697846437 0.6838265504 0.6976133125 0.7111291833 0.7177809096 0.7308633053 0.7436402849 0.7561000746 0.7682321047 0.7800270277 0.7914767278 0.8025743239 0.8133141650 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
0.8236918193 0.8337040556 0.8478939443 0.8571252761 0.8658498100 0.8775732165 0.8860668003 0.8971217868 0.9050326201 0.9138115525 0.9256633342 0.9337193359 0.9434083202 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
0.9519662613 0.9594752857 0.9675154893 0.9732995167 0.9797727554 0.9865173509 0.9907701814 0.9954216476 0.9982894105 0.9999991572 

# FDR-Adjusted p-value.
new_p <- p.adjust(pv, method = "fdr", n = length(pv))
summary(new_p)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1187  0.7549  0.9592  0.8522  1.0000  1.0000 

# Histogram of p-values.
par(mfrow=c(1,2))
hist(pv, nclass=100)
hist(new_p, nclass=100)

# Export Top Hits.
length((1:n)[pv < 0.05])
[1] 1454
length((1:n)[pv < 0.005])
[1] 244
length((1:n)[new_p < 0.10])
[1] 0
length((1:n)[new_p < 0.05])
[1] 0
# Export hits with p-value < 0.005.
infectionallvshlt <- cbind(total_matrix_unk_qn_infectionall[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
writecontent <- cbind(pv[(1:n)[pv < 0.005]],total_matrix_unk_qn_infectionall[(1:n)[pv < 0.005],],total_matrix_unk_qn_hlt[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(infectionallvshlt)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_INFALLvsHLT_p0_005Gr.txt", sep = "\t", ncolumn = 25)
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_INFALLvsHLT_p0_005Gr.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Export hits with p-value < 0.005 (whole time course).
writecontentall <- cbind(pv[(1:n)[pv < 0.005]], total_matrix_unk_qn[(1:n)[pv < 0.005],])
write(c("PROBE_ID","p-value",colnames(total_matrix_unk_qn)), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_INFALLvsHLT_p0_005Gr_ts.txt", sep = "\t", ncolumn = 53)
write.table(writecontentall, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_INFALLvsHLT_p0_005Gr_ts.txt", sep = "\t", row.names = TRUE, col.names = FALSE, append = TRUE)
# Heatmap of Hits -- Full Time Series View.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(writecontentall[,2:(length(writecontentall[1,]))], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_INFALLvsHLT_p0_005Gr_ts_hm_color.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()
total_heatmap <- heatmap.2(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmapc <- heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
png('/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_NoGrAC_SVA_INFALLvsHLT_p0_005Gr_ts_hm_color_ori.png')
heatmap.colorplay(writecontentall[,2:(length(writecontentall[1,]))], Rowv=FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,16), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
dev.off()





############################################################
# APPENDIX.1. Hijacking heatmap.2 to Calculate Row / Column Medians instead of Row / Column Means in Row / Column Scaling
############################################################
heatmap.median <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- apply(x, 1, median, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- apply(x, 2, median, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!missing(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    x.scaled <- x.scaled - x.scaled[,7]
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


############################################################
# APPENDIX.2. Hijacking heatmap.2 to Use Day 255 as Baseline in Row Scaling.
############################################################
heatmap.d255 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- x[,7]
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- apply(x, 2, median, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!missing(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    x.scaled <- x.scaled - x.scaled[,7]
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}


############################################################
# APPENDIX.3. Further Play with heatmap.2 to Display Different Colors
############################################################
heatmap.colorplay <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, perl = TRUE,
    ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorder(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorder(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- apply(x, 1, mean, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- apply(x, 2, mean, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 501
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    
    # Gray12-White-Gold
	my_palette1 <- colorRampPalette(c("gray12", "white", "gold"))(n = 500)
    # Gray12-White-Red
	my_palette2 <- colorRampPalette(c("gray12", "white", "red"))(n = 500)
    # Gray12-White-Blue
	my_palette3 <- colorRampPalette(c("gray12", "white", "blue"))(n = 500)
    # Gray12-White-Purple
	my_palette4 <- colorRampPalette(c("gray12", "white", "turquoise"))(n = 500)

    # # Turquoise-White-Gold
	# my_palette1 <- colorRampPalette(c("turquoise", "white", "gold"))(n = 500)
    # # Blue-White-Red
	# my_palette2 <- colorRampPalette(c("blue", "white", "red"))(n = 500)
    # # Green-White-Orange
	# my_palette3 <- colorRampPalette(c("green", "white", "orange"))(n = 500)
    # # Gray12-White-Purple
	# my_palette4 <- colorRampPalette(c("gray12", "white", "purple"))(n = 500)

    # # White-Turquoise
	# my_palette1 <- colorRampPalette(c("white", "gold"))(n = 500)
    # # White-Red
	# my_palette2 <- colorRampPalette(c("white", "red"))(n = 500)
    # # White-Gold
	# my_palette3 <- colorRampPalette(c("white", "blue"))(n = 500)
    # # White-Purple
	# my_palette4 <- colorRampPalette(c("white", "turquoise"))(n = 500)

	xtr <- t(x)
	rownamex <- row.names(xtr)
	color <-  matrix(nrow = nr, ncol = 500)
	par(mfrow = c(nr,1), mar=c(0,5,0,5))
	for (i in nr:1) {
		if (length(name <- grep("HRV", rownamex[i], ignore.case = TRUE))) {
			color[i,] <- my_palette1
		} else if (length(name <- grep("RSV", rownamex[i], ignore.case = TRUE))) {
			color[i,] <- my_palette2
		} else if (length(name <- grep("ADV", rownamex[i], ignore.case = TRUE))) {
			color[i,] <- my_palette3
		} else {
			color[i,] <- my_palette4
		}
		mat <- matrix(rnorm(nr), nrow=nr, ncol=1)
		mat <- as.matrix(x[,i])
		image(mat, axes = FALSE, col = color[i,], breaks = breaks, ...)
	}

	# image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        # c(0, nr), axes = FALSE, xlab = "", ylab = "", col = color[1:nr,], 
        # ...)
    
    
    
    # image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        # c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        # breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
        cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        par(mar = c(5, 4, 2, 1), cex = 0.75)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        axis(1, at = xv, labels = lv)
        if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
        else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
        else mtext(side = 1, "Value", line = 2)
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
        }
        else title("Color Key")
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    invisible(retval)
}

