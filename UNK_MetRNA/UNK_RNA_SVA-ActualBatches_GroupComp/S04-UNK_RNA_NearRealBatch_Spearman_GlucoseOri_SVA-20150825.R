############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Calculate Correlation between UNK RNA rescaled analytes with Glucose levels (moving average).
# This updated version uses the "cor.test" function instead of "cor" used originally to obtain p-values.





############################################################
# SPEARMAN CORRELATION COEFFICIENT > 99TH PERCENTILE: 94
# SPEARMAN CORRELATION COEFFICIENT < 1ST PERCENTILE: 94
# PEARSON CORRELATION COEFFICIENT > 99TH PERCENTILE: 94
# PEARSON CORRELATION COEFFICIENT < 1ST PERCENTILE: 94
# SPEARMAN + PEARSON R >  0.4: 105
# SPEARMAN + PEARSON R < -0.4: 19
############################################################





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
# Load data with original numbers.
# Entrez Gene Symbol used.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNASeq_NearRealBatch_ComBat_Rescaled_forGluCorr.txt",  header =TRUE, sep = "\t")
cyt <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK-Glucose_Input-20150206_tp.txt",  header =TRUE, sep = "\t")

# Colume Names.
# Column 1: ANALYTE_ID;
# Column 2-43: UNK (1-60) OVERLAPS;

# ANALYTE_ID	U06RNA	U07RNA	U08RNA	U09RNA	U14RNA	U17RNA	U18RNA	U19RNA	U23RNA	U24RNA	U25RNA	U26RNA	U27RNA	U28RNA	U29RNA	U30RNA	U31RNA	U33RNA	U34RNA	U35RNA	U36RNA	U37RNA	U38RNA	U39RNA	U40RNA	U42RNA	U43RNA	U45RNA	U46RNA	U47RNA	U48RNA	U49RNA	U50RNA	U51RNA	U52RNA	U53RNA	U54RNA	U56RNA	U58RNA	U58_5RNA	U59RNA	U60RNA

# ANALYTE_ID	D_185	D_186	D_255	D_289	D_301	D_322	D_329	D_369	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1045	D_1051	D_1060	D_1109	D_1124


# Prepare dat.
dat_data <- as.matrix(dat[,c(2:43)])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
row.names(dat_data) <- dat$GeneID
# Alternative Column Names (by Day):
colID <- c("D_185", "D_186", "D_255", "D_289", "D_301", "D_322", "D_329", "D_369", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
colnames(dat_data) <- colID
# Rescaling between -1 and 1 for each row (as D_186RZ and D_255RZ were also included when doing the original rescaling).
n <- length(dat_data[,1])
for (i in 1:n) {
	dat_data[i,] <- 2 * (dat_data[i,] - min(dat_data[i,])) / (max(dat_data[i,]) - min(dat_data[i,])) - 1
}
dat_data[is.na(dat_data)] <- 0

# Prepare cyt.
cyt_data <- as.matrix(cyt[,c(2:43)])
nrow <- length(cyt_data[,1])
ncol <- length(cyt_data[1,])
row.names(cyt_data) <- cyt$ANALYTE_ID
# Alternative Column Names (by Day):
colID <- c("D_185", "D_186", "D_255", "D_289", "D_301", "D_322", "D_329", "D_369", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124")
colnames(cyt_data) <- colID
cyt_data_mean <- cyt_data[1,]
timeseries <- c(185, 186, 255, 289, 301, 322, 329, 369, 476, 532, 546, 602, 615, 616, 618, 620, 625, 647, 679, 680, 683, 688, 694, 700, 711, 796, 840, 912, 944, 945, 948, 959, 966, 984, 1029, 1030, 1032, 1045, 1051, 1060, 1109, 1124)
par(mfrow = c(2,1))
barplot(cyt_data_mean, las = 2, main = "Glucose", ylab = "Concentration (mg/dL)", ylim = c(80,150), beside=TRUE, xpd = FALSE)
plot(cyt_data_mean~timeseries, type = "o", pch = 16, main = "Glucose", xlab = "Time (Days)", ylab = "Concentration (mg/dL)")
# Plot only the barplot to pair with heatmaps.
par(mfrow = c(1,1))
barplot(cyt_data_mean, las = 2, main = "Glucose", ylab = "Concentration (mg/dL)", ylim = c(80,150), beside=TRUE, xpd = FALSE, space = 0)


############################################################
# C. SPEARMAN / PEARSON CORRELATION WITH GLUCOSE LEVELS
############################################################
############################################################
# C.1. SPEARMAN CORRELATION WITH GLUCOSE LEVELS
############################################################
n <- length(dat_data[,1]) # Length of the combined dataset.
n
[1] 9335
scor <- c()
rho <- c()
pv <- c()
for (i in 1:n) {
	scor <- cor.test(dat_data[i,], cyt_data_mean, alternative = c("two.sided"), method = c("spearman"))
	scor$estimate[is.na(scor$estimate)] <- 0
	rho <- c(rho, scor$estimate)
	pv <- c(pv, scor$p.value)
}
rhos <- rho

# Histogram of distribution of Spearman correlation coefficient.
hist(rho, nclass = 100, main = "Distribution of Spearman Correlation Coefficient")

summary(rho)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.56970 -0.13110  0.01322  0.01950  0.16740  0.64340 

scorquantiles <- quantile(rho, seq(0,1,0.01), na.rm = TRUE)
plot(scorquantiles, main = "Quantiles of Spearman Correlation Coefficient Distribution")
scorquantiles
          0%           1%           2%           3%           4%           5%           6%           7%           8%           9%          10%          11%          12% 
-0.569702228 -0.405608842 -0.370942099 -0.342401618 -0.323955579 -0.308038811 -0.294556351 -0.280193781 -0.267694340 -0.258036253 -0.247755168 -0.238209025 -0.228267760 
         13%          14%          15%          16%          17%          18%          19%          20%          21%          22%          23%          24%          25% 
-0.219684577 -0.210015133 -0.201663070 -0.193744392 -0.186088126 -0.179102115 -0.171405497 -0.164753698 -0.157199120 -0.149695873 -0.142690236 -0.137335308 -0.131100707 
         26%          27%          28%          29%          30%          31%          32%          33%          34%          35%          36%          37%          38% 
-0.123886667 -0.117623271 -0.110685931 -0.105130868 -0.099420860 -0.093276276 -0.088514261 -0.082677171 -0.076936976 -0.070172705 -0.064330368 -0.059372068 -0.054431359 
         39%          40%          41%          42%          43%          44%          45%          46%          47%          48%          49%          50%          51% 
-0.048977585 -0.044047887 -0.038207283 -0.032366679 -0.026339351 -0.021501604 -0.014277031 -0.008897417 -0.003408777  0.002446596  0.007773788  0.013222478  0.017685333 
         52%          53%          54%          55%          56%          57%          58%          59%          60%          61%          62%          63%          64% 
 0.022660545  0.028396691  0.034694086  0.040823302  0.046646956  0.052932094  0.059413489  0.065792779  0.071973301  0.078118076  0.084085129  0.090089373  0.097041632 
         65%          66%          67%          68%          69%          70%          71%          72%          73%          74%          75%          76%          77% 
 0.103022179  0.109754679  0.115062720  0.120557605  0.126107718  0.132627639  0.138655934  0.146438376  0.152896765  0.160885921  0.167356303  0.175380353  0.183375207 
         78%          79%          80%          81%          82%          83%          84%          85%          86%          87%          88%          89%          90% 
 0.190115250  0.197421662  0.206627583  0.213912114  0.222940225  0.230947208  0.239484224  0.248688042  0.257296252  0.265997082  0.276111299  0.287296339  0.298243943 
         91%          92%          93%          94%          95%          96%          97%          98%          99%         100% 
 0.308304320  0.318637384  0.332265460  0.346785850  0.362197005  0.381689947  0.401295559  0.422144851  0.469110807  0.643439850 

# Output scor values.
scorcombined <- cbind(row.names(dat_data), rho, pv)
write.table(scorcombined, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_Spearman_Coefficient.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

# Above 99th Quantile.
SGr99quantile <- (1:n)[rho > 0.469110807]
length(SGr99quantile)
[1] 94

# Heatmap of Above 99th Quantile.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data[SGr99quantile,], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data[SGr99quantile,], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Write out >99% Quantile data.
writecontent <- cbind(scorcombined[SGr99quantile,], dat_data[SGr99quantile,])
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_SpearmanR_Gr99QT.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

# Below 1st Quantile.
SSm01quantile <- (1:n)[rho < -0.405608842]
length(SSm01quantile)
[1] 94

# Heatmap of Below 1st Quantile.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data[SSm01quantile,], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data[SSm01quantile,], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Write out <1% Quantile data.
writecontent <- cbind(scorcombined[SSm01quantile,], dat_data[SSm01quantile,])
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_SpearmanR_Sm01QT.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)


############################################################
# C.2. PEARSON CORRELATION WITH GLUCOSE LEVELS
############################################################
n <- length(dat_data[,1]) # Length of the combined dataset.
n
[1] 9335
pcor <- c()
rho <- c()
pv <- c()
for (i in 1:n) {
	pcor <- cor.test(dat_data[i,], cyt_data_mean, alternative = c("two.sided"), method = c("pearson"))
	pcor$estimate[is.na(pcor$estimate)] <- 0
	rho <- c(rho, pcor$estimate)
	pv <- c(pv, pcor$p.value)
}
rhop <- rho

# Histogram of distribution of Spearman correlation coefficient.
hist(rho, nclass = 100, main = "Distribution of Pearson Correlation Coefficient")

summary(rho)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.49660 -0.11040  0.01801  0.02280  0.15710  0.54260 

pcorquantiles <- quantile(rho, seq(0,1,0.01), na.rm = TRUE)
plot(pcorquantiles, main = "Quantiles of Pearson Correlation Coefficient Distribution")
pcorquantiles
           0%            1%            2%            3%            4%            5%            6%            7%            8%            9%           10%           11%           12% 
-0.4966456572 -0.3637286726 -0.3293484869 -0.3049265241 -0.2862067717 -0.2704597399 -0.2589351350 -0.2450593295 -0.2338501195 -0.2229967705 -0.2112027557 -0.2022718261 -0.1947663003 
          13%           14%           15%           16%           17%           18%           19%           20%           21%           22%           23%           24%           25% 
-0.1870988849 -0.1791570931 -0.1714357063 -0.1645690362 -0.1573904563 -0.1516988638 -0.1455209345 -0.1389547206 -0.1334308720 -0.1261449367 -0.1208268821 -0.1151137610 -0.1104255179 
          26%           27%           28%           29%           30%           31%           32%           33%           34%           35%           36%           37%           38% 
-0.1046949418 -0.0993333794 -0.0936304322 -0.0879584180 -0.0833353481 -0.0784991998 -0.0734568553 -0.0680453462 -0.0620416849 -0.0564258568 -0.0517688610 -0.0456378743 -0.0403556422 
          39%           40%           41%           42%           43%           44%           45%           46%           47%           48%           49%           50%           51% 
-0.0359134628 -0.0315112407 -0.0261859124 -0.0215766795 -0.0162566435 -0.0110841821 -0.0062870643 -0.0006677895  0.0038486281  0.0088120518  0.0137478175  0.0180071532  0.0240033556 
          52%           53%           54%           55%           56%           57%           58%           59%           60%           61%           62%           63%           64% 
 0.0291587395  0.0342935449  0.0387170880  0.0430463890  0.0476404260  0.0531341058  0.0593613016  0.0643304291  0.0703349776  0.0759784954  0.0808078820  0.0865522637  0.0926114608 
          65%           66%           67%           68%           69%           70%           71%           72%           73%           74%           75%           76%           77% 
 0.0978765659  0.1026168846  0.1084868218  0.1146007777  0.1203072045  0.1254204138  0.1312419608  0.1374855291  0.1438366166  0.1504253570  0.1571277018  0.1643742910  0.1702461453 
          78%           79%           80%           81%           82%           83%           84%           85%           86%           87%           88%           89%           90% 
 0.1759417229  0.1822190246  0.1892663691  0.1966866255  0.2039652643  0.2109102129  0.2172792381  0.2248221253  0.2324546232  0.2393706589  0.2476301746  0.2568114019  0.2660904877 
          91%           92%           93%           94%           95%           96%           97%           98%           99%          100% 
 0.2754830295  0.2862484584  0.2955544276  0.3074757754  0.3220912191  0.3371199272  0.3514675181  0.3747912728  0.4140459390  0.5426081618 

# Output pcor values.
pcorcombined <- cbind(row.names(dat_data), rho, pv)
write.table(pcorcombined, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_Pearson_Coefficient.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

# Above 99th Quantile.
PGr99quantile <- (1:n)[rho > 0.4140459390]
length(PGr99quantile)
[1] 94

# Heatmap of Above 99th Quantile.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data[PGr99quantile,], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data[PGr99quantile,], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Write out >99% Quantile data.
writecontent <- cbind(pcorcombined[PGr99quantile,], dat_data[PGr99quantile,])
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_PearsonR_Gr99QT.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

# Below 1st Quantile.
PSm01quantile <- (1:n)[rho < -0.3637286726]
length(PSm01quantile)
[1] 94

# Heatmap of Below 1st Quantile.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data[PSm01quantile,], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data[PSm01quantile,], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Write out <1% Quantile data.
writecontent <- cbind(pcorcombined[PSm01quantile,], dat_data[PSm01quantile,])
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_PearsonR_Sm01QT.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)


############################################################
# C.3. SCOR VS PCOR
############################################################
plot(rhos~rhop, main="Spearman Correlation Coefficient vs Pearson Correlation Coefficient", xlab = "Pearson", ylab = "Spearman", xlim = c(-1,1), ylim = c(-1,1), pch = 16, cex = 0.5)

# Use an arbitrary cutoff of correlation coefficient of 0.4 (Strength of Correlation: moderate).
# rhos above 0.4.
SGr40 <- (1:n)[rhos > 0.4]
length(SGr40)
[1] 288
PGr40 <- (1:n)[rhop > 0.4]
length(PGr40)
[1] 118
SPGr40 <- (1:n)[rhos > 0.4 & rhop > 0.4]
length(SPGr40)
[1] 105

# Heatmap of Analytes with both Spearman and Pearson Correlation Coefficient > 0.4.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data[SPGr40,], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data[SPGr40,], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
# Write out R > 0.4 data.
writecontent <- cbind(scorcombined[SPGr40,1:3], pcorcombined[SPGr40,2:3], dat_data[SPGr40,])
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_Spearman_Pearson_Gr40.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)

# Use an arbitrary cutoff of correlation coefficient of -0.4 (Strength of Correlation: moderate).
# scor below -0.4.
SSm40 <- (1:n)[rhos < -0.4]
length(SSm40)
[1] 101
PSm40 <- (1:n)[rhop < -0.4]
length(PSm40)
[1] 45
SPSm40 <- (1:n)[rhos < -0.4 & rhop < -0.4]
length(SPSm40)
[1] 19

# Heatmap of Analytes with both Spearman and Pearson Correlation Coefficient < -0.4.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(dat_data[SPSm40,], method = "euclidean"), method = "ward.D")
total_heatmap <- heatmap.2(dat_data[SPSm40,], Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,8), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5, cexRow = 0.8)
# Write out R < -0.4 data.
writecontent <- cbind(scorcombined[SPSm40,1:3], pcorcombined[SPSm40,2:3], dat_data[SPSm40,])
write.table(writecontent, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNA_NearRealBatch_SVA_Glc_Spearman_Pearson_Sm-40.txt", sep = "\t", row.names = FALSE, col.names = TRUE, append = FALSE)




