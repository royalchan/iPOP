############################################################
# TO-DO LIST.
############################################################
# 1. This script is to use randomized number series of different number of categories and check about the strange local enrichment in methylome correlation data.





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
library("vioplot")





############################################################
# B. LOAD DATA AND INITIAL PROCESSING
############################################################
# Load data with original numbers.
# Entrez Gene Symbol used.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_ComBat_Rescaled.txt",  header =TRUE, sep = "\t")

# GeneID	D_0	D_4	D_21	D_116	D_185	D_186	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_329	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045	D_1051	D_1060	D_1109	D_1124

# Prepare dat.
dat_data <- as.matrix(dat[,c(2:length(dat[1,]))])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
row.names(dat_data) <- dat$GeneID


############################################################
# C. RANDOM NUMBER SERIES GENERATION
############################################################
# Length of random scalar is the same length as that of each row of dat_data.
n <- length(dat_data[1,])
n
[1] 57

# Random Order.
set.seed(7)
randscalar <- matrix(nrow = 10, ncol = 57)
for (i in 1:10) {
	randscalar[i,] <- sample(1:i, n, replace = TRUE, prob = NULL)
}

# Fixed Order.
fixscalar <- matrix(nrow = 10, ncol = 57)
for (i in 1:10) {
	series <- rep(c(1:i), times = round(n/i) + 1)
	fixscalar[i,] <- series[1:n]
}

# Real Order.
realorder <- c(1, 0.728571429, 0.1, 0, 0, 0, 0, 1, 0.859090909, 0.777272727, 0.695454545, 0.572727273, 0.409090909, 0.163636364, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.84, 0.72, 0.6, 0.3, 0.1, 0, 1, 0.857142857, 0.728571429, 0.514285714, 0.257142857, 0.1, 0, 0, 0, 0, 0, 1, 0.8775, 0.81, 0.5625, 0.405, 0.1, 1, 0.870967742, 0.812903226, 0.638709677, 0.435483871, 0.261290323, 0.1, 0, 0)


############################################################
# C. CORRELATION CALCULATION AND VISUALIZATION
############################################################
############################################################
# C.1. SPEARMAN CORRELATION COEFFICIENT CALCULATION AND VISUALIZATION (RANDOM CATEGORIES)
############################################################
l <- length(dat_data[,1])
rho <- c()
cat <- c()
scor <- c()
for (k in 1:10) {
	for (m in 1:l) {
		scor <- cor.test(dat_data[m,], randscalar[k,], alternative = c("two.sided"), method = c("spearman"))
		scor$estimate[is.na(scor$estimate)] <- 0
		rho <- c(rho, scor$estimate)
		cat <- c(cat, k)
	}
}

# Violin Plots.
# P_Value.
dat_data_rho <- data.frame(SPEARMAN_RHO = rho, CATEGORY = cat)
p <- ggplot(dat_data_rho[9336:93350,], aes(x=as.factor(CATEGORY), y=SPEARMAN_RHO))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:9], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "NUMBER_OF_CATEGORIES", y = "SPEARMAN_RHO")

# Write Output Corr.
write.table(cbind(rho, cat), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_RandomCatCorr_Spearman.txt",  col.name = TRUE, row.name = FALSE, sep = "\t", quote = FALSE)


############################################################
# C.2. PEARSON CORRELATION COEFFICIENT CALCULATION AND VISUALIZATION (RANDOM CATEGORIES)
############################################################
l <- length(dat_data[,1])
rho <- c()
cat <- c()
pcor <- c()
for (k in 1:10) {
	for (m in 1:l) {
		pcor <- cor.test(dat_data[m,], randscalar[k,], alternative = c("two.sided"), method = c("pearson"))
		pcor$estimate[is.na(pcor$estimate)] <- 0
		rho <- c(rho, pcor$estimate)
		cat <- c(cat, k)
	}
}

# Violin Plots.
# P_Value.
dat_data_rho <- data.frame(PEARSON_R = rho, CATEGORY = cat)
p <- ggplot(dat_data_rho[9336:93350,], aes(x=as.factor(CATEGORY), y= PEARSON_R))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:9], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "NUMBER_OF_CATEGORIES", y = "PEARSON_R")

# Write Output Corr.
write.table(cbind(rho, cat), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_RandomCatCorr_Pearson.txt",  col.name = TRUE, row.name = FALSE, sep = "\t", quote = FALSE)


############################################################
# C.3. SPEARMAN CORRELATION COEFFICIENT CALCULATION AND VISUALIZATION (FIXED CATEGORIES)
############################################################
l <- length(dat_data[,1])
rho <- c()
cat <- c()
scor <- c()
for (k in 1:10) {
	for (m in 1:l) {
		scor <- cor.test(dat_data[m,], fixscalar[k,], alternative = c("two.sided"), method = c("spearman"))
		scor$estimate[is.na(scor$estimate)] <- 0
		rho <- c(rho, scor$estimate)
		cat <- c(cat, k)
	}
}

# Violin Plots.
# P_Value.
dat_data_rho <- data.frame(SPEARMAN_RHO = rho, CATEGORY = cat)
p <- ggplot(dat_data_rho[9336:93350,], aes(x=as.factor(CATEGORY), y=SPEARMAN_RHO))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:9], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "NUMBER_OF_CATEGORIES", y = "SPEARMAN_RHO")

# Write Output Corr.
write.table(cbind(rho, cat), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_FixedCatCorr_Spearman.txt",  col.name = TRUE, row.name = FALSE, sep = "\t", quote = FALSE)


############################################################
# C.4. PEARSON CORRELATION COEFFICIENT CALCULATION AND VISUALIZATION (FIXED CATEGORIES)
############################################################
l <- length(dat_data[,1])
rho <- c()
cat <- c()
pcor <- c()
for (k in 1:10) {
	for (m in 1:l) {
		pcor <- cor.test(dat_data[m,], fixscalar[k,], alternative = c("two.sided"), method = c("pearson"))
		pcor$estimate[is.na(pcor$estimate)] <- 0
		rho <- c(rho, pcor$estimate)
		cat <- c(cat, k)
	}
}

# Violin Plots.
# P_Value.
dat_data_rho <- data.frame(PEARSON_R = rho, CATEGORY = cat)
p <- ggplot(dat_data_rho[9336:93350,], aes(x=as.factor(CATEGORY), y= PEARSON_R))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:9], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "NUMBER_OF_CATEGORIES", y = "PEARSON_R")

# Write Output Corr.
write.table(cbind(rho, cat), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_FixedCatCorr_Pearson.txt",  col.name = TRUE, row.name = FALSE, sep = "\t", quote = FALSE)


############################################################
# C.5. SPEARMAN CORRELATION COEFFICIENT CALCULATION AND VISUALIZATION (REAL CATEGORIES)
############################################################
l <- length(dat_data[,1])
rho <- c()
cat <- c()
scor <- c()
for (k in 1:1) {
	for (m in 1:l) {
		scor <- cor.test(dat_data[m,], realorder, alternative = c("two.sided"), method = c("spearman"))
		scor$estimate[is.na(scor$estimate)] <- 0
		rho <- c(rho, scor$estimate)
		cat <- c(cat, k)
	}
}

# Violin Plots.
# P_Value.
dat_data_rho <- data.frame(SPEARMAN_RHO = rho, CATEGORY = cat)
p <- ggplot(dat_data_rho, aes(x=as.factor(CATEGORY), y=SPEARMAN_RHO))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:1], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "REAL_CATEGORIES", y = "SPEARMAN_RHO")

# Write Output Corr.
write.table(cbind(rho, cat), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_RealCatCorr_Spearman.txt",  col.name = TRUE, row.name = FALSE, sep = "\t", quote = FALSE)


############################################################
# C.6. PEARSON CORRELATION COEFFICIENT CALCULATION AND VISUALIZATION (REAL CATEGORIES)
############################################################
l <- length(dat_data[,1])
rho <- c()
cat <- c()
pcor <- c()
for (k in 1:1) {
	for (m in 1:l) {
		pcor <- cor.test(dat_data[m,], realorder, alternative = c("two.sided"), method = c("pearson"))
		pcor$estimate[is.na(pcor$estimate)] <- 0
		rho <- c(rho, pcor$estimate)
		cat <- c(cat, k)
	}
}

# Violin Plots.
# P_Value.
dat_data_rho <- data.frame(PEARSON_R = rho, CATEGORY = cat)
p <- ggplot(dat_data_rho, aes(x=as.factor(CATEGORY), y= PEARSON_R))
p + geom_violin() + geom_boxplot(width=.1, fill=rainbow(150)[1:1], outlier.colour=NA) + stat_summary(fun.y=median, geom="point", fill="white", shape=21, size=2.5) + labs(x = "REAL_CATEGORIES", y = "PEARSON_R")

# Write Output Corr.
write.table(cbind(rho, cat), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150924-Pass15-RandomizedCorr/UNK_RNASeq_NearRealBatch_RealCatCorr_Pearson.txt",  col.name = TRUE, row.name = FALSE, sep = "\t", quote = FALSE)




