############################################################
# TO-DO LIST.
############################################################
DONE. # Compare and plot variant expression by MAF (overlapping with Nimblegen array variants).





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


############################################################
# B. U TEST P-VALUE DISTRIBUTION BY MAF IN DIFFERENT GROUPS.
############################################################
############################################################
# B.1. RAW DATA LOADING.
############################################################
# Load data.

dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150926-Pass16-RareVarExpression/UNK_AlleleExpMeanByAllele_NimblgenVar_MAF.txt",  header =TRUE, sep = "\t")

# Column Names.
# Chr	Position	Maternal/Paternal	Chrom_Pos_Ref_Alt	ESP_EUR_MAF	1000G_EUR_MAF	Ref	Allele_Count_Mean_Ref_hg19	Allele_Count_Mean_Ref_mat	Allele_Count_Mean_Ref_pat	Allele_Count_Mean_Ref_mp	Alt	Allele_Count_Mean_Alt_hg19	Allele_Count_Mean_Alt_mat	Allele_Count_Mean_Alt_pat	Allele_Count_Mean_Alt_mp	HetHom

datrownames <- row.names(dat) <- paste(dat[,1],"_",dat[,2], sep = "")
n <- length(dat[,1])
n
[1] 1346
dim(dat)
[1] 1346   17


############################################################
# B.2. EXPRESSION DISTRIBUTION BY MAF.
############################################################
# COLUMN NAMES FOR dat_data:
# ESP_EUR_MAF	1000G_EUR_MAF	Allele_Count_Mean_Ref_hg19	Allele_Count_Mean_Ref_mat	Allele_Count_Mean_Ref_pat	Allele_Count_Mean_Ref_mp	Allele_Count_Mean_Alt_hg19	Allele_Count_Mean_Alt_mat	Allele_Count_Mean_Alt_pat	Allele_Count_Mean_Alt_mp
dat_data <- as.matrix(dat[,c(5:6, 8:11, 13:16)])
# dat_data[is.na(dat_data)] <- 0
row.names(dat_data) <- datrownames


############################################################
# B.2.1. DIVIDE DATA BY MAF (3 CATEGORIES).
############################################################
############################################################
# B.2.1.1. DIVIDE DATA BY ESP MAF.
############################################################
# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,1] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.01 & dat_data[,1] < 0.05],])
dat_data_c <-  as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.05],])


############################################################
# B.2.1.1.1. REF.
############################################################
# hg19 Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,3]))), log10(as.numeric(as.character(dat_data_b[,3]))), log10(as.numeric(as.character(dat_data_c[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("MAF >= 0.05", length(dat_data_c[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.7885, df = 2, p-value = 0.248
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,4]))), log10(as.numeric(as.character(dat_data_b[,4]))), log10(as.numeric(as.character(dat_data_c[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("MAF >= 0.05", length(dat_data_c[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.3833, df = 2, p-value = 0.3037
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,5]))), log10(as.numeric(as.character(dat_data_b[,5]))), log10(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF >= 0.05", length(dat_data_c[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.7176, df = 2, p-value = 0.257
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,6]))), log10(as.numeric(as.character(dat_data_b[,6]))), log10(as.numeric(as.character(dat_data_c[,6])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,6])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,6])), rep("MAF >= 0.05", length(dat_data_c[,6]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.3313, df = 2, p-value = 0.3117
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.1.1.2. ALT.
############################################################
# hg19 Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,7]))), log10(as.numeric(as.character(dat_data_b[,7]))), log10(as.numeric(as.character(dat_data_c[,7])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,7])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,7])), rep("MAF >= 0.05", length(dat_data_c[,7]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 1.3559, df = 2, p-value = 0.5076
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,8]))), log10(as.numeric(as.character(dat_data_b[,8]))), log10(as.numeric(as.character(dat_data_c[,8])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,8])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,8])), rep("MAF >= 0.05", length(dat_data_c[,8]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 1.8363, df = 2, p-value = 0.3993
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Paternal Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,9]))), log10(as.numeric(as.character(dat_data_b[,9]))), log10(as.numeric(as.character(dat_data_c[,9])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,9])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,9])), rep("MAF >= 0.05", length(dat_data_c[,9]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.0258, df = 2, p-value = 0.3632
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,10]))), log10(as.numeric(as.character(dat_data_b[,10]))), log10(as.numeric(as.character(dat_data_c[,10])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,10])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,10])), rep("MAF >= 0.05", length(dat_data_c[,10]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 1.8086, df = 2, p-value = 0.4048
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.1.1.3. REF+ALT (TOTAL).
############################################################
# Calculate Sum of Ref + Alt.
dat_data_a_tot1 <- dat_data_a[,3] + dat_data_a[,7]
dat_data_b_tot1 <- dat_data_b[,3] + dat_data_b[,7]
dat_data_c_tot1 <- dat_data_c[,3] + dat_data_c[,7]
dat_data_a_tot2 <- dat_data_a[,4] + dat_data_a[,8]
dat_data_b_tot2 <- dat_data_b[,4] + dat_data_b[,8]
dat_data_c_tot2 <- dat_data_c[,4] + dat_data_c[,8]
dat_data_a_tot3 <- dat_data_a[,5] + dat_data_a[,9]
dat_data_b_tot3 <- dat_data_b[,5] + dat_data_b[,9]
dat_data_c_tot3 <- dat_data_c[,5] + dat_data_c[,9]
dat_data_a_tot4 <- dat_data_a[,6] + dat_data_a[,10]
dat_data_b_tot4 <- dat_data_b[,6] + dat_data_b[,10]
dat_data_c_tot4 <- dat_data_c[,6] + dat_data_c[,10]

# hg19 Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot1))), log10(as.numeric(as.character(dat_data_b_tot1))), log10(as.numeric(as.character(dat_data_c_tot1)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot1)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot1)), rep("MAF >= 0.05", length(dat_data_c_tot1))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.0608, df = 2, p-value = 0.3569
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot2))), log10(as.numeric(as.character(dat_data_b_tot2))), log10(as.numeric(as.character(dat_data_c_tot2)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot2)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot2)), rep("MAF >= 0.05", length(dat_data_c_tot2))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.1325, df = 2, p-value = 0.3443
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot3))), log10(as.numeric(as.character(dat_data_b_tot3))), log10(as.numeric(as.character(dat_data_c_tot3)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot3)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot3)), rep("MAF >= 0.05", length(dat_data_c_tot3))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.4456, df = 2, p-value = 0.2944
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot4))), log10(as.numeric(as.character(dat_data_b_tot4))), log10(as.numeric(as.character(dat_data_c_tot4)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot4)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot4)), rep("MAF >= 0.05", length(dat_data_c_tot4))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.0679, df = 2, p-value = 0.3556
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.1.2. DIVIDE DATA BY 1000G MAF.
############################################################
# Divide Peptides by 1000G MAF [0, 0.01), [0.01, 0.05), [0.05, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.01 & dat_data[,2] < 0.05],])
dat_data_c <-  as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.05],])


############################################################
# B.2.1.2.1. REF.
############################################################
# hg19 Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,3]))), log10(as.numeric(as.character(dat_data_b[,3]))), log10(as.numeric(as.character(dat_data_c[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("MAF >= 0.05", length(dat_data_c[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.5257, df = 2, p-value = 0.2829
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,4]))), log10(as.numeric(as.character(dat_data_b[,4]))), log10(as.numeric(as.character(dat_data_c[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("MAF >= 0.05", length(dat_data_c[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.6213, df = 2, p-value = 0.2696
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,5]))), log10(as.numeric(as.character(dat_data_b[,5]))), log10(as.numeric(as.character(dat_data_c[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("MAF >= 0.05", length(dat_data_c[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.468, df = 2, p-value = 0.2911
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,6]))), log10(as.numeric(as.character(dat_data_b[,6]))), log10(as.numeric(as.character(dat_data_c[,6])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,6])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,6])), rep("MAF >= 0.05", length(dat_data_c[,6]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.3111, df = 2, p-value = 0.3149
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.1.2.2. ALT.
############################################################
# hg19 Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,7]))), log10(as.numeric(as.character(dat_data_b[,7]))), log10(as.numeric(as.character(dat_data_c[,7])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,7])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,7])), rep("MAF >= 0.05", length(dat_data_c[,7]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.0624, df = 2, p-value = 0.3566
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,8]))), log10(as.numeric(as.character(dat_data_b[,8]))), log10(as.numeric(as.character(dat_data_c[,8])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,8])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,8])), rep("MAF >= 0.05", length(dat_data_c[,8]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.4873, df = 2, p-value = 0.2883
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Paternal Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,9]))), log10(as.numeric(as.character(dat_data_b[,9]))), log10(as.numeric(as.character(dat_data_c[,9])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,9])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,9])), rep("MAF >= 0.05", length(dat_data_c[,9]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.4199, df = 2, p-value = 0.2982
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,10]))), log10(as.numeric(as.character(dat_data_b[,10]))), log10(as.numeric(as.character(dat_data_c[,10])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,10])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,10])), rep("MAF >= 0.05", length(dat_data_c[,10]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.2717, df = 2, p-value = 0.3212
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.1.2.3. REF+ALT (TOTAL).
############################################################
# Calculate Sum of Ref + Alt.
dat_data_a_tot1 <- dat_data_a[,3] + dat_data_a[,7]
dat_data_b_tot1 <- dat_data_b[,3] + dat_data_b[,7]
dat_data_c_tot1 <- dat_data_c[,3] + dat_data_c[,7]
dat_data_a_tot2 <- dat_data_a[,4] + dat_data_a[,8]
dat_data_b_tot2 <- dat_data_b[,4] + dat_data_b[,8]
dat_data_c_tot2 <- dat_data_c[,4] + dat_data_c[,8]
dat_data_a_tot3 <- dat_data_a[,5] + dat_data_a[,9]
dat_data_b_tot3 <- dat_data_b[,5] + dat_data_b[,9]
dat_data_c_tot3 <- dat_data_c[,5] + dat_data_c[,9]
dat_data_a_tot4 <- dat_data_a[,6] + dat_data_a[,10]
dat_data_b_tot4 <- dat_data_b[,6] + dat_data_b[,10]
dat_data_c_tot4 <- dat_data_c[,6] + dat_data_c[,10]

# hg19 Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot1))), log10(as.numeric(as.character(dat_data_b_tot1))), log10(as.numeric(as.character(dat_data_c_tot1)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot1)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot1)), rep("MAF >= 0.05", length(dat_data_c_tot1))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.2259, df = 2, p-value = 0.3286
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot2))), log10(as.numeric(as.character(dat_data_b_tot2))), log10(as.numeric(as.character(dat_data_c_tot2)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot2)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot2)), rep("MAF >= 0.05", length(dat_data_c_tot2))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.4974, df = 2, p-value = 0.2869
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot3))), log10(as.numeric(as.character(dat_data_b_tot3))), log10(as.numeric(as.character(dat_data_c_tot3)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot3)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot3)), rep("MAF >= 0.05", length(dat_data_c_tot3))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.3942, df = 2, p-value = 0.3021
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot4))), log10(as.numeric(as.character(dat_data_b_tot4))), log10(as.numeric(as.character(dat_data_c_tot4)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot4)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot4)), rep("MAF >= 0.05", length(dat_data_c_tot4))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.2012, df = 2, p-value = 0.3327
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2. DIVIDE DATA BY MAF (5 CATEGORIES).
############################################################
############################################################
# B.2.2.1. DIVIDE DATA BY ESP MAF.
############################################################
# Divide Peptides by ESP MAF [0, 0.01), [0.01, 0.05), [0.05, 0.10), [0.10, 0.50), [0.50, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,1] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.01 & dat_data[,1] < 0.05],])
dat_data_c <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.05 & dat_data[,1] < 0.10],])
dat_data_d <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.10 & dat_data[,1] < 0.50],])
dat_data_e <- as.matrix(dat_data[(1:n)[dat_data[,1] >= 0.50],])


############################################################
# B.2.2.1.1. REF.
############################################################
# hg19 Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,3]))), log10(as.numeric(as.character(dat_data_b[,3]))), log10(as.numeric(as.character(dat_data_c[,3]))), log10(as.numeric(as.character(dat_data_d[,3]))), log10(as.numeric(as.character(dat_data_e[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,3])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,3])), rep("MAF >= 0.50", length(dat_data_e[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.7119, df = 4, p-value = 0.4464
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,4]))), log10(as.numeric(as.character(dat_data_b[,4]))), log10(as.numeric(as.character(dat_data_c[,4]))), log10(as.numeric(as.character(dat_data_d[,4]))), log10(as.numeric(as.character(dat_data_e[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,4])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,4])), rep("MAF >= 0.50", length(dat_data_e[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.1977, df = 4, p-value = 0.5253
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,5]))), log10(as.numeric(as.character(dat_data_b[,5]))), log10(as.numeric(as.character(dat_data_c[,5]))), log10(as.numeric(as.character(dat_data_d[,5]))), log10(as.numeric(as.character(dat_data_e[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,5])), rep("MAF >= 0.50", length(dat_data_e[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.5599, df = 4, p-value = 0.4688
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,6]))), log10(as.numeric(as.character(dat_data_b[,6]))), log10(as.numeric(as.character(dat_data_c[,6]))), log10(as.numeric(as.character(dat_data_d[,6]))), log10(as.numeric(as.character(dat_data_e[,6])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,6])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,6])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,6])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,6])), rep("MAF >= 0.50", length(dat_data_e[,6]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.2604, df = 4, p-value = 0.5152
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2.1.2. ALT.
############################################################
# hg19 Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,7]))), log10(as.numeric(as.character(dat_data_b[,7]))), log10(as.numeric(as.character(dat_data_c[,7]))), log10(as.numeric(as.character(dat_data_d[,7]))), log10(as.numeric(as.character(dat_data_e[,7])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,7])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,7])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,7])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,7])), rep("MAF >= 0.50", length(dat_data_e[,7]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.2399, df = 4, p-value = 0.6917
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,8]))), log10(as.numeric(as.character(dat_data_b[,8]))), log10(as.numeric(as.character(dat_data_c[,8]))), log10(as.numeric(as.character(dat_data_d[,8]))), log10(as.numeric(as.character(dat_data_e[,8])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,8])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,8])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,8])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,8])), rep("MAF >= 0.50", length(dat_data_e[,8]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.5855, df = 4, p-value = 0.6294
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,9]))), log10(as.numeric(as.character(dat_data_b[,9]))), log10(as.numeric(as.character(dat_data_c[,9]))), log10(as.numeric(as.character(dat_data_d[,9]))), log10(as.numeric(as.character(dat_data_e[,9])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,9])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,9])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,9])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,9])), rep("MAF >= 0.50", length(dat_data_e[,9]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.8646, df = 4, p-value = 0.5807
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,10]))), log10(as.numeric(as.character(dat_data_b[,10]))), log10(as.numeric(as.character(dat_data_c[,10]))), log10(as.numeric(as.character(dat_data_d[,10]))), log10(as.numeric(as.character(dat_data_e[,10])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,10])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,10])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,10])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,10])), rep("MAF >= 0.50", length(dat_data_e[,10]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.5866, df = 4, p-value = 0.6292
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2.1.3. REF+ALT (TOTAL).
############################################################
# Calculate Sum of Ref + Alt.
dat_data_a_tot1 <- dat_data_a[,3] + dat_data_a[,7]
dat_data_b_tot1 <- dat_data_b[,3] + dat_data_b[,7]
dat_data_c_tot1 <- dat_data_c[,3] + dat_data_c[,7]
dat_data_d_tot1 <- dat_data_d[,3] + dat_data_d[,7]
dat_data_e_tot1 <- dat_data_e[,3] + dat_data_e[,7]
dat_data_a_tot2 <- dat_data_a[,4] + dat_data_a[,8]
dat_data_b_tot2 <- dat_data_b[,4] + dat_data_b[,8]
dat_data_c_tot2 <- dat_data_c[,4] + dat_data_c[,8]
dat_data_d_tot2 <- dat_data_d[,4] + dat_data_d[,8]
dat_data_e_tot2 <- dat_data_e[,4] + dat_data_e[,8]
dat_data_a_tot3 <- dat_data_a[,5] + dat_data_a[,9]
dat_data_b_tot3 <- dat_data_b[,5] + dat_data_b[,9]
dat_data_c_tot3 <- dat_data_c[,5] + dat_data_c[,9]
dat_data_d_tot3 <- dat_data_d[,5] + dat_data_d[,9]
dat_data_e_tot3 <- dat_data_e[,5] + dat_data_e[,9]
dat_data_a_tot4 <- dat_data_a[,6] + dat_data_a[,10]
dat_data_b_tot4 <- dat_data_b[,6] + dat_data_b[,10]
dat_data_c_tot4 <- dat_data_c[,6] + dat_data_c[,10]
dat_data_d_tot4 <- dat_data_d[,6] + dat_data_d[,10]
dat_data_e_tot4 <- dat_data_e[,6] + dat_data_e[,10]

# hg19 Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot1))), log10(as.numeric(as.character(dat_data_b_tot1))), log10(as.numeric(as.character(dat_data_c_tot1))), log10(as.numeric(as.character(dat_data_d_tot1))), log10(as.numeric(as.character(dat_data_e_tot1)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot1)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot1)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot1)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot1)), rep("MAF >= 0.50", length(dat_data_e_tot1))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.958, df = 4, p-value = 0.5649
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot2))), log10(as.numeric(as.character(dat_data_b_tot2))), log10(as.numeric(as.character(dat_data_c_tot2))), log10(as.numeric(as.character(dat_data_d_tot2))), log10(as.numeric(as.character(dat_data_e_tot2)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot2)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot2)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot2)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot2)), rep("MAF >= 0.50", length(dat_data_e_tot2))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.9123, df = 4, p-value = 0.5726
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot3))), log10(as.numeric(as.character(dat_data_b_tot3))), log10(as.numeric(as.character(dat_data_c_tot3))), log10(as.numeric(as.character(dat_data_d_tot3))), log10(as.numeric(as.character(dat_data_e_tot3)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot3)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot3)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot3)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot3)), rep("MAF >= 0.50", length(dat_data_e_tot3))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.2801, df = 4, p-value = 0.5121
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot4))), log10(as.numeric(as.character(dat_data_b_tot4))), log10(as.numeric(as.character(dat_data_c_tot4))), log10(as.numeric(as.character(dat_data_d_tot4))), log10(as.numeric(as.character(dat_data_e_tot4)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot4)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot4)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot4)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot4)), rep("MAF >= 0.50", length(dat_data_e_tot4))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 2.9613, df = 4, p-value = 0.5643
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2.2. DIVIDE DATA BY 1000G MAF.
############################################################
# Divide Peptides by 1000G MAF [0, 0.01), [0.01, 0.05), [0.05, 0.10), [0.10, 0.50), [0.50, Max].
dat_data_a <- as.matrix(dat_data[(1:n)[dat_data[,2] < 0.01],])
dat_data_b <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.01 & dat_data[,2] < 0.05],])
dat_data_c <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.05 & dat_data[,2] < 0.10],])
dat_data_d <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.10 & dat_data[,2] < 0.50],])
dat_data_e <- as.matrix(dat_data[(1:n)[dat_data[,2] >= 0.50],])


############################################################
# B.2.2.2.1. REF.
############################################################
# hg19 Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,3]))), log10(as.numeric(as.character(dat_data_b[,3]))), log10(as.numeric(as.character(dat_data_c[,3]))), log10(as.numeric(as.character(dat_data_d[,3]))), log10(as.numeric(as.character(dat_data_e[,3])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,3])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,3])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,3])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,3])), rep("MAF >= 0.50", length(dat_data_e[,3]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.6463, df = 4, p-value = 0.456
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,4]))), log10(as.numeric(as.character(dat_data_b[,4]))), log10(as.numeric(as.character(dat_data_c[,4]))), log10(as.numeric(as.character(dat_data_d[,4]))), log10(as.numeric(as.character(dat_data_e[,4])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,4])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,4])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,4])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,4])), rep("MAF >= 0.50", length(dat_data_e[,4]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.5244, df = 4, p-value = 0.4742
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,5]))), log10(as.numeric(as.character(dat_data_b[,5]))), log10(as.numeric(as.character(dat_data_c[,5]))), log10(as.numeric(as.character(dat_data_d[,5]))), log10(as.numeric(as.character(dat_data_e[,5])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,5])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,5])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,5])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,5])), rep("MAF >= 0.50", length(dat_data_e[,5]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.6332, df = 4, p-value = 0.4579
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,6]))), log10(as.numeric(as.character(dat_data_b[,6]))), log10(as.numeric(as.character(dat_data_c[,6]))), log10(as.numeric(as.character(dat_data_d[,6]))), log10(as.numeric(as.character(dat_data_e[,6])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,6])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,6])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,6])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,6])), rep("MAF >= 0.50", length(dat_data_e[,6]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.3749, df = 4, p-value = 0.4972
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Ref Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2.2.2. ALT.
############################################################
# hg19 Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,7]))), log10(as.numeric(as.character(dat_data_b[,7]))), log10(as.numeric(as.character(dat_data_c[,7]))), log10(as.numeric(as.character(dat_data_d[,7]))), log10(as.numeric(as.character(dat_data_e[,7])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,7])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,7])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,7])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,7])), rep("MAF >= 0.50", length(dat_data_e[,7]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.0487, df = 4, p-value = 0.5497
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,8]))), log10(as.numeric(as.character(dat_data_b[,8]))), log10(as.numeric(as.character(dat_data_c[,8]))), log10(as.numeric(as.character(dat_data_d[,8]))), log10(as.numeric(as.character(dat_data_e[,8])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,8])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,8])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,8])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,8])), rep("MAF >= 0.50", length(dat_data_e[,8]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.271, df = 4, p-value = 0.5135
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,9]))), log10(as.numeric(as.character(dat_data_b[,9]))), log10(as.numeric(as.character(dat_data_c[,9]))), log10(as.numeric(as.character(dat_data_d[,9]))), log10(as.numeric(as.character(dat_data_e[,9])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,9])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,9])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,9])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,9])), rep("MAF >= 0.50", length(dat_data_e[,9]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.1615, df = 4, p-value = 0.5312
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a[,10]))), log10(as.numeric(as.character(dat_data_b[,10]))), log10(as.numeric(as.character(dat_data_c[,10]))), log10(as.numeric(as.character(dat_data_d[,10]))), log10(as.numeric(as.character(dat_data_e[,10])))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a[,10])), rep("0.01 <= MAF < 0.05", length(dat_data_b[,10])), rep("0.05 <= MAF < 0.10", length(dat_data_c[,10])), rep("0.10 <= MAF < 0.50", length(dat_data_d[,10])), rep("MAF >= 0.50", length(dat_data_e[,10]))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.0168, df = 4, p-value = 0.555
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Alt Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 


############################################################
# B.2.2.1.3. REF+ALT (TOTAL).
############################################################
# Calculate Sum of Ref + Alt.
dat_data_a_tot1 <- dat_data_a[,3] + dat_data_a[,7]
dat_data_b_tot1 <- dat_data_b[,3] + dat_data_b[,7]
dat_data_c_tot1 <- dat_data_c[,3] + dat_data_c[,7]
dat_data_d_tot1 <- dat_data_d[,3] + dat_data_d[,7]
dat_data_e_tot1 <- dat_data_e[,3] + dat_data_e[,7]
dat_data_a_tot2 <- dat_data_a[,4] + dat_data_a[,8]
dat_data_b_tot2 <- dat_data_b[,4] + dat_data_b[,8]
dat_data_c_tot2 <- dat_data_c[,4] + dat_data_c[,8]
dat_data_d_tot2 <- dat_data_d[,4] + dat_data_d[,8]
dat_data_e_tot2 <- dat_data_e[,4] + dat_data_e[,8]
dat_data_a_tot3 <- dat_data_a[,5] + dat_data_a[,9]
dat_data_b_tot3 <- dat_data_b[,5] + dat_data_b[,9]
dat_data_c_tot3 <- dat_data_c[,5] + dat_data_c[,9]
dat_data_d_tot3 <- dat_data_d[,5] + dat_data_d[,9]
dat_data_e_tot3 <- dat_data_e[,5] + dat_data_e[,9]
dat_data_a_tot4 <- dat_data_a[,6] + dat_data_a[,10]
dat_data_b_tot4 <- dat_data_b[,6] + dat_data_b[,10]
dat_data_c_tot4 <- dat_data_c[,6] + dat_data_c[,10]
dat_data_d_tot4 <- dat_data_d[,6] + dat_data_d[,10]
dat_data_e_tot4 <- dat_data_e[,6] + dat_data_e[,10]

# hg19 Mapping Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot1))), log10(as.numeric(as.character(dat_data_b_tot1))), log10(as.numeric(as.character(dat_data_c_tot1))), log10(as.numeric(as.character(dat_data_d_tot1))), log10(as.numeric(as.character(dat_data_e_tot1)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot1)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot1)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot1)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot1)), rep("MAF >= 0.50", length(dat_data_e_tot1))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.2854, df = 4, p-value = 0.5112
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_hg19)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal Genome Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot2))), log10(as.numeric(as.character(dat_data_b_tot2))), log10(as.numeric(as.character(dat_data_c_tot2))), log10(as.numeric(as.character(dat_data_d_tot2))), log10(as.numeric(as.character(dat_data_e_tot2)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot2)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot2)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot2)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot2)), rep("MAF >= 0.50", length(dat_data_e_tot2))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.3536, df = 4, p-value = 0.5005
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Mat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Paternal Mapping Ref Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot3))), log10(as.numeric(as.character(dat_data_b_tot3))), log10(as.numeric(as.character(dat_data_c_tot3))), log10(as.numeric(as.character(dat_data_d_tot3))), log10(as.numeric(as.character(dat_data_e_tot3)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot3)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot3)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot3)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot3)), rep("MAF >= 0.50", length(dat_data_e_tot3))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.3116, df = 4, p-value = 0.5071
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_Pat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 

# Maternal + Paternal Genome Mapping Ref+Alt Allele Expression Distribution.
datdatapforbp <- data.frame(datmean = c(log10(as.numeric(as.character(dat_data_a_tot4))), log10(as.numeric(as.character(dat_data_b_tot4))), log10(as.numeric(as.character(dat_data_c_tot4))), log10(as.numeric(as.character(dat_data_d_tot4))), log10(as.numeric(as.character(dat_data_e_tot4)))), category = c(rep("0 <= MAF < 0.01", length(dat_data_a_tot4)), rep("0.01 <= MAF < 0.05", length(dat_data_b_tot4)), rep("0.05 <= MAF < 0.10", length(dat_data_c_tot4)), rep("0.10 <= MAF < 0.50", length(dat_data_d_tot4)), rep("MAF >= 0.50", length(dat_data_e_tot4))))
test <- kruskal.test(datdatapforbp[,1]~ datdatapforbp[,2])
test
############
	Kruskal-Wallis rank sum test

data:  datdatapforbp[, 1] by datdatapforbp[, 2]
Kruskal-Wallis chi-squared = 3.1192, df = 4, p-value = 0.5381
############
boxplot(datdatapforbp[,1]~ datdatapforbp[,2], ylab = "Log10(Read_Count_MatPat)", main = paste("UNK Total (Ref+Alt) Allele Expression for Different MAF\n(Kruskal-Wallis Rank Sum Test P = ", sprintf("%1.2E", test$p.value), ")"), col = c("gold", "gray", "turquoise", "red", "darkblue"))
stripchart(as.numeric(datdatapforbp[,1])~ datdatapforbp[,2], vertical = TRUE, method = "jitter", pch = 16, cex = 0.5, col = c("turquoise", "red", "gold", "darkblue", "gray"), bg = "bisque", add = TRUE) 




