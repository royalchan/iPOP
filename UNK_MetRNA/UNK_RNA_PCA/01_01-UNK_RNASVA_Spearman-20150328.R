############################################################
# TO-DO LIST.
############################################################
# Calculate Spearman Correlation Coefficient of Combined UNK RNA Data.





############################################################
# A. Library Used.
############################################################
############################################################
# A.1. Full Libraries.
############################################################
library("gplots")
library("ggplot2")
library("reshape2")
library("preprocessCore")
library("coin")
library("magic")
library("RColorBrewer")
library("limma")
library("LDheatmap")
library("MASS")
library("lomb")
library("cts")
library("Sushi")
library("heatmap3")
library("geneplotter")





############################################################
# B. QUANTILE NORMALIZED LOG2 DATA LOADING AND RESCALING.
############################################################
# Load data.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass11-SampleSpearmanCorrCytoscape/UNK_RNASeq_hg19_QNsep_ComBat_Rescaled.txt",  header =TRUE, sep = "\t")

# attach(dat)
# SAMPLE Numbers.
# GeneID	D_0	D_4	D_21	D_116	D_185	D_186	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_329	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045	D_1051	D_1060	D_1109	D_1124

row.names(dat) <- dat$GeneID
dat_data <- as.matrix(dat[,2:length(dat[1,])])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)
datcolnames <- colnames(dat)[2:length(dat[1,])]


############################################################
# C.1. TIME SERIES.
############################################################
# Day Series.
dayseries <- c(0, 4, 21, 116, 185, 186, 255, 289, 290, 292, 294, 297, 301, 307, 311, 322, 329, 369, 380, 400, 476, 532, 546, 602, 615, 616, 618, 620, 625, 630, 647, 679, 680, 683, 688, 694, 700, 711, 735, 796, 840, 912, 944, 945, 948, 959, 966, 984, 1029, 1030, 1032, 1038, 1045, 1051, 1060, 1109, 1124)


############################################################
# C.2. SPEARMAN CORRELATION
############################################################
# Spearman Correlation.
scor <- cor(dat_data, use="everything", method="spearman")
# Output Spearman Correlation Coefficient.
write.table(scor, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass11-SampleSpearmanCorrCytoscape/UNK_RNASVA_Spearman_Rho.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)

scor2 <- scor^2
LDheatmap(scor2, genetic.distance=NULL, distances="physical", LDmeasure="r", title="UNK RNA Post-SVA-Spearman Correlation (Time Points)", add.map=FALSE, add.key=TRUE, geneMapLocation=0.15, geneMapLabelX=0.5, geneMapLabelY=0.3, SNP.name=colnames(dat_data), color=rainbow(1000), newpage=TRUE, name="ldheatmap", vp.name=NULL, pop=FALSE)
# Output Spearman Correlation Coefficient Squared.
write.table(scor2, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass11-SampleSpearmanCorrCytoscape/UNK_RNASVA_Spearman_RhoSquared.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)




