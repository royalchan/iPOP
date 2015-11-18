############################################################
# TO-DO LIST.
############################################################
# Calculate Spearman Correlation Coefficient of UNK Nimblegen Hits.





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
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass11-SampleSpearmanCorrCytoscape/UNK_RNASeq_hg19_QNSVA_AllHits.txt",  header =TRUE, sep = "\t")

# attach(dat)
# SAMPLE Numbers.
# PROBE_ID	UNK-1	UNK-2	UNK-3	UNK-5	UNK-6	UNK-8	UNK-9	UNK-10	UNK-11	UNK-12	UNK-13	UNK-14	UNK-15	UNK-16	UNK-17	UNK-19	UNK-20	UNK-21	UNK-23	UNK-24	UNK-25	UNK-26	UNK-27	UNK-29	UNK-30	UNK-31	UNK-32	UNK-33	UNK-34	UNK-35	UNK-36	UNK-37	UNK-38	UNK-39	UNK-40	UNK-41	UNK-42	UNK-43	UNK-45	UNK-46	UNK-47	UNK-48	UNK-49	UNK-50	UNK-51	UNK-52	UNK-53	UNK-54	UNK-55	UNK-56
# PROBE_ID	D_0	D_4	D_21	D_116	D_185	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045

row.names(dat) <- dat$PROBE_ID
dat_data <- as.matrix(dat[,2:length(dat[1,])])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
datrownames <- row.names(dat)
datcolnames <- colnames(dat)[2:length(dat[1,])]


############################################################
# C.1. TIME SERIES.
############################################################
# Day Series.
dayseries_ori <- c(0, 4, 21, 116, 185, 255, 289, 290, 292, 294, 297, 301, 307, 311, 322, 369, 380, 400, 476, 532, 546, 602, 615, 618, 620, 625, 630, 647, 679, 680, 683, 688, 694, 700, 711, 735, 796, 840, 912, 944, 945, 948, 959, 966, 984, 1029, 1030, 1032, 1038, 1045)
dayseries <- dayseries_ori + 123


############################################################
# C.2. SPEARMAN CORRELATION
############################################################
# Spearman Correlation.
scor <- cor(dat_data, use="everything", method="spearman")
# Output Spearman Correlation Coefficient.
write.table(scor, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass11-SampleSpearmanCorrCytoscape/UNK_RNASVA_AllHits_Spearman_Rho.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)

scor2 <- scor^2
LDheatmap(scor2, genetic.distance=NULL, distances="physical", LDmeasure="r", title="UNK RNA All Hits Post-SVA-Spearman Correlation (Time Points)", add.map=FALSE, add.key=TRUE, geneMapLocation=0.15, geneMapLabelX=0.5, geneMapLabelY=0.3, SNP.name=colnames(dat_data), color=rainbow(1000), newpage=TRUE, name="ldheatmap", vp.name=NULL, pop=FALSE)
# Output Spearman Correlation Coefficient Squared.
write.table(scor2, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass11-SampleSpearmanCorrCytoscape/UNK_RNASVA_AllHits_Spearman_RhoSquared.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)




