############################################################
# TO-DO LIST.
############################################################
DONE. # 1. Remove potential batch effect with SVA for ACF Hits removed ShapiroWilk data.





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

library("sva")
library("bladderbatch")
data("bladderdata")
library("pamr")
library("limma")
###################################################
# SVA Tutorial
pheno = pData(bladderEset)
edata = exprs(bladderEset)
mod = model.matrix(~as.factor(cancer), data = pheno)
mod0 = model.matrix(~1, data = pheno)
n.sv = num.sv(edata, mod, method = "leek")
n.sv
svobj = sva(edata, mod, mod0, n.sv = n.sv)
pValues = f.pvalue(edata, mod, mod0)
qValues = p.adjust(pValues, method = "BH")
modSv = cbind(mod, svobj$sv)
mod0Sv = cbind(mod0, svobj$sv)
pValuesSv = f.pvalue(edata, modSv, mod0Sv)
qValuesSv = p.adjust(pValuesSv, method = "BH")
fit = lmFit(edata, modSv)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)
eb = eBayes(fitContrasts)
topTableF(eb, adjust = "BH")
batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)
pValuesComBat = f.pvalue(combat_edata,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")
modBatch = model.matrix(~as.factor(cancer) + as.factor(batch),data=pheno)
mod0Batch = model.matrix(~as.factor(batch),data=pheno)
pValuesBatch = f.pvalue(edata,modBatch,mod0Batch)
qValuesBatch = p.adjust(pValuesBatch,method="BH")
n.sv = num.sv(edata,mod,vfilter=2000,method="leek")
svobj = sva(edata,mod,mod0,n.sv=n.sv,vfilter=2000)
set.seed(12354)
trainIndicator = sample(1:57,size=30,replace=F)
testIndicator = (1:57)[-trainIndicator]
trainData = edata[,trainIndicator]
testData = edata[,testIndicator]
trainPheno = pheno[trainIndicator,]
testPheno = pheno[testIndicator,]
mydata = list(x=trainData,y=trainPheno$cancer)
mytrain = pamr.train(mydata)
table(pamr.predict(mytrain,testData,threshold=2),testPheno$cancer)
trainMod = model.matrix(~cancer,data=trainPheno)
trainMod0 = model.matrix(~1,data=trainPheno)
trainSv = sva(trainData,trainMod,trainMod0)
fsvaobj = fsva(trainData,trainMod,trainSv,testData)
mydataSv = list(x=fsvaobj$db,y=trainPheno$cancer)
mytrainSv = pamr.train(mydataSv)
table(pamr.predict(mytrainSv,fsvaobj$new,threshold=1),testPheno$cancer)
library(zebrafishRNASeq)
library(genefilter)
data(zfGenes)
filter = apply(zfGenes, 1, function(x) length(x[x>5])>=2)
filtered = zfGenes[filter,]
genes = rownames(filtered)[grep("^ENS", rownames(filtered))]
controls = grepl("^ERCC", rownames(filtered))
group = as.factor(rep(c("Ctl", "Trt"), each=3))
dat0 = as.matrix(filtered)
mod1 = model.matrix(~group)
mod0 = cbind(mod1[,1])
svseq = svaseq(dat0,mod1,mod0,n.sv=1)$sv
plot(svseq,pch=19,col="blue")
sup_svseq = svaseq(dat0,mod1,mod0,controls=controls,n.sv=1)$sv
plot(sup_svseq, svseq,pch=19,col="blue")
###################################################

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
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_Log2_QN_id_no_grac.txt",  header =TRUE, sep = "\t")

# Colume Names.
# Column 1: PROBE_ID;
# Columns 2-52: UNK (0-56);
# PROBE_ID	UNK-0	UNK-1	UNK-2	UNK-3	UNK-5	UNK-6	UNK-8	UNK-9	UNK-10	UNK-11	UNK-12	UNK-13	UNK-14	UNK-15	UNK-16	UNK-17	UNK-19	UNK-20	UNK-21	UNK-23	UNK-24	UNK-25	UNK-26	UNK-27	UNK-29	UNK-30	UNK-31	UNK-32	UNK-33	UNK-34	UNK-35	UNK-36	UNK-37	UNK-38	UNK-39	UNK-40	UNK-41	UNK-42	UNK-43	UNK-45	UNK-46	UNK-47	UNK-48	UNK-49	UNK-50	UNK-51	UNK-52	UNK-53	UNK-54	UNK-55	UNK-56
# PROBE_ID	D_-123	D_0	D_4	D_21	D_116	D_185	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045

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

dat_data <- as.matrix(dat[,c(2:52)])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
row.names(dat_data) <- dat$PROBE_ID
# Alternative Column Names (by Day):
colID <- c("D_-123", "D_0", "D_4", "D_21", "D_116", "D_185", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045")
colnames(dat_data) <- colID
batch <- c("1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "3", "4", "4", "4", "4", "4", "4", "4", "4", "4", "4", "4", "4", "5", "5", "5", "5")
status <- c("HLT", "HRV", "HRV", "HRV", "HLT", "HLT", "HLT", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HLT", "HLT", "HLT", "HLT", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HRV", "HRV", "HRV", "HRV", "HRV")

my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmap <- heatmap.2(dat_data, Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "none", symkey = TRUE, keysize = 0.5)


############################################################
# C. REMOVING BATCH EFFECT WITH SVA
############################################################
# TO DO THIS I ASSIGNED STATUS TO EACH TIME POINT.
pheno <- as.data.frame(cbind(colID, batch, status))

# Using ComBat.
modsv <- model.matrix(~as.factor(status), data = pheno)
mod0sv <- model.matrix(~1, data = pheno)
combat_dat_data <- ComBat(dat = dat_data, batch = batch, mod = modsv, numCovs = NULL, par.prior = TRUE, prior.plots = FALSE)
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmap <- heatmap.2(combat_dat_data, Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "none", symkey = TRUE, keysize = 0.5)

# Output ComBat-ted Data.
write.table(combat_dat_data, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_Log2_QN_id_NoGrAC_ComBat.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


############################################################
# D. RE-QUANTILE NORMALIZATION AND RESCALING
############################################################
# Quantile Normalization.
total_matrix <- data.matrix(combat_dat_data)
total_matrix_qn <- normalizeQuantiles(total_matrix, tie=TRUE)
row.names(total_matrix_qn) <- row.names(dat_data)
colnames(total_matrix_qn) <- colnames(dat_data)

# Heatmap.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmap <- heatmap.2(total_matrix_qn, Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "none", symkey = TRUE, keysize = 0.5)

# Output Rescaled Data.
write.table(total_matrix_qn, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_Log2_QN_id_NoGrAC_ComBat_QN.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Rescaling between -1 and 1 for each row.
n <- length(total_matrix_qn[,1])
total_matrix_qnrs <- total_matrix_qn
for (i in 1:n) {
	total_matrix_qnrs[i,] <- 2 * (total_matrix_qnrs[i,] - min(total_matrix_qnrs[i,])) / (max(total_matrix_qnrs[i,]) - min(total_matrix_qnrs[i,])) - 1
}
total_matrix_qnrs[is.na(total_matrix_qnrs)] <- 0
total_heatmap <- heatmap.2(total_matrix_qnrs, Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,12), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)

# Output Rescaled Data.
write.table(total_matrix_qnrs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/UNK-NimblegenAll/A_UNK/20150609-Pass22_SVA_PostTSGr_GroupComp/Nimblegen_All_UNK_ShapiroWilk_Log2_QN_id_NoGrAC_ComBat_Rescaled.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


############################################################
# E. COMPARISON OF PRE- VS POST-COMBAT DATA BY HISTOGRAMS
############################################################
# Histogram.
hist(dat_data, nclass = 200, col = rgb(0,0,1,1/2), xlab = "METHOD 5: COMBAT", main = "COMPARISON OF UNK Nimblegen Data (PRE- AND POST-COMBAT)", ylim = c(0, 70000))
hist(total_matrix_qn, nclass = 200, col = rgb(0,1,0,1/4), add = TRUE)
legendtext <- c("Pre-ComBat", "Post-ComBat")
legend("topright", legendtext, fill = c(rgb(0,0,1,1/2), rgb(0,1,0,1/4)), cex =1.2)

# Density Plot.
allsignals <- c(as.vector(dat_data), as.vector(total_matrix_qn))
dendat <- data.frame(INTENSITY = allsignals, CATEGORIES = c(rep("Pre-ComBat", length(as.vector(dat_data))), rep("Post-Combat", length(as.vector(total_matrix_qn)))))
ggplot(dendat, aes(x = INTENSITY, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle("COMPARISON OF UNK Nimblegen Data (PRE- AND POST-COMBAT)")




