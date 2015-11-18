############################################################
# TO-DO LIST.
############################################################
# Using the actual batch information will break ComBat (as all RSVs belong to its own batch which leads to error "Lapack routine dgesv: system is exactly singular: U[15,15] = 0", and singletons will also lead to error "dim(X) must have a positive length").
# As a result, I used near actual batches (by fusing the first HRV and UNK-5,6,7,8 to one batch, and fuse RSV and the next batch of 5 healthy time points in one batch).




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
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNASeq_hg19_norep_QNsep_Unstitched.txt",  header =TRUE, sep = "\t")

# Colume Names.
# Column 1: GeneID;
# Columns 2-58: UNK (1-60);
# Columns 59-60: UNK-07 & 08 RiboZero;
# GeneID	U01RNA	U02RNA	U03RNA	U05RNA	U06RNA	U07RNA	U08RNA	U09RNA	U10RNA	U11RNA	U12RNA	U13RNA	U14RNA	U15RNA	U16RNA	U17RNA	U18RNA	U19RNA	U20RNA	U21RNA	U23RNA	U24RNA	U25RNA	U26RNA	U27RNA	U28RNA	U29RNA	U30RNA	U31RNA	U32RNA	U33RNA	U34RNA	U35RNA	U36RNA	U37RNA	U38RNA	U39RNA	U40RNA	U41RNA	U42RNA	U43RNA	U45RNA	U46RNA	U47RNA	U48RNA	U49RNA	U50RNA	U51RNA	U52RNA	U53RNA	U54RNA	U55RNA	U56RNA	U58RNA	U58_5RNA	U59RNA	U60RNA	U07RNARZ	U08RNARZ

# GeneID	D_0	D_4	D_21	D_116	D_185	D_186	D_255	D_289	D_290	D_292	D_294	D_297	D_301	D_307	D_311	D_322	D_329	D_369	D_380	D_400	D_476	D_532	D_546	D_602	D_615	D_616	D_618	D_620	D_625	D_630	D_647	D_679	D_680	D_683	D_688	D_694	D_700	D_711	D_735	D_796	D_840	D_912	D_944	D_945	D_948	D_959	D_966	D_984	D_1029	D_1030	D_1032	D_1038	D_1045	D_1051	D_1060	D_1109	D_1124	D_186RZ	D_255RZ

dat_data <- as.matrix(dat[,c(2:60)])
nrow <- length(dat_data[,1])
ncol <- length(dat_data[1,])
row.names(dat_data) <- dat$GeneID
# Alternative Column Names (by Day):
colID <- c("D_0", "D_4", "D_21", "D_116", "D_185", "D_186", "D_255", "D_289", "D_290", "D_292", "D_294", "D_297", "D_301", "D_307", "D_311", "D_322", "D_329", "D_369", "D_380", "D_400", "D_476", "D_532", "D_546", "D_602", "D_615", "D_616", "D_618", "D_620", "D_625", "D_630", "D_647", "D_679", "D_680", "D_683", "D_688", "D_694", "D_700", "D_711", "D_735", "D_796", "D_840", "D_912", "D_944", "D_945", "D_948", "D_959", "D_966", "D_984", "D_1029", "D_1030", "D_1032", "D_1038", "D_1045", "D_1051", "D_1060", "D_1109", "D_1124", "D_186RZ", "D_255RZ")
colnames(dat_data) <- colID
method <- c("PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "PA", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ", "RZ")
# Batch info for Poly-A vs Ribozero.
# batch <- c("1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2")
# Actual Batch info.
# batch <- c("1", "1", "2", "3", "4", "4", "5", "6", "6", "6", "6", "6", "6", "6", "6", "7", "7", "7", "7", "7", "8", "8", "8", "8", "9", "9", "9", "9", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "10", "11", "11", "11", "11", "11", "11", "11", "11", "11", "11", "11", "11", "12", "12")
# Simplified Batch Info, as RSV cannot be in its own batch, nor can any singletons.
batch <- c("1", "1", "1", "1", "1", "1", "1", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "2", "3", "3", "3", "3", "4", "4", "4", "4", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "5", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "6", "7", "7")
status <- c("HRV", "HRV", "HRV", "HLT", "HLT", "HLT", "HLT", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HLT", "HLT", "HLT", "HLT", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "HLT", "HLT", "HLT")

# Heatmap.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmap <- heatmap.2(dat_data, Rowv=TRUE, Colv=FALSE, col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = TRUE, keysize = 0.5)


############################################################
# C. STITCHING DATA WITH SVA
############################################################
# TO DO THIS I ASSIGNED STATUS TO EACH TIME POINT.
pheno <- as.data.frame(cbind(colnames(dat_data), batch, status))

modsv <- model.matrix(~as.factor(status), data = pheno)
mod0sv <- model.matrix(~1, data = pheno)
combat_dat_data <- ComBat(dat = dat_data, batch = batch, mod = modsv, numCovs = NULL, par.prior = TRUE, prior.plots = FALSE)
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmap <- heatmap.2(combat_dat_data, Rowv=TRUE, Colv=FALSE, col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = TRUE, keysize = 0.5)


############################################################
# D. RE-QUANTILE NORMALIZATION AND RESCALING
############################################################
# Quantile Normalization.
# No Linear adjusting with the two time points (186  and 255, PA and RZ) was done as it made it worse after ComBat.
total_matrix <- data.matrix(combat_dat_data)
total_matrix_qn <- normalizeQuantiles(total_matrix, tie=TRUE)
row.names(total_matrix_qn) <- row.names(dat_data)
colnames(total_matrix_qn) <- colnames(dat_data)

# Heatmap.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
total_heatmap <- heatmap.2(total_matrix_qn, Rowv=TRUE, Colv=FALSE, col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = TRUE, keysize = 0.5)

# Output Quantile Normalized Data.
write.table(total_matrix_qn, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNASeq_NearRealBatch_ComBat_QN.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)

# Rescaling between -1 and 1 for each row.
n <- length(total_matrix_qn[,1])
total_matrix_qnrs <- total_matrix_qn[,1:57]
for (i in 1:n) {
	total_matrix_qnrs[i,] <- 2 * (total_matrix_qnrs[i,] - min(total_matrix_qnrs[i,])) / (max(total_matrix_qnrs[i,]) - min(total_matrix_qnrs[i,])) - 1
}
total_matrix_qnrs[is.na(total_matrix_qnrs)] <- 0
total_heatmap <- heatmap.2(total_matrix_qnrs, Rowv=TRUE, Colv=FALSE, col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)
total_heatmap <- heatmap.2(total_matrix_qnrs, Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-1,-0.3,length=167), seq(-0.3,0.3,length=167), seq(0.3,1,length=167)), keysize = 0.5)

# Output Rescaled Data.
write.table(total_matrix_qnrs, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150324-Pass08-SVA-ActualBatches/UNK_RNASeq_NearRealBatch_ComBat_Rescaled.txt", sep = "\t", row.names = TRUE, col.names = TRUE, append = FALSE)


############################################################
# E. COMPARISON OF THE 4 STITCHING METHODS BY HISTOGRAMS
############################################################
# Histogram.
hist(dat_data, nclass = 200, col = rgb(0,0,1,1/2), xlab = "METHOD 5: COMBAT", main = "COMPARISON OF UNK RNA-SEQ (PRE- AND POST-COMBAT)")
hist(total_matrix_qn, nclass = 200, col = rgb(0,1,0,1/4), add = TRUE)
legendtext <- c("Pre-ComBat", "Post-ComBat")
legend("topright", legendtext, fill = c(rgb(0,0,1,1/2), rgb(0,1,0,1/4)), cex =1.2)

# Density Plot.
allsignals <- c(as.vector(dat_data), as.vector(total_matrix_qn))
dendat <- data.frame(INTENSITY = allsignals, CATEGORIES = c(rep("Pre-ComBat", length(as.vector(dat_data))), rep("Post-Combat", length(as.vector(total_matrix_qn)))))
ggplot(dendat, aes(x = INTENSITY, fill = CATEGORIES)) + geom_density(alpha = 0.5) + ggtitle("COMPARISON OF UNK RNA-SEQ (PRE- AND POST-COMBAT)")




