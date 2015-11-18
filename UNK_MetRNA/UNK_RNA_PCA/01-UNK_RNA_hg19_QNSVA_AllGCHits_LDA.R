####################################
# A. LIBRARIES USED.
####################################
library(MASS)
library(ggbiplot)
library("RColorBrewer")
library("gplots")
library("ggplot2")
library(scales)
library(gridExtra)
####################################


####################################
# B. RAW DATA LOADING AND PROCESSING.
####################################
# Load data.
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150328-Pass10-LDA_of_SVAHits/UNK_RNASeq_hg19_QNSVA_AllGCHits_PCA_PCx.txt", sep= "\t", header = TRUE, row.names = 1)
datt <- t(dat)
status <- c("HRV", "HRV", "HRV", "HLT", "HLT", "HLT", "HLT", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HLT", "HLT", "HLT", "HLT", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "HLT")
status_detailed <- c("HRV1", "HRV1", "HRV1", "HLT", "HLT", "HLT", "HLT", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HRV2", "HRV2", "HRV2", "HRV2", "HRV2", "HRV2", "HLT", "ADV1", "ADV1", "ADV1", "ADV1", "ADV1", "ADV1", "HLT", "HLT", "HLT", "HLT", "HLT", "ADV2", "ADV2", "ADV2", "ADV2", "ADV2", "ADV2", "HRV3", "HRV3", "HRV3", "HRV3", "HRV3", "HRV3", "HRV3", "HLT", "HLT")
status_dh <- c("INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY")
datts1 <- as.data.frame(cbind(datt, status_dh))
datts2 <- as.data.frame(cbind(datt, status))
datts3 <- as.data.frame(cbind(datt, status_detailed))

# apply LDA. 
datts1.lda <- lda(formula = status_dh ~ ., data = datts1, prior = c(1,1)/2)
datts2.lda <- lda(formula = status ~ ., data = datts2, prior = c(1,1,1,1)/4)
datts3.lda <- lda(formula = status_detailed ~ ., data = datts3, prior = c(1,1,1,1,1,1,1)/7)

prop.lda1 = datts1.lda$svd^2/sum(datts1.lda$svd^2)
prop.lda2 = datts2.lda$svd^2/sum(datts2.lda$svd^2)
prop.lda3 = datts3.lda$svd^2/sum(datts3.lda$svd^2)

plda1 <- predict(object = datts1.lda, newdata = datts1)
plda2 <- predict(object = datts2.lda, newdata = datts2)
plda3 <- predict(object = datts3.lda, newdata = datts3)

# Visualization (The Two Group one will not work in 2D plotting as it has only one dimension).
# Heatmaps.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hm1 <- heatmap.2(t(plda1$posterior), Rowv = FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(0.4,0.49,length=167), seq(0.49,0.51,length=167), seq(0.51,0.6,length=167)), keysize = 0.5, cexRow = 1.2)
hm2 <- heatmap.2(t(plda2$posterior), Rowv = FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(0.2,0.24,length=167), seq(0.24,0.26,length=167), seq(0.26,0.3,length=167)), keysize = 0.5, cexRow = 1.2)
hm3 <- heatmap.2(t(plda3$posterior), Rowv = FALSE, Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,7), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(0.09,0.14,length=167), seq(0.14,0.15,length=167), seq(0.15,0.2,length=167)), keysize = 0.5, cexRow = 1.2)

# Dot Plots.
dataset2 = data.frame(Status = status, lda = plda2$x)
p2 <- ggplot(dataset2) + geom_point(aes(lda.LD1, lda.LD2, colour = Status, shape = Status), size = 2.5) + 
  labs(x = paste("LD1 (", percent(prop.lda2[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda2[2]), ")", sep=""))
grid.arrange(p2)
 
dataset3 = data.frame(Status = status_detailed, lda = plda3$x)
p3 <- ggplot(dataset3) + geom_point(aes(lda.LD1, lda.LD2, colour = Status), size = 2.5) + 
  labs(x = paste("LD1 (", percent(prop.lda3[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda3[2]), ")", sep=""))
grid.arrange(p3)

# grid.arrange(p2, p3, nrow = 1)





####################################
# EXAMPLE PCA VS LDA.
####################################
require(MASS)
require(ggplot2)
require(scales)
require(gridExtra)
 
pca <- prcomp(iris[,-5],
              center = TRUE,
              scale. = TRUE) 
 
prop.pca = pca$sdev^2/sum(pca$sdev^2)
 
r <- lda(Species ~ ., 
           iris, 
           prior = c(1,1,1)/3)
 
prop.lda = r$svd^2/sum(r$svd^2)
 
plda <- predict(object = lda,
                newdata = iris)
 
dataset = data.frame(species = iris[,"Species"],
                     pca = pca$x, lda = plda$x)
 
p1 <- ggplot(dataset) + geom_point(aes(lda.LD1, lda.LD2, colour = species, shape = species), size = 2.5) + 
  labs(x = paste("LD1 (", percent(prop.lda[1]), ")", sep=""),
       y = paste("LD2 (", percent(prop.lda[2]), ")", sep=""))
 
p2 <- ggplot(dataset) + geom_point(aes(pca.PC1, pca.PC2, colour = species, shape = species), size = 2.5) +
  labs(x = paste("PC1 (", percent(prop.pca[1]), ")", sep=""),
       y = paste("PC2 (", percent(prop.pca[2]), ")", sep=""))
 
grid.arrange(p1, p2)




