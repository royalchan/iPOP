####################################
A. LIBRARIES USED.
####################################
library(devtools)
# install_github("ggbiplot", "vqv")
library(ggbiplot)
library("RColorBrewer")
library("gplots")

# ggbiplot manual.
ggbiplot(pcobj, choices = 1:2, scale = 1, pc.biplot =
  TRUE, obs.scale = 1 - scale, var.scale = scale, groups =
  NULL, ellipse = FALSE, ellipse.prob = 0.68, labels =
  NULL, labels.size = 3, alpha = 1, var.axes = TRUE, circle
  = FALSE, circle.prob = 0.69, varname.size = 3,
  varname.adjust = 1.5, varname.abbrev = FALSE, ...)
Arguments
pcobj
an object returned by prcomp() or princomp()
choices
which PCs to plot
scale
covariance biplot (scale = 1), form biplot (scale = 0). When scale = 1, the inner product between the variables approximates the covariance and the distance between the points approximates the Mahalanobis distance.
obs.scale
scale factor to apply to observations
var.scale
scale factor to apply to variables
pc.biplot
for compatibility with biplot.princomp()
groups
optional factor variable indicating the groups that the observations belong to. If provided the points will be colored according to groups
ellipse
draw a normal data ellipse for each group?
ellipse.prob
size of the ellipse in Normal probability
labels
optional vector of labels for the observations
labels.size
size of the text used for the labels
alpha
alpha transparency value for the points (0 = TRUEransparent, 1 = opaque)
circle
draw a correlation circle? (only applies when prcomp was called with scale = TRUE and when var.scale = 1)
var.axes
draw arrows for the variables?
varname.size
size of the text for variable names
varname.adjust
adjustment factor the placement of the variable names, >= 1 means farther from the arrow
varname.abbrev
whether or not to abbreviate the variable names
####################################


####################################
A. RAW DATA LOADING AND PROCESSING.
####################################
# Load data
dat <- read.table("/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150325-Pass09-PCA_of_SVAHits/UNK_RNASeq_hg19_QNSVA_AllGCHits.txt", sep= "\t", header = TRUE, row.names = 1)
datt <- t(dat)
status <- c("HRV", "HRV", "HRV", "HLT", "HLT", "HLT", "HLT", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HLT", "HLT", "HLT", "HLT", "HLT", "ADV", "ADV", "ADV", "ADV", "ADV", "ADV", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HRV", "HLT", "HLT")
status_detailed <- c("HRV1", "HRV1", "HRV1", "HLT", "HLT", "HLT", "HLT", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "RSV", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HLT", "HRV2", "HRV2", "HRV2", "HRV2", "HRV2", "HRV2", "HLT", "ADV1", "ADV1", "ADV1", "ADV1", "ADV1", "ADV1", "HLT", "HLT", "HLT", "HLT", "HLT", "ADV2", "ADV2", "ADV2", "ADV2", "ADV2", "ADV2", "HRV3", "HRV3", "HRV3", "HRV3", "HRV3", "HRV3", "HRV3", "HLT", "HLT")
status_dh <- c("INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "HEALTHY", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "INFECTED", "HEALTHY", "HEALTHY")


# apply PCA - scale. = TRUE is highly advisable, but default is FALSE. 
datt.pca <- prcomp(datt, center = TRUE, scale. = TRUE)

print(datt.pca)        
Standard deviations:
 [1] 1.614200e+01 9.876824e+00 8.914177e+00 8.520986e+00 8.177953e+00 7.060773e+00 6.736652e+00 6.278533e+00 5.868498e+00 5.306742e+00 5.176158e+00 4.772815e+00 4.555560e+00
[14] 4.289933e+00 4.151638e+00 3.909617e+00 3.814209e+00 3.637962e+00 3.443077e+00 3.368894e+00 3.278701e+00 3.277419e+00 3.112629e+00 3.070006e+00 2.921203e+00 2.910103e+00
[27] 2.870969e+00 2.869081e+00 2.815206e+00 2.719669e+00 2.676196e+00 2.655699e+00 2.587029e+00 2.531125e+00 2.448879e+00 2.431376e+00 2.418706e+00 2.333704e+00 2.289876e+00
[40] 2.259164e+00 2.170975e+00 2.158189e+00 2.122049e+00 2.092384e+00 2.037824e+00 2.032658e+00 1.993353e+00 1.896944e+00 1.880970e+00 1.795553e+00 1.712612e+00 1.631103e+00
[53] 1.598226e+00 1.467380e+00 1.442625e+00 1.385014e+00 7.168490e-15

Rotation: ... (omitted)

# plot method
plot(datt.pca, type = "l")

# Output Transposed Scores Matrix.
write.table(t(datt.pca$x), file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150325-Pass09-PCA_of_SVAHits/UNK_RNASeq_hg19_QNSVA_AllGCHits_PCA_PCx.txt", sep = "\t", row.names = TRUE, col.names = TRUE)
# Output Loadings Matrix.
write.table(datt.pca$rotation, file = "/Users/Rui/Rui Folders/Royal-Rock/Work/Snyder Lab Work Selected/Results/20130610-UNK45to60-RNA-Seq/Analyses/20150325-Pass09-PCA_of_SVAHits/UNK_RNASeq_hg19_QNSVA_AllGCHits_PCA_PCrotation.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Plotting Transposed Scores Matrix.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(t(datt.pca$x), method = "euclidean"), method = "ward.D")
# Min is -41.31807; Max is 31.79488.
total_heatmap <- heatmap.2(t(datt.pca$x), Rowv=as.dendrogram(hc_ts), Colv=FALSE, dendrogram = "row", col = my_palette, margins=c(5,5), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-42,-17,length=167), seq(-17,8,length=167), seq(8,33,length=167)), keysize = 0.5)
total_heatmap <- heatmap.2(t(datt.pca$x), Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,5), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-42,-17,length=167), seq(-17,8,length=167), seq(8,33,length=167)), keysize = 0.5)

# Plotting Loadings Matrix.
my_palette <- colorRampPalette(c("turquoise", "gray12", "gold"))(n = 500)
hc_ts <- hclust(dist(t(datt.pca$rotation), method = "euclidean"), method = "ward.D")
# Min is -0.3482751; Max is 0.3775531.
total_heatmap <- heatmap.2(datt.pca$rotation, Colv=as.dendrogram(hc_ts), Rowv=FALSE, dendrogram = "column", col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-0.35,-0.1,length=167), seq(-0.1,0.13,length=167), seq(0.13,0.38,length=167)), keysize = 0.5)
total_heatmap <- heatmap.2(datt.pca$rotation, Rowv=FALSE, Colv=FALSE, col = my_palette, margins=c(5,6), density.info = "none", trace = "none", scale = "none", symkey = FALSE, breaks = c(seq(-0.35,-0.1,length=167), seq(-0.1,0.13,length=167), seq(0.13,0.38,length=167)), keysize = 0.5)


# summary method
summary(datt.pca)
Importance of components:
                           PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8    PC9    PC10    PC11    PC12    PC13    PC14    PC15    PC16    PC17    PC18    PC19
Standard deviation     16.1420 9.87682 8.91418 8.52099 8.17795 7.06077 6.73665 6.27853 5.8685 5.30674 5.17616 4.77282 4.55556 4.28993 4.15164 3.90962 3.81421 3.63796 3.44308
Proportion of Variance  0.2262 0.08468 0.06898 0.06303 0.05805 0.04328 0.03939 0.03422 0.0299 0.02445 0.02326 0.01977 0.01801 0.01598 0.01496 0.01327 0.01263 0.01149 0.01029
Cumulative Proportion   0.2262 0.31086 0.37984 0.44287 0.50092 0.54420 0.58360 0.61781 0.6477 0.67215 0.69541 0.71519 0.73320 0.74918 0.76414 0.77741 0.79004 0.80152 0.81181
                          PC20    PC21    PC22    PC23    PC24    PC25    PC26    PC27    PC28    PC29    PC30    PC31    PC32    PC33    PC34    PC35    PC36    PC37    PC38
Standard deviation     3.36889 3.27870 3.27742 3.11263 3.07001 2.92120 2.91010 2.87097 2.86908 2.81521 2.71967 2.67620 2.65570 2.58703 2.53113 2.44888 2.43138 2.41871 2.33370
Proportion of Variance 0.00985 0.00933 0.00932 0.00841 0.00818 0.00741 0.00735 0.00715 0.00715 0.00688 0.00642 0.00622 0.00612 0.00581 0.00556 0.00521 0.00513 0.00508 0.00473
Cumulative Proportion  0.82167 0.83100 0.84032 0.84873 0.85691 0.86432 0.87167 0.87883 0.88597 0.89285 0.89927 0.90549 0.91161 0.91742 0.92298 0.92819 0.93332 0.93840 0.94313
                          PC39    PC40    PC41    PC42    PC43   PC44   PC45    PC46    PC47    PC48    PC49   PC50    PC51    PC52    PC53    PC54    PC55    PC56      PC57
Standard deviation     2.28988 2.25916 2.17098 2.15819 2.12205 2.0924 2.0378 2.03266 1.99335 1.89694 1.88097 1.7956 1.71261 1.63110 1.59823 1.46738 1.44263 1.38501 7.168e-15
Proportion of Variance 0.00455 0.00443 0.00409 0.00404 0.00391 0.0038 0.0036 0.00359 0.00345 0.00312 0.00307 0.0028 0.00255 0.00231 0.00222 0.00187 0.00181 0.00167 0.000e+00
Cumulative Proportion  0.94768 0.95211 0.95620 0.96024 0.96415 0.9679 0.9716 0.97514 0.97859 0.98172 0.98479 0.9876 0.99013 0.99244 0.99466 0.99653 0.99833 1.00000 1.000e+00

# Plot Cumulative Proportion.
cp <- c("0.2262", "0.31086", "0.37984", "0.44287", "0.50092", "0.54420", "0.58360", "0.61781", "0.6477", "0.67215", "0.69541", "0.71519", "0.73320", "0.74918", "0.76414", "0.77741", "0.79004", "0.80152", "0.81181", "0.82167", "0.83100", "0.84032", "0.84873", "0.85691", "0.86432", "0.87167", "0.87883", "0.88597", "0.89285", "0.89927", "0.90549", "0.91161", "0.91742", "0.92298", "0.92819", "0.93332", "0.93840", "0.94313", "0.94768", "0.95211", "0.95620", "0.96024", "0.96415", "0.9679", "0.9716", "0.97514", "0.97859", "0.98172", "0.98479", "0.9876", "0.99013", "0.99244", "0.99466", "0.99653", "0.99833", "1.00000", "1.000e+00")
pcname <- paste("PC_", c(1:57), sep = "")
plot(cp, axes = F, xlab="", ylab="", type = "b")
axis(2)
axis(1, at = 1:length(pcname), labels=pcname)
title(main = "", xlab="Principal Components", ylab = "Cumulative Proportion of Variance")
# Plot Proportion of Variance.
pv <- c("0.2262", "0.08468", "0.06898", "0.06303", "0.05805", "0.04328", "0.03939", "0.03422", "0.0299", "0.02445", "0.02326", "0.01977", "0.01801", "0.01598", "0.01496", "0.01327", "0.01263", "0.01149", "0.01029", "0.00985", "0.00933", "0.00932", "0.00841", "0.00818", "0.00741", "0.00735", "0.00715", "0.00715", "0.00688", "0.00642", "0.00622", "0.00612", "0.00581", "0.00556", "0.00521", "0.00513", "0.00508", "0.00473", "0.00455", "0.00443", "0.00409", "0.00404", "0.00391", "0.0038", "0.0036", "0.00359", "0.00345", "0.00312", "0.00307", "0.0028", "0.00255", "0.00231", "0.00222", "0.00187", "0.00181", "0.00167", "0.000e+00")
plot(pv, axes = F, xlab="", ylab="", ylim = c(0,0.3), type = "b")
axis(2)
axis(1, at = 1:length(pcname), labels=pcname)
title(main = "", xlab="Principal Components", ylab = "Proportion of Variance")

# Plot the first 2 PCs -- 4 groups.
g <- ggbiplot(datt.pca, obs.scale = 1, var.scale = 1, 
              groups = status, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

# Plot the first 2 PCs -- 7 groups.
g <- ggbiplot(datt.pca, obs.scale = 1, var.scale = 1, 
              groups = status_detailed, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)

# Plot the first 2 PCs -- 2 groups.
g <- ggbiplot(datt.pca, obs.scale = 1, var.scale = 1, 
              groups = status_dh, ellipse = TRUE, 
              circle = TRUE, var.axes = FALSE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g)




