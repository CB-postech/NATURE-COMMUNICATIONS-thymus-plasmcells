##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

library(DropletUtils)
setwd("D:/OneDrive - dgist.ac.kr/BCR_postech/")
dir.H7 <- "D:/OneDrive - dgist.ac.kr/BCR_postech/raw_data/H7/raw_feature_bc_matrix/"
list.files(dir.H7)
set.seed(0)
H7 <- read10xCounts(dir.H7)
class(counts(H7))


br.out <- barcodeRanks(counts(H7))
sort(H7$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.H7 <- emptyDrops(counts(H7))
is.cell.H7 <- e.out.H7$FDR <= 0.05
sum(is.cell.H7, na.rm=TRUE)
colnames(H7) = colData(H7)$Barcode

cd.H7 = counts(H7)[,which(e.out.H7$FDR <= 0.05)]
save(cd.H7,file="counts_cells_dgCmatrix.RData")
load(file="counts_cells_dgCmatrix.RData")
# write.csv(as.matrix(cd.H7), file="counts.csv")


H7 <- H7[,which(e.out.H7$FDR < 0.05)]
library(scater)
library(EnsDb.Mmusculus.v79)

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(H7)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(H7)$CHR <- location
summary(location=="MT")
H7 <- calculateQCMetrics(H7, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
H7 <- runPCA(H7, use_coldata=TRUE)
# reducedDimNames(H7)
plotReducedDim(H7, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(H7$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H7$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(H7$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(H7$pct_counts_Mito <= 10 )] <- 1
filtering[which(H7$pct_counts_Mito > 10 )] <- 0
H7$filtering <- filtering
table(H7$filtering)

ggplot(as.data.frame(reducedDim(H7)), 
       aes(x=PC1, y=PC2, color = as.factor(H7$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(reducedDim(H7)),
       aes(x=PC1, y=PC2, color = H7$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
ggplot(as.data.frame(reducedDim(H7)),
       aes(x=PC1, y=PC2, color = H7$pct_counts_Mito)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
ggplot(as.data.frame(reducedDim(H7)),
       aes(x=PC1, y=PC2, color = H7$log10_total_features_by_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())




filtering <- numeric()
filtering[which( H8$pct_counts_Mito <= 10 )] <- 1
filtering[which( H8$pct_counts_Mito > 10 )] <- 0
H8$filtering <- filtering
table(H8$filtering)

ggplot(as.data.frame(reducedDim(H8)), 
       aes(x=PC1, y=PC2, color = as.factor(H8$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(reducedDim(H8)),
       aes(x=PC1, y=PC2, color = H8$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
ggplot(as.data.frame(reducedDim(H8)),
       aes(x=PC1, y=PC2, color = H8$pct_counts_Mito)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
ggplot(as.data.frame(reducedDim(H8)),
       aes(x=PC1, y=PC2, color = H8$log10_total_features_by_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())



##### Code for
##### Author: Eunmin Lee (elee4472@dgist.ac.kr)
##### Last Update: 2018/09/06

library(DropletUtils)
setwd("D:/OneDrive - dgist.ac.kr/BCR_postech/")
dir.H8 <- "D:/OneDrive - dgist.ac.kr/BCR_postech/raw_data/H8/raw_feature_bc_matrix/"
list.files(dir.H8)
set.seed(0)
H8 <- read10xCounts(dir.H8)
class(counts(H8))


br.out <- barcodeRanks(counts(H8))
sort(H8$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.H8 <- emptyDrops(counts(H8))
is.cell.H8 <- e.out.H8$FDR <= 0.05
sum(is.cell.H8, na.rm=TRUE)
colnames(H8) = colData(H8)$Barcode

cd.H8 = counts(H8)[,which(e.out.H8$FDR <= 0.05)]
save(cd.H8,file="counts_cells_dgCmatrix.RData")
load(file="counts_cells_dgCmatrix.RData")
# write.csv(as.matrix(cd.H8), file="counts.csv")




H8 <- H8[,which(e.out.H8$FDR < 0.05)]
library(scater)
library(EnsDb.Mmusculus.v79)

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(H8)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(H8)$CHR <- location
summary(location=="MT")
H8 <- calculateQCMetrics(H8, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
H8 <- runPCA(H8, use_coldata=TRUE)
# reducedDimNames(H8)
plotReducedDim(H8, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(H8$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(H8$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(H8$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")

filtering <- numeric()
filtering[which(H8$pct_counts_Mito <= 10)] <- 1
filtering[which(H8$pct_counts_Mito > 10)] <- 0
H8$filtering <- filtering
table(H8$filtering)

ggplot(as.data.frame(H8@reducedDims$PCA_coldata), 
       aes(x=PC1, y=PC2, color = as.factor(H8$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(H8@reducedDims$PCA_coldata),
       aes(x=PC1, y=PC2, color = H8$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())



dir.A4 <- "D:/OneDrive - dgist.ac.kr/BCR_postech/raw_data/dataset2/A4/raw_feature_bc_matrix/"
list.files(dir.A4)
set.seed(0)
A4 <- read10xCounts(dir.A4)
class(counts(A4))


br.out <- barcodeRanks(counts(A4))
sort(A4$total_counts,decreasing = TRUE)
# # Making a plot.
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")
abline(h=br.out$knee, col="dodgerblue", lty=2)
abline(h=br.out$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
       legend=c("knee", "inflection"))

set.seed(100)
e.out.A4 <- emptyDrops(counts(A4))
is.cell.A4 <- e.out.A4$FDR <= 0.05
sum(is.cell.A4, na.rm=TRUE)
colnames(A4) = colData(A4)$Barcode

cd.A4 = counts(A4)[,which(e.out.A4$FDR <= 0.05)]
save(cd.A4,file="A4_counts_cells_dgCmatrix.RData")
load(file="counts_cells_dgCmatrix.RData")
# write.csv(as.matrix(cd.A4), file="counts.csv")


A4 <- A4[,which(e.out.A4$FDR < 0.05)]
library(scater)
library(EnsDb.Mmusculus.v79)

location <- mapIds(EnsDb.Mmusculus.v79, keys=rowData(A4)$ID, 
                   column="SEQNAME", keytype="GENEID")

rowData(A4)$CHR <- location
summary(location=="MT")
A4 <- calculateQCMetrics(A4, feature_controls=list(Mito=which(location=="MT")))

# Dimensionality reduction plots
A4 <- runPCA(A4, use_coldata=TRUE)
# reducedDimNames(A4)
plotReducedDim(A4, use_dimred = "PCA_coldata")


# par(mfrow=c(1,3))

dev.off()

hist(A4$log10_total_counts, 
     breaks=100, col="grey80",
     xlab="Log-total UMI count")

hist(A4$log10_total_features_by_counts, 
     breaks=100, col="grey80",
     xlab="Log-total number of expressed features")

hist(A4$pct_counts_Mito, 
     breaks=100, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")


filtering <- numeric()
filtering[which(A4$log10_total_counts > 3.0 & A4$pct_counts_Mito <= 10 )] <- 1
filtering[which(A4$log10_total_counts <= 3.0 | A4$pct_counts_Mito > 10 )] <- 0
A4$filtering <- filtering
table(A4$filtering)

library(RColorBrewer)
ggplot(as.data.frame(reducedDim(A4)), 
       aes(x=PC1, y=PC2, color = as.factor(A4$filtering))) +
  geom_point() + 
  theme_bw() + 
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(reducedDim(A4)),
       aes(x=PC1, y=PC2, color = A4$log10_total_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

ggplot(as.data.frame(reducedDim(A4)),
       aes(x=PC1, y=PC2, color = A4$pct_counts_Mito)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())
ggplot(as.data.frame(reducedDim(A4)),
       aes(x=PC1, y=PC2, color = A4$log10_total_features_by_counts)) +
  geom_point() + scale_colour_gradientn(colours = rev(brewer.pal(11, "RdBu"))) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.title = element_blank(),
        legend.key = element_blank())

colnames(H7) <- paste0(colnames(H7),"_H7")
colnames(H8) <- paste0(colnames(H8),"_H8")
colnames(A4) <- paste0(colnames(A4),"_A4")


saveRDS(A4, file = "D:/OneDrive - dgist.ac.kr/BCR_postech/dropletfiltered_mitox/dropletfilterA4.rds")
saveRDS(H7, file = "D:/OneDrive - dgist.ac.kr/BCR_postech/dropletfiltered_mitox/dropletfilterH7.rds")
saveRDS(H8, file = "D:/OneDrive - dgist.ac.kr/BCR_postech/dropletfiltered_mitox/dropletfilterH8.rds")

