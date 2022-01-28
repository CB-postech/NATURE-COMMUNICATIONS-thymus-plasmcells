##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

library(scater)
library(SingleCellExperiment)
setwd("D:/OneDrive - dgist.ac.kr/BCR_postech/")

### Load DropletUtils Filtered Singlecellexperiment object
H7 <- readRDS(file = "D:/OneDrive - dgist.ac.kr/BCR_postech/dropletfiltered_mitox/dropletfilterH7.rds")
H8 <- readRDS(file = "D:/OneDrive - dgist.ac.kr/BCR_postech/dropletfiltered_mitox/dropletfilterH8.rds")
A4 <- readRDS(file = "D:/OneDrive - dgist.ac.kr/BCR_postech/dropletfiltered_mitox/dropletfilterA4.rds")

### Filtered Singlecellexperiment Object
H7_filtered <- H7[,H7$filtering==1]
H8_filtered <- H8[,H8$filtering==1]
A4_filtered <- A4[,A4$filtering==1]

### Merge after droplet QC
library(Seurat)
sce_merge <- cbind(H7_filtered,H8_filtered,A4_filtered,deparse.level = 1)

### filtering genes
library(Matrix)
keep_feature <- rowSums(counts(sce_merge) != 0) > 3
sce_merge <- sce_merge[keep_feature, ]

### Normalization
library(scran)
library(scater)
set.seed(123)
clusters <- quickCluster(sce_merge, method="igraph")
table(clusters)
sce_merge <- computeSumFactors(sce_merge, cluster=clusters)
sce_merge <- logNormCounts(sce_merge)

### Feature Selection
dec <- modelGeneVar(sce_merge)
hvg <- getTopHVGs(dec, n=1000)


plot(sizeFactors(sce_merge),colSums(logcounts(sce_merge)), log="xy",
     ylab="Library size (kilo)", xlab = "Size factor")

### Make Seurat object
library(Seurat)
seurat <- as.Seurat(sce_merge)
seurat@assays$RNA@var.features <- hvg
seurat <- ScaleData(seurat)
PCA = 20
seurat<- RunPCA(seurat, npcs = PCA, weight.by.var = FALSE)
seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA,  force.recalc = TRUE)
seurat<- FindClusters(seurat, reduction.type="pca", dims= 1:PCA, save.SNN = TRUE, force.recalc = TRUE)
seurat<- RunTSNE(seurat, dims = 1:PCA, do.fast = T,seed.use = 42, perplexity=100)
seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42, n.neighbors = 30)

#Add metadata
seurat@meta.data$sample <- substring(rownames(seurat@meta.data), 20,22)
seurat@meta.data$dataset <- "Dataset1"
seurat@meta.data$dataset[seurat@meta.data$sample=="A4"] <- "Dataset2"

### Convert ensembl gene name to gene symol

load("ensemblGenes2019-07-11.RData")
rownames(seurat@assays$RNA@counts) <- uniquifyFeatureNames(rownames(seurat@assays$RNA@counts),
                                                           ensemblGenes[rownames(seurat@assays$RNA@counts),"external_gene_name"])
rownames(seurat@assays$RNA@data) <- uniquifyFeatureNames(rownames(seurat@assays$RNA@data),
                                                         ensemblGenes[rownames(seurat@assays$RNA@data),"external_gene_name"])
rownames(seurat@assays$RNA@scale.data) <- uniquifyFeatureNames(rownames(seurat@assays$RNA@scale.data),
                                                               ensemblGenes[rownames(seurat@assays$RNA@scale.data),"external_gene_name"])

### Remove T cell and macrophage
rm_T_mac_subset <- subset(first_seurat, ident=c(14,15),invert=T)

### Normalization
library(scran)
library(scater)
library(SingleCellExperiment)
rm_T_mac_sce <- as.SingleCellExperiment(rm_T_mac_subset)
set.seed(123)
clusters <- quickCluster(rm_T_mac_sce, method="igraph")
table(clusters)
rm_T_mac_sce<- computeSumFactors(rm_T_mac_sce, cluster=clusters)
rm_T_mac_sce<- logNormCounts(rm_T_mac_sce)

### Feature Selection
dec <- modelGeneVar(rm_T_mac_sce)
hvg <- getTopHVGs(dec, n=1000)
plot(sizeFactors(rm_T_mac_sce),colSums(logcounts(rm_T_mac_sce)), log="xy",
     ylab="Library size (kilo)", xlab = "Size factor")

### Make Seurat Object
library(Seurat)
seurat <- as.Seurat(rm_T_mac_sce)
seurat@assays$RNA@var.features <- hvg
seurat@assays$RNA@var.features <- seurat@assays$RNA@var.features[
  !grepl("Ig",seurat@assays$RNA@var.features)]
seurat@reductions$PCA <- NULL
seurat@reductions$TSNE <- NULL
seurat@reductions$UMAP <- NULL

seurat <- ScaleData(seurat)

PCA = 15
seurat<- RunPCA(seurat, npcs = PCA, weight.by.var = FALSE)
seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA,  force.recalc = TRUE)
seurat<- FindClusters(seurat, reduction.type="pca", dims= 1:PCA, save.SNN = TRUE, force.recalc = TRUE)
seurat<- RunTSNE(seurat, dims = 1:PCA, do.fast = T,seed.use = 42, perplexity=100)
seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42, n.neighbors = 30)

### Batch Correction
library(harmony)
library(dplyr)
set.seed(123456)
seurat_H <- seurat %>% 
  RunHarmony("dataset", plot_convergence = F, max.iter.harmony = 100)

seurat_H <- seurat_H %>% 
  RunUMAP(reduction = "harmony", dims = 1:PCA,seed.use = 42) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:PCA) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

rm_mito_high_subset <- subset(second_seurat, ident=c(8),invert=T)
library(scran)
library(scater)
library(SingleCellExperiment)
rm_mito_high_sce <- as.SingleCellExperiment(rm_mito_high_subset)
set.seed(123)
clusters <- quickCluster(rm_mito_high_sce, method="igraph")
table(clusters)
rm_mito_high_sce<- computeSumFactors(rm_mito_high_sce, cluster=clusters)
rm_mito_high_sce<- logNormCounts(rm_mito_high_sce)
dec <- modelGeneVar(rm_mito_high_sce)
hvg <- getTopHVGs(dec, n=1000)

plot(sizeFactors(rm_mito_high_sce),colSums(logcounts(rm_mito_high_sce)), log="xy",
     ylab="Library size (kilo)", xlab = "Size factor")

### Make Seurat Object
library(Seurat)
seurat <- as.Seurat(rm_mito_high_sce)
seurat@assays$RNA@var.features <- hvg
seurat@assays$RNA@var.features <- seurat@assays$RNA@var.features[
  !grepl("Ig",seurat@assays$RNA@var.features)]
seurat@reductions$PCA <- NULL
seurat@reductions$TSNE <- NULL
seurat@reductions$UMAP <- NULL

seurat <- ScaleData(seurat)

PCA = 25
seurat<- RunPCA(seurat, npcs = PCA,weight.by.var = F)
seurat<- FindNeighbors(seurat, reduction="pca", dims= 1:PCA,  force.recalc = TRUE)
seurat<- FindClusters(seurat, reduction.type="pca", dims= 1:PCA, save.SNN = TRUE, force.recalc = TRUE)
seurat<- RunTSNE(seurat, dims = 1:PCA, do.fast = T,seed.use = 42, perplexity=100)
seurat<- RunUMAP(seurat, dims = 1:PCA, seed.use = 42, n.neighbors = 30)

