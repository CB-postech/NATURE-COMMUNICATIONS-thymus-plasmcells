##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

library(Seurat)
### Load data
final_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/real_final_BCR_seurat_harmony.rds")

### Subset Plasma cells
plasma_subset <- subset(final_seurat, ident=c(3,5,6,8,7),invert=T)

### Normalization
library(scran)
library(scater)
library(SingleCellExperiment)
plasma_sce <- as.SingleCellExperiment(plasma_subset)
set.seed(123)
clusters <- quickCluster(plasma_sce, method="igraph")
table(clusters)
plasma_sce<- computeSumFactors(plasma_sce, cluster=clusters)
plasma_sce<- logNormCounts(plasma_sce)

### Feature Selection
dec <- modelGeneVar(plasma_sce)
hvg <- getTopHVGs(dec, n=1000)
plot(sizeFactors(plasma_sce),colSums(logcounts(plasma_sce)), log="xy",
     ylab="Library size (kilo)", xlab = "Size factor")

### Make Seurat object
library(Seurat)
seurat <- as.Seurat(plasma_sce)
seurat@assays$RNA@var.features <- hvg
seurat@assays$RNA@var.features <- seurat@assays$RNA@var.features[
  !grepl("Ig",seurat@assays$RNA@var.features)]
seurat@reductions$PCA <- NULL
seurat@reductions$TSNE <- NULL
seurat@reductions$UMAP <- NULL

seurat <- ScaleData(seurat)
i=13
PCA = i
seurat<- RunPCA(seurat, npcs = PCA,weight.by.var = F)
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
  FindClusters(resolution = 1.0) %>% 
  identity()
### Cell type annotation
plasma_seurat_13$Iso_annot <- "IgM"
plasma_seurat_13$Iso_annot[plasma_seurat_13$seurat_clusters %in% c(0,3)] <- "IgE1"
plasma_seurat_13$Iso_annot[plasma_seurat_13$seurat_clusters %in% c(4)] <- "IgE2"
plasma_seurat_13$Iso_annot[plasma_seurat_13$seurat_clusters %in% c(5)] <- "IgG"
plasma_seurat_13$Iso_annot[plasma_seurat_13$seurat_clusters %in% c(2,7)] <- "IgA"
plasma_seurat_13$Iso_annot <- factor(plasma_seurat_13$Iso_annot, levels = c("IgA","IgE1","IgE2","IgG","IgM"))


Iso_seurat <- plasma_seurat_13
Idents(Iso_seurat) <- Iso_seurat$Iso_annot

### Find Isotype specific DEG
Iso_marker <- FindAllMarkers(Iso_seurat)
Iso_seurat <- ScaleData(Iso_seurat,rownames(Iso_seurat))
sig_up_Iso_marker <- Iso_marker[Iso_marker$avg_logFC>0.5 & Iso_marker$p_val_adj<0.05,]
sig_up_Iso_marker$gene[grepl("Rp",sig_up_Iso_marker$gene)]

DoHeatmap(Iso_seurat,label = F,group.colors = c("#00CC00","#FF0000","#993800","#0000FF","#CCCC00"),
          unique(sig_up_Iso_marker$gene[!grepl("Rp",sig_up_Iso_marker$gene)]))+ 
  scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))+
  ggsave("I:/Postech_BCR/plasma_heatmap_version1.png",height = 10)
