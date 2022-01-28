##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

library(scater)
library(SingleCellExperiment)
setwd("D:/OneDrive - dgist.ac.kr/BCR_postech/")

### Load Data
seurat <- readRDS("Paper_final/data/real_final_BCR_seurat_noharmony.rds")

### Batch correction
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

DimPlot(seurat_H,label = T)

### Cell type annotation
seurat_H$figure_annot <- "Trans B"
seurat_H$figure_annot[seurat_H$seurat_clusters %in% c(0,1,4,6,7,9)] <- "PC"
seurat_H$figure_annot[seurat_H$seurat_clusters %in% c(8)] <- "PB"
seurat_H$figure_annot[seurat_H$seurat_clusters %in% c(2,11)] <- "Mem B"
seurat_H$figure_annot[seurat_H$seurat_clusters %in% c(3,10)] <- "Mat B"
seurat_H$figure_annot <- factor(seurat_H$figure_annot, levels = c("Trans B","Mat B","Mem B","PB","PC"))
final_seurat$figure_annot <- factor(final_seurat$figure_annot, levels=c("Trans B","Mat B","Mem B","PB","PC"))


celltype_seurat <- final_seurat
Idents(celltype_seurat) <- celltype_seurat$figure_annot
### Find celltype specific marker
celltype_marker <- FindAllMarkers(celltype_seurat)
upsig_celltype_marker <- celltype_marker[celltype_marker$avg_logFC>1& celltype_marker$p_val_adj<0.05,]

celltype_seurat <- ScaleData(celltype_seurat, rownames(celltype_seurat))

### Heatmap
DoHeatmap(celltype_seurat,label = F,
          unique(upsig_celltype_marker$gene))+ scale_fill_gradientn(colors = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100))+ggsave("D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/figure/final_all_seurat_heatmap_logFC1.png",height = 25)
