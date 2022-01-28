##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

library(Seurat)
### Load data
setwd("D:/OneDrive - dgist.ac.kr/BCR_postech/")
BCR_seurat_H <- readRDS("D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/real_final_BCR_seurat_harmony.rds")
filelist = list.files("D:/OneDrive - dgist.ac.kr/BCR_postech/merge/stcr_source/",full.names = T)
sapply(filelist, source)

### Load BCR Data
ssetA4 = Load10xVDJ('raw_data/dataset2/A4/outs/', filtered = TRUE)
ssetH7 = Load10xVDJ('raw_data/H9/outs', filtered = TRUE)
ssetH8 = Load10xVDJ('raw_data/H10/outs', filtered = TRUE)
ssetA4= SetSample(ssetA4,"A4")
ssetH7= SetSample(ssetH7,"H7")
ssetH8= SetSample(ssetH8,"H8")

sset_B = MergeSamples(ssetA4,ssetH7) 
sset_B = MergeSamples(sset_B,ssetH8) 

### Merge BCR information and Transcriptome
sset_B = ImportSeurat(sset_B,BCR_seurat_H, version="v3")

library(reshape2)
B_cdr3_mat = SetCDR3Mat(sset_B, cell_ind = colnames(BCR_seurat_H),
                        chain1 = c("IGH"), chain2 = c("IGK"))
B_cdr3_mat2 = SetCDR3Mat(sset_B, cell_ind = colnames(BCR_seurat_H),
                         chain1 = c("IGH"), chain2 = c("IGL"))

B_cdr3_pair_mat1 = SetCDR3PairMat(sset_B, B_cdr3_mat, chain1 = c("IGH"), chain2 = c("IGK"))
B_cdr3_pair_mat2 = SetCDR3PairMat(sset_B, B_cdr3_mat2, chain1 = c("IGH"), chain2 = c("IGL"))

### IGH_IGK pair
mat_cdr3_pair1 = B_cdr3_pair_mat1[[1]]
df_clonotype1 = B_cdr3_pair_mat1[[2]]
rownames(mat_cdr3_pair1)[rowSums(mat_cdr3_pair1)>0]
tail(df_clonotype1)

### IGH_IGL pair
mat_cdr3_pair2 = B_cdr3_pair_mat2[[1]]
colnames(mat_cdr3_pair2) = substring(colnames(mat_cdr3_pair2),10)
colnames(mat_cdr3_pair2) = as.numeric(colnames(mat_cdr3_pair2)) + 1187
colnames(mat_cdr3_pair2) = paste0("clonotype",colnames(mat_cdr3_pair2))
colnames(mat_cdr3_pair2) = as.factor(colnames(mat_cdr3_pair2))
df_clonotype2 = B_cdr3_pair_mat2[[2]]
df_clonotype2$clonotype = as.numeric(df_clonotype2$clonotype)
df_clonotype2$clonotype = df_clonotype2$clonotype + 1187
df_clonotype2$clonotype = paste0("clonotype",df_clonotype2$clonotype)
df_clonotype2$clonotype = as.factor(df_clonotype2$clonotype)


mat_cdr3_pair = cbind(mat_cdr3_pair1,mat_cdr3_pair2)
df_clonotype = rbind(df_clonotype1,df_clonotype2)
B_cdr3_pair_mat = list()
B_cdr3_pair_mat[[1]] =mat_cdr3_pair
B_cdr3_pair_mat[[2]] =df_clonotype

### Consider only heavy chain
heavy_pass <- subset(sset_B@contig, sset_B@contig$chain=="IGH")
### Filter cells without clonotype information
heavy_only_pass <- subset(heavy_pass,heavy_pass$cellName %in% rownames(mat_cdr3_pair)[rowSums(mat_cdr3_pair)==1])

heavy_only_pass<- subset(heavy_only_pass,heavy_only_pass$cellName %in% colnames(BCR_seurat_H)[BCR_seurat_H$figure_annot %in% c("PB","PC")])
plamamat_cdr3_pair <- mat_cdr3_pair[rownames(mat_cdr3_pair) %in% colnames(BCR_seurat_H)[BCR_seurat_H$figure_annot %in% c("PB","PC")],]


write.csv(table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHA",heavy_only_pass$c_gene)],]))/table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHA",heavy_only_pass$c_gene)],])>0)[2],
          "D:/OneDrive - dgist.ac.kr/BCR_postech/IGHA_plasma_clonesize_ratio.csv")
write.csv(table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHE",heavy_only_pass$c_gene)],]))/table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHE",heavy_only_pass$c_gene)],])>0)[2],
          "D:/OneDrive - dgist.ac.kr/BCR_postech/IGHE_plasma_clonesize_ratio.csv")
write.csv(table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHM",heavy_only_pass$c_gene)],]))/table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHM",heavy_only_pass$c_gene)],])>0)[2],
          "D:/OneDrive - dgist.ac.kr/BCR_postech/IGHM_plasma_clonesize_ratio.csv")
write.csv(table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHG",heavy_only_pass$c_gene)],]))/table(colSums(plamamat_cdr3_pair[heavy_only_pass$cellName[grepl("IGHG",heavy_only_pass$c_gene)],])>0)[2],
          "D:/OneDrive - dgist.ac.kr/BCR_postech/IGHG_plasma_clonesize_ratio.csv")
