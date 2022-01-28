##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

###Kidera Factor
#install.packages("Peptides")
### Load Data
library(Peptides)
final_BCR_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/real_final_BCR_seurat_harmony.rds")
filelist = list.files("D:/OneDrive - dgist.ac.kr/BCR_postech/merge/stcr_source/",full.names = T)
sapply(filelist, source)
ssetA4 = Load10xVDJ('raw_data/dataset2/A4/outs/', filtered = TRUE)
ssetH7 = Load10xVDJ('raw_data/H9/outs', filtered = TRUE)
ssetH8 = Load10xVDJ('raw_data/H10/outs', filtered = TRUE)
ssetA4= SetSample(ssetA4,"A4")
ssetH7= SetSample(ssetH7,"H7")
ssetH8= SetSample(ssetH8,"H8")

sset_B = MergeSamples(ssetA4,ssetH7) 
sset_B = MergeSamples(sset_B,ssetH8) 

### Merge BCR information and Transcriptome
sset_B = ImportSeurat(sset_B,final_BCR_seurat, version="v3")

Heavy_B <- subset(sset_B@contig, sset_B@contig$cellName %in% colnames(final_BCR_seurat)[grep("PB|PC",final_BCR_seurat$figure_annot)]&
                    sset_B@contig$chain=="IGH")
Heavy_B$isotype <- substring(Heavy_B$c_gene,1,4) 

#Kidera Factors
Kidera <- kideraFactors(Heavy_B$cdr3)
Kidera_df <- data.frame(matrix(unlist(Kidera), nrow=length(Kidera), byrow=T))
Minkowski_Kidera <- dist(Kidera_df, method = "minkowski",p=4)

df <- prcomp(Minkowski_Kidera)
df$cellname <- Heavy_B$cellName
df$isotype <- Heavy_B$Heavy_B
PC_df <- as.data.frame(df$x) 
PC_df$isotype<- Heavy_B$isotype
library(ggplot2)
ggplot(PC_df,
       aes(x=PC1, y=PC2, color = isotype)) +scale_colour_manual(values = c("#00CC00","#FF0000","#0000FF","#CCCC00","#B3AEB2"))+
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

write.csv(PC_df,"D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data//Kidera_factor_PCA.csv")
