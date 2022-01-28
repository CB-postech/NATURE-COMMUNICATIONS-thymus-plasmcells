##### Code for
##### Author: Eun Seo Park (evergreen@dgist.ac.kr)

### Data Load
heavy_pass = read.table("D:/OneDrive - dgist.ac.kr/BCR_postech/changeO/0806/heavy_merge_seurat_prod_parse-select.tsv",sep = "\t",header = T)
merge_seurat <- readRDS("D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/real_final_BCR_seurat_harmony.rds")

### Split B cell, Plasma cell in Seurat
B_cell_id <- colnames(merge_seurat)[merge_seurat$figure_annot %in% c("Pre-B","Mat B")]
memB_cell_id <- colnames(merge_seurat)[merge_seurat$figure_annot %in% c("Mem B")]
### Pre_B_cell_id <- colnames(merge_seurat)[merge_seurat$simple_annot %in% c("Pre_B")]
Plasma_cell_id <- colnames(merge_seurat)[merge_seurat$figure_annot %in% c("PB","PC")]

### Split B cell, Plasma cell in BCR Dataset
heavy_pass$cell_id =paste0(substring(heavy_pass$cell_id,20),"-",substring(heavy_pass$cell_id,1,16))

heavy_pass$cell_type <- "Plasma"
heavy_pass$cell_type[heavy_pass$cell_id %in% B_cell_id] <- "B"
heavy_pass$cell_type[heavy_pass$cell_id %in% memB_cell_id] <- "memB"
### heavy_pass$cell_type[heavy_pass$cell_id %in% Pre_B_cell_id] <- "Pre_B"
heavy_pass$vcall_celltype <- paste0(heavy_pass$v_call_10x,"_",heavy_pass$cell_type)

heavy_pass$sample <- substring(heavy_pass$cell_id,1,2)
heavy_pass$vcall_sample <- paste0(heavy_pass$v_call_10x,"_",heavy_pass$sample)

B_pass = subset(heavy_pass, as.character(heavy_pass$cell_id) %in% B_cell_id)
"%ni%" <- Negate("%in%")
rm_matureB_pass = subset(B_pass,B_pass$c_call %ni% c("IGHA*03","IGHE*01","IGHG1*01","IGHG2A*01","IGHG2B*02","IGHG3*01"))

Plasma_pass = subset(heavy_pass, as.character(heavy_pass$cell_id) %in% Plasma_cell_id)
filtered_Plasma_pass = subset(Plasma_pass, Plasma_pass$vcall_sample %in% rm_matureB_pass$vcall_sample)

filtered_Plasma_pass$cdr <- paste0(filtered_Plasma_pass$cdr1,filtered_Plasma_pass$cdr2)
write.csv(gsub("\\.","-",filtered_Plasma_pass[filtered_Plasma_pass$vcall_sample=="IGHV2-9*02_A4","cdr"]),"D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/weblogo_SHM_cdr12.csv")

test_pass$cdr <- paste0(test_pass$cdr1,test_pass$cdr2)
######subset_allpass 
filtered_pass = rbind(rm_matureB_pass,filtered_Plasma_pass)

filtered_pass$simple_celltype <- filtered_pass$cell_type
# filtered_pass$simple_celltype[filtered_pass$simple_celltype =="Pre_B"] <- "B"
filtered_pass$check <- paste0(filtered_pass$v_call_10x,"_",filtered_pass$sample,"_",filtered_pass$simple_celltype)


filtered_pass$cdr_SHM <-"SHM"
filtered_pass$cdr1_SHM <-"SHM"
filtered_pass$cdr2_SHM <-"SHM"
check_pass <- data.frame()
filtered_pass$cdr <- paste0(filtered_pass$cdr1,filtered_pass$cdr2)
for (i in unique(filtered_pass$vcall_sample)){
  test_pass <- subset(filtered_pass, filtered_pass$vcall_sample==i)
  test_pass$cdr_SHM[test_pass$cell_type=="B"] <- "B"
  test_pass$cdr1_SHM[test_pass$cell_type=="B"] <- "B"
  test_pass$cdr2_SHM[test_pass$cell_type=="B"] <- "B"
  B_pass <- subset(test_pass, test_pass$cell_type=="B")
  test_pass$cdr_SHM[test_pass$cdr_SHM=="SHM"&test_pass$cdr %in% as.character(B_pass$cdr)] <-"noSHM"
  test_pass$cdr1_SHM[test_pass$cdr1_SHM=="SHM"&test_pass$cdr1 %in% as.character(B_pass$cdr1)] <-"noSHM"
  test_pass$cdr2_SHM[test_pass$cdr2_SHM=="SHM"&test_pass$cdr2 %in% as.character(B_pass$cdr2)] <-"noSHM"
  
  check_pass <- rbind(test_pass, check_pass)
}
B_pass <- subset(check_pass, check_pass$cell_type=="B")
B_pass$cdr <- paste0(B_pass$cdr1,test_pass$cdr2)
for (i in unique(B_pass$vcall_sample)){
  tmp_pass <- subset(B_pass,B_pass$vcall_sample==i)
  if(length(unique(tmp_pass$cdr))>1){
    print(i)
    print(tmp_pass$cdr)
    print(table(tmp_pass$cdr))
  }
}
library(dplyr)
library(alakazam)

check_pass3 <- data.frame()
for (i in unique(filtered_pass$vcall_sample)){
  test_pass <- subset(filtered_pass, filtered_pass$vcall_sample==i)
  test_pass$cdr_SHM[test_pass$cell_type=="B"] <- "B"
  test_pass$cdr1_SHM[test_pass$cell_type=="B"] <- "B"
  test_pass$cdr2_SHM[test_pass$cell_type=="B"] <- "B"
  test_pass$cdr <- paste0(test_pass$cdr1,test_pass$cdr2)
  B_pass <- subset(test_pass, test_pass$cell_type=="B")
  for (j in 1:length(rownames(test_pass))){
    
    if(length(unique(B_pass$cdr))==1){
      test_pass$cdr1_SHM_f[j]<- seqDist(as.character(test_pass$cdr1)[j], unique(as.character(B_pass$cdr1)))
      test_pass$cdr2_SHM_f[j] <- seqDist(as.character(test_pass$cdr2)[j], unique(as.character(B_pass$cdr2)))
      test_pass$cdr_SHM_f[j] <- seqDist(as.character(test_pass$cdr)[j], unique(as.character(B_pass$cdr)))
    }else if(length(unique(B_pass$cdr))==0){
      print("No CDR")
    }else{
      print("Over one CDR")
      test_pass$cdr1_SHM_f[j] <- "No"
      test_pass$cdr2_SHM_f[j] <- "No"
      test_pass$cdr_SHM_f[j] <- "No"
    }
  }
  check_pass3 <- rbind(test_pass, check_pass3)
}
Plasma_pass3 <- subset(check_pass3, check_pass3$cell_type=="Plasma")
Plasma_pass2 <- subset(Plasma_pass3, Plasma_pass3$cdr_SHM_f !="No")



num_CDR1 <-Plasma_pass2 %>% group_by(cdr1_SHM_f) %>% summarise(count=n()) %>% mutate(perc=count/sum(count))
num_CDR2 <-Plasma_pass2 %>% group_by(cdr2_SHM_f) %>% summarise(count=n()) %>% mutate(perc=count/sum(count))
num_CDR <-Plasma_pass2 %>% group_by(cdr_SHM_f) %>% summarise(count=n()) %>% mutate(perc=count/sum(count))

write.csv(num_CDR1,"D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/all_SHM_CDR1_raw_data.csv")
write.csv(num_CDR2,"D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/all_SHM_CDR2_raw_data.csv")
write.csv(num_CDR,"D:/OneDrive - dgist.ac.kr/BCR_postech/Paper_final/data/all_SHM_CDR_raw_data.csv")



