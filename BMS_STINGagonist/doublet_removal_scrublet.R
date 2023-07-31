##### preparation of matrix in R #####
################################
library(dplyr)
library(Seurat)
setwd('~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/rds/')

############# data import ##################
obj.srt= readRDS('NK_Marco.MN1_8.23.01.14.rds')
obj.srt@meta.data[1:3,]
###### prepare the matrix #######
obj.srt@assays$RNA@data[1:3,1:3]
df = obj.srt@assays$RNA@data %>% data.frame(check.names = F) %>% t() 
rownames(df) = gsub(pattern = '-', replacement = '_' ,rownames(df))
df[1:3,1:3]
df %>% write.csv(file = '../data/MN1_MN8.count.matrix.csv')

###### import the scrublet result #######
scr=read.csv('../data/MN1_MN8.scr.csv', row.names = 1)
scr[1:3,1:2]
# correct the rownames of scr to match seurat object
rownames(scr) = gsub(pattern = '_1', replacement = '-1' ,rownames(scr))
scr[1:3,]
scr$predicted_doublet %>% table()
obj.srt[['doublet_score']] = scr$doublet_score
obj.srt[['predicted_doublet']] =scr$predicted_doublet

# save rds
obj.srt %>% saveRDS('NK_Marco.MN1_8.23.01.14.rds')

# subset 500-20000, 25, singlet only
# save rds
obj.tmp = subset(obj.srt, nCount_RNA >=500 & nCount_RNA <= 20000 & percent.mt <= 25 & predicted_doublet == 'False')
obj.tmp %>% saveRDS('NK_Marco.MN1_8.500_20000_25_singlet.23.07.27.rds')

# subset each comparison set
# function
subset_and_save <- function(obj.tmp, groups, file_name) {
  # Subset the Seurat object based on the orig.ident column
  obj.tmp2 <- subset(obj.tmp, orig.ident %in% groups)
  # Save the subsetted Seurat object as RDS file
  saveRDS(obj.tmp2, file = file_name)
}

# save rds:
grs_mn1_mn2 <- c('MN1', 'MN2')
subset_and_save(obj.tmp, grs_mn1_mn2, 'NK_MN1.MN2.23.07.27.rds')

grs_mn3_mn4 <- c('MN3', 'MN4')
subset_and_save(obj.tmp, grs_mn3_mn4, 'NK_MN3.MN4.23.07.27.rds')

grs_mn3_mn4 <- c('MN5', 'MN6')
subset_and_save(obj.tmp, grs_mn3_mn4, 'NK_MN5.MN6.23.07.27.rds')
