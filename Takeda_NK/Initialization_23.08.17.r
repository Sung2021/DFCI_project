## 2023 Takeda NK project
## Patient sample
## 30293
## 30298
## 2023.08.17
setwd('~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/')

## load required packages
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

patient = '30293'
patient = '30298'
dir= paste0('~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/raw_data/',patient,'/filtered_feature_bc_matrix')
list.files(paste0(dir))

## Create seurat object
obj.raw = Read10X(dir)
obj.srt = CreateSeuratObject(counts = obj.raw, project = patient)
obj.srt #36601 features across 23169 samples within 1 assay 

obj.srt@meta.data[1:3,]
obj.srt@meta.data %>% tail()

## add mitochondrial content
obj.srt[["percent.mt"]] <- PercentageFeatureSet(obj.srt, pattern = "^MT-")

## add sample info
add_sample_info <- function(meta_data, sample_names) {
  meta_data$sample = ''
  for (i in 1:length(sample_names)) {
    meta_data[grepl(paste0('-', i), rownames(meta_data)),]$sample = sample_names[i]
  }
  meta_data$sample = factor(meta_data$sample, levels = sample_names)
  return(meta_data)
}

# add sample info using function:
sample_names <- c('CTL', 'NKC', 'TAK', 'COM')
obj.srt@meta.data <- add_sample_info(obj.srt@meta.data, sample_names)

obj.srt@meta.data[1:3,]
obj.srt@meta.data$sample %>% table()

dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/'
obj.srt %>% saveRDS(paste0(dir,'rds/','Takeda.NK.', patient, '.23.08.17','.rds'))
obj.srt@meta.data %>% write.csv(paste0(dir,'/data/','Takeda.NK.', patient, '.meta.23.08.17','.csv'))


## scrublet
## prepare count mtx
# import rds
dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/rds'
list.files(dir)
obj.srt1= readRDS(paste0(dir,list.files(dir)[1]))
obj.srt2= readRDS(paste0(dir,list.files(dir)[2]))

## prepare count mtx using function
write_count_matrix_csv <- function(obj.srt, output_dir, file_name_prefix) {
  data <- obj.srt@assays$RNA@data
  df <- data %>% data.frame() %>% t()
  rownames(df) <- gsub(pattern = '-', replacement = '_', rownames(df))
  output_file <- file.path(output_dir, paste0(file_name_prefix, '.count.mtx.csv'))
  df %>% write.csv(file = output_file)
}

# write_count_matrix_csv
dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/raw_data'
output_directory <- dir
file_prefix1 <- 'Takeda.NK.30293'
write_count_matrix_csv(obj.srt=obj.srt1, output_directory, file_prefix1)
file_prefix2 <- 'Takeda.NK.30298'
write_count_matrix_csv(obj.srt=obj.srt2, output_directory, file_prefix2)

####### add scrublet info
update_seurat_object <- function(dir, patient, obj.srt) {
  df <- read.csv(paste0(dir, 'Takeda.NK.', patient, '.scr.csv'), row.names = 1)
  rownames(df) <- gsub(pattern = '\\.', replacement = '-', rownames(df))
  
  obj.srt[['doublet_score']] <- df$doublet_score
  obj.srt[['predicted_doublet']] <- df$predicted_doublet
  
  return(obj.srt)
}

# Usage
dir <- '~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/raw_data/'
obj.srt1 <- update_seurat_object(dir, patient='30293', obj.srt = obj.srt1)
obj.srt2 <- update_seurat_object(dir, patient='30298', obj.srt = obj.srt2)

dir <- '~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/rds/'
obj.srt1 %>% saveRDS(paste0(dir,'Takeda.NK.30293.23.08.17.rds'))
obj.srt2 %>% saveRDS(paste0(dir,'Takeda.NK.30298.23.08.17.rds'))

## Merge data
## in request form: analysis of MT30269 and MT30271
dir <- '~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/rds/'
obj.srt1= readRDS(paste0(dir,'Takeda.NK.30293.23.08.17.rds'))
obj.srt2= readRDS(paste0(dir,'Takeda.NK.30298.23.08.17.rds'))
## merge 
obj.srt = merge(obj.srt1, y = c(obj.srt2), add.cell.ids = c('TN30293','TN30298'))
## save merged object
dir <- '~/Desktop/DF/DFCI_Paweletz/2023_Takeda_NK/rds/'
obj.srt %>% saveRDS(paste0(dir,'Takeda.NK.30293.30298.23.08.17.rds'))
obj.srt@meta.data %>% write.csv(paste0(dir,'Takeda.NK.30293.30298.meta','.csv'))

