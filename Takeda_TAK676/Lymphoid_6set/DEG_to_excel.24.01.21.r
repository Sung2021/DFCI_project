library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

dir <- "~/Desktop/DF/DFCI_Paweletz/2023_Takeda_PD1_NK/"
obj.srt = readRDS(paste0(dir,"rds/obj.cond6_merged.lymphoid.rds"))

## CD4>0,CD8A=0,CD8B=0 population    
# Find T cell populations
g1= "CD8A"
g2= "CD8B"
g3= "CD4"

obj = obj.srt

obj@meta.data[,g1] = obj@assays$RNA@data[g1,]
obj@meta.data[,g2] = obj@assays$RNA@data[g2,]
obj@meta.data[,g3] = obj@assays$RNA@data[g3,]

# CD4>0,CD8A=0,CD8B=0 
rs = obj@meta.data %>% filter(!!(sym(g1)) == 0 &!!(sym(g2)) == 0 & !!(sym(g3)) > 0) %>% rownames()
# number of cells table
obj@meta.data[rs,] %>% dplyr::select(orig.ident, sample) %>% table() %>% data.frame() %>% 
  tidyr::spread(sample, Freq) %>% DT::datatable()


# DEG:   
#   (212C+271C) vs (269C+273C+293C+298C)   
# and    
# (212T+271T) vs (269T+273T+293T+298T).   
# Could you please, provide excel files for the significant genes.   

# CD4 population 
obj.T = subset(obj.srt, cells = rs)
obj.T@meta.data$T_id = paste0(obj.T@meta.data$orig.ident, "_", obj.T@meta.data$sample)
obj.T@meta.data$T_id = factor(obj.T@meta.data$T_id)
# Rename levels 
levels(obj.T@meta.data$T_id) = c("293_C","293_T","298_C","298_T","212_C","212_T","269_C","269_T","271_C","271_T","273_C","273_T") 

# Prepare DEG analysis
group1 = c("212_C","271_C") # (212C+271C)
group2 = c("269_C","273_C","293_C","298_C") # (269C+273C+293C+298C)
obj.T@meta.data = obj.T@meta.data %>% mutate(compare = ifelse(T_id %in% group1, "212C+271C", ifelse(T_id %in% group2, "269C+273C+293C+298C", "other")))
# 
# # For FindMarkers with DESeq2 option, this is required 
obj.T_CTL = subset(obj.T, compare != "other") 

g1 = "212C+271C"
g2 = "269C+273C+293C+298C"
Idents(obj.T_CTL) = 'compare'
logfc=log2(1)

mks.CTL =FindMarkers(obj.T_CTL, ident.1 = g1, ident.2 = g2, 
                     logfc.threshold = logfc, test.use = "DESeq2", slot = "counts")

# which(is.na(mks))
mks = mks.CTL
mks = mks %>% filter(!is.na(p_val))
pval=0.05
fc=1.2
# mks$avg_log2FC = mks$avg_log2FC -log2(fc)

mks = mks %>% mutate(DE=ifelse(avg_log2FC >= log2(fc) & p_val_adj < pval, 'UP',
                               ifelse(avg_log2FC < -log2(fc) & p_val_adj < pval, 'DN','no_sig')))
# mks$DE %>% table()
# mks %>% filter(p_val == "NA")
mks$DE = factor(mks$DE, levels = c('UP','DN','no_sig'))

mks  %>% ggplot(aes(avg_log2FC, -log10(p_val_adj), color = DE)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = c('red', 'blue', 'grey')) +
  theme_classic() +
  geom_vline(xintercept = c(-log2(fc), log2(fc)), color = 'grey') +
  geom_hline(yintercept = -log10(0.05), color = 'grey') +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggeasy::easy_center_title() ## to center title

mks.CTL.DE = mks

mks.CTL.DE %>% filter(DE != ("no_sig")) %>% write.csv(paste0(dir,"data/Lymphoid_CD4_CTL_DEG_significant_genes.24.01.21.csv"))
mks.CTL.DE %>% write.csv(paste0(dir,"data/Lymphoid_CD4_CTL_DEG_all_genes.24.01.21.csv"))

###############################
## TAK ##
###############################
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

dir <- "~/Desktop/DF/DFCI_Paweletz/2023_Takeda_PD1_NK/"
obj.srt = readRDS(paste0(dir,"rds/obj.cond6_merged.lymphoid.rds"))

## CD4>0,CD8A=0,CD8B=0 population    
# Find T cell populations
g1= "CD8A"
g2= "CD8B"
g3= "CD4"

obj = obj.srt

obj@meta.data[,g1] = obj@assays$RNA@data[g1,]
obj@meta.data[,g2] = obj@assays$RNA@data[g2,]
obj@meta.data[,g3] = obj@assays$RNA@data[g3,]

# CD4>0,CD8A=0,CD8B=0 
rs = obj@meta.data %>% filter(!!(sym(g1)) == 0 &!!(sym(g2)) == 0 & !!(sym(g3)) > 0) %>% rownames()
# number of cells table
obj@meta.data[rs,] %>% dplyr::select(orig.ident, sample) %>% table() %>% data.frame() %>% 
  tidyr::spread(sample, Freq) %>% DT::datatable()

# DEG:   
#   (212C+271C) vs (269C+273C+293C+298C)   
# and    
# (212T+271T) vs (269T+273T+293T+298T).   
# Could you please, provide excel files for the significant genes.   

# CD4 population 
obj.T = subset(obj.srt, cells = rs)
obj.T@meta.data$T_id = paste0(obj.T@meta.data$orig.ident, "_", obj.T@meta.data$sample)
obj.T@meta.data$T_id = factor(obj.T@meta.data$T_id)
# Rename levels 
levels(obj.T@meta.data$T_id) = c("293_C","293_T","298_C","298_T","212_C","212_T","269_C","269_T","271_C","271_T","273_C","273_T") 


# Prepare DEGs 
group1 = c("212_T","271_T") # (212T+271T)
group2 = c("269_T","273_T","293_T","298_T") # (269T+273T+293T+298T)
obj.T@meta.data = obj.T@meta.data %>% mutate(compare = ifelse(T_id %in% group1, "212T+271T", ifelse(T_id %in% group2, "269T+273T+293T+298T", "other")))
# 
# # For FindMarkers with DESeq2 option, this is required 
obj.T_TAK = subset(obj.T, compare != "other") # 622 cells

g1 = "212T+271T"
g2 = "269T+273T+293T+298T"
Idents(obj.T_TAK) = 'compare'
logfc=log2(1)

mks.TAK =FindMarkers(obj.T_TAK, ident.1 = g1, ident.2 = g2, 
                     logfc.threshold = logfc, test.use = "DESeq2", slot = "counts")

mks = mks.TAK
mks = mks %>% filter(!is.na(p_val))
pval=0.05
fc=1.2
# mks$avg_log2FC = mks$avg_log2FC -log2(fc)

mks = mks %>% mutate(DE=ifelse(avg_log2FC >= log2(fc) & p_val_adj < pval, 'UP',
                               ifelse(avg_log2FC < -log2(fc) & p_val_adj < pval, 'DN','no_sig')))
# mks$DE %>% table()
# mks %>% filter(p_val == "NA")
mks$DE = factor(mks$DE, levels = c('UP','DN','no_sig'))
mks  %>% ggplot(aes(avg_log2FC, -log10(p_val_adj), color = DE)) +
  geom_point(size = 1, alpha = 0.5) +
  scale_color_manual(values = c('red', 'blue', 'grey')) +
  theme_classic() +
  geom_vline(xintercept = c(-log2(fc), log2(fc)), color = 'grey') +
  geom_hline(yintercept = -log10(0.05), color = 'grey') +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  ggeasy::easy_center_title() ## to center title

mks.TAK.DE = mks

mks.TAK.DE %>% filter(DE != ("no_sig")) %>% write.csv(paste0(dir,"data/Lymphoid_CD4_TAK_DEG_significant_genes.24.01.21.csv"))
mks.TAK.DE %>% write.csv(paste0(dir,"data/Lymphoid_CD4_TAK_DEG_all_genes.24.01.21.csv"))
