---
title: "Takeda PD1 + NK merged: Lymphoid T cell request"
author: "Sung Rye Park"
date: "`r format(Sys.Date())`"
output:  
  rmdformats::robobook: 
    code_folding: hide 
    number_sections: FALSE
    toc_depth: 6
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=F, warning=FALSE, message=FALSE, results = 'asis')
options(warn = F)
```

```{r, echo=FALSE}
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
```

# Lymphoid Data   

This dataset constitutes a refined subset derived from lymphoid annotation performed by Elena.   
```{r}
dir <- "~/Desktop/DF/DFCI_Paweletz/2023_Takeda_PD1_NK/"
obj.srt = readRDS(paste0(dir,"rds/obj.cond6_merged.lymphoid.rds"))
```

## Number of Lymphoid cells  
```{r}
df1 =obj.srt@meta.data %>% select(orig.ident) %>% table() %>% data.frame()
df2 =obj.srt@meta.data %>% select(orig.ident, sample) %>% table() %>% 
  data.frame() %>% tidyr::spread(sample, Freq)
df = cbind(df1, df2[,2:3])
colnames(df)[2] = "Total"
df %>% DT::datatable()
```

## Number of CD8+T cells    

**T cell definition : CD8 > 0 & CD4=0**   

A T cell is defined by the presence of CD8 and the absence of CD4.     

**CD8A > 0 or CD8B > 0 & CD4 =0**     
```{r}
## CD8A+ population ratio    
# Find T cell populations
g1= "CD8A"
g2= "CD8B"
g3= "CD4"

obj = obj.srt

obj@meta.data[,g1] = obj@assays$RNA@data[g1,]
obj@meta.data[,g2] = obj@assays$RNA@data[g2,]
obj@meta.data[,g3] = obj@assays$RNA@data[g3,]

# CD8A > 0 or CD8B > 0 & CD4 =0 
rs = obj@meta.data %>% filter(!!(sym(g1)) > 0 |!!(sym(g2)) > 0 & !!(sym(g3)) == 0) %>% rownames()
# number of cells table
obj@meta.data[rs,] %>% dplyr::select(orig.ident, sample) %>% table() %>% data.frame() %>% 
  tidyr::spread(sample, Freq) %>% DT::datatable()
```


**CD8A > 0 & CD8B > 0 & CD4 =0**    
```{r}
## CD8A+ population ratio    
# Find T cell populations
g1= "CD8A"
g2= "CD8B"
g3= "CD4"

obj = obj.srt

obj@meta.data[,g1] = obj@assays$RNA@data[g1,]
obj@meta.data[,g2] = obj@assays$RNA@data[g2,]
obj@meta.data[,g3] = obj@assays$RNA@data[g3,]

# CD8A > 0 & CD8B >0 & CD4 =0 
rs = obj@meta.data %>% filter(!!(sym(g1)) > 0 & !!(sym(g2)) > 0 & !!(sym(g3)) == 0) %>% rownames()
# number of cells table
obj@meta.data[rs,] %>% dplyr::select(orig.ident, sample) %>% table() %>% data.frame() %>% 
  tidyr::spread(sample, Freq) %>% DT::datatable()
```

# Requested table    


```{r}
# Subset T cells only
rs = obj@meta.data %>% filter(!!(sym(g1)) > 0 |!!(sym(g2)) > 0 & !!(sym(g3)) == 0) %>% rownames()
obj.T = subset(obj.srt, cells = rs)
obj.T@meta.data$T_id = paste0(obj.T@meta.data$orig.ident, "_", obj.T@meta.data$sample)
obj.T@meta.data$T_id = factor(obj.T@meta.data$T_id)
levels(obj.T@meta.data$T_id) = c("293_C","293_T","298_C","298_T","212_C","212_T",
                                 "269_C","269_T","271_C","271_T","273_C","273_T")
```

```{r}
# Find genes 
genes = read.csv(paste0(dir, "info/CD8A_B_Features_all genes_20231222.csv"), row.names = 1, check.names = F)
genes[,"298_C Average"] = as.numeric(genes[,"298_C Average"])
gs=rownames(genes) # 708 genes   
```


**Requested order of samples**    
   
"212_C","271_C","269_C","273_C","293_C","298_C",   
"212_T","271_T","269_T","273_T","293_T","298_T"    

t-test :   
     1. between R vs NR in CTL     
     2. between R vs NR in TAK     
```{r}
### average expression of selected genes
Idents(obj.T) = 'T_id'
df = AverageExpression(obj.T, features = gs) 
df.rna = df$RNA

cols =c("212_C","271_C","269_C","273_C","293_C","298_C",
        "212_T","271_T","269_T","273_T","293_T","298_T")
df.rna =df.rna[,cols] %>% data.frame(check.names = F)

# t-테스트를 수행하는 함수
perform_t_test <- function(row_values) {
  group1 <- row_values[1:2]  # 212_C, 271_C
  group2 <- row_values[3:6]  # 269_C, 273_C, 293_C, 298_C
  result <- t.test(group1, group2)
  return(result$p.value)
}

# 새로운 열을 담을 벡터 초기화
t_test_results <- as.numeric(nrow(df.rna))

# 각 행에 대해 t-테스트를 수행하고 결과를 벡터에 저장
for (i in seq(nrow(df.rna))) {
  t_test_results[i] <- perform_t_test(df.rna[i, c("212_C", "271_C", "269_C", "273_C", "293_C", "298_C")])
}

# 새로운 열 추가
df.rna$t_test_C <- t_test_results


# 새로운 열을 담을 벡터 초기화
t_test_results <- as.numeric(nrow(df.rna))

# 각 행에 대해 t-테스트를 수행하고 결과를 벡터에 저장
for (i in seq(nrow(df.rna))) {
  t_test_results[i] <- perform_t_test(df.rna[i, c("212_T", "271_T", "269_T", "273_T", "293_T", "298_T")])
}

# 새로운 열 추가
df.rna$t_test_T <- t_test_results
df.rna_backup = df.rna
```



```{r}
# log2FC columns
# DEG first
# 1. R/NR in CTL
# 2. R/NR in TAK

#obj.T@meta.data$T_id %>% unique()
R_C = paste0(c("212","271"),"_C")
R_T = paste0(c("212","271"),"_T")
NR_C = paste0(c("273","269","298","293"),"_C")
NR_T = paste0(c("273","269","298","293"),"_T")
obj.T@meta.data = obj.T@meta.data %>% mutate(DE_id = ifelse(T_id %in% R_C, "R_C", ifelse(T_id %in% R_T, "R_T", ifelse(T_id %in% NR_C, "NR_C","NR_T"))))
#obj.T@meta.data$DE_id %>% table()

obj = subset(obj.T,features= rownames(genes)) # downsizing features with provided geneset

Idents(obj) = "DE_id"
logfc=log2(1)
mks_C =FindMarkers(obj, ident.1 = 'R_C', ident.2 = 'NR_C', logfc.threshold = logfc, 
                 test.use="wilcox", features = rownames(genes))

mks_T =FindMarkers(obj, ident.1 = 'R_T', ident.2 = 'NR_T', logfc.threshold = logfc, 
                   test.use="wilcox", features = rownames(genes))

mks_C$Comparison = "RoverNR_CTL"
mks_T$Comparison = "RoverNR_TAK"

# Add DEG information 
mks_C = mks_C %>% mutate(DE = ifelse(avg_log2FC > 0 & p_val <= 0.05, "UP", ifelse(avg_log2FC < 0 & p_val <= 0.05, "DN", "no_sig")))
mks_C$DE %>% table()
mks_T = mks_T %>% mutate(DE = ifelse(avg_log2FC > 0 & p_val <= 0.05, "UP", ifelse(avg_log2FC < 0 & p_val <= 0.05, "DN", "no_sig")))
rownames(genes) %>% length()

rs= intersect(rownames(mks_C), intersect(rownames(mks_T), rownames(genes)))
#rs %>% length() # 692 
```


```{r}
df.rna= df.rna_backup
rs= intersect(rownames(mks_C), intersect(rownames(mks_T), rownames(genes)))
df.rna = df.rna[rs,]
df.rna$log2FC_CTL = mks_C[rownames(df.rna),"avg_log2FC"]
df.rna$log2FC_TAK = mks_T[rownames(df.rna),"avg_log2FC"]
#df.rna$DE_info_CTL = mks_C[rownames(df.rna),"DE_info"]
#df.rna$DE_info_TAK = mks_T[rownames(df.rna),"DE_info"]
```





## Average expression and t-test between R vs NR    

```{r}
df.rna %>% DT::datatable(extensions = "Buttons", options = list(dom="Bfrtip", buttons=c("csv","excel"), pageLength=10))
```

## Differentially Expressed genes   


```{r}
cat("UP regulated genes in CTL form R/NR comparison" , "\n")
cat("\n")
mks_C[rownames(df.rna),] %>% filter(DE== "UP") %>% rownames()
cat("\n")
cat("DOWN regulated genes in CTL form R/NR comparison" , "\n")
cat("\n")
mks_C[rownames(df.rna),] %>% filter(DE== "DN") %>% rownames()
cat("\n")
cat("\n")
cat("UP regulated genes in TAK form R/NR comparison" , "\n")
cat("\n")
mks_T[rownames(df.rna),] %>% filter(DE== "UP") %>% rownames()
cat("\n")
cat("DOWN regulated genes in TAK form R/NR comparison" , "\n")
cat("\n")
mks_T[rownames(df.rna),] %>% filter(DE== "DN") %>% rownames()
cat("\n")
```

```{r}
up1=mks_C[rownames(df.rna),] %>% filter(DE== "UP") %>% rownames()
dn1=mks_C[rownames(df.rna),] %>% filter(DE== "DN") %>% rownames()
up2=mks_T[rownames(df.rna),] %>% filter(DE== "UP") %>% rownames()
dn2=mks_T[rownames(df.rna),] %>% filter(DE== "DN") %>% rownames()

ven_list = list(UP_CTL=up1, UP_TAK=up2)
ven_out <- VennDetail::venndetail(ven_list)

plot(ven_out, type = "upset")

ven_list = list(DOWN_CTL=dn1, DOWN_TAK=dn2)
ven_out <- VennDetail::venndetail(ven_list)

plot(ven_out, type = "upset")

```


```{r}
df.genes = genes[rownames(df.rna),]
df.tmp = cbind(df.genes[,1:12], df.rna[,1:12])
```


```{r}
plot_scatter_correlation <- function(df=df.tmp, cols_wanted) {
  # 원하는 열 선택
  cols_selected <- colnames(df)[grep(cols_wanted, colnames(df))]
  df_selected <- df[, cols_selected]
  
  # 간단한 산점도 및 회귀선
  ggplot(df_selected, aes(x = !!sym(names(df_selected)[1]), y = !!sym(names(df_selected)[2]))) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "grey") +
    annotate("text", x = max(df_selected[, 1]), y = max(df_selected[, 2]),
             label = paste("Cor:", round(cor(df_selected[, 1], df_selected[, 2]), 3))) +
    ggtitle(cols_wanted)
}

# 함수 호출
#colnames(df.tmp)

plot_scatter_correlation(df.tmp, "212_C")
plot_scatter_correlation(df.tmp, "271_C")
plot_scatter_correlation(df.tmp, "273_C")
plot_scatter_correlation(df.tmp, "269_C")
plot_scatter_correlation(df.tmp, "293_C")
#plot_scatter_correlation(df.tmp, "298_C")

cols_wanted = "212_C"
df= df.tmp
cols_selected <- colnames(df)[grep(cols_wanted, colnames(df))]
df_selected <- df[, cols_selected]

# 간단한 산점도 및 회귀선
ggplot(df_selected, aes(x = !!sym(names(df_selected)[1]), y = !!sym(names(df_selected)[2]))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  annotate("text", x = max(df_selected[, 1]), y = max(df_selected[, 2]),
           label = paste("Cor:", round(cor(df_selected[, 1], df_selected[, 2]), 3))) +
    ggtitle(cols_wanted)

plot_scatter_correlation(df.tmp, "212_T")
plot_scatter_correlation(df.tmp, "271_T")
plot_scatter_correlation(df.tmp, "273_T")
plot_scatter_correlation(df.tmp, "269_T")
plot_scatter_correlation(df.tmp, "293_T")
plot_scatter_correlation(df.tmp, "298_T")
```

