---
title: "scGSVA workflow"
author: "Sung Rye Park"
date: "`r format(Sys.Date())`"
output:  
  rmdformats::readthedown: 
    code_folding: show  
    number_sections: TRUE
    toc_depth: 6
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, results = 'asis')
options(warn = F)

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(DT)
```

# scGSVA template

This is a guideline to use scGSVA analysis

R code that generated the plot.


# Import data   
```{r}
# dir='~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/rds/'
# obj.srt = readRDS(paste0(dir,'NK_MN3.MN4.23.07.27.rds'))
dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda/'
obj.srt = readRDS(paste0(dir,'rds/MT30271.30273.23.07.31.rds'))
```

## Overview of data
```{r}
# obj.srt@meta.data %>% head() %>% DT::datatable(options = list(scrollX = T))
obj.srt@meta.data %>% head() %>% 
  gt::gt() 
```

# scGSVA 

scGSVA: GSVA for single cell RNA seq analysis.  
https://rpubs.com/bioguo/803521.  
url : devtools::install_github("guokai8/scGSVA").  


## Import reference gene set

### Azizi data
```{r}
# Import geneset as a long format
# dir='~/Desktop/DF/DFCI_Paweletz/gmt/'
# celltypes= read.csv(paste0(dir,'azizi.gmt.long_format.csv'), row.names = 1)
```

### HALLMARK data
```{r}
# Import geneset as a long format
hallmark <- msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, gene_symbol)
colnames(hallmark) = c('term','gene')
celltypes = hallmark
```

```{r}
# Define the class "Annot"
setClass("Annot", representation(
  species = "character",
  anntype = "character",
  keytype = "character",
  annot = "data.frame"
))

# Create an instance of the "Annot" class
df <- data.frame(
  GeneID = celltypes$gene,
  GOALL = celltypes$term,
  Annot = celltypes$term
)

species <- "human"
anntype <- "GO"
keytype <- "SYMBOL"

result <- new("Annot",
              species = species,
              anntype = anntype,
              keytype = keytype,
              annot = df)
```


## Run scGSVA 
```{r}
set.seed(123)   
library(scGSVA)

# run scGSVA
gsva.out <- scgsva(obj.srt, result, verbose=FALSE, cores = 8)

# save gsva.out if necessary
# dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda/'
# gsva.out %>% saveRDS(paste0(dir, 'rds/MT30271.30273.scGSVA.HALLMARK.23.08.21.rds'))
```

## Heatmap
```{r, fig.width=20, fig.height=10, eval=F}
# transpose output
df =t(gsva.out@gsva)

# generate annotation
df.anno=obj.srt@meta.data[,c('res','sample','main_cell_type')] %>% 
  arrange(main_cell_type, res,sample, .by_group=TRUE)

# set color
my.color=c(colorRampPalette(colors = c("#2874A6","white"))(20), 
           colorRampPalette(colors = c("white","#D35400"))(80))

# pheatmap::pheatmap(gsva.out@gsva[rownames(df.anno),], cluster_rows = F, 
#                    cluster_cols = F,
#                    show_rownames = F,
#                    annotation_row = df.anno,
#                    col = my.color)

# draw heatmap
pheatmap::pheatmap(t(gsva.out@gsva[rownames(df.anno),]), cluster_rows = F, 
                   cluster_cols = F,
                   show_rownames = T, show_colnames = F,
                   annotation_col = df.anno,
                   col = my.color)
```


## Heatmap with column gaps
```{r,fig.width=20, fig.height=10}

# generate annotation
df.anno=obj.srt@meta.data[,c('res','sample','main_cell_type')] %>% 
  arrange(main_cell_type, res,sample, .by_group=TRUE)

# set color
my.color=c(colorRampPalette(colors = c("#2874A6","white"))(20), 
           colorRampPalette(colors = c("white","#D35400"))(80))


# Calculate the gap widths based on the annotation
unique_cells <- unique(df.anno$main_cell_type)
main_cell_type_freq <- df.anno$main_cell_type %>% table() %>% data.frame() %>% pull(Freq)
accumulated_sum <- cumsum(main_cell_type_freq)

# Draw heatmap with gaps
pheatmap::pheatmap(
  t(gsva.out@gsva[rownames(df.anno),]), 
  cluster_rows = FALSE, 
  cluster_cols = FALSE,
  show_rownames = TRUE, 
  show_colnames = FALSE,
  annotation_col = df.anno,
  col = my.color,
  gaps_col = accumulated_sum
)
```
