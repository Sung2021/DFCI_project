---
## title: "NK project: MN1 vs MN2"
## output: html_document
---
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo=F, fig.align = "center", message=F, warning=F, fig.height = 5, fig.width = 6)
```
## R Markdown
This is an R Markdown document.

The details concerning MN1 and MN2 are currently confidential.


```{r libraries, echo=FALSE}
# Load required libraries
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
```

## Import pre-processed data
```{r rds_data, echo=T}
## import data
dir= '~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/rds/'
obj.srt = readRDS(paste0(dir,'NK_MN1.MN2.23.07.27.rds'))
```
  
  
#### Display the number of cells in each sample 
```{r, echo=F}
obj.srt@meta.data %>% select(orig.ident) %>% table() %>% data.frame() %>% DT::datatable()
```


## Visualize UMAP results

Visualize samples in UMAP

Visualize clusters of res 0.2 in UMAP

Visualize cell_type in UMAP

Visualize main cell_type in UMAP

```{r visualization}
DimPlot(obj.srt, group.by = 'orig.ident') + theme(plot.title = element_blank())
DimPlot(obj.srt, group.by = 'RNA_snn_res.0.2',) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'RNA_snn_res.0.2', label = T, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'cell_type', label = F, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'cell_type', label = T, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'main_cell_type', label = T, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'main_cell_type', label = F, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 

```

# cell type marker genes (manually chosen)
```{r genes, echo=F}
gene_list <- c("CTSW", "S100A4", "GZMB", "XCL2", "CD52", "NKG7", "VWF", "CLEC14A", "PECAM1", "MCAM", "CD93", "COL1A1", "COL1A2", "COL3A1", "DCN", "COL6A3", "EPCAM", "ASCL1", "CDKN2A")

```

```{r , echo=T}
print(gene_list)
```

```{r, echo=F}
# Reverse the order of 'ident' levels
obj.srt@meta.data$rev.cluster <- factor(obj.srt@meta.data$RNA_snn_res.0.2, 
                                        levels = rev(levels(obj.srt@meta.data$RNA_snn_res.0.2)))
```

```{r, fig.height=5, fig.width=10}
# Create the violin plot with the reversed order of identities
VlnPlot(obj.srt, features = gene_list, group.by = "rev.cluster",
        stack = T, flip = F) + NoLegend()
```

## DEG analysis
```{r, echo=FALSE}
# Import mks files
dir='~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/data/mks/re_analysis/MN1.MN2/'
mks.nk= read.csv(paste0(dir,'NK','.mn1.2.deg.csv'), row.names = 1)
mks.huvec=read.csv(paste0(dir,'HUVEC','.mn1.2.deg.csv'), row.names = 1)
mks.fibro=read.csv(paste0(dir,'Fibro','.mn1.2.deg.csv'), row.names = 1)
mks.corl=read.csv(paste0(dir,'CORL','.mn1.2.deg.csv'), row.names = 1)

# plot volcano function
plot_volcano_text=function(mks=mks, celltype='NK'){
  mks$DE = factor(mks$DE, levels = c('UP','DN','no_sig'))
  mks$gene = rownames(mks)
  fc=1.2
  p=mks %>% 
    ggplot(aes(avg_log2FC, -log10(p_val_adj), color=DE)) + 
    geom_point(size=1, alpha=0.5) + 
    scale_color_manual(values = c('red','blue','grey')) +
    theme_classic() +
    geom_vline(xintercept = c(-log2(fc),log2(fc)), color='grey') +
    geom_hline(yintercept = -log10(0.05),color='grey') +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    ggtitle(paste0(celltype,' population')) +
    # Use geom_text with hjust argument to adjust label position
    geom_text(aes(label = labels), size = 2.5, show.legend = FALSE, hjust = 0, nudge_x = 0.01) +
    ggeasy::easy_center_title() ## to center title
  
  print(p)
}
```

### Volcano plot of Differentially Expressed Genes of MN2 vs MN1
```{r, echo=T, fig.height=5, fig.width=7}
plot_volcano_text(mks = mks.nk,celltype = 'NK')
plot_volcano_text(mks = mks.huvec,celltype = 'HUVEC')
plot_volcano_text(mks = mks.fibro,celltype = 'fibroblast')
plot_volcano_text(mks = mks.corl,celltype = 'CORL47')

```

  
UP regulated genes in NK from MN2 
```{r, echo=F}
mks.nk %>% filter(DE == 'UP') %>% select(gene) %>% pull()
```
  
DN regulated genes in NK from MN2 
```{r, echo=F}
mks.nk %>% filter(DE == 'DN') %>% select(gene) %>% pull()
```
  
UP regulated genes in HUVEC from MN2 
```{r, echo=F}
mks.huvec %>% filter(DE == 'UP') %>% select(gene) %>% pull()
```
  
DN regulated genes in HUVEC from MN2 
```{r, echo=F}
mks.huvec %>% filter(DE == 'DN') %>% select(gene) %>% pull()
```
  
UP regulated genes in Fibroblast from MN2 
```{r, echo=F}
mks.fibro %>% filter(DE == 'UP') %>% select(gene) %>% pull()
```
  
DN regulated genes in Fibroblast from MN2 
```{r, echo=F}
mks.fibro %>% filter(DE == 'DN') %>% select(gene) %>% pull()
```
  
UP regulated genes in CORL47 from MN2 
```{r, echo=F}
mks.corl %>% filter(DE == 'UP') %>% select(gene) %>% pull()
```
  
DN regulated genes in CORL47 from MN2 
```{r, echo=F}
mks.corl %>% filter(DE == 'DN') %>% select(gene) %>% pull()
```
  

# Violin plot of UP/DN regulated genes 
```{r, echo=F}
# Function to create Violin Plots by compare levels
plot_vln_by_compare <- function(main_cell_type, seurat_obj, gene_vector) {
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(compare = ifelse(main_cell_type == main_cell_type & orig.ident == 'MN1', paste0(main_cell_type, ' MN1'),
                            ifelse(main_cell_type == main_cell_type & orig.ident == 'MN2', paste0(main_cell_type, ' MN2'), 'other')))
  
  compare_levels <- c(paste0(main_cell_type, ' MN1'), paste0(main_cell_type, ' MN2'))
  
  Idents(seurat_obj) <- 'compare'
  VlnPlot(seurat_obj, features = gene_vector, idents = compare_levels, stack = TRUE, flip = TRUE) + NoLegend()
}

# Define the main_cell_type and gene_vector for each cell type
cell_types <- c('CORL47', 'HUVEC', 'Fibroblast', 'NK')

# Upregulated genes in MN2
gene_vectors_up <- list(
  c("ISG15", "IFI6", "EIF2AK2", "STAT1", "IFI27", "SAMD9", "PLSCR1", "HES4", "TIMP1"),
  c("ISG15", "IL11", "IFI6", "IFI27", "AREG", "STC1", "CD55", "IFIT2", "IFIT1",  "TFPI2"),
  c("ISG15", "IFI27", "IFI6", "MX1", "IFIT1", "IFI44L", "BST2", "EIF2AK2", "OAS3", "SAMHD1", "OAS2"),
  c("ISG15", "IFI6", "MX1", "ISG20", "IFIT1", "RSAD2", "IFIT3", "HERC5", "PLSCR1",  "TRIM22")
)
# Downregulated genes in MN2
gene_vectors_down <- list(
  c("PHGDH", "HIST1H1B", "HIST1H1E", "WARS", "EIF4EBP1", "PMF1", "SHMT2", "HIST1H1C", "MTHFD2"),
  c("IL32", "TAGLN", "TUBA1B", "PTX3", "TUBA1A", "NES", "TUBB", "FDPS", "ANXA2", "NUPR1"),
  c("FABP4", "WSB1", "SESN3", "GZMB", "XCL2", "SLC5A3", "ANGPT2", "CTSW", "CCL14", "NOP53", "IL32"),
  c("LTB", "CTSW", "CSF2", "NCR3", "GZMA", "S100A4", "ITGB2", "HSP90AB1", "LGALS1", "TNFRSF18")
)
```


```{r, echo=F}
# Create and arrange the plots using cowplot
# up
plots1 <- lapply(1:length(cell_types), function(i) {
  main_cell_type <- cell_types[i]
  gene_vector <- gene_vectors_up[[i]]
  plot_vln_by_compare(main_cell_type, obj.srt, gene_vector)
})

# Create and arrange the plots using cowplot
# down
plots2 <- lapply(1:length(cell_types), function(i) {
  main_cell_type <- cell_types[i]
  gene_vector <- gene_vectors_down[[i]]
  plot_vln_by_compare(main_cell_type, obj.srt, gene_vector)
})


```

UP regulated genes in MN2 compared to MN1
```{r, echo=TRUE, fig.height= 20, fig.width= 20}
# combine plots
cowplot::plot_grid(plotlist = plots1, nrow = 1)
```

DOWN regulated genes in MN3 compared to MN1
```{r, echo=TRUE, fig.height= 20, fig.width= 20}
# combine plots
cowplot::plot_grid(plotlist = plots2, nrow = 1)
```


```{r, echo=FALSE}
## up/dn regulated genes in each cell type in MN2
dir= '~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/data/mks/re_analysis/MN1.MN2/'
DEG.combined= read.csv(paste0(dir, 'DEG.combined.csv'), row.names = 1)
celltypes= levels(factor(DEG.combined$cell_type))

DE_condition= 'UP'
for(i in 1:length(celltypes)){assign(paste0(celltypes[i], '_vector'), DEG.combined %>% filter(cell_type== celltypes[i], DE== DE_condition) %>% pull(gene))}

ven_list = list(CORL47 = CORL47_vector,
                Fibroblast = Fibroblast_vector,
                HUVEC = HUVEC_vector,
                NK = NK_vector)

## draw venn diagram
ven_out <- VennDetail::venndetail(ven_list)
```

# UP regulated shared genes
```{r,fig.height= 5, fig.width= 5}
plot(ven_out)
plot(ven_out, type = "upset")
```

```{r}
## generate the shared gene list matrix
df = ven_out@result$Subset %>% table() %>% data.frame()
ven.info = df$. %>% factor() %>% levels() ## types of sharing among groups
df$gene =''
for(i in 1:nrow(df)){
  df[i, ]$gene = ven_out@result[ven_out@result$Subset == ven.info[i], ] %>% 
    dplyr::select(Detail) %>% 
    pull() %>% paste0(collapse = ',')
}
```

Information of shared genes
```{r, echo=TRUE}
DT::datatable(df)
```


```{r, echo=FALSE}
DE_condition= 'DN'
for(i in 1:length(celltypes)){assign(paste0(celltypes[i], '_vector'), DEG.combined %>% filter(cell_type== celltypes[i], DE== DE_condition) %>% pull(gene))}

ven_list = list(CORL47 = CORL47_vector,
                Fibroblast = Fibroblast_vector,
                HUVEC = HUVEC_vector,
                NK = NK_vector)

## draw venn diagram
ven_out <- VennDetail::venndetail(ven_list)
```

# DOWN regulated shared genes
```{r,fig.height= 5, fig.width= 5}
plot(ven_out)
plot(ven_out, type = "upset")

```

```{r}
## generate the shared gene list matrix
df = ven_out@result$Subset %>% table() %>% data.frame()
ven.info = df$. %>% factor() %>% levels() ## types of sharing among groups
df$gene =''
for(i in 1:nrow(df)){
  df[i, ]$gene = ven_out@result[ven_out@result$Subset == ven.info[i], ] %>% 
    dplyr::select(Detail) %>% 
    pull() %>% paste0(collapse = ',')
}
```

Information of shared genes
```{r, echo=TRUE}
DT::datatable(df)
```

### GSEA 
```{r}
dir= '~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/data/'
gsea.all = read.csv(paste0(dir,'mks/re_analysis/MN1.MN2/GESA.HALLMARK.all.csv'), row.names = 1)

```


```{r}
# function to draw GESA plot
plot_gsea_nes <- function(data, title) {
  data <- gsea.all[gsea.all$celltype == title,]
  p=ggplot(data, aes(reorder(ID, NES), NES)) +
    geom_col(aes(fill=p.adjust)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(title)) + 
    theme_classic() +
    scale_fill_gradient(low = 'red', high = '#E5E7E9') +
    theme(axis.text.x= element_text(size=5, face = 'bold'),
          axis.text.y= element_text(size=6, face = 'bold'), 
          axis.title =element_text(size=10)) +
    ggeasy::easy_center_title()
  print(p)
}
```


#### Individual GSEA plot
```{r, echo=TRUE, fig.height=6,fig.width=6}
# Usage
plot_gsea_nes(gsea.all, title= 'CORL47')
plot_gsea_nes(gsea.all, title= 'NK')
plot_gsea_nes(gsea.all, title= 'Fibroblast')
plot_gsea_nes(gsea.all, title= 'HUVEC')

```


#### Combined GSEA plot

```{r}
# descending order of NES
df =gsea.all %>% select(celltype,ID, NES) %>% tidyr::spread(ID, NES)
colsm = colSums(df[,2:ncol(df)])
colsm <- colsm[order(-colsm)]
# names(colsm)

# Get the order of the names in colsm
order_x <- names(colsm)

# Convert ID column to a factor with the desired order
gsea.all$ID <- factor(gsea.all$ID, levels = order_x)

# color palette for plot
custom_palette <- colorRampPalette(c("#2874A6", "white", "#FC4119"))(1000)

```

```{r, fig.height=7, fig.width=15}
# Plot the ggplot with the reordered x-axis
gsea.all %>%
  ggplot(aes(x = ID, y = celltype)) +
  geom_point(aes(size = -(p.adjust), color = NES)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_gradientn(colors = custom_palette) +
  labs(color = "NES", size = "FDR")
```

