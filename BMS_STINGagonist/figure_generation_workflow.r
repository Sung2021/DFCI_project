##### NK project reanalysis #####
################################
library(dplyr)
library(Seurat)
library(ggplot2)
setwd('~/Desktop/DF/DFCI_Barbie/DFCI_Barbie_NK_Marco/rds/')

## MN3.MN4 dataset
# obj.srt %>% saveRDS('NK_MN3.MN4.23.07.27.rds')
obj.srt = readRDS('NK_MN3.MN4.23.07.27.rds')
## perform default analysis
perform_default_analysis <- function(obj.srt, n_features = 3000, n_pcs = 30, 
                                     dims_for_neighbors = 1:20, 
                                     resolutions = c(0.1, 0.2, 0.4, 0.8), 
                                     umap_dims = 1:10) {
  # Step 1: Find variable features
  obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst', nfeatures = n_features)
  
  # Step 2: Scale and normalize data
  all_genes <- rownames(obj.srt)
  obj.srt <- ScaleData(obj.srt, features = all_genes)
  obj.srt <- NormalizeData(obj.srt)
  
  # Step 3: Run PCA
  obj.srt <- RunPCA(obj.srt, features = VariableFeatures(object = obj.srt), npcs = n_pcs)
  
  # Step 4: Find neighbors
  obj.srt <- FindNeighbors(obj.srt, dims = dims_for_neighbors)
  
  # Step 5: Find clusters
  obj.srt <- FindClusters(obj.srt, resolution = resolutions)
  
  # Step 6: Run UMAP
  obj.srt <- RunUMAP(obj.srt, dims = umap_dims)
  
  # Return the Seurat object with analysis results
  return(obj.srt)
}

# apply
obj.srt <- perform_default_analysis(obj.srt)

# Visualize PCA results
DimPlot(obj.srt, reduction = 'pca')
ElbowPlot(obj.srt) ## decide number of PC to use
# Visualize UMAP results
DimPlot(obj.srt)
DimPlot(obj.srt, group.by = 'orig.ident') + theme(plot.title = element_blank())
DimPlot(obj.srt, group.by = 'RNA_snn_res.0.2',) + theme(plot.title = element_blank()) 
# save meta data
obj.srt@meta.data %>% write.csv('../data/mks/MN3.MN4.meta.23.07.27.csv')

find_and_save_markers <- function(obj.srt, cluster_id, logfc_threshold = 1.2, 
                                  test_method = 'wilcox', min_percent = 0.25) {
  Idents(obj.srt) = cluster_id
  all.markers = FindAllMarkers(obj.srt, logfc.threshold = log2(logfc_threshold), 
                               only.pos = TRUE, 
                               test.use = test_method, min.pct = min_percent)
  return(all.markers)
}
# apply:
all.markers= find_and_save_markers(obj.srt= obj.srt, cluster_id = 'RNA_snn_res.0.2', 
                                   logfc_threshold = 1.2, test_method = 'wilcox', min_percent = 0.25)
output_file = '../data/mks/re_analysis/MN3.MN4.res.0.2.all_markers.23.07.27.csv'
write.csv(all.markers, file = output_file, row.names = T)


## generate top 50 genes table 
df = all.markers %>% dplyr::select(cluster) %>% table() %>% data.frame()
df$genes=''
for(i in 1:nrow(df)){
  cl=df[i,]$cluster 
  df[i,]$genes = all.markers %>% filter(cluster== cl) %>% arrange(desc(avg_log2FC)) %>%
    head(50) %>% dplyr::select(gene) %>% pull() %>% paste0(collapse = ',')
}
df %>% write.csv(file = '../data/mks/re_analysis/MN3.MN4.res.0.2.all_markers.top50.23.07.27.csv')


## Visualize data

# annotation was done manually
DimPlot(obj.srt, group.by = 'RNA_snn_res.0.2', label = T, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'RNA_snn_res.0.2', label = F, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'orig.ident', label = F, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 

obj.srt@meta.data$cell_type = factor(obj.srt@meta.data$RNA_snn_res.0.2)
levels(obj.srt@meta.data$cell_type) = c("NK", "HUVEC_1", "fibroblast1", "Monocyte", "fibroblast2", 
                                        "fibroblast3", "CORL47_1", "CORL47_2", "HUVEC_2", "UKN_1", "UKN_2")
DimPlot(obj.srt, group.by = 'cell_type', label = T, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'cell_type', label = F, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
obj.srt@meta.data$main_cell_type = factor(obj.srt@meta.data$RNA_snn_res.0.2)
levels(obj.srt@meta.data$main_cell_type) = c("NK", "HUVEC", "fibroblast", "Monocyte", "fibroblast", 
                                        "fibroblast", "CORL47", "CORL47", "HUVEC", "UKN", "UKN")
DimPlot(obj.srt, group.by = 'main_cell_type', label = T, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 
DimPlot(obj.srt, group.by = 'main_cell_type', label = F, label.box = T, label.size = 3) + theme(plot.title = element_blank()) 

# cell type marker genes (manually chosen)
gene_list <- c("GNLY", "GZMB", "XCL2", "NKG7", "PECAM1", "CDH5", "ANGPT2", "LYVE1", "CCL14", "CD93", "MCAM", "LAMA4", "CLEC14A", "COL3A1", "DCN", "COL1A2", "COL1A1", "CD68", "LYZ", "CCL2", "CCL3", "APOBEC3A", "EPCAM", "ASCL1", "BIRC5")
# VlnPlot(obj.srt, features = gene_list, group.by = 'RNA_snn_res.0.2',
#         stack = T, flip = F) + NoLegend()

# Reverse the order of 'ident' levels
obj.srt@meta.data$rev.cluster <- factor(obj.srt@meta.data$RNA_snn_res.0.2, 
                                        levels = rev(levels(obj.srt@meta.data$RNA_snn_res.0.2)))

# Create the violin plot with the reversed order of identities
VlnPlot(obj.srt, features = gene_list, group.by = "rev.cluster",
        stack = T, flip = F) + NoLegend()

# heatmap of top 5 cluster defining genes 
gs=all.markers %>% group_by(cluster) %>% 
  top_n(5, avg_log2FC) %>% select(gene) %>% pull()
DoHeatmap(obj.srt, features = gs, group.by = 'RNA_snn_res.0.2')

## distribution of cells 
obj.srt@meta.data[1:3,]
colors <- viridis::viridis_pal(option = "D")(11)
obj.srt@meta.data %>% ggplot(aes(orig.ident, fill=rev.cluster)) + 
  geom_bar(position = 'fill') +coord_flip() +theme_bw() +scale_fill_manual(values = colors)
colors <- viridis::viridis_pal(option = "D")(6)
obj.srt@meta.data %>% ggplot(aes(orig.ident, fill=main_cell_type)) + 
  geom_bar(position = 'fill') +coord_flip() +theme_bw() +scale_fill_manual(values = colors)
ggplot(obj.srt@meta.data, aes(orig.ident, fill = factor(main_cell_type, levels = rev(levels(main_cell_type))))) + 
  geom_bar(position = 'fill') +coord_flip() +theme_bw() +scale_fill_manual(values = colors) +
  theme(legend.title =element_blank())


## DEG gene analysis
mks.function= function(obj.srt, celltype='NK'){
  obj.srt@meta.data =obj.srt@meta.data %>% 
    mutate(compare=ifelse(main_cell_type == celltype & orig.ident == 'MN3','g1',
                          ifelse(main_cell_type == celltype & orig.ident == 'MN4','g2','other')))
  # add compare to the ident
  Idents(obj.srt) = 'compare'
  logfc=log2(1)
  mks =FindMarkers(obj.srt, ident.1 = 'g2', ident.2 = 'g1', 
                   logfc.threshold = logfc)
  pval=0.05
  fc=1.2
  mks = mks %>% mutate(DE=ifelse(avg_log2FC >= log2(fc) & p_val_adj < pval, 'UP',
                                 ifelse(avg_log2FC <= -log2(fc) & p_val_adj < pval, 'DN','no_sig')))
  mks$DE = factor(mks$DE, levels = c('UP','DN','no_sig'))
  mks$gene = rownames(mks)
  mks =mks %>% mutate(labels= ifelse(DE == 'UP', gene, ifelse(DE=='DN',gene,'other')))
  mks =mks %>% arrange(desc(avg_log2FC))
  return(mks)
}

mks.nk=mks.function(obj.srt = obj.srt, celltype = 'NK')
mks.huvec=mks.function(obj.srt = obj.srt, celltype = 'HUVEC')
mks.fibro=mks.function(obj.srt = obj.srt, celltype = 'fibroblast')
mks.mono=mks.function(obj.srt = obj.srt, celltype = 'Monocyte')
mks.corl=mks.function(obj.srt = obj.srt, celltype = 'CORL47')

# save DEG files 
l <- list(NK = mks.nk, 
          HUVEC = mks.huvec, 
          Fibro = mks.fibro, 
          Mono = mks.mono, 
          CORL = mks.corl)

lapply(names(l), function(name) {
  write.csv(l[[name]], paste0(name, '.mn3.4.deg.csv'), row.names = T)
})

# plot volcano
p=mks.huvec %>% 
  ggplot(aes(avg_log2FC, -log10(p_val_adj), color=DE)) + 
  geom_point(size=1, alpha=0.5) + 
  scale_color_manual(values = c('red','blue','grey')) +
  theme_classic() +
  geom_vline(xintercept = c(-log2(fc),log2(fc)), color='grey') +
  geom_hline(yintercept = -log10(0.05),color='grey') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggtitle(paste0(celltype,' population')) +
  ggrepel::geom_text_repel(aes(label=labels), size = 2.5, 
                           min.segment.length = Inf, show.legend = F) +
  ggeasy::easy_center_title() ## to center title
p



## plot label looks somewhat weird for some labels
plot_volcano=function(mks=mks, celltype='NK'){
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
    ggrepel::geom_text_repel(aes(label=labels), size = 2.5,show.legend = F) +
    ggeasy::easy_center_title() ## to center title
  print(p)
}
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

plot_volcano(mks = mks.nk,celltype = 'NK')
plot_volcano(mks = mks.huvec,celltype = 'HUVEC')
plot_volcano(mks = mks.fibro,celltype = 'fibroblast')
plot_volcano(mks = mks.mono,celltype = 'Monocyte')
plot_volcano(mks = mks.corl,celltype = 'CORL47')

## other version
mks= read.csv('../data/mks/re_analysis/HUVECmn3.4.deg.csv', row.names = 1)
plot_volcano(mks = mks,celltype = 'HUVEC')

## GSEA HALLMARK

library(clusterProfiler)

perform_GSEA <- function(res, ref, title, pvalueCutoff = 1) {
  ranking <- function(res) {
    df <- res$avg_log2FC
    names(df) <- res$gene
    df <- sort(df, decreasing = TRUE)
    return(df)
  }
  
  ranked.res <- ranking(res)
  
  x <- clusterProfiler::GSEA(geneList = ranked.res,
                             TERM2GENE = ref,
                             pvalueCutoff = pvalueCutoff,
                             pAdjustMethod = "BH",
                             verbose = TRUE,
                             seed = TRUE)
  
  result <- x@result %>% arrange(desc(NES))
  result <- result[, c('NES', 'pvalue', 'p.adjust', 'core_enrichment', 'ID')]
  result$celltype = title
  return(result)
}

# Replace 'res', 'hallmark', and 'title' with your actual data
# perform_GSEA(res, hallmark, title)

title='MONOCYTE'
gsea.mono = perform_GSEA(res = mks.mono, ref = hallmark, title = title)
write.csv(gsea.mono, file = paste0(title, ' ', 'GSEA.HALLMARK.csv'))
gsea.mono %>% nrow()

title='HUVEC'
gsea.huvec = perform_GSEA(res = mks.huvec, ref = hallmark, title = title)
write.csv(gsea.huvec, file = paste0(title, ' ', 'GSEA.HALLMARK.csv'))
gsea.huvec %>% nrow()

title='Fibroblast'
gsea.fibro = perform_GSEA(res = mks.fibro, ref = hallmark, title = title)
write.csv(gsea.fibro, file = paste0(title, ' ', 'GSEA.HALLMARK.csv'))
gsea.fibro %>% nrow()

title='NK'
gsea.nk = perform_GSEA(res = mks.nk, ref = hallmark, title = title)
write.csv(gsea.nk, file = paste0(title, ' ', 'GSEA.HALLMARK.csv'))
gsea.nk %>% nrow()

title='CORL47'
gsea.corl = perform_GSEA(res = mks.corl, ref = hallmark, title = title)
write.csv(gsea.corl, file = paste0(title, ' ', 'GSEA.HALLMARK.csv'))
gsea.corl %>% nrow()


gsea.corl[1:3,]

## find the shared GSEA IDs
gsea.list = list(NK=gsea.nk$ID,
                 MONO=gsea.mono$ID,
                 HUVEC=gsea.huvec$ID,
                 CORL=gsea.corl$ID,
                 FIBRO=gsea.fibro$ID)
## generate venndiagram data
## generate venn diagram 
ven_gsea <- VennDetail::venndetail(gsea.list)
ven_gsea@result[1:3,]
IDs=ven_gsea@result[ven_gsea@result$Subset =='Shared',]$Detail

gsea.all=rbind(gsea.corl[IDs,c('ID','celltype','NES','pvalue','p.adjust')],
               gsea.huvec[IDs,c('ID','celltype','NES','pvalue','p.adjust')],
               gsea.mono[IDs,c('ID','celltype','NES','pvalue','p.adjust')],
               gsea.nk[IDs,c('ID','celltype','NES','pvalue','p.adjust')],
               gsea.fibro[IDs,c('ID','celltype','NES','pvalue','p.adjust')])
rownames(gsea.all) = paste0(gsea.all$ID,'_',gsea.all$celltype)
gsea.all %>% write.csv('GESA.HALLMARK.all.csv')

# read csv 
gsea.all = read.csv('../data/mks/re_analysis/GESA.HALLMARK.all.csv', row.names = 1)

## GSEA NES plot
ggplot(gsea.huvec, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill=p.adjust)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0(title)) + 
  theme_classic() +
  scale_fill_gradient(low = 'red', high = '#E5E7E9') +
  theme(axis.text.x= element_text(size=5, face = 'bold'),
        axis.text.y= element_text(size=6, face = 'bold'), 
        axis.title =element_text(size=10)) +
  ggeasy::easy_center_title() ## to center title

nk.gsea=read.csv('../data/mks/re_analysis/NK GSEA.HALLMARK.csv')
huvec.gsea=read.csv('../data/mks/re_analysis/HUVEC GSEA.HALLMARK.csv')
fibro.gsea=read.csv('../data/mks/re_analysis/Fibroblast GSEA.HALLMARK.csv')
mono.gsea=read.csv('../data/mks/re_analysis/MONOCYTE GSEA.HALLMARK.csv')
corl.gsea=read.csv('../data/mks/re_analysis/CORL47 GSEA.HALLMARK.csv')

## GSEA plot by function
plot_gsea <- function(data, title) {
  p <- ggplot(data, aes(reorder(ID, NES), NES)) +
    geom_col(aes(fill=p.adjust)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(title)) + 
    theme_classic() +
    scale_fill_gradient(low = 'red', high = '#E5E7E9') +
    theme(axis.text.x= element_text(size=5, face = 'bold'),
          axis.text.y= element_text(size=8, face = 'bold'), 
          axis.title =element_text(size=10))
  
  print(p)
}
plot_gsea(nk.gsea, ' HALLMARK pathways')
plot_gsea(huvec.gsea, ' HALLMARK pathways')
plot_gsea(fibro.gsea, ' HALLMARK pathways')
plot_gsea(mono.gsea, ' HALLMARK pathways')
plot_gsea(corl.gsea, ' HALLMARK pathways')

## GSEA plot version 2
plot_gsea <- function(data, title) {
  p <- ggplot(data, aes(reorder(ID, NES), NES)) +
    geom_col(aes(fill = ifelse(NES > 0, 'positively enriched', 'negatively enriched'))) +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste0(title),
         fill = "Enrichment") +  # Add legend title here
    theme_classic() +
    scale_fill_manual(values = c('blue', 'red')) +
    theme(axis.text.x = element_text(size = 5, face = 'bold'),
          axis.text.y = element_text(size = 8, face = 'bold'), 
          axis.title = element_text(size = 10))
  
  print(p)
}
plot_gsea(data, ' HALLMARK pathways')
p1=plot_gsea(nk.gsea, ' NK DEG MN4 vs MN3')
p2=plot_gsea(huvec.gsea, ' HUVEC DEG MN4 vs MN3')
p3=plot_gsea(fibro.gsea, ' FIBROBLAST DEG MN4 vs MN3')
p4=plot_gsea(mono.gsea, ' MONOCYTE DEG MN4 vs MN3')
p5=plot_gsea(corl.gsea, ' CORL47 DEG MN4 vs MN3')

# Assuming you have already created p1, p2, p3, p4, and p5 plots

# List of plots
plot_list <- list(p1, p2, p3, p4, p5)

# Output directory
output_dir <- '../data/mks/re_analysis/'

# Iterate over the list of plots and save each plot as a PNG file
for (i in seq_along(plot_list)) {
  filename <- paste0(output_dir, 'GSEAplot_', i, '.png')
  png(filename, width = 10, height = 8, units = 'in', res = 300)
  print(plot_list[[i]])
  dev.off()
}


## sankey plot 
library(ggsankey)
library(ggplot2)
library(dplyr)

data =obj.srt@meta.data[, c('orig.ident','RNA_snn_res.0.2','cell_type','main_cell_type')]

# Function to create a Sankey plot
data =data[,c(1,3,4)]
create_sankey_plot <- function(data, column_names, title, show_labels = FALSE){
  df <- data %>%
    make_long(!!column_names[1], !!column_names[2], !!column_names[3])
  
  dagg <- df %>%
    group_by(node) %>%
    tally()
  
  df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
  
  pl <- ggplot(df2, aes(x = x,
                        next_x = next_x,
                        node = node,
                        next_node = next_node,
                        fill = factor(node),
                        label = paste0(node," n=", n))
  )
  
  pl <- pl + geom_sankey(flow.alpha = 0.5, color = "gray40", show.legend = show_labels)
  if (show_labels) {
    pl <- pl + geom_sankey_label(size = 3, color = "white", fill = "gray40", hjust = 1)  # Change hjust value to 1 (right-aligned)
  }
  
  pl <- pl + theme_bw()
  pl <- pl + theme(legend.position = "none")
  pl <- pl + theme(axis.title = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   panel.grid = element_blank())
  pl <- pl + scale_fill_viridis_d(option = "inferno")
  pl <- pl + labs(title = "Sankey diagram")
  pl <- pl + labs(subtitle = "cell type distribution by MN3 MN4 identity")
  
  pl <- pl + labs(fill = 'Nodes')
  
  return(pl)
}

# Assuming obj.srt@meta.data is your data and you want to use specific columns for the Sankey plot
column_names <- c('orig.ident', 'RNA_snn_res.0.2', 'cell_type', 'main_cell_type')
column_names <- c('orig.ident', 'main_cell_type','cell_type')

# Create the Sankey plot using the specified columns
sankey_plot <- create_sankey_plot(data = obj.srt@meta.data, column_names, title = "My Sankey Plot", show_labels = TRUE)

# Display the plot
print(sankey_plot)

## GSEA dotplot
gsea.all = read.csv('../data/mks/re_analysis/GESA.HALLMARK.all.csv', row.names = 1)

gsea.all[1:3,]
custom_palette <- colorRampPalette(c("#2874A6", "white", "#FC4119"))(1000)

gsea.all %>%
  ggplot(aes(x = ID, y = celltype)) +
  geom_point(aes(size = NES, color = -(p.adjust))) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_gradientn(colors = custom_palette) 

# corrected 23.08.02

# descending order of NES
df =gsea.all %>% select(celltype,ID, NES) %>% tidyr::spread(ID, NES)
colsm = colSums(df[,2:ncol(df)])
colsm <- colsm[order(-colsm)]
names(colsm)

# Get the order of the names in colsm
order_x <- names(colsm)

# Convert ID column to a factor with the desired order
gsea.all$ID <- factor(gsea.all$ID, levels = order_x)

# Plot the ggplot with the reordered x-axis
gsea.all %>%
  ggplot(aes(x = ID, y = celltype)) +
  geom_point(aes(size = -(p.adjust), color = NES)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_gradientn(colors = custom_palette) +
  labs(color = "NES", size = "FDR")

#################################################
## generate the combined DGE file
library(dplyr)

# Define the cell types in the desired order
cell_types <- c('CORL47', 'Fibroblast', 'HUVEC', 'Monocyte', 'NK')

# Get the list of files matching the pattern 'deg.csv'
fs <- list.files(path = '.', pattern = 'deg.csv', full.names = TRUE)

# Function to read and preprocess each CSV file
read_and_preprocess_csv <- function(file_path, cell_type) {
  df <- read.csv(file_path)  # Read the CSV file
  df$cell_type <- cell_type  # Add the 'cell_type' column and fill with the specified value
  return(df)
}

# Read and preprocess all CSV files, and combine them using rbind
combined_data <- do.call(rbind, mapply(read_and_preprocess_csv, fs, cell_types, SIMPLIFY = FALSE))
combined_data = combined_data[,-1]
combined_data$cell_type %>% table()
combined_data %>% write.csv('DEG.combined.csv')

#######################################
## VlnPlot 
# Function to create Violin Plots by compare levels
plot_vln_by_compare <- function(main_cell_type, seurat_obj, gene_vector) {
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    mutate(compare = ifelse(main_cell_type == main_cell_type & orig.ident == 'MN3', paste0(main_cell_type, ' MN3'),
                            ifelse(main_cell_type == main_cell_type & orig.ident == 'MN4', paste0(main_cell_type, ' MN4'), 'other')))
  
  compare_levels <- c(paste0(main_cell_type, ' MN3'), paste0(main_cell_type, ' MN4'))
  
  Idents(seurat_obj) <- 'compare'
  VlnPlot(seurat_obj, features = gene_vector, idents = compare_levels, stack = TRUE, flip = TRUE) + NoLegend()
}

# Define the main_cell_type and gene_vector for each cell type
cell_types <- c('CORL47', 'HUVEC', 'fibroblast', 'NK', 'Monocyte')

# Upregulated genes in MN4
gene_vectors_up <- list(
  c("ISG15", "MT2A", "STAT1", "IRF1", "B2M", "GBP1", "PARP14", "IFI6", "TUBA1A", "MLLT11", "IFIT3", "IRF9", "DLX6-AS1", "TCF4", "IFI27"),
  c("CXCL10", "CXCL11", "CXCL9", "IFI27", "ISG15", "GBP1", "IFIT3", "IFI6", "MX1", "IFI44L", "IFIT2", "BST2", "RNF213", "IFIT1", "WARS"),
  c("ISG15", "IFI27", "CXCL10", "IFIT3", "IFI6", "GBP1", "WARS", "IFIT2", "BST2", "MX1", "IL32", "PMAIP1", "LAP3", "CXCL11", "IDO1"),
  c("ISG15", "IFIT2", "IFIT3", "IFI6", "MX1", "IFIT1", "ISG20", "RSAD2", "HERC5", "OASL", "IFI44L", "XAF1", "STAT1", "PMAIP1", "LY6E"),
  c("ISG15", "IFIT2", "CCL8", "RSAD2", "APOBEC3A", "IFIT3", "IFITM3", "IFI27", "TNFSF10", "IFI6", "CXCL10", "LY6E", "MX1", "ISG20", "IFIT1"))

# Downregulated genes in MN4
gene_vectors_down <- list(
  c("MYDGF", "PPIB", "SERP1", "SSR3", "SEC11C", "AL627171.2", "SURF2", "EDF1", "PRDX4", "OSTC", "SSR2", "SIVA1", "SEC61B", "TCEA1", "NOL7"),
  c("LYVE1", "IGFBP4", "EGR1", "CD34", "CLDN5", "SPRY1", "GCHFR", "CCL14", "SESN3", "FABP5", "PCDH17", "RAMP2", "AHR", "TM4SF18", "MGST2"),
  c("MGP", "CLU", "COL3A1", "DCN", "AREG", "LUM", "ELN", "A2M", "GPX3", "CEBPD", "FGF7", "ADH1B", "FOS", "C11orf96", "APOE"),
  c("HIST1H1D", "DUSP2", "KLRB1", "NOP53", "EEF1B2", "ZFP36L2", "EEF1G", "LTB", "MATK", "TPT1", "PNRC1", "EEF1D", "EIF4B", "IKZF1", "RACK1"),
  c("CYP27A1", "G0S2", "S100A4", "MRC1", "VCAN", "FBP1", "NME2", "APRT", "LYZ", "NACA", "SNHG29", "LTA4H", "EEF1B2", "CD14", "IL1B")
)

# Combine both lists into one
# gene_vectors <- c(gene_vectors_up, gene_vectors_down)
# names(gene_vectors) <- cell_types

# Create and arrange the plots using cowplot
# up
# down
plots <- lapply(1:length(cell_types), function(i) {
  main_cell_type <- cell_types[i]
  gene_vector <- gene_vectors_up[[i]]
  plot_vln_by_compare(main_cell_type, obj.srt, gene_vector)
})
plots <- lapply(1:length(cell_types), function(i) {
  main_cell_type <- cell_types[i]
  gene_vector <- gene_vectors[[i]]
  plot_vln_by_compare(main_cell_type, obj.srt, gene_vector)
})

# combine plots
cowplot::plot_grid(plotlist = plots, nrow = 1)
