# https://satijalab.org/seurat/articles/integration_introduction.html#perform-integration

obj.srt@meta.data[1:3,]
obj.srt$intgr_set= paste0(obj.srt$orig.ident,"_",obj.srt$sample)

# Integration of 6 samples + 2 conditions (12 sets)   
obj.list <- SplitObject(obj.srt, split.by = "intgr_set")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

#features 
# Perform integration (two steps)
obj.anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                      anchor.features = features)
# this command creates an 'integrated' data assay
obj.intgr <- IntegrateData(anchorset = obj.anchors)

# Perform an integrated analysis

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(obj.intgr) <- "integrated"

# Run the standard workflow for visualization and clustering
obj.intgr <- ScaleData(obj.intgr, verbose = FALSE)
obj.intgr <- RunPCA(obj.intgr, npcs = 20, verbose = FALSE)
obj.intgr <- RunUMAP(obj.intgr, reduction = "pca", dims = 1:20)
obj.intgr <- FindNeighbors(obj.intgr, reduction = "pca", dims = 1:20)
obj.intgr <- FindClusters(obj.intgr, resolution = c(0.1,0.2,0.5,0.8))

#obj.intgr %>% saveRDS(paste0(dir,"rds/sung.lymphoid_intgr.rds"))

DimPlot(obj.intgr, group.by="intgr_set")

DimPlot(obj.intgr, group.by="sample")
DimPlot(obj.intgr, group.by="orig.ident")
DimPlot(obj.intgr, group.by="integrated_snn_res.0.1", label = T)




find_and_save_markers <- function(obj.srt, cluster_id, logfc_threshold = 1.2, 
                                  test_method = 'wilcox', min_percent = 0.25) {
  Idents(obj.srt) = cluster_id
  all.markers = FindAllMarkers(obj.srt, logfc.threshold = log2(logfc_threshold), 
                               only.pos = TRUE, 
                               test.use = test_method, min.pct = min_percent)
  return(all.markers)
}


res="integrated_snn_res.0.1"
mks= find_and_save_markers(obj.srt= obj.intgr, cluster_id = res, 
                                   logfc_threshold = 1.2, test_method = 'wilcox', min_percent = 0.25)

genes = c("CD3D","CD8A","CD4",
          "CD14","APOE",
          "CD79A",
          "KRT7",
          "NKG7")
VlnPlot(obj.intgr, features = genes, group.by = res, stack = T, flip = T)

obj.srt@meta.data[1:3,] 
obj.srt@meta.data[rownames(obj.intgr@meta.data),]$intgr_cluster = obj.intgr@meta.data$integrated_snn_res.0.1

obj.srt@meta.data[rownames(obj.intgr@meta.data), "intgr_cluster"] = obj.intgr@meta.data$integrated_snn_res.0.1
DimPlot(obj.srt, group.by="intgr_cluster", label = T, split.by = "intgr_cluster", ncol = 5)


res ="intgr_cluster"
mks.obj.srt= find_and_save_markers(obj.srt= obj.srt, cluster_id = res, 
                           logfc_threshold = 1.2, test_method = 'wilcox', min_percent = 0.25)
