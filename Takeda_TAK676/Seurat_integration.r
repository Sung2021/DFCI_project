obj.list <- SplitObject(obj.srt, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)

#features 
# Perform integration
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
obj.intgr <- RunPCA(obj.intgr, npcs = 30, verbose = FALSE)
obj.intgr <- RunUMAP(obj.intgr, reduction = "pca", dims = 1:30)
obj.intgr <- FindNeighbors(obj.intgr, reduction = "pca", dims = 1:30)
obj.intgr <- FindClusters(obj.intgr, resolution = c(0.1,0.2,0.5,0.8))
