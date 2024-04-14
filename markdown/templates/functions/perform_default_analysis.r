## perform default analysis
perform_default_analysis <- function(obj.srt, n_features = 3000, n_pcs = 30, 
                                     dims_for_neighbors = 1:20, 
                                     resolutions = c(0.1, 0.2, 0.4, 0.8,1), 
                                     umap_dims = 1:20) {
  # Step 1: Find variable features
  obj.srt <- FindVariableFeatures(obj.srt, selection.method = 'vst', nfeatures = n_features)
  
  # Step 2: Scale and normalize data
  all_genes <- rownames(obj.srt)
  obj.srt <- NormalizeData(obj.srt)
  obj.srt <- ScaleData(obj.srt, features = all_genes)
  
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
