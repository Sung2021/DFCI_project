dir='~/Desktop/DF/DFCI_Paweletz/2023_Takeda_PD1_NK/'
# Find Markers based on cluster
# Define function
# Apply function
find_and_save_markers <- function(obj.srt, cluster_id, logfc_threshold = 1.2, 
                                  test_method = 'wilcox', min_percent = 0.25) {
  Idents(obj.srt) = cluster_id
  all.markers = FindAllMarkers(obj.srt, logfc.threshold = log2(logfc_threshold), 
                               only.pos = TRUE, 
                               test.use = test_method, min.pct = min_percent)
  return(all.markers)
}
# apply:

resolution_values <- c(0.2, 0.4, 0.8, 1)
# Loop through each resolution value
for (resolution_number in resolution_values) {
  # Construct the resolution name
  res <- paste0("RNA_snn_res.", resolution_number)
  all.markers= find_and_save_markers(obj.srt= obj.srt, cluster_id = res, 
                                   logfc_threshold = 1.2, test_method = 'wilcox', min_percent = 0.25)
  output_file = paste0(dir,"data/mks/Merge6_lymphoid_", res, ".markers.csv")
  write.csv(all.markers, file = output_file, row.names = T)
}
