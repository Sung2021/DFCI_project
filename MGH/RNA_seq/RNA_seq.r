# Load required packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggrepel)
library(DEGreport)
library(purrr)

# Function to run DESeq2 and get results for all clusters comparing group A vs. group B
get_dds_resultsAvsB <- function(x, A, B){
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  counts <- pb[[clusters[x]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

  #all(rownames(cluster_metadata) == colnames(cluster_counts))

  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)

  # Transform counts for data visualization
  rld <- rlog(dds, blind=TRUE)

  # Plot PCA
  DESeq2::plotPCA(rld, intgroup = "group_id")
  ggsave(paste0("results/", clusters[x], "_specific_PCAplot.png"))

  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)

  # Plot heatmap
  png(paste0("results/", clusters[x], "_specific_heatmap.png"))
  pheatmap(rld_cor, annotation = cluster_metadata[, c("group_id"), drop = F])
  dev.off()

  # Run DESeq2 differential expression analysis
  dds <- DESeq(dds)

  # Plot dispersion estimates
  png(paste0("results/", clusters[x], "_dispersion_plot.png"))
  plotDispEsts(dds)
  dev.off()

  # Output results of Wald test for contrast for A vs B
  contrast <- c("group_id", levels(cluster_metadata$group_id)[A], levels(cluster_metadata$group_id)[B])

  # resultsNames(dds)
  res <- results(dds, 
                 contrast = contrast,
                 alpha = 0.05)

  res <- lfcShrink(dds, 
                   contrast =  contrast,
                   res = res)
  # Set thresholds
  padj_cutoff <- 0.05

  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()

  write.csv(res_tbl,
            paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)

  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)

  write.csv(sig_res,
            paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)

  ## ggplot of top genes
  normalized_counts <- counts(dds, 
                              normalized = TRUE)

  ## Order results by padj values
  top20_sig_genes <- sig_res %>%
    dplyr::arrange(padj) %>%
    dplyr::pull(gene) %>%
    head(n = 20)


  top20_sig_norm <- data.frame(normalized_counts) %>%
    rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% top20_sig_genes)

  gathered_top20_sig <- top20_sig_norm %>%
    gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")

  gathered_top20_sig <- inner_join(ei[, c("sample_id", "group_id" )], gathered_top20_sig, by = c("sample_id" = "samplename"))

  ## plot using ggplot2
  ggplot(gathered_top20_sig) +
    geom_point(aes(x = gene, 
                   y = normalized_counts, 
                   color = group_id), 
               position = position_jitter(w = 0.1, h = 0)) +
    scale_y_log10() +
    xlab("Genes") +
    ylab("log10 Normalized Counts") +
    ggtitle("Top 20 Significant DE Genes") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("DESeq2/pairwise/", clusters[x], "_", levels(cluster_metadata$group_id)[A], "_vs_", levels(cluster_metadata$group_id)[B], "_top20_DE_genes.png"))
}

# Function to run likelihood ratio test (LRT) on all clusters
get_dds_LRTresults <- function(x){
  cluster_metadata <- metadata[which(metadata$cluster_id == clusters[x]), ]
  rownames(cluster_metadata) <- cluster_metadata$sample_id
  counts <- pb[[clusters[x]]]
  cluster_counts <- data.frame(counts[, which(colnames(counts) %in% rownames(cluster_metadata))])

  #all(rownames(cluster_metadata) == colnames(cluster_counts))        

  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)

  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)

  # Extract results
  res_LRT <- results(dds_lrt)

  # Create a tibble for LRT results
  res_LRT_tb <- res_LRT %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>% 
    as_tibble()

  # Save all results
  write.csv(res_LRT_tb,
            paste0("DESeq2/lrt/", clusters[x], "_LRT_all_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)

  # Subset to return genes with padj < 0.05
  sigLRT_genes <- res_LRT_tb %>% 
    filter(padj < 0.05)

  # Save sig results
  write.csv(sigLRT_genes,
            paste0("DESeq2/lrt/", clusters[x], "_LRT_sig_genes.csv"),
            quote = FALSE, 
            row.names = FALSE)

  # Transform counts for data visualization
  rld <- rlog(dds_lrt, blind=TRUE)

  # Extract the rlog matrix from the object and compute pairwise correlation

 values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)


  # Obtain rlog values for those significant genes
  cluster_rlog <- rld_mat[sigLRT_genes$gene, ]

  cluster_meta_sig <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]

  # # Remove samples without replicates
  # cluster_rlog <- cluster_rlog[, -1]
  # cluster_metadata <- cluster_metadata[which(rownames(cluster_metadata) %in% colnames(cluster_rlog)), ]


  # Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
  cluster_groups <- degPatterns(cluster_rlog, metadata = cluster_meta_sig, time = "group_id", col=NULL)
  ggsave(paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.png"))

  # Let's see what is stored in the `df` component
  write.csv(cluster_groups$df,
            paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.csv"),
            quote = FALSE, 
            row.names = FALSE)

  saveRDS(cluster_groups, paste0("DESeq2/lrt/", clusters[x], "_LRT_DEgene_groups.rds"))
  save(dds_lrt, cluster_groups, res_LRT, sigLRT_genes, file = paste0("DESeq2/lrt/", clusters[x], "_all_LRTresults.Rdata"))
}

# Assuming 'metadata', 'clusters', 'pb', and 'gg_df' are the relevant data and variables.
# Run the script on all clusters comparing stim condition relative to control condition
dir.create("DESeq2/pairwise")
map(1:length(clusters), get_dds_resultsAvsB, A = 2, B = 1)

# Run LRT on all clusters
dir.create("DESeq2/lrt")
map(1:length(clusters), get_dds_LRTresults)

# Rest of the code for the heatmap, volcano plot, and other visualizations.
# ... (The previous code for the heatmap, volcano plot, and other visualization parts. Use the provided code.)
