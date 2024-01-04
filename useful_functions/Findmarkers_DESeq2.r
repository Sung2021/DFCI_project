obj.lymphoid@meta.data[1:3,]
obj.lymphoid@meta.data$response_treatment = paste0(obj.lymphoid@meta.data$response,"_",obj.lymphoid@meta.data$sample)
obj.lymphoid@meta.data$response_treatment = factor(obj.lymphoid@meta.data$response_treatment)
levels(obj.lymphoid@meta.data$response_treatment)

obj.lymphoid@meta.data %>% dplyr::select(response_treatment) %>% table() %>%  View()

obj.lymphoid@meta.data %>% dplyr::select(response_treatment,lymphoid_celltype ) %>% table() %>%  View()


# Effector T  

update_gene_Status <- function(obj, g) {
  exp <- obj@assays$RNA@data[g, ]
  
  # Assuming that 'CD8A' is already a column in 'obj@meta.data'
  obj@meta.data[, paste0("exp_",g)] <- exp

    obj@meta.data <- obj@meta.data %>%
    mutate(!!paste0(g, "_pos") := ifelse(.data[[paste0("exp_",g)]] > 0, paste0(g, "+"), paste0(g, "-")))
  
  return(obj)
}

obj.lymphoid <- update_gene_Status(obj.lymphoid, 'CD8A')
obj.lymphoid <- update_gene_Status(obj.lymphoid, 'CD8B')


obj.tmp = subset(obj.lymphoid, response_treatment %in% c("NR_CTL","R_CTL") & 
                   lymphoid_celltype == "CD8_effector")
obj.tmp2 = subset(obj.tmp, CD8A_pos == "CD8A+" | CD8B_pos == "CD8B+" )

# To run FindMarkers with DESeq2, the next step is required. 
# Or subset CD8A+ | CD8B+, there will be no 0 in "counts", so no need to add 1.
# obj.tmp[["RNA"]]@counts<-as.matrix(obj.tmp[["RNA"]]@counts)+1

obj.tmp2 <- NormalizeData(obj.tmp2)  
obj.tmp2 <- ScaleData(object = obj.tmp2, features = VariableFeatures(object = obj.tmp2))

obj.tmp2@meta.data =obj.tmp2@meta.data %>% 
  mutate(compare=ifelse(response_treatment == "R_CTL", "R_CTL", "NR_CTL"))
Idents(obj.tmp2) = 'compare'
logfc=log2(1)
g1="NR_CTL"
g2="R_CTL"
mks =FindMarkers(obj.tmp2, ident.1 = g2, ident.2 = g1, 
                 logfc.threshold = logfc, test.use = "DESeq2")
# mks =FindMarkers(obj.tmp, ident.1 = g2, ident.2 = g1, 
#                  logfc.threshold = logfc, test.use = "DESeq2", slot = "counts")
pval=0.05
fc=1.2
mks = mks %>% mutate(DE=ifelse(avg_log2FC >= log2(fc) & p_val_adj < pval, 'UP',
                               ifelse(avg_log2FC <= -log2(fc) & p_val_adj < pval, 'DN','no_sig')))
mks$DE = factor(mks$DE, levels = c('UP','DN','no_sig'))
mks$gene = rownames(mks)
mks =mks %>% mutate(labels= ifelse(DE == 'UP', gene, ifelse(DE=='DN',gene,'other')))
mks =mks %>% arrange(desc(avg_log2FC))


mks %>% 
  ggplot(aes(avg_log2FC, -log10(p_val_adj), color=DE)) + 
  geom_point(size=1, alpha=0.5) + 
  scale_color_manual(values = c('red','blue','grey')) +
  theme_classic() +
  geom_vline(xintercept = c(-log2(fc),log2(fc)), color='grey') +
  geom_hline(yintercept = -log10(0.05),color='grey') +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  ggtitle("Resist/Naive") +
  ggeasy::easy_center_title() ## to center title
