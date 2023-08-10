## color the population (one)
meta = obj.srt@meta.data
cgs= meta$cty %>% unique()
cgs %>% class()
cgs[2]
ps = list()
for(i in 1:length(cgs)){
  rs= meta %>% filter(cty == cgs[i]) %>% rownames()
  meta$color = 'other'
  meta[rs,]$color = as.character(cgs[i])
  
  p=obj.srt@reductions$umap@cell.embeddings %>% data.frame() %>% 
    ggplot(aes(UMAP_1, UMAP_2, color=meta$color)) + geom_point(size=0.1, alpha=0.5) +
    theme_classic() +scale_color_manual(values = c('red', 'grey')) +
    ggtitle(paste0('cell type: ', cgs[i])) +
    theme(legend.title = element_blank()) + 
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(axis.title=element_text(size=20),plot.title = element_text(size=20), legend.text = element_text(size = 15)) 
  ps[[i]] = p
}
ps[[10]]
cgs[1]
meta$color %>% table()


for(i in 1:length(cgs)){
  print(cgs[i])
}
