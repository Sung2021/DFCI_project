library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

dir = "~/Desktop/DF/ref/"

gtf_path <- paste0(dir,"hg38_gtf_ucsc/hg38.refGene.gtf")
gtf <- rtracklayer::import(gtf_path)

txdb <- makeTxDbFromGFF(gtf_path, format="gtf")

# 유전자별로 시작점과 끝점 가져오기
genes <- genes(txdb)

# 유전자 ID와 해당 유전자의 길이를 포함하는 데이터 프레임 생성
gene_lengths <- with(as.data.frame(genes), {
  gene_id = mcols(genes)$gene_id
  start = start
  end = end
  width = end - start + 1
  data.frame(gene_id, start, end, width)
})

head(gene_lengths)

gene_lengths %>% dim()

# Check 
gs = c("TACSTD2","EGFR") 
gene_lengths[gene_lengths$gene_id %in% gs,]

# Save output 
gene_lengths %>% write.csv(paste0(dir,"hg38_refGene_geneLength.csv"))
