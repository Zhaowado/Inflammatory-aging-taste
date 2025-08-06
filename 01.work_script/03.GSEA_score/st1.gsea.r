## Get the parameters
parser = argparse::ArgumentParser(description = "rscript for GSVA")

parser$add_argument('-i', dest = 'input', help = 'input DEG files')
parser$add_argument('-c', dest = 'celltype', help = 'select celltype for scoring')

opts = parser$parse_args()

## lib
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

markers <- read.csv(opts$input)
tastes <- markers[markers$celltype %in% opts$celltype & markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25,]
genes <- tastes$gene
gene=bitr(genes,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Mm.eg.db")
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
tastes$SYMBOL <- tastes$gene
tastes <- tastes %>%
  inner_join(gene,by="SYMBOL")

tastes <- tastes %>%
  arrange(desc(avg_log2FC))

geneList = tastes$avg_log2FC
names(geneList) <- tastes$ENTREZID

gse.KEGG <- gseKEGG(geneList,
                    organism = "mmu", # 人 hsa
                    pvalueCutoff = 1,
                    pAdjustMethod = "BH")
write.csv(gse.KEGG,"M24_others_gsea_kegg.csv")
