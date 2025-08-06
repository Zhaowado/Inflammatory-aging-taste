library(Seurat)
library(dplyr)

parser = argparse::ArgumentParser(description = "geneset score")

parser$add_argument('-f', dest = 'file', help = 'input gene set file')
parser$add_argument('-r', dest = 'rds', help = 'input rds')
parser$add_argument('-o', dest = 'output', help = 'output rds')

opts = parser$parse_args()

obj <- readRDS(opts$rds)
DefaultAssay(obj) <- "RNA"

gene_set <- read.table(opts$file, sep = "\t",quote = F)
genes <- gene_set$V1

############
present_genes <- unique(genes[genes %in% rownames(obj)])
if (length(present_genes) < nrow(gene_set)) {   
  missing_genes <- setdiff(genes, rownames(obj))  
  warning("The following features are not present in the object: ", paste(missing_genes, collapse = ", "))  
}  
###########
obj <- AddModuleScore(object = obj, features = list(present_genes), ctrl = 100, name = "Score")
scoreColname <- "Score1"

saveRDS(obj,opts$output)
