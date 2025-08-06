library(dplyr)
library(Seurat)
library(patchwork) 
library(monocle)
library(ggplot2)

obj <- readRDS("../../taste.final.filter.rds")
taste <- subset(obj,figname == "MTC")

DefaultAssay(taste) <- "SCT"
data <- as.matrix(taste@assays$RNA@data)
pd <- new('AnnotatedDataFrame', data = taste@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
my_cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(my_cds),
                                    num_cells_expressed >= 10))
diff <- differentialGeneTest(my_cds[expressed_genes,],fullModelFormulaStr="~timepoint",cores=1)

deg <- subset(diff, qval < 0.05) 
deg <- deg[order(deg$qval,decreasing=F),]
write.table(deg,"degs.xls",sep = "\t",quote = F)

ordergene <- rownames(deg) 
my_cds <- setOrderingFilter(my_cds, ordergene)  

my_cds <- reduceDimension(my_cds, max_components = 2,
                       method = 'DDRTree')
my_cds <- orderCells(my_cds) 
my_cds <- orderCells(my_cds)
  GM_state <- function(my_cds){
    if (length(unique(pData(my_cds)$State)) > 1){
      T0_counts <- table(pData(my_cds)$State, pData(my_cds)$timepoint)[,"P0"]
      return(as.numeric(names(T0_counts)[which
                                         (T0_counts == max(T0_counts))]))
    } else {
      return (1)
    }
  }
my_cds <- orderCells(my_cds,root_state = GM_state((my_cds)))

saveRDS(my_cds,"order.timepoint.cds")
