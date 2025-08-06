#### load packages ---- 
library(Seurat)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(GSVA)
library(tidyverse)
library(msigdbr)
library(org.Mm.eg.db)


####  Load data----

obj <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/taste.integrated.rds")
meta <- read.table("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/01.bigclass/bigclass_annotation.txt")
seurat.obj <- AddMetaData(object = obj, metadata = meta)


Idents(seurat.obj) <- "timepoint"


# 计算各组别间细胞的表达量均值
expr <- AverageExpression(seurat.obj, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)

Mus_GO <- msigdbr(species = "Mus musculus", 
                  category = "C5",
                  subcategory = "GO:BP") %>%  
  dplyr::select(gs_name,gene_symbol)
Mus_GO_Set = Mus_GO %>% split(x = .$gene_symbol, f = .$gs_name)

gsva.go <- gsva(expr, 
                gset.idx.list = Mus_GO_Set, 
                kcdf="Gaussian",
                method = "gsva",
                parallel.sz=1)

#保存GSVA结果
gsva.go2 <- as.data.frame(gsva.go)
gsva.go2 <- rownames_to_column(gsva.go2, var = "way") 
gsva.go2$way <- gsva.go2$way |>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

new_order <- c("way","P0","P14","M2","M12","M24")
gsva.go2 <- gsva.go2[, new_order]
write.csv(gsva.go2, "All_gsva_GO.csv", quote = F, row.names = F)
