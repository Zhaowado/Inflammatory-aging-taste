## Get the parameters
parser = argparse::ArgumentParser(description = "rscript for scRNAseq data integrated")

parser$add_argument('-i', dest = 'input', help = 'input directory')
parser$add_argument('-o', dest = 'out', help = 'out directory')
parser$add_argument('-n', dest = 'name', help = 'sample name')

opts = parser$parse_args()

## lib
library(Seurat)
library(dplyr)

## Merge
setwd(opts$input)
obj.list <- list.files(opts$input, pattern = '\\.rds$')

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})


obj.anchors <- FindIntegrationAnchors(object.list = obj.list,
                                      anchor.features = features,reduction = "rpca")

obj.combined <- IntegrateData(anchorset = obj.anchors)

DefaultAssay(obj.combined) <- "integrated"

obj.combined <- ScaleData(obj.combined)
obj.combined <- RunPCA(obj.combined, npcs = 30)
obj.combined <- RunUMAP(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- FindNeighbors(obj.combined, reduction = "pca", dims = 1:30)
obj.combined <- FindClusters(obj.combined, resolution = 0.5)

saveRDS(obj.combined,paste0(opts$out, '/', opts$name, "_integrated.rds"))
