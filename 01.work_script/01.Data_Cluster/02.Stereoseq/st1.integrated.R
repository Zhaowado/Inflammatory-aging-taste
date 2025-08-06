library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

parser = argparse::ArgumentParser(description = "rscript for spatial data and corresponding single cell data integrated")

parser$add_argument('-i', dest = 'input', help = 'input directory')
parser$add_argument('-o', dest = 'out', help = 'out directory')
parser$add_argument('-s', dest = 'sample', help = 'sample name')

opts = parser$parse_args()

# 提取参数
input_dir <- opts$input_dir
output_dir <- opts$output_dir
sample <- opts$sample_name

### 动态生成时间点列表
timepoint <- c(
  paste0("P0_", sample),
  paste0("P14_", sample),
  paste0("M2_", sample),
  paste0("M12_", sample),
  paste0("M24_", sample)
)

### load scRNAseq data and annotation
obj <- readRDS(file.path(input_dir, "taste.integrated.rds"))
cell_anno <- read.csv(file.path(input_dir, "01.anno/anno_7.29.csv"))

cell_annos <- cell_anno$celltype
names(cell_annos) <- as.character(cell_anno$cluster)
obj$celltype <- cell_annos[as.character(obj$integrated_snn_res.0.5)]

### 子集数据
sc <- subset(obj, orig.ident %in% timepoint)

### load corresponding spatial data
obj2 <- readRDS(file.path(input_dir)
obj3 <- merge(sc, obj2)

### 数据整合流程
DefaultAssay(obj3) <- "RNA"
cv.list <- SplitObject(obj3, split.by = "orig.ident")

# 标准化和特征选择
cv.list <- lapply(cv.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = cv.list)
cv.anchors <- FindIntegrationAnchors(object.list = cv.list, anchor.features = features)
cv.combined <- IntegrateData(anchorset = cv.anchors)

### 下游分析
DefaultAssay(cv.combined) <- "integrated"
cv.combined <- ScaleData(cv.combined, verbose = TRUE)
cv.combined <- RunPCA(cv.combined, npcs = 30, verbose = TRUE)
cv.combined <- RunUMAP(cv.combined, reduction = "pca", dims = 1:30)
cv.combined <- FindNeighbors(cv.combined, reduction = "pca", dims = 1:30)
cv.combined <- FindClusters(cv.combined, resolution = 0.5)

### 保存结果（路径基于输出目录和样本名称）
saveRDS(cv.combined, file.path(output_dir, paste0(sample, ".integrated.rds")))
write.table(
  cv.combined@meta.data,
  file.path(output_dir, paste0(sample, ".meta.txt")),
  sep = "\t",
  quote = FALSE
)

### 绘图参数
cluster_number <- length(unique(cv.combined@meta.data$seurat_clusters))
cluster_Palette <- if (cluster_number <= 25) {
  c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
    'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
    'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
    'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
} else {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
}

### 绘制UMAP
pdf(file.path(output_dir, "umap.pdf"), width = 16, height = 12)
print(DimPlot(cv.combined, reduction = 'umap', label = TRUE, cols = cluster_Palette, repel = TRUE, raster = FALSE))
dev.off()

pdf(file.path(output_dir, "umap_anno_sc.pdf"), width = 16, height = 12)
print(DimPlot(cv.combined, reduction = 'umap', group.by = "celltype", label = TRUE, cols = cluster_Palette, repel = TRUE, raster = FALSE))
dev.off()

### 寻找标记基因
DefaultAssay(cv.combined) <- 'RNA'
markers <- FindAllMarkers(cv.combined, min.pct = 0.1, logfc.threshold = 0.25)
markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]

write.table(
  markers,
  file.path(output_dir, paste0(sample, ".AllMarkers.xls")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(
  topn,
  file.path(output_dir, paste0(sample, ".AllMarkers.top30.xls")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

### 绘制批次效应图
pdf(file.path(output_dir, "umap_timepoint.pdf"), width = 60, height = 5)
print(DimPlot(cv.combined, reduction = 'umap', group.by = 'seurat_clusters', split.by = 'orig.ident', cols = cluster_Palette, ncol = 12, raster = FALSE))
dev.off()
