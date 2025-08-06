library(Seurat)
library(ggplot2)
library(harmony)
library(dplyr)
parser = argparse::ArgumentParser(description = 'Script to recluster integrated sample obj')
parser$add_argument('-d', dest = 'dims', default = 20, help = 'dims for umap, default 20')
parser$add_argument('-r', dest = 'resolution', default = 0.5, help = 'cluster resolution, default 0.5')
parser$add_argument('-rds', dest = 'rds', help = 'merged rds file name')
parser$add_argument('-c', dest = 'cluster', help = 'the cluster you want to subset')
parser$add_argument('-b', dest = 'batch', help = 'which column used to correct batch effect')
parser$add_argument('-n', dest = 'normalize', default = 'sct', choices = c('sct', 'log'), help = 'normalization method for seurat, sct for sctransform, log for LogNormalize, default sct')
opts = parser$parse_args()

cluster <- as.character(opts$cluster)
####cluster 
obj <- readRDS(opts$rds)
robj <- subset(obj,seurat_clusters == cluster)
DefaultAssay(robj) <- "RNA"

obj.list <- SplitObject(robj, split.by = opts$batch)
##############################################################
# now let normalize the data one by one
##############################################################
normalize <- function(obj, method = 'sct', nfeatures = 2000){
  assay = 'RNA'
  if (method == 'sct'){
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- SCTransform(obj, method = "glmGamPoi", variable.features.n = nfeatures, 
                       assay = assay, verbose = FALSE)
  }else if (method == 'log'){
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- ScaleData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  }
  return(obj)
}


obj.list <- lapply(X = obj.list, FUN = function(x){
  x <- normalize(x, method = opts$normalize)
#  x <- RunPCA(x)
#  x <- RunUMAP(x, dims = 1:as.numeric(opts$dims))
})
print("sucessfully normalize the data one by one")
##############################################################
# now let intergrate all slices onto one
##############################################################
features <- SelectIntegrationFeatures(object.list = obj.list)
obj.integrated <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)], merge.data = TRUE)

DefaultAssay(object = obj.integrated) <- "SCT"
VariableFeatures(obj.integrated) <- features

##############################################################
# now let do cluster
##############################################################
obj.integrated <- RunPCA(object = obj.integrated, assay = 'SCT', verbose = FALSE)
obj.integrated <- RunHarmony(obj.integrated, group.by.vars = opts$batch, 
                             reduction = 'pca', assay.use = 'SCT', reduction.save = 'harmony'
                             )
obj.integrated <- RunUMAP(object = obj.integrated, reduction = "harmony", assay = 'SCT', dims = 1:as.numeric(opts$dims))

obj.integrated <- FindNeighbors(object = obj.integrated,dims = 1:as.numeric(opts$dims), reduction = "harmony")
obj.integrated <- FindClusters(object = obj.integrated, resolution = as.numeric(opts$resolution))

###################
# save metadata and rds
robj <- obj.integrated

cluster_Palette <- c('#33a02c','#ff69b4','#1f78b4','#ff4500','#b2df8a','#a6cee3','#a020f0','#ff7f00','#66cdaa','#db7093','#1e90ff',
                     '#fdbf6f','#6a3d9a','#ffd700','#b15928','#8dd3c7','#ffffb3','#bebada','#ffc1c1','#80b1d3','#fdb462','#b3de69',
                     '#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf',
                     '#999999','#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2','#8c510a',
                     '#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#ffff33','#a65628','#f781bf',
                     '#999999','#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#fddaec','#f2f2f2','#8c510a','#bf812d',
                     '#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#c51b7d')
DimPlot(robj,group.by=c(opts$batch,"seurat_clusters"),cols = cluster_Palette,label = T,repel = T,raster=FALSE)
ggsave("harmony_batch.pdf",width = 12,height = 6)
saveRDS(robj,paste0("Recluster",cluster,".rds"))

# find markers
DefaultAssay(robj) <- 'RNA'
markers <- FindAllMarkers(robj, min.pct = 0.1, logfc.threshold = 0.25)
markers <- markers[with(markers, order(cluster, -avg_log2FC)), ]
write.table(markers, paste0('recluster',cluster,'.AllMarkers.xls'), sep = '\t', quote = FALSE, row.names = F)

topn <- markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.table(topn, paste0('recluster',cluster,'.AllMarkers.top30.xls'), sep='\t', quote = FALSE, row.names = F)


###转移回去label
robj$seurat_clusters <- as.character(robj$seurat_clusters)

for (i in 1:length(unique(robj$seurat_clusters))-1) {
  robj$seurat_clusters[which(robj$seurat_clusters %in% c(as.character(i)))] <-paste0(cluster,"_",as.character(i))
}

obj$recluster <- obj$seurat_clusters
obj$recluster <- as.character(obj$recluster)
obj$recluster[which(obj$recluster == cluster)] <- robj$seurat_clusters

#######Todo: order how to custom
#obj$recluster <- factor(as.character(obj$recluster), levels = c(paste0("0_", 0:3), 1:44))

DimPlot(obj,group.by = "recluster",cols = cluster_Palette,label = T,repel = T,raster=FALSE)
ggsave("umap_ori.pdf")
