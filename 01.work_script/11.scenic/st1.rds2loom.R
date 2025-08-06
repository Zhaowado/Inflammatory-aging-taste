library(Seurat)
library(SeuratDisk)

obj <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/02.sub_celltype/01.taste.cells/recluster_tastecells.rds")
obj <- subset(obj,timepoint == "M24")
obj$celltype2 <- gsub(" ",".",obj$celltype2)
obj.list <- SplitObject(obj, split.by = 'celltype2')
for (cell in unique(obj$celltype2)) {
  obj.list[[cell]] <- subset(obj.list[[cell]],celltype2 == cell)
  obj.loom <- as.loom(obj.list[[cell]], assay = 'RNA', filename = paste0(cell,'.loom'), verbose = TRUE,overwrite = T)
  obj.loom$close_all()
  saveRDS(obj.list[[cell]],paste0(cell,'.rds'))
}

library(reticulate)
use_condaenv('/dellfsqd2/ST_OCEAN/USER/guoxinyu1/bin/miniconda3/envs/st/')
sc <- import("scanpy")

for (time in unique(obj$celltype2)) {
  dat <- sc$read_loom(paste0(time,'.loom'))
  dat$X <- dat$layers['counts']
  dat$write_loom(paste0(time,'.loom'))
}
