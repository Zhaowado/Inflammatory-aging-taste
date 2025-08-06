args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]

library(Seurat)
spa <- readRDS(paste0("../00.data/",sample,".sp.rds"))
DefaultAssay(spa) <- "RNA"

sc <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/scRNA/04.cluster/02.integrate_all/taste.integrated.rds")
sc <- subset(sc, orig.ident %in% c(paste0("P0_",sample),paste0("P14_",sample),paste0("M2_",sample),paste0("M12_",sample),paste0("M24_",sample)))

i2 <- FindTransferAnchors(
  reference = sc,
  query = spa,
  features = rownames(spa),
  reduction = 'cca',
  reference.assay = 'RNA',
  query.assay = 'RNA'
)

saveRDS(i2,paste0(sample,"_anchors.rds"))

refdata <- GetAssayData(
  object = sc,
  assay = 'RNA',
  slot = 'data'
)

imputation <- TransferData(
  anchorset = i2,
  query = spa,
  refdata = refdata,
  weight.reduction = 'pca'
)

saveRDS(imputation,paste0(sample,"_gene.pre.rds"))
