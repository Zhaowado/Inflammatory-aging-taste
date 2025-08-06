library(CellChat)
library(patchwork)
library(Seurat)

#### Load data ---- 
output.dir <- ("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/04.cellchat/06.M24-M2_MF-TC-TPC-MTC")
obj1 <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/taste.integrated.rds")
meta <- read.table("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/01.bigclass/bigclass_annotation.txt")
obj2 <- AddMetaData(object = obj1, metadata = meta)
DefaultAssay(obj2) <- "RNA"

Idents(obj2) <- "timepoint"
obj3 <- subset(obj2, idents=c("M2","M24"))

Idents(obj3) <- "figname"
seurat.obj <- subset(obj3, idents=c("MF","TC","TPC","MTC"))

names(table(Idents(seurat.obj)))
print(table(Idents(seurat.obj)))

for (kk in unique(seurat.obj$timepoint)){
  obj <- subset(seurat.obj,timepoint %in% kk)
  batch <- unique(obj@meta.data$timepoint)
  expr <- obj@assays$RNA@data
  data.input <- expr
  dim(data.input)
  data.input[1:4, 1:4]
  meta <- as.data.frame(Idents(obj))
  colnames(meta) <- "labels"
  
  unique(meta$labels) # check the cell labels
  cellchat <-
    createCellChat(object = data.input,
                   meta = meta,
                   group.by = "labels")
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <-
    setIdent(cellchat, ident.use = "labels") 
  levels(cellchat@idents) # show factor levels of the cell labels
  groupSize <- 
    as.numeric(table(cellchat@idents)) # number of cells in each cell group
  
  CellChatDB <-
    CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)       #展示配受体饼图
  
  dplyr::glimpse(CellChatDB$interaction)
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <-
    subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  future::plan("multisession", workers = 1) # do parallel,或者multisession
  cellchat <- updateCellChat(cellchat)
  cellchat <-
    identifyOverExpressedGenes(cellchat) # take a short time
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse) # take a short time,PPI.human
  
  #### Compute the communication probability and infer cellular communication network ----
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #### Extract the inferred cellular communication network as a data frame ----
  #pass
  #### Infer the cell-cell communication at a signaling pathway level ----
  cellchat <- computeCommunProbPathway(cellchat)
  #### Calculate the aggregated cell-cell communication network ----
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents)) #每个细胞类型的数量
  saveRDS(cellchat,file.path(output.dir,file=paste0(batch,"_cellchat.rds")))
