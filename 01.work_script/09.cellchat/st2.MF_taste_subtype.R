library(CellChat)
library(patchwork)
library(Seurat)
library(ggpubr)
library(cowplot)
options(stringsAsFactors = FALSE)

sample <- "All"
obj <- readRDS("../../01.annotation/taste.integrated.rds")
meta <- read.table("../../01.annotation/sub_anno_all_9.9.txt",sep = "\t")

obj@meta.data <- meta
obj <- subset(obj,figname %in% c("MTC","TPC","MF"))

for (time in unique(obj$timepoint)){
robj <- subset(obj, timepoint == time)
##########prepare input
data.input <- GetAssayData(robj, assay = "RNA", slot = "data") 
meta <- data.frame(labels = robj@meta.data$celltype_sub,row.names = colnames(data.input))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
########Add cell information into meta slot of the object (Optional)
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

######Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

#######Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
####Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 5)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat,paste0(time,"_cellchat.rds"))
