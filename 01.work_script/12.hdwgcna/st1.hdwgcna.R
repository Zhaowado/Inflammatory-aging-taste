# single-cell analysis package
library(Seurat)
# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
#enableWGCNAThreads(nThreads = 3)

robj <- readRDS("../../02.sub_celltype/01.taste.cells/taste.final.filter.rds")
DefaultAssay(robj) <- "RNA"

robj <- subset(robj,timepoint == "M24")

#Set up Seurat object for WGCNA  (after this step, do not support subsetting obj)
robj <- SetupForWGCNA(
  robj,
  gene_select = "fraction", # the gene selection approach,还可以是variable（hvg）和custom
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  group.by = "celltype2",
  wgcna_name = "aging_taste"
)

#Construct metacells in each group 
obj <- MetacellsByGroups(
  robj,
  group.by = "celltype2",
  reduction = "pca",
  k = 10, #The number of nearest neighbors cells to be aggregated 
  max_shared = 10, # maximum number of shared cells between two metacells
  min_cells = 30, # exclude groups that are smaller than a specified number of cells. this cell mean single cell!! Errors are likely to arise if the selected value for min_cells is too low. Default 100 
  ident.group = 'celltype2'
)

#normalize metacell expression matrix
seurat_obj <- NormalizeMetacells(obj)

################################
#Co-expression network analysis
################################

#Set up the expression matrix, 通过group_name选择感兴趣的细胞类型（们）
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Type I MTC","Type II MTC","Type III MTC"), # the name of the group of interest in the group.by column
  group.by = "celltype2",
  assay = "RNA",
  slot = "data"
)

#Select soft-power threshold
#TestSoftPowers 函数模拟了共表达网络在不同软实力阈值下与无标度图的相似程度
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid" This should be consistent with the network chosen for ConstructNetwork
)

#Co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  overwrite_tom = TRUE,
  tom_name = 'MTC_TPC' # name of the topoligical overlap matrix written to disk
)

pdf("module.pdf",width = 6,height = 3)
PlotDendrogram(seurat_obj, main='MTC and TPC hdWGCNA Dendrogram')
dev.off()

####计算模块特征基因,group.by.vars对 ME 应用 Harmony 批量校正，从而产生协调的模块特征基因 (hME)
seurat_obj <- ModuleEigengenes(
  seurat_obj
)

#对每个模块的特征基因进行排名，基于每个模块的特征基因计算模块连接性（kME）最好与SetDatExpr一致
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'celltype2', group_name = c("Type I MTC","Type II MTC","Type III MTC")
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "MTC and TPC"
)

#可视化特征基因的排名
pdf("module_gene.pdf",width = 12,height = 8)
PlotKMEs(seurat_obj, ncol=5)
dev.off()

pdf("module_umap.pdf",width = 14,height = 8)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)
wrap_plots(plot_list, ncol=4)
dev.off()

pdf("module_corr.pdf")
ModuleCorrelogram(seurat_obj)
dev.off()

saveRDS(seurat_obj,"M24_MTC_TPC.rds")


ModuleNetworkPlot(
  seurat_obj, 
  outdir='ModuleNetworks', # new folder name
  n_inner = 15, # number of genes in inner ring
  n_outer = 20, # number of genes in outer ring
  n_conns = 500,     # show all of the connections:Inf
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)
#############
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
pdf("merge_select_network.pdf")
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods[c(1,8,16)] # only select 5 modules
)
dev.off()

pdf("merge_all_network.pdf")
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 5, n_other=8, #展示前5个枢纽基因和9哥其他基因
  edge_prop = 0.75,
  mods = 'all'
)
dev.off()
