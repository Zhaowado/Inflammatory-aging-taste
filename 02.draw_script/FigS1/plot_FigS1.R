########FigS1a########
library(Seurat)
library(ggplot2)
colors <- read.csv("../02.spatialplot/colors.csv")

for (sample in c("CV","FL","FF")) {
obj <- readRDS(paste0(sample,".integrated.rds"))
meta <- read.table(paste0(sample,".anno.txt"),sep = "\t")
obj@meta.data <- meta
obj$celltype[which(obj$celltype == "Immature taste cells")] <- "Mature taste cells"
figname <- colors$figname
names(figname) <- colors$celltype
obj$figname_sc <- figname[obj$celltype]


color <- colors$color
names(color) <- colors$figname

DimPlot(obj,group.by = "figname_sc",label = F,cols = color,raster = F)
ggsave(paste0(sample,"_scUMAP.pdf"),width = 8,height = 6)

obj$figname_all <- figname[obj$celltype2]
DimPlot(obj,group.by = "figname_all",label = F,cols = color,raster = F)
ggsave(paste0(sample,"_allUMAP.pdf"),width = 8,height = 6)
DimPlot(obj,group.by = c("figname_sc","figname_all"),label = F,cols = color,raster = F)
ggsave(paste0(sample,"_all.split.UMAP.pdf"),width = 18,height = 6)
}

########FigS1bcd########
library(Seurat)
library(cowplot)
library(ggplotify)

SpatialFeaturePlot_gene1 <- function(x,features,name.assay,slice,ppt=0.5){
  library(ggplot2)
  if (name.assay == "SCT"){
    title <- "True_exp"
  }else{
    title <- "Pre_exp"
  }
  da1 <- x@meta.data
  #  da2 <- x@assays$SCT@data[features,]
  da2 <- x@assays[[name.assay]]@data[features,]
  da1$expression <- da2
    da3 <- subset(da1,orig.ident == slice)
    p1 <- ggplot(da3,aes(x,y,color=expression))+
      geom_point(size=ppt)+
      scale_color_gradient2(low = "#6baed6",mid = "#ffeda0",high = "#e31a1c",midpoint = sum(range(da2))/2)+
      theme_bw()+ coord_fixed(ratio = 1)+labs(x = "",y = "")+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank())
    p <- as.ggplot(p1) +
      theme(legend.position=c(0.8,1.03), legend.direction="horizontal")+
      labs(x=unique(da3$stage))+
      ggtitle(title)+
      theme(legend.title=element_blank(),
            legend.key.size = unit(50, "pt"), #图例大小
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.x=element_text(size = 30),
            axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.y=element_blank())
    print(p)
}


