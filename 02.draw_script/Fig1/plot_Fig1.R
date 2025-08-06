########Fig1b########
library(Seurat)
library(ggplot2)
library(tidydr)
library(dplyr)
obj <- readRDS("../taste.integrated.rds")
meta <- read.table("bigclass_annotation.txt",sep = "\t")
obj@meta.data <- meta
obj$CellType <- paste0(obj$celltype," (",obj$figname, ")")

colors <- read.csv("../colors.csv")
colors$CellType <- paste0(colors$celltype," (",colors$figname, ")")

color.palette <- colors$color
names(color.palette) <- colors$CellType
obj$CellType <- factor(obj$CellType,levels = c("Adipose cells (AC)","Basal cells (BC)","Epithelial barrier cells (EBC)","Endothelial cells (EC)","Filiform (FLF)","Mucous acinar cells (MAC)","Macrophages (MF)","Myofibroblasts (MFB)","Mature taste cells (MTC)","Pericyte cells (PC)","Proliferating cells (PFC)","Schwann cells (SC)","Serous demilune cells (SDC)","Strat
ified epithelial cells (SEC)","Serous glandular cells (SGC)","T cells (TC)","Tongue muscle cells (TMC)","Taste progenitor cells (TPC)"))

DimPlot(obj,group.by = "CellType",label = F,cols = color.palette,raster=T,pt.size = 1) +
  tidydr::theme_dr(xlength = 0.2,
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        legend.key.height = unit(0.4,'cm'),
        legend.key.width = unit(0.5,'cm'),
        axis.title = element_text(face = 2,hjust = 0.03))
ggsave("umap.pdf",width = 10,height = 8)

nums <- data.frame(table(meta$figname))
colnames(nums)[1:2] <- c("cell_type","num")
nums <- nums %>% mutate(log = log(num))

fills <- cols[names(cols) %in% unique(meta$figname)]
p2 <- ggplot(nums,aes(x = cell_type, y = log,fill = cell_type)) +
  geom_bar(stat = "identity") + labs(y = "log cell numbers") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, colour = "black",hjust = 1),
        panel.grid=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(size=1,colour="black"))
ggsave(paste0("allcelltype_cell.num.pdf"),p2,width = 8,height = 4)
########Fig1c########
library(Seurat)
library(tidyverse)
library(viridisLite)
library(ggplot2)
obj <- readRDS("taste.integrated.rds")
DefaultAssay(obj) <- "integrated"
meta <- read.table("sub_anno_all_9.9.txt",sep = "\t")
obj@meta.data <- meta
meta$UMAP1 <- Embeddings(obj,reduction = "umap")[,1]
meta$UMAP2 <- Embeddings(obj,reduction = "umap")[,2]

meta$timepoint <- factor(meta$timepoint,levels = c("P0","P14","M2","M12","M24"))
dat_bg <- meta[,-(which(colnames(meta)=="timepoint"))]

my_color_palette <- colorRampPalette(c("#FFFFFF", "#e0533f","black"))

meta$sample <- factor(meta$sample,levels = c("CV","FL","FF"))
dat_bg <- meta[,-(which(colnames(meta)=="sample"))]
my_color_palette <- colorRampPalette(c("#FFFFFF", "#e0533f","black"))
density_plot <- ggplot(meta, aes(x=UMAP1, y=UMAP2)) +
  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=F) + #ndensity calculates the normalized density for each sample--otherwise density would be affected by the number of cells for each sample, which is variable
  geom_point(data=dat_bg, shape=16, size=0.1, alpha=0.2, color="white") +
  #scale_fill_gradientn(colours=my_color_palette(100), name="Density") +
  scale_fill_gradientn(colours=viridisLite::mako(100), name="Density") +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_wrap(~sample, ncol=3) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12, color="black"),
        axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"))
ggsave("density_sample.png",width = 15,height = 5,density_plot)

########Fig1d########
library(ggplot2)
meta <- read.table("bigclass_annotation.txt",sep = "\t")
df <- data.frame(table(meta$figname,meta$orig.ident))

colnames(df)[1:2] <- c("cell_type","Tissue1")
df$Tissue1 <- factor(df$Tissue1,levels = c("P0_CV","P14_CV","M2_CV","M12_CV","M24_CV",
                                           "P0_FL","P14_FL","M2_FL","M12_FL","M24_FL",
                                           "P0_FF","P14_FF","M2_FF","M12_FF","M24_FF"))
colors <- read.csv("../colors.csv")
mycolors1 <- colors$color
names(mycolors1) <- colors$figname
p1 <- ggplot(df, aes( x = Tissue1, y = Freq, fill = cell_type,))+
  geom_bar( position = "fill",stat = "identity",width = 0.7)+theme_bw()+    #position "fill"元素堆叠且高度标准化为1,"stack"为真实高度
  scale_fill_manual(values=mycolors1) +
  labs(x = "",y = "proportion") +
  theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45,size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
       # panel.border = element_blank(), # 去除绘图面板的边框  
       # axis.line = element_blank(), # 去除坐标轴线  
       # axis.ticks = element_blank() # 去除坐标轴刻度线
       )


ggsave(paste0("celltype_rio_barplot_all.pdf"),p1,width = 18,height = 6)

########Fig1e########
library(ggplot2)
library(Seurat)

obj <- readRDS("../../taste.integrated.rds")
meta <- read.table("../bigclass_annotation.txt",sep = "\t")
obj@meta.data <- meta
DefaultAssay(obj) <- "RNA"
markers <- read.csv("markers_celltype_2.14.csv")
markers <- markers[rev(seq_len(nrow(markers))), ]

for (celltype in unique(markers$celltype)) {
  gene<-intersect(rownames(obj@assays$RNA@data),markers[which(markers$celltype==celltype),"marker"])
  print(celltype,length(gene))
  for (j in gene) {
    pdf(paste(celltype,j,"featureplot.pdf",sep = "_"))
    p1<-FeaturePlot(obj,features = j,pt.size=0.3, cols = c("grey77","red"), order = T,raster=FALSE,slot = "data")
    print(p1)
    dev.off()
  }
}

########Fig1f########
colors <- read.csv("../colors.csv")

mapping <- colors$figname
names(mapping) <- colors$celltype
txt$figname <- mapping[txt$celltype2]

cluster_Palette <- colors$color
names(cluster_Palette) <- as.character(colors$figname)

pall <- list()
for (id in sampleid) {
  da1 <- txt[txt$orig.ident == id,]
  p <- ggplot(da1,aes(x,y,color = figname))+
    geom_point(size=1)+
    scale_color_manual(values = cluster_Palette)+
    theme_bw()+
    theme(legend.title=element_blank(),
          # legend.key.size = unit(50, "pt"), #图例大小
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.x=element_text(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.y=element_text()) + coord_fixed(ratio = 1) + labs(x=id, y = unique(da1$stage))
  pall[[id]] <- as.ggplot(p)
 # + theme(legend.position.inside=c(0.8,1.03), legend.direction="horizontal",legend.key.size = unit(100, "pt"))
}
pdf(paste0(sample,"_SpatialDimplot_allcelltype.pdf"),width = 8*length(unique(txt$orig.ident)),height = 12)
print(plot_grid(plotlist = pall,ncol = 4))
dev.off()

########Fig1g########
#### load packages ---- 
library(Seurat)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
output.dir <- paste0("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/09.Fig1/")  

#### Load data ----
obj <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/taste.integrated.rds")
meta <- read.table("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/01.bigclass/bigclass_annotation.txt",sep = "\t")
obj@meta.data <- meta
DefaultAssay(obj) <- "RNA"
Idents(obj) <- "figname"

markers <- read.table("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/01.bigclass/04.DEGs/taste.celltype.AllMarkers.xls",sep = "\t")
colnames(markers) <- markers[1,]
markers <- markers[-c(1),]

markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC)

markers$cluster<-factor(markers$cluster,levels = c("AC","BC","EBC","EC","FLF","MAC","MF","MFB",
                                                   "MTC","PC","PFC","SC","SDC","SEC","SGC","TC",
                                                   "TMC","TPC"))
markers <- markers[order(markers$cluster), ]

genes <- unique(markers$gene) 
aver_dt <- AverageExpression(obj,
                             features= genes,
                             group.by = 'figname',
                             slot= 'data') 
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt[1:6,1:6]
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))
aver_dtt <- t(scale(t(aver_dt)))

celltype_info <- sort(obj$figname)
table(celltype_info)
celltype_info<-factor(celltype_info,levels = c("AC","BC","EBC","EC","FLF","MAC","MF","MFB",
                                               "MTC","PC","PFC","SC","SDC","SEC","SGC","TC",
                                               "TMC","TPC"))

colors <- read.csv("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/colors.csv")
mycolors1 <- colors$color
names(mycolors1) <- colors$figname

cell_type <- data.frame(colnames(aver_dtt))
colnames(cell_type) <- 'CellType'

col_anno <- HeatmapAnnotation(df = cell_type,
                              show_annotation_name = F,
                              show_legend = FALSE, #不显示legend
                              gp = gpar(col = 'white', lwd = 1),
                              simple_anno_size = unit(0.3, "cm"),
                              col = list(CellType = mycolors1))

markers_show <- c("Cd36","Aqp1",#"Egfl7", AC
                  "Krt5","Krt14",#"Krt15",  BC
                  "Sprr1a","Krt75",#"Cnfn","Krt6b",            EBC
                  "Lyve1", "Mmrn1",#   "Fgl2","Prox1",         EC
                  "Krt36","Krt76",#"Krt24",                     FLF
                  "Tff2","Muc5b",# "Smgc","Muc19",            MAC
                  "C1qa","C1qb",#MF"Fcer1g","Cd74",
                  "Col1a1","Bgn",#MFB"Thbs4","Mfap4",
                  "Krt8","Msln",#MTC,"Otop1""Gng13","Hes6", 
                  "Acta2","Rgs5", #PC"Myl9","Tagln",
                  "Mki67","Ube2c",#PFC"Birc5", "Top2a",
                  "Mpz","Plp1", #SC,"Gatm","Fgl2"
                  "Dcpp1","Dmbt1", #SDC"Ltf", "Dcpp3",
                  "Krt4","Krt13",#SEC"Mt4",  
                  "Bpifb1","Lipf",# "Bpifa2","Sval2",             SGC
                  "Il17a","Ptprc",        #"Icos","Srgn", "Rgs1", TC
                  "Tnni2","Myog",#"Acta1","Tnnt3",                TMC
                  "Ptch1","Fst"#,"Lgr5","Krt24","Hes6","Bdnf" ,"Stmn1"          TPC
)
gene_pos <- which(rownames(aver_dtt)%in%markers_show)


labs <- rownames(aver_dt)[gene_pos]

row_anno1 <- rowAnnotation(
  gene = anno_mark(at = gene_pos, 
                   labels = labs,
                   labels_gp = gpar(cex = 0.9, col = "black", fontface = 'italic',fontsize = 8))) # 

ht <- Heatmap(aver_dtt,name = 'expression',
              row_order = 1:nrow(aver_dtt),
              col = mycol2,
              cluster_columns = F,
              cluster_rows = F,
              show_row_names = F,
              show_column_names = T, 
              column_names_side = c('top'), 
              column_names_rot = 45, 
              column_names_gp = gpar(fontsize = 7), 
              right_annotation = row_anno1,
              top_annotation = col_anno,
              heatmap_legend_param = list(title = "Z-score", 
                                          title_gp=gpar(fontsize=8),
                                          labels_gp=gpar(fontsize=8),
                                          legend_direction="vertical")) 
pdf(file.path(output.dir,"Averageheatmap_expression-V2.pdf"),width=3.8,height=4.7)
print(ht)
dev.off()

########Fig1H########
#### load packages ---- 
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
output.dir <- paste0("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/09.Fig1/")  

#### Load data ----
markers <- read.table("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/01.bigclass/04.DEGs/taste.celltype.AllMarkers.xls",sep = "\t",header=T)

for (i in unique(markers$cluster)) {
  col<-c("SYMBOL","ENTREZID")
  marker <- markers[markers$cluster == i,] %>% .[.$p_val_adj < 0.05 & .$avg_log2FC >0.25,]
  da1 <- AnnotationDbi::select(org.Mm.eg.db,columns=col,keytype="SYMBOL",keys=as.character(marker$gene))
  ego<-enrichGO(OrgDb = "org.Mm.eg.db",gene=da1$ENTREZID,ont="BP",pvalueCutoff=0.05,readable=TRUE)
  write.csv(ego,file.path(output.dir,paste0("GO_",i,".csv")),row.names = F)
}
go_all <- NULL
GO2 <- NULL
for (i in unique(markers$cluster)) {
  GO1 <- read.csv(paste0("GO_",i,".csv"))
  GO1$celltype <- i
  go_all <- rbind(GO1,go_all)
  GO2 <- rbind(GO1[1:2,],GO2)
}
write.csv(go_all,file.path(output.dir,"allcelltype_GO.csv"))

#### Load data ----
colors <- read.csv("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/colors.csv")
mycolors1 <- colors$color
names(mycolors1) <- colors$figname

GO2 <- read.csv("allcelltype_GO_select-V2.csv",row.names = "X")
GO2$lgP <- -log10(GO2$pvalue)
GO2$Description <- (GO2$Description) |> str_to_sentence()

p <- ggdotchart(GO2,x = "Description",y = "lgP",color="celltype",
                palette = mycolors1,
                sorting = "asc", sort.by.groups = TRUE,
                add = "segments",add.params = list(color = "grey", size = 1), 
                group = "celltype",dot.size = 3,
                ggtheme = theme_pubclean()) + 
  font("x.text", size = 8, vjust = 0.5) + ylab("-Log10(P value)") + xlab("")+
  theme(axis.text.y = element_text(colour='black',size = 8,angle=90,hjust=0.5),
        axis.text.x = element_text(colour='black',size = 8),
        axis.title=element_text(size=8),
        # axis.line=element_line(size=0.2),
        axis.ticks=element_line(size=0.4),
        axis.ticks.length=unit(0.1,"cm"),
        legend.position = "none")

pdf(file.path(output.dir,"enrichment_way_bangplot-V3.pdf"),width = 5.4,height = 4.5)
print(p)
dev.off()
