####Fig2A####
library(ggplot2)
library(Seurat)
obj <- readRDS("../taste.integrated.rds")
meta <- read.table("sub_anno_all_9.29.txt",sep = "\t")
obj@meta.data <- meta

DimPlot(obj,group.by = "class",cols=c("#BBDD78","#caa2f4","#5F97C6"),label = F,raster=FALSE)
ggsave("umap_class_all.pdf",width = 8,height = 6)

####Fig2B-D####
colors <- read.csv("../colors.csv")
color <- colors$color
names(color) <- colors$figname
Idents(obj) <- "figname"
Idents(obj) <- factor(Idents(obj),levels = c("SC","PC","AC","MFB","MF","TMC","TC","EC","TPC","SDC","MTC","MAC","BC","SEC","EBC","FLF","PFC","SGC"))
VlnPlot(obj,features = c("Epcam","Krt13","Vim"),cols = color,pt.size = 0,ncol = 3)
ggsave("vln.pdf",width = 14,height = 2)

FeaturePlot(obj,features= c("Epcam","Krt13","Vim"),raster = T,min.cutoff="q1",ncol = 1)
ggsave("FeaturePlot.pdf",width = 5,height = 12)

####Figure2E####
library(pheatmap)
library(readxl)

tis <- "E"
sgo <- read_xlsx(paste0(tis,"-GO-ren-select-20250411.xlsx"),sheet = "Sheet1")
  sgo <- as.data.frame(sgo)
  rownames(sgo) <- sgo$way
  sgo <- sgo[,-1]

  anno_col <- data.frame(Timepoint = c("P0","P14","M2","M12","M24"))
  anno_col$Timepoint <- factor(anno_col$Timepoint)
  rownames(anno_col) <-  colnames(sgo)

  anno_color <- list(Timepoint = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                                   M12='#ff69b4',M24="#E66101"))

  if (tis == "E"){
    gap <- 16
  }else if (tis == "K"){
    gap <- 16
  }else if (tis == "M"){
    gap <- 13
  }


  pdf(paste0(tis,"_select_GSVA_go.score_timepoint.pdf"),width = 10, height = 14)
  print(pheatmap(sgo,cluster_rows = F,cluster_cols = F,scale = "row",
                 #treeheight_row = 5, #聚类树高度
                 fontsize=14,fontsize_row = 14, fontsize_col = 14,
                 cellwidth = 25,cellheight = 25,  #单元格长宽
                 annotation_colors = anno_color,
                 annotation_col = anno_col,
                 gaps_row = gap,
                 cutree_rows = 2,
                 colorRampPalette(c("#053061", "#2166AC", "#4393C3","#92C5DE","#D1E5F0","#FFFFFF","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(100)))
  dev.off()

####Figure2F####
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)

#### Load data ----
obj <- readRDS("../../../02.sub_celltype/01.taste.cells/recluster_tastecells.rds")
obj <- subset(obj,figname == "MTC" & timepoint %in% c("M2","M24"))
Idents(obj) <- "timepoint"
Idents(obj) <- factor(Idents(obj), levels = c("M2","M24"))
genes <- c("Edn1","Hnrnpa1","Hsf1", #salt
            "Azgp1","Car6","Pip", #bitter
            "Asic3","Scnn1a","Pkd1l3", #sour
            "Reep2","Tas1r3","Itpr3", #sweet
            "Tas1r1","Gnat1","Calhm1")#umami
p<-DotPlot(obj,unique(genes),assay="RNA",group.by = "timepoint")
data<-p$data

plot4 <- ggplot(data,aes(features.plot, id, fill=avg.exp.scaled, size=0))+
  geom_point(shape = 21,colour="black", size=6)+
  scale_fill_gradient2(low = "#4b84b3",mid = "white",high = "#E09293",midpoint = 0)+theme_bw()+
  theme(axis.text.y = element_text(size=8,color="black"),
        axis.text.x = element_text(size=8,color="black",face="italic",angle=90,vjust=0.5,hjust = 1),
        strip.background = element_rect(color="black",size=1,linetype="solid"),
        axis.line=element_line(color="black",size=0),
        strip.text.y.left=element_text(angle=0,face="bold"),
        legend.position="top",
        legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        panel.spacing.x = unit(0.5, "cm"),
        panel.spacing.y = unit(0.2, "cm"))+
  labs(x = "", y = "")+
  scale_y_discrete(position = "left")
pdf('Dotplot_selectmarker_long.pdf',width = 6, height = 4)
plot4
dev.off()

####Figure2G####
library(Seurat)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(stringr)

#### Load data ----
setwd("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/02.GSVA/06.all_compare_time_M2-M24/")
output.dir <- paste0("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/02.GSVA/06.all_compare_time_M2-M24")  

selected_way.go2 <- c(
  "Sensory perception of bitter taste",
  "Sensory perception of sour taste",
  "Cellular response to salt",
  "Sensory perception of sweet taste",
  "Sensory perception of umami taste")

gsva.go2 <- read.csv("gsva_GO.csv")
sub_data <- gsva.go2 %>% 
  filter(way %in% selected_way.go2) 


library(reshape2)
data <- melt(sub_data, id.vars = "way")

p <- ggplot(data,aes(x=way,y=value,fill=variable))+
  geom_col(position = position_dodge(width = 0.8), 
           width = 0.8, size = 0.3, colour = 'black') +
  # geom_hline(yintercept=c(log10(0.05),-log10(0.05)),linetype="dashed",color="black")+
  scale_fill_manual(values=c("M24"="#f89f68","M2"="#4b84b3"))+
  labs(x="",y="score")+
  theme_bw()+
  theme(
    panel.grid=element_blank(),
    # panel.border=element_rect(linewidth=0.2),
    axis.text.x=element_text(size=8,color="black",angle=0,vjust=0.5,hjust = 0.5),
    axis.text.y=element_text(size=8,color="black"),
    axis.title=element_text(size=8),
    # aspect.ratio=0.75,
    # plot.margin=margin(t=1,r=1,b=1,l=1,unit="cm"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(color = "black",size = 8),
    legend.key.size = unit(0.4, 'cm')
  )+
  coord_flip()

pdf(file.path(output.dir, paste0("all_barplot_tasteway-wd.pdf")), width = 4, height = 2.8)
print(p)
dev.off()

####Figure 2H####
library('plyr')
library('dplyr')
library('ggalluvial')
library('reshape2')
library('tidyverse')
library('data.table')
library('ggsci')

arg <- commandArgs(trailingOnly = TRUE)
prefix=arg[1]
meta <- read.table("../../sub_anno_all_9.29.txt",sep = "\t")
epi <- meta[meta$class == "Keratinized epithelium",]

Sample_cluster_num<-acast(epi,figname~timepoint,length)
d=apply(Sample_cluster_num,2,sum) #列求和
Sample_cluster_percent=sweep(Sample_cluster_num,2,d,'/') #每列÷列求和

meta2=as.data.frame(as.table(Sample_cluster_percent))
colnames(meta2)=c("Cluster","Sample","Freq")
meta2$Sample=factor(meta2$Sample,levels=c("P0","P14","M2","M12","M24"))

#label标签
add_label=function(x,y,z,tag){#数据框，分组值，显示阈值,显示标签
  data=get(x)
  Taxonomies=ddply(data,y, transform,value= rev(Freq))
  Taxonomies=ddply(Taxonomies,y, transform,label_y=cumsum(value) - 0.5*value)
  Taxonomies$label=Taxonomies[,{{tag}}]  #数据框列名为参数，使用{{}}
  for(i in 1:nrow(Taxonomies)){  if(Taxonomies[i,"Freq"] > z){  #只标记宽柱子
    Taxonomies[i,"label"] = Taxonomies[i,"label"]
  }else{
    Taxonomies[i,"label"] = NA
  }
  }
  Taxonomies=ddply(Taxonomies,y, transform,label=rev(label))
  return(Taxonomies)
}

colors <- read.csv("/dellfsqd2/ST_LBI/USER/zhaowandong/projects/taste/01.chongci/01.scRNA/01.annotation/colors.csv")
cluster_Palette <- colors$color
names(cluster_Palette) <- as.character(colors$figname)

pdf(paste0(prefix,"_Sample_Cluster_barplot.pdf",sep=''),width=10,height=6)

p1=meta2 %>% ggplot(aes(x=Sample, y=Freq, fill=Cluster, stratum=Cluster, alluvium=Cluster)) +
  geom_flow(width = 0.35, curve_type = "linear", alpha=0.5) +
  geom_stratum(width = 0.33, alpha=0.9, color=NA) +
  scale_fill_manual(values = cluster_Palette) +
  geom_alluvium(width = 0.33,curve_type = "linear", fill=NA, color="#f4e2de") +
  # scale_fill_nejm()+
  scale_y_continuous(expand = c(0, 0),
                     name = NULL) +
  labs(x="",y="") +
  theme_classic() +
  theme(plot.background = element_rect(fill='white', color='white'),
        panel.grid = element_blank(),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
print(p1)
dev.off()

####Figure2I####
library(ggplot2)
colors <- read.csv("../../colors.csv")
mapping <- colors$figname
names(mapping) <- colors$celltype
txt$figname <- mapping[txt$celltype2]
cluster_Palette <- colors$color
names(cluster_Palette) <- as.character(colors$figname)
p <- ggplot(txt,aes(x,y,color = figname))+
  geom_point(size=1)+
  scale_color_manual(values = cluster_Palette)+
  theme_bw()+
  theme(legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.x=element_text(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.y=element_text()) + coord_fixed(ratio = 1) + labs(x=NULL,y=NULL)
####Figure2J####
cutoff <- quantile(distance_df$distance, 0.95)
distance_df$distance_truncated <- pmin(distance_df$distance, cutoff)  # 截断数据
ggplot(distance_df, aes(x = timepoint, y = distance_truncated, fill = timepoint)) +
  geom_boxplot() +
  labs(
    x = "timepoint",
    y = "Distance to nearest BC (µm)"
  ) +
  theme_bw()+
  scale_fill_manual(values = c("M2" = "#E7B800", "M24" = "#E66101")) +
  theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 14)) +
  stat_compare_means(comparisons = list(c("M2", "M24")), method = "t.test",label = "p.signif")
ggsave("Distance_MTC2BC.pdf")

####Figure2K,L####
library(Seurat)
library(ggsignif)
library(ggpubr)
library(reshape2)
library(ggplot2)

p <- ggplot(data = obj@meta.data, aes(x=timepoint,y=Score1,fill=timepoint)) +
  geom_violin(trim = T, scale = "width",lwd=0.2,alpha=0.7)+
  geom_boxplot(width=.2,lwd=0.3,alpha=0.7,aes(fill = timepoint))+
  scale_fill_manual(values = c('#66cdaa','#1f78b4',"#E7B800",'#ff69b4',"#E66101"))+
  theme_bw()+
  theme(axis.text = element_text(colour="black",size=8),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(colour="black",size=0.2),
        panel.border = element_blank(), #去掉边框
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="",y="positive_regulate_epi_migration")+
  geom_signif(comparisons = list(c("M2","M24")),
              map_signif_level = T,
              test = "t.test"
  )+NoLegend()

####Figure2M####
shh <- as.data.frame(obj@assays$RNA@data["Shh",])
colnames(shh)[1] <- "exp"
shh <- subset(shh,exp > 0)

obj$cell <- rownames(obj@meta.data)
robj <- subset(obj, cell %in% rownames(shh))

df.im <- as.data.frame(table(robj$timepoint))
df.all <- as.data.frame(table(obj$timepoint))
df <- merge(df.im,df.all, by = "Var1",suffixes = c(".im", ".all"))
df <- df[,c("Var1","Freq.im","Freq.m")]
df <- as.data.frame(t(df[,-1]))
df <- df[,c("P0","P14","M2","M12","M24")]

rownames(df)[1] <- "immature"
rownames(df)[2] <- "mature"

df$state <- rownames(df)
df_long <- melt(df,id.vars="state")

ggplot(df_long, aes(x = variable, y = value, fill = state)) +
  geom_col(position = "fill",stat = "identity") +
  labs(x = "timepoint", y = "cell number") +
  theme_minimal() +
  scale_fill_manual(values = c("immature" = "#fdded7","mature" = "#c1e0db")) +
  theme(legend.title=element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        panel.border = element_blank(),
        )

ggsave("MTC_im_m.barplot.pdf")
