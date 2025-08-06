####Figure3D,E####
#### load packages ---- 
library(Seurat)
library(ggplot2)
library(ggVennDiagram)
output.dir <- ("~/practice6/test/Vein_wand")

#### 1. Load data ----
deg <- read.csv("~/practice6/test/Vein_wand/deg_updown.csv")

#### 2. 四组比较UP plot ----
celltypes <- unique(deg$celltype)
all_data <- data.frame() 
for (i in celltypes) {
  
  deg_up <- deg %>% filter(change == "UP") %>% filter(avg_log2FC > 0.5)
  deg_up_celltype <- deg_up %>% filter(celltype %in% i)
  
  A  <-  deg_up_celltype %>% filter(comparison == "M24_vs_P0") %>% pull(gene)
  B  <-  deg_up_celltype %>% filter(comparison == "M24_vs_P14") %>% pull(gene)
  C  <-  deg_up_celltype %>% filter(comparison == "M24_vs_M2") %>% pull(gene)
  D  <-  deg_up_celltype %>% filter(comparison == "M24_vs_M12") %>% pull(gene)
  
  x<-list(`M24 vs P0` = A,
          `M24 vs P14` = B,
          `M24 vs M2` = C,
          `M24 vs M12` = D)
  
  library(ggvenn)
  data <- list_to_data_frame(x)
  data$celltype <- i 
  p1 <- ggvenn(data = data, 
               show_percentage = F, 
               stroke_color ="black",
               stroke_size =0.4,
               set_name_color ="black",
               set_name_size =2,
               text_color = "black",
               text_size =3)+
    scale_fill_brewer(palette = "Set3")+
    labs(title = i) 
  
  pdf(file.path(output.dir,paste0("Vein_uplogfc0.5_",i,".pdf")),width=4,height=3)
  print(p1)
  dev.off()
  
  write.csv(data,file.path(output.dir,paste0("Vein_uplogfc0.5_",i,".csv")))
  
  all_data <- rbind(all_data, data)
}

write.csv(all_data, file.path(output.dir, "Vein_uplogfc0.5.csv"))


#### 3. 四组比较DOWN plot ----
celltypes <- unique(deg$celltype)
all_data <- data.frame() 
for (i in celltypes) {
  
  deg_down <- deg %>% filter(change == "DOWN") %>% filter(avg_log2FC < -0.5)
  deg_down_celltype <- deg_down %>% filter(celltype %in% i)
  
  A  <-  deg_down_celltype %>% filter(comparison == "M24_vs_P0") %>% pull(gene)
  B  <-  deg_down_celltype %>% filter(comparison == "M24_vs_P14") %>% pull(gene)
  C  <-  deg_down_celltype %>% filter(comparison == "M24_vs_M2") %>% pull(gene)
  D  <-  deg_down_celltype %>% filter(comparison == "M24_vs_M12") %>% pull(gene)
  
  x<-list(`M24 vs P0` = A,
          `M24 vs P14` = B,
          `M24 vs M2` = C,
          `M24 vs M12` = D)
  
  library(ggvenn)
  data <- list_to_data_frame(x)
  data$celltype <- i
  p1 <- ggvenn(data = data, 
               show_percentage = F, 
               stroke_color ="black",
               stroke_size =0.4,
               set_name_color ="black",
               set_name_size =2,
               text_color = "black",
               text_size =3)+
    scale_fill_brewer(palette = "Set3")+
    labs(title = i) 
  
  pdf(file.path(output.dir,paste0("Vein_downlogfc0.5_",i,".pdf")),width=4,height=3)
  print(p1)
  dev.off()
  
  write.csv(data,file.path(output.dir,paste0("Vein_downlogfc0.5_",i,".csv")))
  
  all_data <- rbind(all_data, data)
}

write.csv(all_data, file.path(output.dir, "Vein_downlogfc0.5.csv"))


#### 4. 特定细胞类型 down plot ----
mycolor_down <- c("#E0F2E6","#8BC1AF","#398E88","#005D67")

i <- "MTC"
i <- "TPC"

deg_down <- deg %>% filter(change == "DOWN") %>% filter(avg_log2FC < -0.5)
deg_down_celltype <- deg_down %>% filter(celltype %in% i)

A  <-  deg_down_celltype %>% filter(comparison == "M24_vs_P0") %>% pull(gene)
B  <-  deg_down_celltype %>% filter(comparison == "M24_vs_P14") %>% pull(gene)
C  <-  deg_down_celltype %>% filter(comparison == "M24_vs_M2") %>% pull(gene)
D  <-  deg_down_celltype %>% filter(comparison == "M24_vs_M12") %>% pull(gene)

x<-list(`M24 vs P0` = A,
        `M24 vs P14` = B,
        `M24 vs M2` = C,
        `M24 vs M12` = D)

library(ggvenn)
data <- list_to_data_frame(x)
data$celltype <- i 
p1 <- ggvenn(data = data, 
             show_percentage = F, 
             stroke_color ="black",
             stroke_size =0.4,
             set_name_color ="black",
             set_name_size =3,
             text_color = "black",
             text_size =3)+
  scale_fill_manual(values = mycolor_down)+
  labs(title = i) 

pdf(file.path(output.dir,paste0("Vein_downlogfc0.5_",i,".pdf")),width=2.5,height=2.5)
print(p1)
dev.off()

#### 4.特定细胞类型 up plot ----
mycolor_up <- c("#FADDC3", "#F8B78E", "#F3885B", "#EA4C3B")

i <- "MTC"
i <- "TPC"

deg_up <- deg %>% filter(change == "UP") %>% filter(avg_log2FC > 0.5)
deg_up_celltype <- deg_up %>% filter(celltype %in% i)

A  <-  deg_up_celltype %>% filter(comparison == "M24_vs_P0") %>% pull(gene)
B  <-  deg_up_celltype %>% filter(comparison == "M24_vs_P14") %>% pull(gene)
C  <-  deg_up_celltype %>% filter(comparison == "M24_vs_M2") %>% pull(gene)
D  <-  deg_up_celltype %>% filter(comparison == "M24_vs_M12") %>% pull(gene)

x<-list(`M24 vs P0` = A,
        `M24 vs P14` = B,
        `M24 vs M2` = C,
        `M24 vs M12` = D)

library(ggvenn)
data <- list_to_data_frame(x)
data$celltype <- i 
p1 <- ggvenn(data = data, 
             show_percentage = F, 
             stroke_color ="black",
             stroke_size =0.4,
             set_name_color ="black",
             set_name_size =3,
             text_color = "black",
             text_size =3)+
  scale_fill_manual(values = mycolor_up)+
  labs(title = i) 

pdf(file.path(output.dir,paste0("Vein_uplogfc0.5_",i,".pdf")),width=2.5,height=2.5)
print(p1)
dev.off()

setwd("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/01.deg/03.MTC_singletime/")
output.dir <- paste0("/dellfsqd2/ST_LBI/USER/xingwenlu/Taste/09.Fig2/")  

## function ----
GO <- read.csv("updownlogfc0.5.enrichmentGO.csv")

select_up_way <- c("oxidative phosphorylation",
                   "generation of precursor metabolites and energy",
                   "response to oxidative stress",
                   "regulation of inflammatory response",
                   "regulation of interleukin-8 production",
                   "apoptotic mitochondrial changes"
)
GO_up <- GO %>% filter(change == "UP",
                       Description %in% select_up_way)

select_down_way <- c(
  "macroautophagy",
  "protein stabilization",
  "tissue homeostasis",
  "sensory organ morphogenesis",
  "sensory system development",
  "Wnt signaling pathway"
)
GO_down <- GO %>% filter(change == "DOWN",
                         Description %in% select_down_way)

all_results <- rbind(GO_up, GO_down)

all_results$Description <- factor(all_results$Description,levels=all_results$Description)

mytheme<-theme(
  axis.title.x=element_text(colour='black',size=8),
  axis.title.y=element_text(colour='black',size=8),
  axis.line=element_line(colour='black',size=0.3),
  axis.text=element_text(colour='black',size=8),
  axis.text.y=element_blank(),
  axis.ticks=element_line(colour='black',size=0.4),
  axis.ticks.length.x=unit(0.1,"cm"),
  axis.ticks.length.y=unit(0,"cm"),
  plot.title=element_text(colour='black',size=8,hjust=0.5),
  legend.title=element_text(colour='black',size=8),
  legend.text=element_text(colour='black',size=8)
)

all_results$gene  <- sapply(strsplit(all_results$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))

p <- ggplot(data=all_results,aes(x=-log10(pvalue),y=rev(Description),fill=change))+
  geom_bar(stat="identity",width=0.6,alpha=1)+
  geom_text(size=2.8,aes(x=0.1,y=rev(Description),label=Description),hjust=0)+  
  scale_x_continuous(expand=c(0,0))+  
  scale_fill_manual(values=c("#64A79A","#F3885B"))+
  geom_text(data = all_results,            
            aes(x = 0.1, y = rev(Description), label = gene, color = change),            
            size = 2.8,fontface = 'italic', hjust = 0, vjust = 2.3)+
  scale_color_manual(values=c("#64A79A","#F3885B"))+
  scale_y_discrete(expand = c(0.07,0))+
  theme_classic()+
  labs(x="-Log10(P value)",
       y="",
       title="Pathway enrichment")+
  mytheme+
  NoLegend()

pdf(file.path(output.dir,"MTC_updown_way_barplot.pdf"),width=3.5,height=4.2)
print(p)
dev.off()

GO <- read.csv("updownlogfc0.5.enrichmentGO.csv")
select_up_way <- c(
  "oxidative phosphorylation",
  "respiratory electron transport chain",
  "response to oxidative stress",
  "regulation of leukocyte cell-cell adhesion",
  "regulation of leukocyte proliferation",
  "intrinsic apoptotic signaling pathway")
GO_up <- GO %>% filter(change == "UP",
                       Description %in% select_up_way)
select_down_way <- c(
  "cell-cell signaling by wnt",
  "stem cell differentiation",
  "epithelial cell development",
  "protein stabilization",
  "morphogenesis of a branching structure",
  "BMP signaling pathway"
)
GO_down <- GO %>% filter(change == "DOWN",
                         Description %in% select_down_way)

all_results <- rbind(GO_up, GO_down)

all_results$Description <- factor(all_results$Description,levels=all_results$Description)

mytheme<-theme(
  axis.title.x=element_text(colour='black',size=8),
  axis.title.y=element_text(colour='black',size=8),
  axis.line=element_line(colour='black',size=0.3),
  axis.text=element_text(colour='black',size=8),
  axis.text.y=element_blank(),
  axis.ticks=element_line(colour='black',size=0.4),
  axis.ticks.length.x=unit(0.1,"cm"),
  axis.ticks.length.y=unit(0,"cm"),
  plot.title=element_text(colour='black',size=8,hjust=0.5),
  legend.title=element_text(colour='black',size=8),
  legend.text=element_text(colour='black',size=8)
)

all_results$gene  <- sapply(strsplit(all_results$geneID , "/"), function(x) paste(x[1:5], collapse = "/"))

p <- ggplot(data=all_results,aes(x=-log10(pvalue),y=rev(Description),fill=change))+
  geom_bar(stat="identity",width=0.6,alpha=1)+
  geom_text(size=2.8,aes(x=0.1,y=rev(Description),label=Description),hjust=0)+  
  scale_x_continuous(expand=c(0,0))+  
  scale_fill_manual(values=c("#64A79A","#F3885B"))+
  geom_text(data = all_results,            
            aes(x = 0.1, y = rev(Description), label = gene, color = change),            
            size = 2.8,fontface = 'italic', hjust = 0, vjust = 2.3)+ 
  scale_color_manual(values=c("#64A79A","#F3885B"))+
  scale_y_discrete(expand = c(0.07,0))+
  theme_classic()+
  labs(x="-Log10(P value)",
       y="",
       title="Pathway enrichment")+
  mytheme+
  NoLegend()

pdf(file.path(output.dir,"TPC_updown_way_barplot.pdf"),width=3.5,height=4.2)
print(p)
dev.off()
####Figure3F####
DimPlot(obj,group.by="celltype2",label = F,cols= color,pt.size = 0.5) +
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
####Figure3G####
obj <- readRDS("../recluster_tastecells.rds")
colors <- read.csv("../color.csv")
color <- colors$color
names(color) <- colors$celltype
df <- table(obj$celltype2,obj$timepoint)
df_sub <- df[,c("P0","P14","M2","M12","M24")]

Cellratio <- prop.table(as.matrix(df_sub), margin = 1) #margin 1 = row, 2 = column, default is NULL
  tmp <- melt(Cellratio)
  colnames(tmp)[1:2] <- c("celltype","timepoint")
  ggplot(tmp,aes(x = timepoint, y = celltype, color = celltype,size = value)) +
    geom_point() +theme_bw()+
    scale_color_manual(values = color)+
    theme(plot.background = element_rect(fill="white"),   #背景（不包括绘图区域）
          panel.grid.major = element_blank(),     #网格线空白
          panel.grid.minor = element_blank(),     #网格线空白
          axis.title.x = element_blank(),  # 取消x轴标题  
          axis.title.y = element_blank(),  # 取消y轴标题 
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(hjust = 1,size = 10),
          axis.text.y = element_text(size = 10),
          panel.background = element_rect(fill="white"),
          legend.key = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white"), # 设置图例背景为黑色  
          legend.text = element_text(color = "black"))
  ggsave("all_cellrio_dotplot.pdf",width = 5,height = 5)
####Figure3H####
pdf("umap_trackline.pdf")
plot(obj@reductions$umap@cell.embeddings,pch=16,asp=1,col = color_cell,cex = 0.5)
lines(SlingshotDataSet(crv1), lwd = 3, col = c('black'))
dev.off()

