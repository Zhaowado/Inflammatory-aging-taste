####Figure5A####
colors <- read.csv("colors.csv")
color <- colors$color

color <- colors$color
names(color) <- colors$celltype

DimPlot(obj,group.by="celltype2",cols = color,label = F) +
  tidydr::theme_dr(xlength = 0.2,
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        legend.text = element_text(size = 14),
        legend.title = element_blank())

####Figure5B####
FeaturePlot(obj,features = c("Bgn","Cd207","Cd209a","Cx3cr1","Folr2","Mki67","S100a9"))

####Figure5C####
library(scRNAtoolVis)
p1 <- jjVolcano(diffData = markers,
                pvalue.cutoff = 0.01,legend.position = c(1,1.2),
                tile.col = color)

####Figure5D####
p <- pheatmap(as.matrix(aver_expr),scale = "row",
               fontsize=14,fontsize_row=14,fontsize_col=14,
               cluster_rows = F,cluster_cols = FALSE, 
               show_colnames = F,show_rownames =T, 
               annotation_col = anno_col, #列注释
               annotation_colors = anno_color, #注释配色
               gaps_col = seq(5,5*(length(unique(obj$celltype2))),by = 5), #在第几列后分开
               gaps_row = cumsum(sapply(gene.list,length)), #在第几行后分开
               cutree_rows = length(gene.list), #行聚类分成几类
               color = colorRampPalette(c("#5E3C99","white","firebrick3"))(100),
               border_color = NA) #描边颜色

####Figure5E####
plot_cg <- function(celltype,genes){
    cell_index <- WhichCells(obj, expression = celltype2 == as.character(celltype))
    if (length(cell_index) == 0) {
        stop(paste("Cell type", celltype, "not found in the object."))
    }

    expr_data <- data.frame()
    for (gene in genes) {
        gene_expr <- FetchData(obj, vars = gene, cells = cell_index,slot = "scale.data")
        gene_expr$timepoint <- obj$timepoint[cell_index]
        gene_expr$gene <- gene
        colnames(gene_expr)[1] <- "expression"
        expr_data <- rbind(expr_data, gene_expr)
  }
  expr_data <- expr_data %>%
  mutate(
    timepoint = factor(timepoint, levels = c("P0","P14","M2","M12","M24")) 
  )
  
  # 每个基因在每个时间点的平均表达（用于绘制独立折线）
  gene_summary <- expr_data %>%
    group_by(gene, timepoint) %>%
    summarise(
      mean_expr = mean(expression),
      .groups = "drop"
    )
  colors <- read.csv("../../colors.csv")
  color <- colors$color
  names(color) <- colors$celltype
  
  p <- ggplot(gene_summary, aes(x = timepoint, y = mean_expr, color = gene, group = gene)) +
  geom_smooth(
    aes(group = 1, color = NULL, fill = "Overall Trend"),  # 强制所有基因合并为一组
    method = "loess",
    se = TRUE,                 # 显示置信区间（阴影）
    size = 1,                  # 中心线粗细
    color = "black",           # 中心线颜色
    linetype = "dashed",       # 中心线为虚线
    fill = color[celltype],    # 阴影颜色
    alpha = 0.2,               # 阴影透明度
    span = 0.8                 # 控制平滑度（0-1，值越大越平滑）
  ) +
  # 美化图形
  scale_x_discrete(limits = c("P0","P14","M2","M12","M24")) +  # 确保 x 轴按时间顺序排列
  #scale_color_manual(values = c( "#C49A6C","#BFA698","#A65A47","#8A9A65","#9A7B8C"))+ #基因颜色
  theme_bw() +
  theme(
      plot.title = element_text(size = 30,hjust = 0.5),  # 增大图标题字号
      axis.title = element_text(size = 28,hjust = 0.5),  # 增大坐标轴标签字号
      axis.text = element_text(size = 26),   # 增大坐标轴刻度文字字号
      legend.title = element_text(size = 28),# 增大图例标题字号
      legend.text = element_text(size = 26)  # 增大图例文本字号
  ) +
  labs(
    title = celltype,
    y = "Mean Expression",
    x = "",
    color = "Gene",
    ) 
   # + coord_cartesian(ylim = c(-1, 8)) 

  return(p)
}
####Figure5F####
p1 <- ggboxplot(obj@meta.data, x = "celltype2",y = "score",width = 0.8,
                fill = "celltype2",palette = "npg",
                xlab = F,
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.2, #误差条大小
                size=0.5, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right") +
                theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c(`Cd209a+ MF` = '#fc9a9a', `Cx3cr1+ MF` = '#b37557', `Folr2+ MF` = '#b49d99',
                               `Mki67+ MF` = '#96873b', `S100a9+ MF` = '#b24646',`Bgn+ MF` = '#fe9e37', `Cd207+ MF` = '#ed4437'))
####Figure5G####
netVisual_individual(cellchat, signaling="MIF", pairLR.use= "MIF_CD74_CD44" , layout = "chord" )
####Figure5H####
interaction_names <- c("Thbs1_Sdc4","Thbs1_Sdc1","Thbs1_Cd47","Grn_Sort1")
interaction_names <- toupper(interaction_names)
pairLR <- data.frame(interaction_name = interaction_names)
pdf(paste0(time,"cellchat_ccc_M2T.pdf"),width = 8,height = 4)
print(netVisual_bubble(cellchat, sources.use = c("Folr2+ MF","Cx3cr1+ MF","S100a9+ MF"),
                       targets.use = c("Type II MTC","Type III MTC","Lgr5+ TPC","Mki67+ TPC","Shh+ TPC","Krt5+ TPC"),
                       pairLR.use = pairLR,remove.isolate = T,font.size = 15) + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1)))
dev.off()
####Figure5I####
VlnPlot(robj,features = "Grn",group.by = "celltype_sub",split.by = "timepoint",cols = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                            M12='#ff69b4',M24="#E66101"),split.plot = TRUE)
VlnPlot(robj,features = "Thbs1",group.by = "celltype_sub",split.by = "timepoint",cols = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                            M12='#ff69b4',M24="#E66101"),split.plot = TRUE)
VlnPlot(robj,features = "Sort1",group.by = "celltype_sub",split.by = "timepoint",cols = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                            M12='#ff69b4',M24="#E66101"),split.plot = TRUE)
VlnPlot(robj,features = "Sdc1",group.by = "celltype_sub",split.by = "timepoint",cols = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                            M12='#ff69b4',M24="#E66101"),split.plot = TRUE)
