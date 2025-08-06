####FigureS2B-D####
# 定义绘图函数：输入Seurat对象和基因列表，输出小提琴图
plot_gene_vln <- function(obj, gene_list) {
  # 检查输入的基因是否在Seurat对象中存在
  valid_genes <- intersect(gene_list, rownames(obj))
  if (length(valid_genes) == 0) {
    stop("基因列表中的基因在Seurat对象中均不存在，请检查基因名！")
  } else if (length(valid_genes) < length(gene_list)) {
    warning(paste("以下基因未在Seurat对象中找到，已自动过滤：", 
                  paste(setdiff(gene_list, valid_genes), collapse = ", ")))
  }
  
  p <- VlnPlot(
    object = obj,
    features = valid_genes,  
    group.by = "orig.ident",
    assay = "RNA",  
    pt.size = 0,  
    combine = TRUE,  
    same.y.lims = TRUE  
  ) + 
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)  
  )
  
  return(p)
}

####FigureS2E,F####
for (tis in c("E","K","M")){
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
}
