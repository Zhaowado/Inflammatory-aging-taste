####Figure6C####
DoHeatmap(obj, features = genes, slot = 'data', raster = F,size = 3,group.colors = color) + scale_fill_gradient2(
  low = rev(c('#d1e5f0', '#67a9cf', '#2166ac')),
  mid = "white",
  high = rev(c('#b2182b', '#ef8a62', '#fddbc7')),
  midpoint = 2,
  guide = "colourbar",
  aesthetics = "fill"
)
####Figure6E####
ggplot(degs,aes(x = celltype,y = gene,color = comline)) +
  geom_point(size = 6)+
 #scale_size_area(max_size = 100) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),       # 移除x轴标题
    axis.title.y = element_blank(),       # 移除y轴标题
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12),# x轴标签倾斜45度，并调整水平对齐
    axis.text.y = element_text(size = 12)
  )
####Figure6F####
library(hdWGCNA)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods[c(1,8,16)] # only select 5 modules
)
