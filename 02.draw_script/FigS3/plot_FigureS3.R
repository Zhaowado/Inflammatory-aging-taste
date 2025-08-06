####FigureS3A,B####
pdf("heatmap.pdf",width = 3,height = 5)
p <-plot_pseudotime_heatmap(my_cds[unique(genes),],
                            num_clusters = 6,
                            cores = 1,
                            show_rownames = F,return_heatmap=F)
print(p)
dev.off()

####FigureS3C####
FeaturePlot(obj,features ="Lgr5",pt.size=0.3, cols = c(lighten("#ea5c6f",0.6),darken("#ea5c6f", 0.4)),raster=FALSE,slot = "data")
ggsave("Lgr5_FeaturePlot.pdf",width = 4,height = 3)

FeaturePlot(obj,features ="Mki67",pt.size=0.3, cols = c(lighten("#f7905a",0.6),darken("#f7905a",0.4)),raster=FALSE,slot = "data")
ggsave("Mki67_FeaturePlot.pdf",width = 4,height = 3)

FeaturePlot(obj,features ="Krt5",pt.size=0.3, cols = c(lighten("#e187cb",0.6),darken("#e187cb",0.4)),raster=FALSE,slot = "data")
ggsave("Krt5_FeaturePlot.pdf",width = 4,height = 3)

FeaturePlot(obj,features ="Shh",pt.size=0.3, cols = c(lighten("#e2b159",0.6),darken("#e2b159",0.4)),raster=FALSE,slot = "data")
ggsave("Shh_FeaturePlot.pdf",width = 4,height = 3)

FeaturePlot(obj,features ="Nupr1",pt.size=0.3, cols = c(lighten("#fb948d",0.6),darken("#fb948d",0.4)),raster=FALSE,slot = "data")
ggsave("Nupr1_FeaturePlot.pdf",width = 4,height = 3)

FeaturePlot(obj,features ="Ovol3",pt.size=0.3, cols = c(lighten("#ebed6f",0.6),darken("#ebed6f",0.4)),raster=FALSE,slot = "data")
ggsave("Ovol3_FeaturePlot.pdf",width = 4,height = 3)

FeaturePlot(obj,features ="Pkd1l3",pt.size=0.3, cols = c(lighten("#b2db87",0.6),darken("#b2db87",0.4)),raster=FALSE,slot = "data")
ggsave("Pkd1l3_FeaturePlot.pdf",width = 4,height = 3)

####FigureS3D####
p<-DotPlot(obj, features = unique(markers$marker)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33'))

p$data$id <- factor(p$data$id, levels = unique(markers$celltype))
p$data$features.plot <- factor(p$data$features.plot, levels = unique(markers$marker))

pdf("marker_dotplot.pdf", width = 4, height=4)
print(p)
dev.off()
