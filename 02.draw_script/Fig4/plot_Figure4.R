####Figure4A,B####
obj$timepoint <- factor(obj$timepoint,levels = c("M2","M24"))
my_comparisons <- list(c("M2","M24"))

p1 <- ggboxplot(obj@meta.data, x = "timepoint",y = "Score1",width = 0.9,
                fill = "timepoint",palette = "npg",
                xlab = F,
                ylab = F,
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.4, #误差条大小
                size=0.5, #箱型图边线的粗细
                outlier.shape=NA, #不显示outlier
                legend = "right") + facet_grid(.~figname,scales = "free",space="free_x") +
                theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                               M12='#ff69b4',M24="#E66101")) +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",label = "p.signif",label.y = c(max(obj$Score1),max(obj$Score1) +0.01,max(obj$Score1) +0.02,max(obj$Score1)+0.03))
####Figure4C####
anno_color <- list(Timepoint = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                                 M12='#ff69b4',M24="#E66101"))

pdf("influmm.genes_all_exp.pdf",width = 6, height = 8)
print(pheatmap(hp,cluster_rows = F,cluster_cols = F,scale = "row",
               treeheight_row = 0, #聚类树高度
               fontsize=10,fontsize_row = 10, fontsize_col = 10,
               cellwidth = 16,cellheight = 16,  #单元格长宽
               annotation_colors = anno_color,
               annotation_col = anno_col,
               annotation_row = gene_list,
               colorRampPalette(c("#5E3C99","white","firebrick3"))(100)))
dev.off()

####Figure4E####
cluster_Palette <- c("#00BA38","#D1392B","#0F253A","#815e99","grey90")
names(cluster_Palette)[1] <- "TPC"
names(cluster_Palette)[2] <- "MTC"
names(cluster_Palette)[3] <- "TC"
names(cluster_Palette)[4] <- "MF"
names(cluster_Palette)[5] <- "Others"

pall <- list()
for (id in sampleid) {
  da1 <- txt[txt$orig.ident == id,]
  p <- ggplot(da1,aes(x,y,color = feature))+
    geom_point(size=0.001,position = position_jitter(width = 0.001, height = 0.001))+
    scale_color_manual(values = cluster_Palette)+
    theme_bw()+
    theme(legend.title=element_blank(),
           legend.key.size = unit(20, "pt"), #图例大小
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.title.x=element_text(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),axis.title.y=element_text()) + coord_fixed(ratio = 1) + labs(x=unique(da1$stage), y = "")
  pall[[id]] <- as.ggplot(p)
  # + theme(legend.position.inside=c(0.8,1.03), legend.direction="horizontal",legend.key.size = unit(100, "pt"))
}
####Figure4F####
plot_result <- function(distance_df){
    distance_df <- do.call(rbind, distance_df)
#    distance_df$timepoint <- factor(distance_df$timepoint, levels = c("P14","M2","M12","M24"))
    distance_df$timepoint <- factor(distance_df$timepoint, levels = c("M2","M24"))
    p <- ggplot(distance_df, aes(x = timepoint, y = distance, fill = timepoint)) +
            geom_boxplot() +
            labs(
                x = "timepoint",
                y = "Distance (µm)"
            ) +
            scale_fill_manual(values = c(P14='#1f78b4',M2="#E7B800",
                                 M12='#ff69b4',M24="#E66101")) +
            stat_compare_means(comparisons = list(c("M2", "M24")), method = "t.test",label = "p.signif") +
            theme_bw()
}
####Figure4I####
netVisual_diffInteraction(cellchat, comparison = c(1, 2), weight.scale = T,
                          color.use = c("#00BA38","#D1392B","#0F253A","#815e99"))
####Figure4J####
compareInteractions(cellchat, show.legend = F, group = c(1,2),
                           color.use = c("#E7B800","#E66101"))
####Figure4K####
pdf("Inflammation.compare.pdf")
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c('#DEB216','#DA600F'),
        signaling = c("MHC-II","CD45","CD86","TNF","CCL","CXCL","ICAM"))
dev.off()

pdf("Growth_factors.compare.pdf")
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c('#DEB216','#DA600F'),
        signaling = c("EGF","NOTCH","IGF","BMP","WNT","TGFb","KIT"))
dev.off()

pdf("ECM.compare.pdf")
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c('#DEB216','#DA600F'),
        signaling = c("THBS","TENASCIN","LAMININ","DESMOSOME","NCAM","CADM","GALECTIN"))
dev.off()

####Figure4L,M####
pdf(file.path(output.dir, paste0("heatmap_incoming-", names(object.list)[1], "_select.pdf")),height=8,width=6)
netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], height = 6, width = 3.5,
                                  color.heatmap = "GnBu",color.use = c("#00BA38","#D1392B","#0F253A","#815e99"))
dev.off()

pdf(file.path(output.dir, paste0("heatmap_incoming-", names(object.list)[2], "_select.pdf")),height=8,width=6)
netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], height = 6, width = 3.5,
                                  color.heatmap = "GnBu",color.use = c("#00BA38","#D1392B","#0F253A","#815e99"))
dev.off()  
