library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(patchwork)

genes <- read.csv("TPC.heatmap.clustergenes.csv")
col<-c("SYMBOL","ENTREZID")

df <- c()
for (cluster in unique(genes$Gene_clusters)) {
  gene_select <- subset(genes,Gene_clusters == cluster)
  gene_select <- gene_select$X
  data <- AnnotationDbi::select(org.Mm.eg.db,columns=col,keytype="SYMBOL",keys=as.character(unique(gene_select)))
  go <- enrichGO(OrgDb = "org.Mm.eg.db",gene=data$ENTREZID,ont="BP",pvalueCutoff=0.05,readable=TRUE)
  df_go <- data.frame(ID = go$ID,Description = go$Description,p.adjust = go$p.adjust,cluster = cluster,count = go$Count,GeneRatio = go$GeneRatio,geneID = go$geneID)
  df[[cluster]] <- df_go
  }

dfs <- df[1:length(unique(genes$Gene_clusters))]
df_all <- do.call(rbind, dfs)
write.csv(df_all,"TPC.timepoint.hp_GO.csv")
