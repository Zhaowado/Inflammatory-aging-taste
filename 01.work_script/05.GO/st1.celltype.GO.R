library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)

parser = argparse::ArgumentParser(description = "GO")

parser$add_argument('-i', dest = 'input', help = 'input DEG file')

opts = parser$parse_args()

markers <- read.table(opts$input,sep = "\t",header=T)


for (i in unique(markers$cluster)) {
  col<-c("SYMBOL","ENTREZID")
  marker <- markers[markers$cluster == i,] %>% .[.$p_val_adj < 0.05 & .$avg_log2FC >0.25,]
  da1 <- AnnotationDbi::select(org.Mm.eg.db,columns=col,keytype="SYMBOL",keys=as.character(marker$gene))
  ego<-enrichGO(OrgDb = "org.Mm.eg.db",gene=da1$ENTREZID,ont="BP",pvalueCutoff=0.05,readable=TRUE)
  write.csv(ego,paste0("./01.GO_result/",i,"_GO.csv"),row.names = F)
}
go_all <- NULL
GO2 <- NULL
for (i in unique(markers$cluster)) {
  GO1 <- read.csv(paste0("./01.GO_result/",i,"_GO.csv"))
  GO1$celltype <- i
  go_all <- rbind(GO1,go_all)
  GO2 <- rbind(GO1[1:2,],GO2)
}
write.csv(go_all,"GO.csv")
