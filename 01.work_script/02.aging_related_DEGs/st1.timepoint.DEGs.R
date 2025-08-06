library(Seurat)
library(dplyr)
library(tibble)

obj <- readRDS('../../../taste.integrated.rds')
meta <- read.table("../../../sub_anno_all_9.9.txt",sep = "\t",header = T)
obj@meta.data <- meta

DefaultAssay(obj) <- "RNA"

timepoints <- c("P0","P14","M2","M12","M24")
markers.list <- list()
for (i in unique(obj$figname)) {
  if (sum(obj$figname == i) >20){
      rrobj <- subset(obj, figname == i)
      for (seq in 1:4) {
        st1 <- timepoints[seq]
        st2 <- timepoints[5]
        robj <- subset(rrobj,timepoint %in% c(st1,st2))
        comline <- paste(st2,st1,sep = "-")
        Idents(robj) <- "timepoint"   ##rely on what to find deg
        count_st1 <- sum(robj$timepoint == st1)
        count_st2 <- sum(robj$timepoint == st2)
        if (count_st1 > 10 & count_st2 >10){

            subset.markers <- FindMarkers(robj,logfc.threshold = 0.25,ident.1 = st2,ident.2 = st1,slot = "data") 
            subset.markers=rownames_to_column(subset.markers, var = "gene")
            subset.markers$celltype <- i
            subset.markers$comline <- comline
            count <- paste(i,comline,sep = "_")
            print(paste0('sucessfully compare ',count))
            markers.list[[count]] <- subset.markers
   }
  }
 }
}

markers <- do.call(rbind,markers.list)
markers <- markers[markers$p_val_adj < 0.05,]
write.table(markers, paste0('All.timepoint.Markers.xls'), sep = '\t', quote = FALSE, row.names = F)
