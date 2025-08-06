library(nichenetr)
library(Seurat)
library(tidyverse)
library(ggplot2)

### Get the parameters
parser = argparse::ArgumentParser(description = 'Script to infer ligand and target network')
parser$add_argument('-s', dest = 'sender', help = 'the sender cell type in ccc')
parser$add_argument('-r', dest = 'receiver', help = 'the receiver cell type in ccc')
parser$add_argument('-l', dest = 'ligand', help = 'the interest ligand in ccc')
opts = parser$parse_args()

sender_celltypes <- opts$sender
receiver <- opts$receiver
ligand <- opts$ligand
############load data#############
obj <- readRDS("../../../../01.annotation/taste.integrated.rds")
meta <- read.table("../../../../01.annotation/sub_anno_all_9.9.txt",sep = "\t")
obj@meta.data <- meta

obj <- subset(obj,celltype_sub %in% c(sender_celltypes,receiver))
obj$timepoint_select <- "Others"
obj$timepoint_select[which(obj$timepoint == "M24")] <- "M24"

############load lr database#############
ligand_target_matrix <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/Tools/nichenet_lrdb/ligand_target_matrix_nsga2r_final_mouse.rds")
lr_network <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/Tools/nichenet_lrdb/lr_network_mouse_21122021.rds")
weighted_networks <- readRDS("/dellfsqd2/ST_LBI/USER/zhaowandong/Tools/nichenet_lrdb/weighted_networks_nsga2r_final_mouse.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

###########perform nichenet analysis#########
Idents(obj) <- "celltype_sub"
expressed_genes_receiver = get_expressed_genes(receiver, obj, pct = 0.05,assay_oi = "RNA")  

expressed_receptors <- intersect(unique(lr_network$to), expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

####Define sender cells and select potential lr again

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, obj, 0.05,assay_oi = "RNA")
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)

####Define interest gene set, (genes of receiver cells which are influenced by ccc envent)
condition_oi <-  "M24"
condition_reference <- "Others"

seurat_obj_receiver <- subset(obj, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "timepoint_select",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% pull(gene)  
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

print(paste0("genset_oi number is ",length(geneset_oi)))
####Define background genes 
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

###########check ligand in cellchat###############
best_upstream_ligands <- c(ligand)

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 500) %>%
  bind_rows() %>% drop_na()

#This step is for better plot, change the target score to 0 if min cutoff  
active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25) 



score <- as.data.frame(active_ligand_target_links)
#score$log <- -log(score$Grn)
score[,"log"] <- -log(score[,ligand])
score$gene <- rownames(score)
write.csv(score,paste(ligand,sender_celltypes,receiver,"lt_infer.csv",sep = "_"))
ggplot(score,aes(gene,ligand)) +geom_bar(stat = "identity", fill = "skyblue") +
  theme_minimal() 
ggsave(paste(ligand,sender_celltypes,receiver,"target_infer.pdf",sep = "_"),width = 3,height = 4)
