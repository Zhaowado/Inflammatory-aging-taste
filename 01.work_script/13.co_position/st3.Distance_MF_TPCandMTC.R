library(ggplot2)
library(cowplot)
txt <- read.table("../../01.annotation/CV.anno.txt",sep = "\t")
txt$stage <- txt$orig.ident
txt$stage[txt$stage %in% c("B01320A6_7","B01320A6_8")] <- "P14"
txt$stage[txt$stage %in% c("B01320A2_3")] <- "M2"
txt$stage[txt$stage %in% c("B01320A3_4","B01320A3_5")] <- "M12"
txt$stage[txt$stage %in% c("B01320A4_1","B01320A4_3")] <- "M24"

txt <- subset(txt, orig.ident %in% c("B01320A6_7","B01320A2_3","B01320A3_5","B01320A4_3") & celltype2 %in% c("Taste progenitor cells","Mature taste cells","Macrophages"))
txt$celltype2[which(txt$celltype2 == "Mature taste cells")] <- "MTC"
txt$celltype2[which(txt$celltype2 == "Taste progenitor cells")] <- "TPC"
txt$celltype2[which(txt$celltype2 == "Macrophages")] <- "MF"

min_dist <- function(coords_former, coords_latter){
    tmp_fromer <- numeric(nrow(coords_former))
    names(tmp_fromer) <- rownames(coords_former)  # 保留行名
    for (i in 1:nrow(coords_former)){
        tmp_fromer[i] <- min(sqrt((coords_former[i,1] - coords_latter[,1])^2 + (coords_former[i,2] - coords_latter[,2])^2))
    }
    return(tmp_fromer)
}

dis_mtc <- list()
dis_tpc <- list()
dis_mf2mtc <- list()
dis_mf2tpc <- list()
for (time in c("P14","M2","M12","M24")){
    data <- txt[txt$stage == time,]
    mtc_coords <- data[data$celltype2 == "MTC",c("x","y")]
    tpc_coords <- data[data$celltype2 == "TPC",c("x","y")]
    mf_coords <- data[data$celltype2 == "MF",c("x","y")]
    #分别计算每个MTC和TPC到MF的距离和MF分别到他俩的距离
    mtc_mf_dist <- min_dist(mtc_coords, mf_coords)
    dis_mtc[[time]] <- data.frame(distance = mtc_mf_dist, timepoint = time)
    
    tpc_mf_dist <- min_dist(tpc_coords, mf_coords)
    dis_tpc[[time]] <- data.frame(distance = tpc_mf_dist, timepoint = time)
    
    mf2mtc_dis <- min_dist(mf_coords, mtc_coords)
    dis_mf2mtc[[time]] <- data.frame(distance = mf2mtc_dis, timepoint = time)
    
    mf2tpc_dis <- min_dist(mf_coords, tpc_coords)
    dis_mf2tpc[[time]] <- data.frame(distance = mf2tpc_dis, timepoint = time)
}


plot_result <- function(distance_df){
    distance_df <- do.call(rbind, distance_df)
    distance_df$timepoint <- factor(distance_df$timepoint, levels = c("P14","M2","M12","M24"))
    p <- ggplot(distance_df, aes(x = timepoint, y = log(distance), fill = timepoint)) +
            geom_boxplot() +
            labs(
                x = "timepoint",
                y = "Log(Distance)"
            ) +
            scale_fill_manual(values = c(P0='#66cdaa',P14='#1f78b4',M2="#E7B800",
                                 M12='#ff69b4',M24="#E66101")) +
            theme_bw()
}

ggsave(plot_result(dis_mtc), file = "mtc_mf_distance.pdf", width = 4, height = 3)
ggsave(plot_result(dis_tpc), file = "tpc_mf_distance.pdf", width = 4, height = 3)
