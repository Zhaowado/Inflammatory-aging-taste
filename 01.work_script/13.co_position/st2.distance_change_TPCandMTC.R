library(ggplot2)
library(ggpubr)

txt <- read.table("../../01.annotation/CV.anno.txt",sep = "\t")
txt$stage <- txt$orig.ident
txt$stage[txt$stage %in% c("B01320A6_7","B01320A6_8")] <- "D14"
txt$stage[txt$stage %in% c("B01320A2_3")] <- "M2"
txt$stage[txt$stage %in% c("B01320A3_4","B01320A3_5")] <- "M12"
txt$stage[txt$stage %in% c("B01320A4_1","B01320A4_3")] <- "M24"

txt <- subset(txt, orig.ident %in% c("B01320A6_7","B01320A2_3","B01320A3_5","B01320A4_3") & celltype2 %in% c("Taste progenitor cells","Mature taste cells"))
txt$celltype2[which(txt$celltype2 == "Mature taste cells")] <- "MTC"
txt$celltype2[which(txt$celltype2 == "Taste progenitor cells")] <- "TPC"

distance_results <- list()
for (time in c("M2","M24")) {
    data <- txt[txt$stage == time,]
    mtc_coords <- data[data$celltype2 == "MTC",c("x","y")]
    tpc_coords <- data[data$celltype2 == "TPC",c("x","y")]
    # 计算每个 MTC 到最近 TPC 的距离
    mtc_distances <- numeric(nrow(mtc_coords))
    for (i in 1:nrow(mtc_coords)) {
    mtc_point <- mtc_coords[i, ]
    dist_vec <- sqrt((tpc_coords$x - mtc_point$x)^2 + (tpc_coords$y - mtc_point$y)^2)
    mtc_distances[i] <- min(dist_vec)
  }
  # 存储结果
  distance_results[[time]] <- data.frame(
    distance = mtc_distances,
    timepoint = time
  )
}
distance_df <- do.call(rbind, distance_results)
distance_df$timepoint <- factor(distance_df$timepoint, levels = c("M2","M24"))

cutoff <- quantile(distance_df$distance, 0.95)
distance_df$distance_truncated <- pmin(distance_df$distance, cutoff)  # 截断数据
ggplot(distance_df, aes(x = timepoint, y = distance_truncated, fill = timepoint)) +
  geom_boxplot() +
  labs(
    x = "timepoint",
    y = "Distance to nearest TPC"
  ) +
  theme_bw()+
  scale_fill_manual(values = c("M2" = "#5983af", "M24" = "#eba372")) +
  theme(
  axis.title = element_text(size = 14),
  axis.text = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 14)) + 
  stat_compare_means(comparisons = list(c("M2", "M24")), method = "t.test",label = "p.signif")
ggsave("Distance_MTC2TPC.pdf")
