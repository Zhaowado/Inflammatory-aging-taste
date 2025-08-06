####FigureS5A####
p <- netVisual_bubble(
  cellchat,
  sources.use = NULL, 
  targets.use = NULL, 
  remove.isolate = FALSE,
  angle.x = 45,
  font.size = 8, 
  title.name = "All Ligand-Receptor Pairs Communication" )

####FigureS5B####
mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))
pdf("M24_split.circle.pdf")
par(mfrow = c(1,2), xpd=TRUE)
for (i in c(3,4)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i],
                   
label.edge= T)
}
