####FigureS6A####
PlotKMEs(seurat_obj, ncol=5)

####FigureS6B####
PlotKMEs(seurat_obj, ncol=5)

plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

