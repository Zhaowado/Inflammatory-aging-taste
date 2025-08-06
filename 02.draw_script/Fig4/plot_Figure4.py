####Figure4D####
celltypes = ["TPC","MTC","TC","MF"]
mycolor = ListedColormap(["#815e99","#D1392B","#0F253A","#00BA38"])
sq.pl.co_occurrence(adata, cluster_key="celltype", clusters=celltypes, palette = mycolor, save=outpre+".co_occurrence.nsplits"+".pdf",figsize=(6*nums,5))
