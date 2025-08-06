import sys,os
import squidpy as sq
import scanpy as sc
import numpy as np
from matplotlib.colors import ListedColormap

adata= sc.read_h5ad(sys.argv[1])
outpre = os.path.basename(sys.argv[1]).split(".h5ad")[0]
type_key = "figname"

#celltypes = ["TPC","MTC","AC","BC","EC","EBC","FLF","MF","MAC","MFB","PC","PFC","SC","SDC","SGC","SEC","TC","TMC"]
#mycolor = ListedColormap(["#00BA38","#D1392B","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0","#F0F0F0"])
celltypes = ["TPC","MTC","TC","MF"]
#mycolor = ListedColormap(["#00BA38","#D1392B","#0F253A","#815e99"])
mycolor = ListedColormap(["#815e99","#D1392B","#0F253A","#00BA38"])
adata = adata[np.isin(adata.obs["figname"], celltypes)]

adata.obsm['X_spatial'] = adata.obs[["x","y"]].values
adata.obsm['spatial'] = adata.obs[["x","y"]].values
#celltypes = adata.obs[type_key].tolist()
adata.obs["celltype"]=adata.obs[type_key].astype('category')

nums = len(celltypes)

#sc.pl.scatter(adata, color="celltype", basis="spatial", save="."+outpre+".pdf")

sq.gr.co_occurrence(adata, cluster_key="celltype",spatial_key="spatial", interval=200,n_splits = 20)
sq.pl.co_occurrence(adata, cluster_key="celltype", clusters=celltypes, palette = mycolor, save=outpre+".co_occurrence.nsplits"+".pdf",figsize=(6*nums,5))

sq.gr.spatial_neighbors(adata, coord_type="generic")
sq.gr.nhood_enrichment(adata, cluster_key="celltype")
sq.pl.nhood_enrichment(adata, cluster_key="celltype", method="ward",save=outpre+".nhood_enrichment"+".pdf",figsize=(6*nums,5))


