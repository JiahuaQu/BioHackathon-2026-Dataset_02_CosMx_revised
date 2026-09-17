#!/usr/bin/env python
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from src.qc import add_qc_metrics,recommend_thresholds,apply_qc,threshold_sweep
from src.clustering import preprocess,resolution_sweep,suggest_resolution
from src.plotting import spatial_scatter
from config.markers import MARKER_PANEL
IN=ROOT/"data/processed/GSM9046088_CosMx_raw_FOV46.h5ad"; OUT=ROOT/"results/FOV46"; FIG=OUT/"figures"; TAB=OUT/"tables"; FIG.mkdir(parents=True,exist_ok=True); TAB.mkdir(parents=True,exist_ok=True)
a=sc.read_h5ad(IN); add_qc_metrics(a)
threshold_sweep(a).to_csv(TAB/"qc_threshold_sweep.csv",index=False)
thr=recommend_thresholds(a,nmads=3.0); (TAB/"recommended_qc_thresholds.json").write_text(json.dumps(thr,indent=2)); apply_qc(a,thr,use_per_fov_flags=False)
bio=~a.var.is_control.to_numpy(); adata=a[a.obs.qc_pass,bio].copy(); sc.pp.filter_genes(adata,min_cells=1); adata.layers["counts"]=adata.X.copy(); preprocess(adata)
res=resolution_sweep(adata); res.to_csv(TAB/"leiden_resolution_sweep.csv",index=False); selected=suggest_resolution(res); key=f"leiden_r{selected:g}"; adata.obs["leiden"]=adata.obs[key].copy(); adata.uns["selected_leiden_resolution"]=selected
sc.pl.umap(adata,color=res.key.tolist(),ncols=3,show=False); plt.gcf().savefig(FIG/"umap_resolution_sweep.png",dpi=250,bbox_inches="tight"); plt.close("all")
spatial_scatter(adata,"leiden",FIG/"spatial_leiden.png",size=8); plt.close("all")
sc.tl.rank_genes_groups(adata,"leiden",method="wilcoxon",use_raw=True,pts=True); sc.get.rank_genes_groups_df(adata,group=None).to_csv(TAB/"leiden_markers.csv",index=False)
present={k:[g for g in v if g in adata.raw.var_names] for k,v in MARKER_PANEL.items()}; present={k:v for k,v in present.items() if v}
if present:
    dp=sc.pl.dotplot(adata,var_names=present,groupby="leiden",use_raw=True,standard_scale="var",dendrogram=True,show=False,return_fig=True); dp.savefig(FIG/"marker_dotplot.png",dpi=250); plt.close("all")
# Edit this mapping after reviewing markers/dotplot/spatial localization.
CLUSTER_ANNOTATIONS={}
adata.obs["cell_type_scanpy"]=adata.obs.leiden.astype(str).map(CLUSTER_ANNOTATIONS).fillna("Unassigned").astype("category")
adata.write_h5ad(OUT/"GSM9046088_FOV46_approach_i.h5ad",compression="gzip"); a.write_h5ad(OUT/"GSM9046088_FOV46_raw_with_qc_flags.h5ad",compression="gzip")
print("Selected resolution heuristic:",selected); print("Edit CLUSTER_ANNOTATIONS after marker review, then rerun.")
