#!/usr/bin/env python
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
import scanpy as sc
from src.plotting import spatial_scatter
from config.markers import MARKER_PANEL
IN=ROOT/'results/FOV46/GSM9046088_FOV46_leiden_sweep.h5ad'; OUT=ROOT/'results/FOV46'; FIG=OUT/'figures/approach_i'; TAB=OUT/'tables'; CFG=ROOT/'config/fov46_annotation.json'
FIG.mkdir(parents=True,exist_ok=True); TAB.mkdir(parents=True,exist_ok=True)
a=sc.read_h5ad(IN); cfg=json.loads(CFG.read_text()); r=float(cfg['resolution']); key=f'leiden_r{r:g}'
if key not in a.obs: raise KeyError(f'{key} not found. Available: {[c for c in a.obs if c.startswith("leiden_r")]}')
a.obs['leiden']=a.obs[key].copy(); a.uns['selected_leiden_resolution']=r
sc.tl.rank_genes_groups(a,'leiden',method='wilcoxon',use_raw=True,pts=True)
sc.get.rank_genes_groups_df(a,group=None).to_csv(TAB/f'leiden_markers_r{r:g}.csv',index=False)
present={ct:[g for g in genes if g in a.raw.var_names] for ct,genes in MARKER_PANEL.items()}; present={ct:g for ct,g in present.items() if g}
if present:
    dp=sc.pl.dotplot(a,var_names=present,groupby='leiden',use_raw=True,standard_scale='var',dendrogram=True,show=False,return_fig=True); dp.savefig(FIG/f'marker_dotplot_selected_r{r:g}.png',dpi=250); plt.close('all')
spatial_scatter(a,'leiden',FIG/f'spatial_selected_r{r:g}.png',size=8); plt.close('all')
ann={str(k):v for k,v in cfg.get('cluster_annotations',{}).items()}
clusters=set(a.obs['leiden'].astype(str).unique()); missing=sorted(clusters-set(ann))
if missing: raise ValueError(f'Missing manual annotations for clusters {missing}. Review dotplot/markers, edit {CFG}, rerun.')
a.obs['cell_type_scanpy']=a.obs['leiden'].astype(str).map(ann).astype('category')
spatial_scatter(a,'cell_type_scanpy',FIG/'spatial_cell_type_scanpy.png',size=8); plt.close('all')
a.write_h5ad(OUT/'GSM9046088_FOV46_approach_i.h5ad',compression='gzip')
print('Approach I complete. Resolution:',r); print(a.obs.cell_type_scanpy.value_counts())
