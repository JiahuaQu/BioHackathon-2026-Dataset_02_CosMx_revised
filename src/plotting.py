from __future__ import annotations
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def spatial_scatter(adata, color, path=None, size=3, title=None):
    if "spatial" not in adata.obsm: raise KeyError("adata.obsm['spatial'] missing")
    xy=adata.obsm["spatial"]; v=adata.obs[color]
    fig,ax=plt.subplots(figsize=(10,8))
    if pd.api.types.is_numeric_dtype(v):
        p=ax.scatter(xy[:,0],xy[:,1],c=v,s=size,linewidths=0,rasterized=True); fig.colorbar(p,ax=ax,label=color)
    else:
        cat=v.astype("category"); cmap=plt.get_cmap("tab20",max(1,len(cat.cat.categories)))
        ax.scatter(xy[:,0],xy[:,1],c=cat.cat.codes,cmap=cmap,s=size,linewidths=0,rasterized=True)
        if len(cat.cat.categories)<=30:
            handles=[plt.Line2D([0],[0],marker="o",linestyle="",color=cmap(i),label=str(c)) for i,c in enumerate(cat.cat.categories)]
            ax.legend(handles=handles,bbox_to_anchor=(1.02,1),loc="upper left",frameon=False)
    ax.set_aspect("equal"); ax.set_xlabel("Global X (px)"); ax.set_ylabel("Global Y (px)"); ax.set_title(title or color); ax.grid(False); fig.tight_layout()
    if path: fig.savefig(path,dpi=300,bbox_inches="tight")
    return fig

def qc_by_fov_plots(adata, out):
    cols=[c for c in ["total_counts","n_genes_by_counts","control_fraction","qc_cell_area"] if c in adata.obs]
    for c in cols:
        fig,ax=plt.subplots(figsize=(max(12,adata.obs["fov"].nunique()*.25),5)); sns.boxplot(data=adata.obs,x="fov",y=c,showfliers=False,ax=ax); ax.tick_params(axis="x",rotation=90); fig.tight_layout(); fig.savefig(out/f"qc_by_fov_{c}.png",dpi=250,bbox_inches="tight"); plt.close(fig)
