#!/usr/bin/env python
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from src.plotting import spatial_scatter
IN=ROOT/"results/FOV46/GSM9046088_FOV46_approach_i_ii.h5ad"; OUT=ROOT/"results/FOV46"; FIG=OUT/"figures"; TAB=OUT/"tables"
a=sc.read_h5ad(IN)
ct=pd.crosstab(a.obs.cell_type_scanpy,a.obs.cell_type_rule); ct.to_csv(TAB/"approach_i_vs_ii_counts.csv")
ctn=pd.crosstab(a.obs.cell_type_scanpy,a.obs.cell_type_rule,normalize="index"); ctn.to_csv(TAB/"approach_i_vs_ii_row_fraction.csv")
a.obs["typing_agreement"]=(a.obs.cell_type_scanpy.astype(str)==a.obs.cell_type_rule.astype(str))
pd.Series({"agreement_fraction":a.obs.typing_agreement.mean(),"rule_unassigned_fraction":a.obs.cell_type_rule.astype(str).eq("Unassigned").mean(),"rule_ambiguous_fraction":a.obs.cell_type_rule.astype(str).eq("Ambiguous").mean()}).to_csv(TAB/"typing_summary.csv")
fig,ax=plt.subplots(figsize=(12,7)); sns.heatmap(ctn,annot=True,fmt=".2f",cmap="viridis",ax=ax); fig.tight_layout(); fig.savefig(FIG/"approach_i_vs_ii_heatmap.png",dpi=250); plt.close(fig)
for col in ["cell_type_scanpy","cell_type_rule","typing_agreement"]:
    spatial_scatter(a,col,FIG/f"spatial_{col}.png",size=8); plt.close("all")
if "final_CT" in a.obs:
    pd.crosstab(a.obs.final_CT,a.obs.cell_type_scanpy).to_csv(TAB/"reference_vs_approach_i.csv")
    pd.crosstab(a.obs.final_CT,a.obs.cell_type_rule).to_csv(TAB/"reference_vs_approach_ii.csv")
a.write_h5ad(OUT/"GSM9046088_FOV46_typing_compared.h5ad",compression="gzip")
