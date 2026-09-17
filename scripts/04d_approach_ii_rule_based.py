#!/usr/bin/env python
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from config.markers import MARKER_PANEL
from src.typing import detection_rate
IN=ROOT/'results/FOV46/GSM9046088_FOV46_approach_i.h5ad'; OUT=ROOT/'results/FOV46'; TAB=OUT/'tables'; TAB.mkdir(parents=True,exist_ok=True)
a=sc.read_h5ad(IN); layer='counts'; Xall=a.layers[layer]
MIN_POSITIVE_GENES=3
rates=[]; availability=[]; positives={}
for ct,genes in MARKER_PANEL.items():
    r=detection_rate(a,genes,layer=layer); rates.extend({'cell_type':ct,'gene':g,'detection_rate':v} for g,v in r.items())
    present=[g for g in genes if g in a.var_names]; missing=[g for g in genes if g not in a.var_names]
    availability.append({'cell_type':ct,'n_requested':len(genes),'n_present':len(present),'min_positive_genes':MIN_POSITIVE_GENES,'available':len(present)>=MIN_POSITIVE_GENES,'present':','.join(present),'missing':','.join(missing)})
    if len(present)<MIN_POSITIVE_GENES:
        positives[ct]=np.zeros(a.n_obs,dtype=bool); continue
    idx=[a.var_names.get_loc(g) for g in present]; x=Xall[:,idx]; x=x.toarray() if sparse.issparse(x) else np.asarray(x)
    npos=(x>0).sum(axis=1); positives[ct]=npos>=MIN_POSITIVE_GENES
    a.obs[f'rule_{ct}_n_positive_genes']=npos
    a.obs[f'rule_{ct}_positive']=positives[ct]
pd.DataFrame(rates).to_csv(TAB/'rule_marker_detection_rates.csv',index=False)
pd.DataFrame(availability).to_csv(TAB/'rule_marker_availability.csv',index=False)
mat=np.column_stack([positives[ct] for ct in MARKER_PANEL]); n=mat.sum(axis=1); names=np.array(list(MARKER_PANEL),dtype=object)
labels=np.full(a.n_obs,'Unassigned',dtype=object); one=n==1; labels[one]=names[mat[one].argmax(axis=1)]; labels[n>1]='Ambiguous'
a.obs['cell_type_rule']=pd.Categorical(labels); a.obs['rule_n_positive_panels']=n
a.obs.cell_type_rule.value_counts().to_csv(TAB/'rule_label_counts.csv')
a.write_h5ad(OUT/'GSM9046088_FOV46_approach_i_ii.h5ad',compression='gzip')
print(a.obs.cell_type_rule.value_counts()); print(pd.DataFrame(availability).to_string(index=False))
