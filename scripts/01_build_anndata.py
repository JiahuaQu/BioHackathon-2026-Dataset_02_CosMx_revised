#!/usr/bin/env python
from pathlib import Path
import sys, numpy as np, scanpy as sc
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
from src.cosmx_io import build_raw_anndata
raw=build_raw_anndata(ROOT/"data/raw")
sc.pp.calculate_qc_metrics(raw,percent_top=None,log1p=False,inplace=True)
out=ROOT/"data/processed/GSM9046088_CosMx_raw.h5ad"; out.parent.mkdir(parents=True,exist_ok=True)
raw.write_h5ad(out,compression="gzip")
print(raw); print("saved",out)
