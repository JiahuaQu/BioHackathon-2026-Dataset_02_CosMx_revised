#!/usr/bin/env python
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import spatialdata as sd
import spatialdata_plot  # register .pl
from src.spatial_plot_legacy import plot_gene_seg_transcripts
S=sd.read_zarr(ROOT/"data/spatial/GSM9046088_CosMx_FOV46.zarr"); OUT=ROOT/"results/FOV46/figures"
# GEO's five flat files do not include morphology images. These overlays therefore show
# cell-boundary segmentation + cell-level expression + individual transcript molecules.
for gene in ["EPCAM","SFTPC","CD3D","MS4A1","C1QA","PECAM1","COL1A1","RGS5"]:
    if gene not in S["table"].var_names: continue
    try:
        ax=plot_gene_seg_transcripts(S,gene,point_size=1.2,point_alpha=.65,figsize=(8,8))
        # spatialdata-plot controls display; save via current matplotlib figure if available.
        import matplotlib.pyplot as plt
        plt.gcf().savefig(OUT/f"sdata_overlay_{gene}.png",dpi=250,bbox_inches="tight"); plt.close("all")
    except Exception as e: print(f"overlay {gene} failed: {e}")
