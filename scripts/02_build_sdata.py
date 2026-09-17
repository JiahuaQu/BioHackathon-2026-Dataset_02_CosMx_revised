#!/usr/bin/env python
"""Build SpatialData directly from the five GEO CosMx flat files.
No morphology image is available in these five GEO files; this object contains table,
cell-boundary polygons, and transcript points. Global pixel coordinates are used.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import dask.dataframe as dd
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from spatialdata import SpatialData
from spatialdata.models import PointsModel, ShapesModel, TableModel
from spatialdata.transformations import Identity
from src.cosmx_io import FILES, build_raw_anndata, unique_cell_id
RAW=ROOT/"data/raw"
adata=build_raw_anndata(RAW)
# polygons: one Polygon per cell, indexed identically to the table instance key
poly=pd.read_csv(RAW/FILES["polygons"])
cell_col="cell_ID" if "cell_ID" in poly.columns else "cellID"
need={"fov",cell_col,"x_global_px","y_global_px"}
if not need.issubset(poly.columns): raise ValueError(f"polygon file missing {need-set(poly.columns)}")
poly["unique_cell_id"]=unique_cell_id(pd.to_numeric(poly.fov).astype(int),pd.to_numeric(poly[cell_col]).astype(int))
records=[]
for uid,g in poly.groupby("unique_cell_id",sort=False):
    xy=g[["x_global_px","y_global_px"]].to_numpy(float)
    if len(xy)>=3:
        geom=Polygon(xy)
        if not geom.is_valid: geom=geom.buffer(0)
        if not geom.is_empty: records.append((uid,geom))
shapes=gpd.GeoDataFrame({"geometry":[x[1] for x in records]},index=pd.Index([x[0] for x in records],name="instance_id"),crs=None)
shapes_model=ShapesModel.parse(shapes,transformations={"global":Identity()})
# TableModel links table rows to shape index.
adata.obs["region"]="cell_boundaries"; adata.obs["region"]=adata.obs["region"].astype("category")
adata.obs["instance_id"]=adata.obs_names.astype(str)
# Keep only cells for which a polygon exists to maintain a valid region/instance mapping.
keep=adata.obs_names.intersection(shapes.index); adata=adata[keep].copy()
table=TableModel.parse(adata,region="cell_boundaries",region_key="region",instance_key="instance_id")
# transcript points remain lazy with dask
TX=RAW/FILES["transcripts"]
head=pd.read_csv(TX,nrows=3); cols=set(head.columns)
x="x_global_px"; y="y_global_px"; feature="target"
if not {x,y,feature}.issubset(cols): raise ValueError(f"transcript file missing required columns: {x,y,feature}")
tx=dd.read_csv(TX,blocksize="64MB",assume_missing=True)
if "cell_ID" in tx.columns and "fov" in tx.columns:
    # SpatialData PointsModel instance key can reference this study-wide ID.
    tx["unique_cell_id"]="fov_"+tx["fov"].astype("Int64").astype(str)+"_cell_"+tx["cell_ID"].astype("Int64").astype(str)
points=PointsModel.parse(tx,coordinates={"x":x,"y":y},feature_key=feature,instance_key="unique_cell_id" if "unique_cell_id" in tx.columns else None,transformations={"global":Identity()})
sdata=SpatialData(points={"transcripts":points},shapes={"cell_boundaries":shapes_model},tables={"table":table})
out=ROOT/"data/spatial/GSM9046088_CosMx.zarr"; out.parent.mkdir(parents=True,exist_ok=True)
sdata.write(out,overwrite=True)
print(sdata); print("saved",out)
