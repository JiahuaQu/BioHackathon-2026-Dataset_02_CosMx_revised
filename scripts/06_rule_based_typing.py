#!/usr/bin/env python
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT))
import pandas as pd
import scanpy as sc
from config.markers import MARKER_PANEL
from src.typing import detection_rate,rule_based_all_markers
IN=ROOT/"results/FOV46/GSM9046088_FOV46_approach_i.h5ad"; OUT=ROOT/"results/FOV46"; TAB=OUT/"tables"
a=sc.read_h5ad(IN)
# Preserve the original requested logic: ALL markers in a candidate panel must have raw count > 0.
# We use layers['counts']; normalization/log1p therefore cannot alter the rule.
rates=[]
for ct,genes in MARKER_PANEL.items():
    r=detection_rate(a,genes); rates.extend({"cell_type":ct,"gene":g,"detection_rate":v} for g,v in r.items())
pd.DataFrame(rates).to_csv(TAB/"rule_marker_detection_rates.csv",index=False)
availability=rule_based_all_markers(a,MARKER_PANEL,layer="counts"); availability.to_csv(TAB/"rule_marker_availability.csv",index=False)
a.obs.cell_type_rule.value_counts().to_csv(TAB/"rule_label_counts.csv")
a.write_h5ad(OUT/"GSM9046088_FOV46_approach_i_ii.h5ad",compression="gzip")
print(a.obs.cell_type_rule.value_counts()); print(availability.to_string(index=False))
