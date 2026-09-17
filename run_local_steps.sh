#!/bin/bash
set -euo pipefail
PYTHON="${PYTHON:-/home/jqu/.conda/envs/spatialdata/bin/python}"
for s in 00_download_cosmx.py 01_build_anndata.py 02_build_sdata.py 03_explore_whole_slide.py 04_extract_fov46.py 05_analyze_fov46.py 06_rule_based_typing.py 07_compare_typing.py 08_spatial_overlay.py; do
  echo "=== $s ==="
  "$PYTHON" "scripts/$s"
done
