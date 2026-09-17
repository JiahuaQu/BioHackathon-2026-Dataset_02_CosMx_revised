#!/bin/bash
#BSUB -P spatial
#BSUB -J cosmx_f46
#BSUB -q superdome
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=6000]"
#BSUB -oo logs/fov46_%J.log
#BSUB -eo logs/fov46_%J.err
set -euo pipefail
module load conda3/202402
ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"
cd "$ROOT"
for s in 04_extract_fov46.py 05_analyze_fov46.py 06_rule_based_typing.py 07_compare_typing.py 08_spatial_overlay.py; do
  "$PYTHON" "scripts/$s"
done
