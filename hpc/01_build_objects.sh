#!/bin/bash
#BSUB -P spatial
#BSUB -J cosmx_build
#BSUB -q superdome
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=8000]"
#BSUB -oo logs/build_%J.log
#BSUB -eo logs/build_%J.err
set -euo pipefail
module load conda3/202402
ROOT="/research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02_CosMx_revised/Dataset_02_CosMx_revised"
PYTHON="/home/jqu/.conda/envs/spatialdata/bin/python"
cd "$ROOT"
echo "Python: $PYTHON"; "$PYTHON" --version
"$PYTHON" scripts/01_build_anndata.py
"$PYTHON" scripts/02_build_sdata.py
