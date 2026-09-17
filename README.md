# Dataset_02 — GSM9046088 CosMx revised workflow

This folder replaces manual GEO downloading with a reproducible code-first pipeline and keeps **raw AnnData**, **processed AnnData**, and **SpatialData** separate.

## Environment

```bash
module load conda3/202402
conda env create -f env/environment.yml
conda activate spatialdata
```

On this HPC, batch scripts deliberately call the environment's Python explicitly:
`/home/jqu/.conda/envs/spatialdata/bin/python`.
Change `ROOT` and `PYTHON` in `hpc/*.sh` if your copied project path differs.

## Run from scratch

```bash
cd /research/rgs01/home/clusterHome/jqu/activities/learning/BioHackathon/BioHackathon-2026/Dataset_02
mkdir -p logs
bsub < hpc/00_download_raw.sh
# after download finishes:
bsub < hpc/01_build_objects.sh
# after objects finish:
bsub < hpc/02_explore_whole_slide.sh
# review QC + resolution figures/tables, then:
bsub < hpc/03_fov46_pipeline.sh
```

Important: `scripts/05_analyze_fov46.py` intentionally leaves `CLUSTER_ANNOTATIONS = {}`. Review `results/FOV46/tables/leiden_markers.csv`, `marker_dotplot.png`, UMAP and spatial plots, fill the mapping, then rerun scripts 05–08.

## Directory roles

- `data/raw/`: five GEO CosMx flat files only; never modified.
- `data/processed/`: raw AnnData and FOV46 raw AnnData.
- `data/spatial/`: whole-slide and FOV46 SpatialData Zarr stores.
- `results/whole_slide/`: QC, FOV-dominance diagnostics, multi-resolution Leiden.
- `results/FOV46/`: Approach I, Approach II, comparisons and spatial overlays.
- `src/`: reusable code; notebooks/scripts should import from here.
- `config/markers.py`: one source of truth for marker panels.

## Analysis logic

1. Download all five GEO files programmatically.
2. Build raw AnnData from expression + metadata + FOV positions.
3. Build SpatialData from AnnData table + polygons + transcript points.
4. Whole-slide QC: global and per-FOV distributions, robust MAD thresholds, threshold grid.
5. Whole-slide expression graph and Leiden resolution sweep. Spatial distance is **not** used in the Scanpy neighbor graph; cluster/FOV association therefore diagnoses expression differences correlated with tissue/FOV.
6. Extract FOV46 from the **raw** AnnData. Extract FOV46 SpatialData with padding checks and strict fallback.
7. Approach I: within-FOV Scanpy clustering + marker-based manual annotation.
8. Approach II: rule-based typing using the user's original `all markers > 0` rule across all candidate panels.
9. Compare Approach I vs II vs `final_CT` if available, and validate spatially.

## QC threshold policy

The code does not claim a universal CosMx cutoff. It uses:
- hard information floors of 20 total transcripts and 10 detected genes;
- robust 3-MAD lower bounds on log1p(counts/genes), with per-FOV flags;
- a 3-MAD upper bound for control-probe fraction;
- robust cell-area bounds when cell area exists;
- a grid of alternative thresholds written to `qc_threshold_sweep.csv` for sensitivity analysis.

Review retention by FOV before accepting thresholds. The goal is to remove low-information/segmentation outliers without selectively deleting entire tissue regions.

## Leiden resolution policy

`resolution_sweep()` runs 0.1–1.2, reports cluster number, tiny-cluster burden, ARI and NMI versus the preceding resolution, and produces multi-resolution UMAPs. The automated selection is only a **heuristic shortlist**: the first stable solution (ARI >= 0.8) with <=5% of cells in tiny clusters. Final selection should also be supported by marker separation and spatial/biological interpretability, analogous in spirit to inspecting a Seurat clustree.
