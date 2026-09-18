# Interactive Jupyter workflow — Dataset_02 CosMx

This directory is the interactive/exploratory counterpart to `scripts/` + `hpc/`. It deliberately avoids `%run ../scripts/...` for the analysis stages so parameters, intermediate objects, tables, and figures remain visible in the notebook.

Recommended order:

1. `00_setup_data_and_objects_interactive.ipynb`
2. `01_whole_slide_exploration_interactive.ipynb`
3. `02_fov46_qc_interactive.ipynb`
4. `03_fov46_leiden_sweep_interactive.ipynb`
5. `04_approach_i_interactive.ipynb`
6. `05_approach_ii_rule_based_interactive.ipynb`
7. `06_compare_approach_i_vs_ii_interactive.ipynb`
8. `07_spatial_visualization_interactive.ipynb`

## Approach-II parameter experiments

In notebook 05 change `MIN_POSITIVE_GENES`. The default tag is automatically derived from the parameter (`min2`, `min3`, `min4`, ...). Saved variants are isolated under:

```text
results/FOV46/approach_ii_variants/<RUN_TAG>/
```

For example:

```text
approach_ii_variants/min2/
approach_ii_variants/min3/
approach_ii_variants/min4/
```

Notebook 06 uses `RULE_RUN_TAG` to choose which Approach-II variant is compared with Approach I and `final_CT`. Notebook 07 uses the same tag for spatial visualization. This prevents experimental results from overwriting each other.

## Production vs exploration

Use these notebooks to inspect and tune parameters. Once a final choice is made, synchronize the parameter/configuration back to `scripts/` and use the HPC jobs for the final reproducible run.





# CosMx interactive Jupyter workflow v3

This revision standardizes all notebooks 02–07 around the same robust project-root initialization and versioned output pattern.

## Key rules

1. `PROJECT_ROOT` is defined exactly once in the first code cell by searching parent directories for both `src/` and `data/`. Folder names such as `jupyter`, `jupyter-2`, or `notebooks` do not matter.
2. Interactive objects remain in memory unless the corresponding `SAVE_*` flag is set to `True`.
3. Every experimental save uses a tag and prints the planned directory, filename, full path, and whether the file already exists before writing.
4. Approach-II sensitivity runs are separated by `RUN_TAG` (`min2`, `min3`, `min4`, `min3_testA`, etc.).
5. Comparison notebook 06 selects a saved Approach-II variant by `RULE_RUN_TAG`, allowing Approach I to be compared independently with each rule threshold.

## Suggested order

- `02_fov46_qc_interactive.ipynb`: tune QC; save under `qc_variants/<QC_TAG>/` only when satisfied.
- `03_fov46_leiden_sweep_interactive.ipynb`: select a QC input, tune PCA/neighbors/resolutions; optionally save under `leiden_variants/<SWEEP_TAG>/`.
- `04_approach_i_interactive.ipynb`: choose resolution, inspect markers, manually annotate; optionally save under `approach_i_variants/<APPROACH_I_TAG>/`.
- `05_approach_ii_rule_based_interactive.ipynb`: change `MIN_POSITIVE_GENES`; each rule is saved independently under `approach_ii_variants/<RUN_TAG>/`.
- `06_compare_approach_i_vs_ii_interactive.ipynb`: set `RULE_RUN_TAG` to compare Approach I/reference with a selected rule variant.
- `07_spatial_visualization_interactive.ipynb`: select the same rule variant, tune spatial/gene-overlay parameters inline, and optionally export tagged figures.

Production `scripts/` and `hpc/` remain the reproducible final pipeline; notebooks are the exploration and parameter-debugging layer.
