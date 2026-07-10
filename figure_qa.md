# Figure QA Record

## Scope and provenance

- Figures: 1-4 main manuscript figures.
- Backend: Python/matplotlib only.
- Plotting script: `generate_submission_figures.py`.
- Source-data preparation script: `prepare_figure_source_data.py`.
- Quantitative source data: 13 CSV files in `source_data/`, indexed by `source_data_manifest.csv`.
- No raw patient-level data, microscopy, blots, or other image panels are included.

## Export verification

| Figure | SVG | PDF | TIFF | PNG | TIFF resolution | PNG resolution |
|---|---|---|---|---|---|---|
| 1 | Pass | Pass | Pass | Pass | 600 dpi | 300 dpi |
| 2 | Pass | Pass | Pass | Pass | 600 dpi | 300 dpi |
| 3 | Pass | Pass | Pass | Pass | 600 dpi | 300 dpi |
| 4 | Pass | Pass | Pass | Pass | 600 dpi | 300 dpi |

SVG files retain text nodes and PDFs use embedded TrueType-compatible text settings. All figures use a white background, consistent sans-serif typography, and a colorblind-considered blue/orange/teal/gray palette. Panel labels are lowercase bold letters. The exported canvases are double-column-width layouts.

## Frontiers raster check

The four TIFF files copied to the Frontiers local submission package are RGB, white-background, 600 dpi rasters. This exceeds the journal checklist's 300 dpi minimum for separately uploaded TIFF/JPEG figures. The SVG and PDF files remain in the primary figure package as editable archival versions.

## Visual QA

- Figure 1: fixed score, cohort roles, and endpoint separation are legible; GSE68465 and HPA are explicitly non-primary evidence layers.
- Figure 2: stage, descriptive Kaplan-Meier, and continuous/adjusted Cox evidence are visually separated.
- Figure 3: external OS evidence is foregrounded; forest labels, confidence intervals, HR text, and C-index comparison do not overlap.
- Figure 4: endpoint is RFS; the gender-stratified sensitivity estimate is labelled as a sensitivity analysis; forest labels, confidence intervals, HR text, and C-index comparison do not overlap.

## Statistical traceability

| Figure | Cohort and endpoint | N / events | Primary display | Supporting display |
|---|---|---:|---|---|
| 2 | TCGA-LUAD OS | 502 / 182 | Continuous Cox axis HR 1.46 (95% CI 1.24-1.73), p=6.91e-06 | Adjusted complete-case model: 484 / 177, HR 1.37 (1.15-1.64), p=4.01e-04 |
| 3 | GSE72094 OS | 398 / 113 | Continuous exact-score Cox HR 1.54 (1.27-1.87), p=1.23e-05 | Adjusted model: 393 / 111, HR 1.51 (1.23-1.85), p=8.68e-05; optimism-adjusted C-index 0.648 to 0.676; LRT p=1.14e-04 |
| 4 | GSE31210 RFS | 204 / 54 | Continuous exact-score Cox HR 1.54 (1.11-2.12), p=0.0091 | Adjusted HR 1.40 (1.01-1.93), p=0.0436; gender-stratified sensitivity HR 1.41 (1.02-1.95), p=0.0389; optimism-adjusted C-index 0.674 to 0.681; LRT p=0.0458 |

The median-split Kaplan-Meier panels are descriptive visualizations only. No biological or technical replicate framework applies because these are retrospective patient-level public cohorts. No multiple-comparison-adjusted claim is made in the figure package.

## Integrity statement

The figures contain no microscopy, blots, gels, or photographs. No image-specific contrast adjustment, cropping, stitching, pseudo-coloring, or reuse assessment is applicable. All plotted quantitative values can be traced to the source-data CSVs and analysis outputs.
