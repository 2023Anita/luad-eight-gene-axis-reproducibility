# Release Manifest

## Included Derived Materials

| Directory | Contents | Public-data status |
|---|---|---|
| `scripts/` | GSE72094/GSE31210 extraction and locked-score analysis scripts; gender-stratified sensitivity script; figure source-data and rendering scripts | Code only |
| `results/` | Derived missingness, Cox, KM, C-index, likelihood-ratio, PH, and sensitivity-result CSVs | Derived from public cohorts |
| `figure_source_data/` | CSVs supporting the planned four figures and their source manifest | Derived/summary data |
| `figures/` | Four rendered figure bundles in SVG, PDF, TIFF, and PNG | Derived visualizations |
| `figure_qa.md` and `rendering_environment.md` | Export, visual-QA, provenance, and Python environment records | Documentation only |
| `zenodo_metadata_draft.json` | Draft metadata for a future persistent record | Not uploaded; no DOI |
| `LICENSE-CODE-MIT.txt` and `LICENSE-DATA-CC-BY-4.0.txt` | Code and derived-material license notices | Reuse terms |

## Explicitly Excluded

- Raw GEO series matrices, raw CEL archives, and platform downloads.
- TCGA raw/processed source matrices that are recoverable from their public resources.
- Submission-system material, correspondence, credentials, browser data, and non-public author information.

## Release Gate

This package is ready for public deposition only after the authors select a repository and license, approve the content, and add a persistent identifier.
