# Release Manifest

## Included Derived Materials

| Directory | Contents | Public-data status |
|---|---|---|
| `scripts/` | Complete TCGA discovery workflow; GSE72094/GSE31210 extraction and locked-score analysis scripts; gender-stratified sensitivity; clean-run verification; figure source-data and rendering scripts | Code only |
| `results/tcga/` | Complete 28-gene screen, retained/not-retained decisions, formal modeling dataset, Cox, log-rank, KM, and stage-summary outputs | Derived from public TCGA Xena data |
| `results/` | External-cohort missingness, Cox, KM, C-index, likelihood-ratio, PH, and sensitivity-result CSVs | Derived from public GEO cohorts |
| `figure_source_data/` | CSVs supporting the planned four figures and their source manifest | Derived/summary data |
| `figures/` | Four rendered figure bundles in SVG, PDF, TIFF, and PNG | Derived visualizations |
| `figure_qa.md` and `rendering_environment.md` | Export, visual-QA, provenance, and Python environment records | Documentation only |
| `zenodo_metadata_draft.json` | Draft metadata for a future persistent record | Not uploaded; no DOI |
| `LICENSE-CODE-MIT.txt` and `LICENSE-DATA-CC-BY-4.0.txt` | Code and derived-material license notices | Reuse terms |
| `CITATION.cff` | Citation metadata for the public GitHub repository | Citation support |
| `TCGA_REPRODUCIBILITY.md`, `environment/`, and `checksums.sha256` | Clean-run instructions, software versions, and release-file integrity hashes | Reproduction support |

## Explicitly Excluded

- Raw GEO series matrices, raw CEL archives, and platform downloads.
- TCGA Xena source matrices and generated `data/` working directories that are recoverable from public resources and the released scripts.
- Submission-system material, correspondence, credentials, browser data, and non-public author information.

## Release Gate

The next release is intended as `v1.1.0`. It closes the TCGA discovery-code gap while retaining the existing repository and licenses. A GitHub release tag will provide the versioned URL used in the manuscript.
