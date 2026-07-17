# Reproducibility Materials for the LUAD Eight-Gene Axis Study

This versioned package contains reproducibility materials for a public multi-cohort lung adenocarcinoma transcriptomic study. Public source matrices are downloaded from their repositories and are not duplicated in Git. Derived result tables, source data, figures, and executable analysis code are included.

## Included Now

- Complete TCGA-LUAD discovery workflow, including the original 28-gene panel, exploratory screen, retained/not-retained decisions, equal-weight axis construction, and formal survival models.
- External-cohort extraction and locked-score analysis scripts.
- Derived probe maps, missingness tables, Cox results, model-comparison results, and proportional-hazards diagnostics.
- Figure source-data CSVs, Python-only figure-generation scripts, rendered SVG/PDF/TIFF/PNG figure bundles, and the figure QA record.
- A Zenodo metadata draft for the persistent record; it is not itself a DOI.
- Clean-run verification, environment records, checksums, and the release manifest.

## Reproduce the TCGA Discovery Phase

Requirements are Python 3.10+ and R with the `survival` package. From the repository root:

```bash
bash scripts/run_tcga_discovery.sh
```

The command downloads the two public TCGA Xena inputs, runs the complete discovery workflow, and compares the generated tables with the versioned expected outputs. See `TCGA_REPRODUCIBILITY.md` for details.

## Release Identifiers

- GitHub repository URL: https://github.com/2023Anita/luad-eight-gene-axis-reproducibility
- Zenodo version DOI: not assigned; GitHub release is the current versioned record.
- License: code is MIT; derived tables, figure source data, figures, and documentation are CC BY 4.0.

## Exclude from Release

- Raw TCGA Xena matrices, raw GEO series-matrix downloads, raw CEL archives, and all large public-source files that can be recovered from their official repository pages.
- Any credentials, browser data, correspondence, submission-system exports, or author personal contact material beyond what appears on the manuscript title page.

## Planned Deposit Metadata

- Title: Reproducibility materials for the LUAD eight-gene hypoxia/proliferation-immune activation axis study.
- Access: public GitHub repository and versioned GitHub release.
- Persistent identifier: GitHub release tag; Zenodo DOI not yet assigned.
