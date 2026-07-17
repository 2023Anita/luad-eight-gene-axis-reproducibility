# TCGA-LUAD Discovery Reproduction

## Scope

This workflow reproduces the exploratory TCGA-LUAD discovery phase reported in
the manuscript. It downloads the public TCGA Xena expression and clinical
matrices, extracts the manually curated 28-gene feasibility panel, derives the
clinical endpoint fields, runs the exploratory gene-level screen, constructs
the retained equal-weight eight-gene axis, and runs the formal survival models.

The eight-gene composition was selected with analyst judgment during an
outcome-informed exploratory phase. It was not selected by a prespecified or
deterministic threshold. The complete retained/not-retained record is stored in
`results/tcga/tcga_candidate_screen_decisions.csv`.

## Requirements

- Python 3.10 or later; no third-party Python packages.
- R 4.3 or later with the `survival` package.
- Internet access for the initial TCGA Xena download.

## Clean Run

From the repository root:

```bash
bash scripts/run_tcga_discovery.sh
```

The source matrices are downloaded to `data/raw/`, and generated files are
written to `data/processed/` and `reports/`. These directories are ignored by
Git because the source data remain available from TCGA Xena and the expected
versioned outputs are already stored under `results/tcga/`.

For a rerun with existing source matrices:

```bash
bash scripts/run_tcga_discovery.sh --skip-download
```

The final command compares every versioned TCGA result table with the clean-run
outputs using exact string comparison or floating-point tolerance and prints:

```text
TCGA clean-run verification: PASS
```

## Data Sources

- TCGA Xena LUAD gene expression: `TCGA.LUAD.sampleMap/HiSeqV2.gz`
- TCGA Xena LUAD clinical matrix: `TCGA.LUAD.sampleMap/LUAD_clinicalMatrix`

The exact public URLs are defined in `scripts/download_luad_starter.py` and used
by `scripts/download_tcga_xena_only.py`.
