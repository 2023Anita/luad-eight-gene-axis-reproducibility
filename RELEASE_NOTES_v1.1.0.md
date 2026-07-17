# v1.1.0: Complete TCGA Discovery Reproduction

This release closes the remaining TCGA-LUAD discovery-code gap in the public
reproducibility package.

## Added

- Complete TCGA Xena download, extraction, clinical harmonization, exploratory
  gene screening, axis construction, formal survival modeling, and verification
  workflow.
- Full 28-gene candidate decision table, including retained and not-retained
  genes; the selection is explicitly documented as analyst-guided and
  outcome-informed.
- Versioned TCGA result tables supporting the discovery analysis.
- One-command clean reproduction through `scripts/run_tcga_discovery.sh`.
- Automated output verification, software-version record, release checksums,
  and a documented successful clean run.

## Reproduce

```bash
bash scripts/run_tcga_discovery.sh
```

The workflow downloads public TCGA-LUAD source matrices from UCSC Xena. Raw
source matrices and regenerated working directories are not included in the
release.
