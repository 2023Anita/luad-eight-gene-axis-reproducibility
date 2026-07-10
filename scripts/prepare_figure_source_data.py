#!/usr/bin/env python3
"""Create figure-source CSVs from the locked analysis outputs without plotting."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
FIGURES = Path(__file__).resolve().parent
SOURCE = FIGURES / "source_data"
SOURCE.mkdir(parents=True, exist_ok=True)

FILES = {
    "figure2_tcga_stage_summary.csv": ROOT / "data" / "processed" / "axis_score_stage_summary_r.csv",
    "figure2_tcga_km_curve.csv": ROOT / "data" / "processed" / "km_curve_axis_group_median.csv",
    "figure2_tcga_cox_results.csv": ROOT / "data" / "processed" / "survival_cox_results.csv",
    "figure3_gse72094_km_curve.csv": ROOT / "rejection_revision_validation" / "results" / "gse72094_km_curve_data.csv",
    "figure3_gse72094_cox_results.csv": ROOT / "rejection_revision_validation" / "results" / "gse72094_cox_results.csv",
    "figure3_gse72094_cindex.csv": ROOT / "rejection_revision_validation" / "results" / "gse72094_bootstrap_cindex.csv",
    "figure3_gse72094_lrt.csv": ROOT / "rejection_revision_validation" / "results" / "gse72094_likelihood_ratio_test.csv",
    "figure4_gse31210_km_curve.csv": ROOT / "rejection_revision_validation" / "results" / "gse31210_km_curve_data.csv",
    "figure4_gse31210_cox_results.csv": ROOT / "rejection_revision_validation" / "results" / "gse31210_cox_results.csv",
    "figure4_gse31210_cindex.csv": ROOT / "rejection_revision_validation" / "results" / "gse31210_bootstrap_cindex.csv",
    "figure4_gse31210_lrt.csv": ROOT / "rejection_revision_validation" / "results" / "gse31210_likelihood_ratio_test.csv",
    "figure4_gse31210_gender_stratified.csv": ROOT / "rejection_revision_validation" / "results" / "gse31210_gender_stratified_sensitivity.csv",
}


def validate_csv(path: Path) -> None:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header:
            raise ValueError(f"Missing CSV header: {path}")


def main() -> None:
    manifest_rows: list[dict[str, str]] = []
    figure1 = SOURCE / "figure1_cohort_hierarchy.csv"
    figure1.write_text(
        "cohort,role,endpoint,eligible_n,events_or_annotation,exact_score\n"
        "TCGA-LUAD,discovery,OS,502,182,yes\n"
        "GSE72094,primary_external,OS,398,113,yes\n"
        "GSE31210,secondary_external,RFS,204,54,yes\n"
        "GSE68465,sensitivity,OS,442,236,no_TIGIT\n"
        "HPA,annotation,not_applicable,0,annotation_only,not_applicable\n",
        encoding="utf-8",
    )
    manifest_rows.append({"source_data_file": figure1.name, "origin": "S2/S3 locked cohort hierarchy"})
    for target_name, source_path in FILES.items():
        validate_csv(source_path)
        target = SOURCE / target_name
        shutil.copyfile(source_path, target)
        manifest_rows.append({"source_data_file": target_name, "origin": str(source_path.relative_to(ROOT))})

    with (SOURCE / "source_data_manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["source_data_file", "origin"])
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {len(manifest_rows)} source-data files")


if __name__ == "__main__":
    main()
