#!/usr/bin/env python3
"""Compare a clean TCGA discovery run with the versioned expected outputs."""

from __future__ import annotations

import csv
import math
import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
GENERATED = ROOT / "data" / "processed"
EXPECTED = ROOT / "results" / "tcga"

FILE_MAP = {
    "km_logrank_gene_screen.csv": "km_logrank_gene_screen.csv",
    "exploratory_gene_effects.csv": "exploratory_gene_effects.csv",
    "survival_cox_results.csv": "survival_cox_results.csv",
    "survival_logrank_results.csv": "survival_logrank_results.csv",
    "axis_score_stage_summary_r.csv": "axis_score_stage_summary.csv",
    "km_curve_axis_group_median.csv": "km_curve_axis_group_median.csv",
    "formal_modeling_dataset.csv": "formal_modeling_dataset.csv",
    "tcga_candidate_screen_decisions.csv": "tcga_candidate_screen_decisions.csv",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def equal_value(left: str, right: str) -> bool:
    if left == right:
        return True
    try:
        a = float(left)
        b = float(right)
    except ValueError:
        return False
    return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)


def compare(generated: Path, expected: Path) -> None:
    actual_rows = rows(generated)
    expected_rows = rows(expected)
    if len(actual_rows) != len(expected_rows):
        raise AssertionError(f"Row-count mismatch: {generated} vs {expected}")
    for index, (actual, reference) in enumerate(zip(actual_rows, expected_rows), start=2):
        if actual.keys() != reference.keys():
            raise AssertionError(f"Column mismatch in {generated}")
        for field in actual:
            if not equal_value(actual[field], reference[field]):
                raise AssertionError(
                    f"Mismatch in {generated.name}, row {index}, field {field}: "
                    f"{actual[field]!r} != {reference[field]!r}"
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare generated TCGA discovery outputs with versioned expected tables."
    )
    parser.add_argument(
        "--generated-dir",
        type=Path,
        default=GENERATED,
        help=f"Generated output directory (default: {GENERATED})",
    )
    parser.add_argument(
        "--expected-dir",
        type=Path,
        default=EXPECTED,
        help=f"Expected output directory (default: {EXPECTED})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for generated_name, expected_name in FILE_MAP.items():
        compare(args.generated_dir / generated_name, args.expected_dir / expected_name)

    decisions = rows(args.generated_dir / "tcga_candidate_screen_decisions.csv")
    retained = {row["gene"] for row in decisions if row["axis_decision"] != "Not selected"}
    expected_retained = {"CA9", "MKI67", "STAT1", "HIF1A", "CD4", "CD8A", "CD8B", "TIGIT"}
    if len(decisions) != 28 or retained != expected_retained:
        raise AssertionError("The 28-gene provenance or retained eight-gene set changed")
    print("TCGA clean-run verification: PASS")


if __name__ == "__main__":
    main()
