#!/usr/bin/env python3
"""Merge the exploratory TCGA screens with the analyst-retention record."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


RISK_GENES = {"CA9", "MKI67", "STAT1", "HIF1A"}
IMMUNE_GENES = {"CD4", "CD8A", "CD8B", "TIGIT"}

RETAINED_BASIS = {
    "MKI67": "Retained as an adverse proliferation marker with the strongest descriptive survival pattern.",
    "CA9": "Retained as an adverse hypoxia marker with strong descriptive survival and stage patterns.",
    "STAT1": "Retained as an adverse interferon-related signal in the exploratory screen.",
    "HIF1A": "Retained for hypoxia-domain coverage despite weak individual survival separation.",
    "CD4": "Retained in the interpretable lymphocyte-related immune component.",
    "CD8A": "Retained in the interpretable cytotoxic lymphocyte-related immune component.",
    "CD8B": "Retained in the interpretable cytotoxic lymphocyte-related immune component.",
    "TIGIT": "Retained in the interpretable lymphocyte/checkpoint immune component.",
}


def read_index(path: Path, key: str) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed", type=Path, default=Path("data/processed"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/processed/tcga_candidate_screen_decisions.csv"),
    )
    args = parser.parse_args()

    survival = read_index(args.processed / "km_logrank_gene_screen.csv", "gene")
    effects = read_index(args.processed / "exploratory_gene_effects.csv", "gene")
    if set(survival) != set(effects):
        raise SystemExit("Candidate genes differ between the survival and clinical-effect screens")

    rows = []
    for gene, screen in sorted(survival.items(), key=lambda item: float(item[1]["logrank_p_approx"])):
        if gene in RISK_GENES:
            decision = "Risk"
        elif gene in IMMUNE_GENES:
            decision = "Immune"
        else:
            decision = "Not selected"
        rows.append(
            {
                "gene": gene,
                "logrank_p_approx": screen["logrank_p_approx"],
                "deceased_minus_living_mean": effects[gene]["deceased_minus_living_mean"],
                "advanced_minus_early_mean": effects[gene]["advanced_minus_early_mean"],
                "median_split_direction": screen["direction"],
                "axis_decision": decision,
                "decision_basis": RETAINED_BASIS.get(
                    gene,
                    "Not retained because the intended model was a compact four-versus-four biological contrast; no deterministic selection threshold was used.",
                ),
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.output} with {len(rows)} candidate decisions")


if __name__ == "__main__":
    main()
