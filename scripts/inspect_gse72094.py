#!/usr/bin/env python3
"""Extract locked-score inputs and clinical metadata from GSE72094."""

from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "raw"
OUT = ROOT / "processed"
GENES = ("CA9", "MKI67", "STAT1", "HIF1A", "CD4", "CD8A", "CD8B", "TIGIT")


def tsv(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n\r")], delimiter="\t", quotechar='"'))


def platform_probes(path: Path) -> dict[str, list[str]]:
    probes = {gene: [] for gene in GENES}
    in_table = False
    header: list[str] = []
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            line = line.rstrip("\n\r")
            if line == "!platform_table_begin":
                in_table = True
                header = []
                continue
            if line == "!platform_table_end":
                break
            if not in_table:
                continue
            row = line.split("\t")
            if not header:
                header = row
                continue
            if len(row) < 4:
                continue
            values = dict(zip(header, row))
            probe = values.get("ID", "").strip()
            symbols = values.get("GeneSymbol", "").upper().replace("///", ";")
            if not probe:
                continue
            tokens = {item.strip() for item in symbols.split(";") if item.strip()}
            for gene in GENES:
                if gene in tokens:
                    probes[gene].append(probe)
    return probes


def parse_series(path: Path, probe_map: dict[str, list[str]]) -> tuple[list[str], dict[str, list[str]], dict[str, list[float]]]:
    sample_ids: list[str] = []
    metadata: dict[str, list[str]] = {}
    expression: dict[str, list[float]] = {}
    target_probes = {probe for probes in probe_map.values() for probe in probes}
    in_matrix = False

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                row = tsv(line)
                label, values = row[0], row[1:]
                if label == "!Sample_geo_accession":
                    sample_ids = values
                elif label == "!Sample_title":
                    metadata["sample_title"] = values
                elif label == "!Sample_characteristics_ch1":
                    for index, value in enumerate(values):
                        if ":" not in value:
                            continue
                        key, item = value.split(":", 1)
                        key = key.strip().lower().replace(" ", "_")
                        column = metadata.setdefault(key, [""] * len(values))
                        column[index] = item.strip()
                continue
            if line.startswith("!series_matrix_table_begin"):
                in_matrix = True
                header = tsv(next(handle))
                if not sample_ids:
                    sample_ids = header[1:]
                continue
            if not in_matrix:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            row = tsv(line)
            if not row or row[0] not in target_probes:
                continue
            try:
                expression[row[0]] = [float(value) for value in row[1:]]
            except ValueError:
                continue
    return sample_ids, metadata, expression


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * fraction
    lower, upper = int(position), min(int(position) + 1, len(ordered) - 1)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    probe_map = platform_probes(RAW / "GPL15048_full.txt")
    sample_ids, metadata, expression = parse_series(RAW / "GSE72094_series_matrix.txt.gz", probe_map)

    selected: dict[str, str] = {}
    mapping_rows: list[dict[str, object]] = []
    for gene, probes in probe_map.items():
        best_probe, best_iqr = "", float("-inf")
        for probe in probes:
            values = expression.get(probe, [])
            score = quantile(values, 0.75) - quantile(values, 0.25) if values else float("nan")
            mapping_rows.append({"gene": gene, "probe": probe, "found_in_matrix": bool(values), "iqr": score})
            if values and score > best_iqr:
                best_probe, best_iqr = probe, score
        if best_probe:
            selected[gene] = best_probe

    with (OUT / "gse72094_probe_mapping.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene", "probe", "found_in_matrix", "iqr"])
        writer.writeheader()
        writer.writerows(mapping_rows)

    clinical_fields = ["sample_id", *sorted(metadata)]
    with (OUT / "gse72094_clinical_raw.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=clinical_fields)
        writer.writeheader()
        for index, sample_id in enumerate(sample_ids):
            writer.writerow({field: sample_id if field == "sample_id" else metadata[field][index] for field in clinical_fields})

    with (OUT / "gse72094_locked_score_expression.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", *GENES])
        writer.writeheader()
        for index, sample_id in enumerate(sample_ids):
            row = {"sample_id": sample_id}
            for gene in GENES:
                probe = selected.get(gene)
                row[gene] = expression[probe][index] if probe else ""
            writer.writerow(row)

    summary = {
        "series": "GSE72094",
        "platform": "GPL15048",
        "sample_count": len(sample_ids),
        "required_genes": list(GENES),
        "mapped_genes": sorted(selected),
        "missing_genes": [gene for gene in GENES if gene not in selected],
        "clinical_fields": sorted(metadata),
        "exact_score_eligible": len(selected) == len(GENES),
    }
    (OUT / "gse72094_feasibility.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
