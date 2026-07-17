#!/usr/bin/env python3
"""Prepare formal modeling inputs for LUAD axis survival analysis."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"

CLINICAL = PROCESSED / "clinical_os_derived.csv"
EXPR = PROCESSED / "luad_tme_immune_gene_expression_subset.csv"

RISK_GENES = ["CA9", "MKI67", "STAT1", "HIF1A"]
IMMUNE_GENES = ["CD4", "CD8A", "CD8B", "TIGIT"]
SENSITIVITY_NO_STAT1 = ["CA9", "MKI67", "HIF1A"]
SENSITIVITY_ALT_IMMUNE = ["CD4", "CD8A", "CD8B", "GZMB", "PRF1"]


def parse_float(value: object) -> float | None:
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def zscores(values: dict[str, float]) -> dict[str, float]:
    vals = list(values.values())
    avg = mean(vals)
    sd = math.sqrt(sum((x - avg) ** 2 for x in vals) / (len(vals) - 1)) if len(vals) > 1 else 1.0
    if sd == 0:
        sd = 1.0
    return {sample: (value - avg) / sd for sample, value in values.items()}


def load_expression() -> dict[str, dict[str, float]]:
    with EXPR.open(encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        samples = header[1:]
        genes = {}
        for row in reader:
            genes[row[0]] = {
                sample: float(value)
                for sample, value in zip(samples, row[1:])
                if value not in ("", "NA", "nan")
            }
    return genes


def load_clinical() -> dict[str, dict[str, object]]:
    out = {}
    with CLINICAL.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row.get("sampleID", "")
            if not sample.endswith("-01"):
                continue
            time = parse_float(row.get("OS_time_days"))
            event = parse_float(row.get("OS_event"))
            if time is None or event is None or time <= 0:
                continue
            age = parse_float(row.get("age"))
            stage = str(row.get("stage_group", "")).strip()
            gender = str(row.get("gender", "")).strip()
            out[sample] = {
                "OS_time_days": time,
                "OS_event": int(event),
                "stage_group": stage,
                "age": age,
                "gender": gender,
            }
    return out


def mean_score(z: dict[str, dict[str, float]], genes: list[str], sample: str) -> float | None:
    values = [z[gene][sample] for gene in genes if sample in z.get(gene, {})]
    if len(values) != len(genes):
        return None
    return sum(values) / len(values)


def quantile_groups(values: list[float], cuts: int) -> list[float]:
    ordered = sorted(values)
    return [ordered[int(len(ordered) * i / cuts)] for i in range(1, cuts)]


def group_by_cut(value: float, cutoffs: list[float], labels: list[str]) -> str:
    for cutoff, label in zip(cutoffs, labels):
        if value < cutoff:
            return label
    return labels[-1]


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    expr = load_expression()
    clinical = load_clinical()
    z = {gene: zscores(expr[gene]) for gene in set(RISK_GENES + IMMUNE_GENES + SENSITIVITY_ALT_IMMUNE)}
    needed = set(RISK_GENES + IMMUNE_GENES)
    samples = sorted(set(clinical).intersection(*(set(z[gene]) for gene in needed)))
    rows = []
    for sample in samples:
        risk = mean_score(z, RISK_GENES, sample)
        immune = mean_score(z, IMMUNE_GENES, sample)
        risk_no_stat1 = mean_score(z, SENSITIVITY_NO_STAT1, sample)
        immune_alt = mean_score(z, SENSITIVITY_ALT_IMMUNE, sample)
        if risk is None or immune is None or risk_no_stat1 is None or immune_alt is None:
            continue
        axis = risk - immune
        axis_no_stat1 = risk_no_stat1 - immune
        axis_alt_immune = risk - immune_alt
        c = clinical[sample]
        rows.append(
            {
                "sampleID": sample,
                "OS_time_days": c["OS_time_days"],
                "OS_event": c["OS_event"],
                "age": "" if c["age"] is None else c["age"],
                "gender": c["gender"],
                "stage_group": c["stage_group"],
                "risk_score": risk,
                "immune_score": immune,
                "axis_score": axis,
                "axis_score_no_STAT1": axis_no_stat1,
                "axis_score_alt_immune": axis_alt_immune,
            }
        )

    axis_values = [float(row["axis_score"]) for row in rows]
    median_cut = sorted(axis_values)[len(axis_values) // 2]
    tertile_cuts = quantile_groups(axis_values, 3)
    quartile_cuts = quantile_groups(axis_values, 4)
    for row in rows:
        axis = float(row["axis_score"])
        risk = float(row["risk_score"])
        immune = float(row["immune_score"])
        row["axis_group_median"] = "high" if axis >= median_cut else "low"
        row["axis_group_tertile"] = group_by_cut(axis, tertile_cuts, ["low", "mid", "high"])
        row["axis_group_quartile"] = group_by_cut(axis, quartile_cuts, ["Q1", "Q2", "Q3", "Q4"])
        row["risk_group_median"] = "high" if risk >= sorted([float(r["risk_score"]) for r in rows])[len(rows) // 2] else "low"
        row["immune_group_median"] = (
            "high" if immune >= sorted([float(r["immune_score"]) for r in rows])[len(rows) // 2] else "low"
        )
        row["four_quadrant"] = f"risk_{row['risk_group_median']}_immune_{row['immune_group_median']}"

    write_csv(
        PROCESSED / "formal_modeling_dataset.csv",
        rows,
        [
            "sampleID",
            "OS_time_days",
            "OS_event",
            "age",
            "gender",
            "stage_group",
            "risk_score",
            "immune_score",
            "axis_score",
            "axis_score_no_STAT1",
            "axis_score_alt_immune",
            "axis_group_median",
            "axis_group_tertile",
            "axis_group_quartile",
            "risk_group_median",
            "immune_group_median",
            "four_quadrant",
        ],
    )
    print(f"Wrote {PROCESSED / 'formal_modeling_dataset.csv'} with {len(rows)} rows")


if __name__ == "__main__":
    main()
