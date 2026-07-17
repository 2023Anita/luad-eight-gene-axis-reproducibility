#!/usr/bin/env python3
"""Local exploratory summaries for LUAD TME/immune starter data."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path
from statistics import mean, median


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"
TODAY = date.today().isoformat()

CLINICAL_PATH = RAW / "TCGA.LUAD.clinicalMatrix.tsv"
EXPR_PATH = PROCESSED / "luad_tme_immune_gene_expression_subset.csv"

KEY_CLINICAL_FIELDS = [
    "sampleID",
    "_PATIENT",
    "sample_type",
    "pathologic_stage",
    "pathologic_T",
    "pathologic_N",
    "pathologic_M",
    "vital_status",
    "days_to_death",
    "days_to_last_followup",
    "age_at_initial_pathologic_diagnosis",
    "gender",
    "tobacco_smoking_history",
]


def parse_float(value: str) -> float | None:
    value = (value or "").strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def stage_group(stage: str) -> str:
    stage = (stage or "").strip().upper()
    if "STAGE IV" in stage:
        return "IV"
    if "STAGE III" in stage:
        return "III"
    if "STAGE II" in stage:
        return "II"
    if "STAGE I" in stage:
        return "I"
    return ""


def derive_os(row: dict[str, str]) -> tuple[int | None, float | None]:
    vital = (row.get("vital_status") or "").strip().upper()
    death_days = parse_float(row.get("days_to_death", ""))
    followup_days = parse_float(row.get("days_to_last_followup", ""))
    if vital == "DECEASED":
        if death_days is not None and death_days > 0:
            return 1, death_days
        return 1, followup_days
    if vital == "LIVING":
        return 0, followup_days
    return None, None


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def load_clinical() -> tuple[list[dict[str, str]], list[str]]:
    with CLINICAL_PATH.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader), list(reader.fieldnames or [])


def load_expression() -> tuple[list[str], dict[str, dict[str, float]]]:
    with EXPR_PATH.open(encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        samples = header[1:]
        genes: dict[str, dict[str, float]] = {}
        for row in reader:
            gene = row[0]
            genes[gene] = {
                sample: float(value)
                for sample, value in zip(samples, row[1:])
                if value not in ("", "NA", "nan")
            }
    return samples, genes


def missingness(rows: list[dict[str, str]], fields: list[str]) -> list[dict[str, object]]:
    out = []
    total = len(rows)
    for field in fields:
        missing = sum(1 for row in rows if not (row.get(field) or "").strip())
        out.append(
            {
                "field": field,
                "rows": total,
                "missing": missing,
                "missing_pct": round(missing / total * 100, 2) if total else "",
            }
        )
    return out


def summarize(values: list[float]) -> dict[str, object]:
    if not values:
        return {"n": 0, "mean": "", "median": "", "sd": ""}
    avg = mean(values)
    sd = math.sqrt(sum((x - avg) ** 2 for x in values) / (len(values) - 1)) if len(values) > 1 else 0.0
    return {"n": len(values), "mean": round(avg, 4), "median": round(median(values), 4), "sd": round(sd, 4)}


def make_clinical_derived(clinical_rows: list[dict[str, str]]) -> list[dict[str, object]]:
    rows = []
    for row in clinical_rows:
        os_event, os_time = derive_os(row)
        rows.append(
            {
                "sampleID": row.get("sampleID", ""),
                "patient_id": row.get("_PATIENT", ""),
                "sample_type": row.get("sample_type", ""),
                "stage_group": stage_group(row.get("pathologic_stage", "")),
                "pathologic_stage": row.get("pathologic_stage", ""),
                "pathologic_T": row.get("pathologic_T", ""),
                "pathologic_N": row.get("pathologic_N", ""),
                "pathologic_M": row.get("pathologic_M", ""),
                "vital_status": row.get("vital_status", ""),
                "OS_event": "" if os_event is None else os_event,
                "OS_time_days": "" if os_time is None else os_time,
                "age": row.get("age_at_initial_pathologic_diagnosis", ""),
                "gender": row.get("gender", ""),
                "tobacco_smoking_history": row.get("tobacco_smoking_history", ""),
            }
        )
    return rows


def gene_summary(genes: dict[str, dict[str, float]], samples: set[str]) -> list[dict[str, object]]:
    rows = []
    for gene, sample_values in sorted(genes.items()):
        values = [value for sample, value in sample_values.items() if sample in samples]
        row = {"gene": gene}
        row.update(summarize(values))
        rows.append(row)
    return rows


def grouped_gene_summary(
    genes: dict[str, dict[str, float]],
    clinical_by_sample: dict[str, dict[str, object]],
    group_field: str,
    allowed_samples: set[str],
) -> list[dict[str, object]]:
    grouped_samples: dict[str, list[str]] = defaultdict(list)
    for sample in allowed_samples:
        group = str(clinical_by_sample.get(sample, {}).get(group_field, ""))
        if group != "":
            grouped_samples[group].append(sample)
    rows = []
    for gene, sample_values in sorted(genes.items()):
        for group, group_samples in sorted(grouped_samples.items()):
            values = [sample_values[sample] for sample in group_samples if sample in sample_values]
            row = {"gene": gene, group_field: group}
            row.update(summarize(values))
            rows.append(row)
    return rows


def exploratory_effects(
    genes: dict[str, dict[str, float]],
    clinical_by_sample: dict[str, dict[str, object]],
    primary_samples: set[str],
) -> list[dict[str, object]]:
    rows = []
    for gene, sample_values in sorted(genes.items()):
        living = []
        deceased = []
        early = []
        advanced = []
        for sample in primary_samples:
            value = sample_values.get(sample)
            if value is None:
                continue
            clinical = clinical_by_sample.get(sample, {})
            if clinical.get("OS_event") == 0:
                living.append(value)
            elif clinical.get("OS_event") == 1:
                deceased.append(value)
            if clinical.get("stage_group") in ("I", "II"):
                early.append(value)
            elif clinical.get("stage_group") in ("III", "IV"):
                advanced.append(value)
        living_summary = summarize(living)
        deceased_summary = summarize(deceased)
        early_summary = summarize(early)
        advanced_summary = summarize(advanced)
        rows.append(
            {
                "gene": gene,
                "living_n": living_summary["n"],
                "deceased_n": deceased_summary["n"],
                "deceased_minus_living_mean": diff(deceased_summary["mean"], living_summary["mean"]),
                "early_stage_n": early_summary["n"],
                "advanced_stage_n": advanced_summary["n"],
                "advanced_minus_early_mean": diff(advanced_summary["mean"], early_summary["mean"]),
            }
        )
    rows.sort(
        key=lambda row: max(
            abs(float(row["deceased_minus_living_mean"] or 0)),
            abs(float(row["advanced_minus_early_mean"] or 0)),
        ),
        reverse=True,
    )
    return rows


def diff(left: object, right: object) -> object:
    if left == "" or right == "":
        return ""
    return round(float(left) - float(right), 4)


def write_report(
    clinical_rows: list[dict[str, object]],
    samples: list[str],
    primary_samples: set[str],
    matched_primary: set[str],
    effects: list[dict[str, object]],
) -> None:
    valid_os = [row for row in clinical_rows if row["sampleID"] in matched_primary and row["OS_event"] != "" and row["OS_time_days"] != ""]
    stage_counts: dict[str, int] = defaultdict(int)
    event_counts: dict[str, int] = defaultdict(int)
    for row in clinical_rows:
        if row["sampleID"] not in matched_primary:
            continue
        stage_counts[str(row["stage_group"] or "missing")] += 1
        event_counts[str(row["OS_event"] if row["OS_event"] != "" else "missing")] += 1
    top_event = sorted(effects, key=lambda row: abs(float(row["deceased_minus_living_mean"] or 0)), reverse=True)[:8]
    top_stage = sorted(effects, key=lambda row: abs(float(row["advanced_minus_early_mean"] or 0)), reverse=True)[:8]
    lines = [
        "# LUAD TME/immune 临床探索摘要",
        "",
        f"- 生成日期：{TODAY}",
        "- 数据来源：已下载的 TCGA Xena LUAD 表达矩阵、临床矩阵，以及本地 28 个 TME/immune 候选基因子集。",
        "- 边界：本报告为探索性描述，不包含正式统计检验、模型训练或 manuscript-ready 结论；本轮未额外下载任何数据。",
        "",
        "## 样本匹配",
        "",
        f"- 表达矩阵样本列：{len(samples)}",
        f"- 表达矩阵 primary tumor 样本列：{len(primary_samples)}",
        f"- 临床可匹配 primary tumor 样本：{len(matched_primary)}",
        f"- 可构建 OS_time / OS_event 的 matched primary tumor 样本：{len(valid_os)}",
        f"- stage_group 分布：{json.dumps(dict(sorted(stage_counts.items())), ensure_ascii=False)}",
        f"- OS_event 分布：{json.dumps(dict(sorted(event_counts.items())), ensure_ascii=False)}",
        "",
        "## 候选信号：按生存状态的均值差",
        "",
        "| 基因 | living_n | deceased_n | deceased - living mean |",
        "|---|---:|---:|---:|",
    ]
    for row in top_event:
        lines.append(
            f"| {row['gene']} | {row['living_n']} | {row['deceased_n']} | {row['deceased_minus_living_mean']} |"
        )
    lines.extend(
        [
            "",
            "## 候选信号：按分期早晚的均值差",
            "",
            "| 基因 | early_n | advanced_n | advanced - early mean |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in top_stage:
        lines.append(
            f"| {row['gene']} | {row['early_stage_n']} | {row['advanced_stage_n']} | {row['advanced_minus_early_mean']} |"
        )
    lines.extend(
        [
            "",
            "## 下一步",
            "",
            "1. 若继续省数据路线：先做正式缺失率表、OS 构建规则、单因素 Cox/KM 的可行性检查。",
            "2. 若要做免疫浸润/TME 评分：需要补充算法或基因集定义，不一定需要下载大文件。",
            "3. 若要做外部验证：需要从 GEO 候选 series 中人工筛 1-2 个真正可用的 LUAD 队列，届时可能需要额外下载几十 MB。",
            "",
            "## 输出文件",
            "",
            "- `data/processed/clinical_os_derived.csv`",
            "- `data/processed/clinical_missingness.csv`",
            "- `data/processed/gene_expression_summary_primary_tumor.csv`",
            "- `data/processed/gene_by_stage_group_summary.csv`",
            "- `data/processed/gene_by_os_event_summary.csv`",
            "- `data/processed/exploratory_gene_effects.csv`",
        ]
    )
    (REPORTS / "luad_tme_clinical_exploration.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    raw_clinical, clinical_fields = load_clinical()
    samples, genes = load_expression()
    clinical_rows = make_clinical_derived(raw_clinical)
    clinical_by_sample = {str(row["sampleID"]): row for row in clinical_rows}
    expr_samples = set(samples)
    primary_samples = {sample for sample in samples if sample.endswith("-01")}
    matched_primary = primary_samples & set(clinical_by_sample)

    write_csv(
        PROCESSED / "clinical_os_derived.csv",
        clinical_rows,
        [
            "sampleID",
            "patient_id",
            "sample_type",
            "stage_group",
            "pathologic_stage",
            "pathologic_T",
            "pathologic_N",
            "pathologic_M",
            "vital_status",
            "OS_event",
            "OS_time_days",
            "age",
            "gender",
            "tobacco_smoking_history",
        ],
    )
    write_csv(
        PROCESSED / "clinical_missingness.csv",
        missingness(raw_clinical, [field for field in KEY_CLINICAL_FIELDS if field in clinical_fields]),
        ["field", "rows", "missing", "missing_pct"],
    )
    write_csv(
        PROCESSED / "gene_expression_summary_primary_tumor.csv",
        gene_summary(genes, matched_primary),
        ["gene", "n", "mean", "median", "sd"],
    )
    write_csv(
        PROCESSED / "gene_by_stage_group_summary.csv",
        grouped_gene_summary(genes, clinical_by_sample, "stage_group", matched_primary),
        ["gene", "stage_group", "n", "mean", "median", "sd"],
    )
    write_csv(
        PROCESSED / "gene_by_os_event_summary.csv",
        grouped_gene_summary(genes, clinical_by_sample, "OS_event", matched_primary),
        ["gene", "OS_event", "n", "mean", "median", "sd"],
    )
    effects = exploratory_effects(genes, clinical_by_sample, matched_primary)
    write_csv(
        PROCESSED / "exploratory_gene_effects.csv",
        effects,
        [
            "gene",
            "living_n",
            "deceased_n",
            "deceased_minus_living_mean",
            "early_stage_n",
            "advanced_stage_n",
            "advanced_minus_early_mean",
        ],
    )
    summary = {
        "expression_samples": len(samples),
        "expression_primary_tumor_samples": len(primary_samples),
        "clinical_rows": len(clinical_rows),
        "matched_primary_tumor_samples": len(matched_primary),
        "matched_any_expression_samples": len(expr_samples & set(clinical_by_sample)),
        "genes": sorted(genes),
        "extra_download_performed": False,
    }
    (PROCESSED / "luad_tme_clinical_exploration_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    write_report(clinical_rows, samples, primary_samples, matched_primary, effects)
    print(f"Wrote report: {REPORTS / 'luad_tme_clinical_exploration.md'}")


if __name__ == "__main__":
    main()
