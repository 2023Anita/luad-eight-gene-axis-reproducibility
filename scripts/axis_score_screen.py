#!/usr/bin/env python3
"""Screen a hypoxia/proliferation versus immune-activation axis."""

from __future__ import annotations

import csv
import math
from datetime import date
from pathlib import Path
from statistics import mean, median


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"
TODAY = date.today().isoformat()

CLINICAL = PROCESSED / "clinical_os_derived.csv"
EXPR = PROCESSED / "luad_tme_immune_gene_expression_subset.csv"

RISK_GENES = ["CA9", "MKI67", "STAT1", "HIF1A"]
IMMUNE_GENES = ["CD4", "CD8A", "CD8B", "TIGIT"]


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
            out[sample] = {
                "OS_time_days": time,
                "OS_event": int(event),
                "stage_group": row.get("stage_group", ""),
                "age": row.get("age", ""),
                "gender": row.get("gender", ""),
            }
    return out


def summarize(values: list[float]) -> dict[str, object]:
    if not values:
        return {"n": 0, "mean": "", "median": ""}
    return {"n": len(values), "mean": round(mean(values), 4), "median": round(median(values), 4)}


def normal_sf_chisq1(stat: float) -> float:
    return math.erfc(math.sqrt(max(stat, 0.0) / 2.0))


def logrank(group_a: list[tuple[float, int]], group_b: list[tuple[float, int]]) -> tuple[float, float]:
    event_times = sorted({time for time, event in group_a + group_b if event == 1})
    observed_a = expected_a = variance_a = 0.0
    for event_time in event_times:
        risk_a = sum(1 for time, _ in group_a if time >= event_time)
        risk_b = sum(1 for time, _ in group_b if time >= event_time)
        events_a = sum(1 for time, event in group_a if time == event_time and event == 1)
        events_b = sum(1 for time, event in group_b if time == event_time and event == 1)
        risk_total = risk_a + risk_b
        events_total = events_a + events_b
        if risk_total <= 1 or events_total == 0:
            continue
        observed_a += events_a
        expected_a += events_total * risk_a / risk_total
        variance_a += (
            risk_a
            * risk_b
            * events_total
            * (risk_total - events_total)
            / (risk_total * risk_total * (risk_total - 1))
        )
    if variance_a <= 0:
        return 0.0, 1.0
    stat = (observed_a - expected_a) ** 2 / variance_a
    return stat, normal_sf_chisq1(stat)


def km_survival_at(records: list[tuple[float, int]], day: float) -> float:
    survival = 1.0
    for event_time in sorted({time for time, event in records if event == 1 and time <= day}):
        at_risk = sum(1 for time, _ in records if time >= event_time)
        events = sum(1 for time, event in records if time == event_time and event == 1)
        if at_risk > 0:
            survival *= 1.0 - events / at_risk
    return survival


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    expr = load_expression()
    clinical = load_clinical()
    z = {gene: zscores(expr[gene]) for gene in RISK_GENES + IMMUNE_GENES}
    samples = sorted(set(clinical).intersection(*(set(z[gene]) for gene in RISK_GENES + IMMUNE_GENES)))
    rows = []
    for sample in samples:
        risk_score = sum(z[gene][sample] for gene in RISK_GENES) / len(RISK_GENES)
        immune_score = sum(z[gene][sample] for gene in IMMUNE_GENES) / len(IMMUNE_GENES)
        axis_score = risk_score - immune_score
        rows.append(
            {
                "sampleID": sample,
                "risk_score": round(risk_score, 6),
                "immune_score": round(immune_score, 6),
                "axis_score": round(axis_score, 6),
                "OS_time_days": clinical[sample]["OS_time_days"],
                "OS_event": clinical[sample]["OS_event"],
                "stage_group": clinical[sample]["stage_group"],
                "age": clinical[sample]["age"],
                "gender": clinical[sample]["gender"],
            }
        )
    cutoff = median([float(row["axis_score"]) for row in rows])
    for row in rows:
        row["axis_group"] = "high" if float(row["axis_score"]) >= cutoff else "low"
    write_csv(
        PROCESSED / "luad_axis_scores.csv",
        rows,
        [
            "sampleID",
            "risk_score",
            "immune_score",
            "axis_score",
            "axis_group",
            "OS_time_days",
            "OS_event",
            "stage_group",
            "age",
            "gender",
        ],
    )

    high = [(float(r["OS_time_days"]), int(r["OS_event"])) for r in rows if r["axis_group"] == "high"]
    low = [(float(r["OS_time_days"]), int(r["OS_event"])) for r in rows if r["axis_group"] == "low"]
    stat, p_value = logrank(high, low)
    stage_rows = []
    for group in ["I", "II", "III", "IV"]:
        values = [float(row["axis_score"]) for row in rows if row["stage_group"] == group]
        item = {"stage_group": group}
        item.update(summarize(values))
        stage_rows.append(item)
    write_csv(PROCESSED / "luad_axis_score_by_stage.csv", stage_rows, ["stage_group", "n", "mean", "median"])
    event_rows = []
    for event in [0, 1]:
        values = [float(row["axis_score"]) for row in rows if int(row["OS_event"]) == event]
        item = {"OS_event": event}
        item.update(summarize(values))
        event_rows.append(item)
    write_csv(PROCESSED / "luad_axis_score_by_os_event.csv", event_rows, ["OS_event", "n", "mean", "median"])

    lines = [
        "# LUAD 最佳路线轴评分验证",
        "",
        f"- 生成日期：{TODAY}",
        "- 路线：hypoxia/proliferation risk axis 对比 immune activation axis。",
        "- 边界：本轮未下载数据；这是路线验证，不是正式模型。",
        "",
        "## 轴定义",
        "",
        f"- risk score = mean z-score({', '.join(RISK_GENES)})",
        f"- immune score = mean z-score({', '.join(IMMUNE_GENES)})",
        "- axis score = risk score - immune score",
        "",
        "## 生存初筛",
        "",
        f"- 可用样本：{len(rows)}",
        f"- axis score 中位数 cutoff：{cutoff:.4f}",
        f"- high axis events：{sum(e for _, e in high)} / {len(high)}",
        f"- low axis events：{sum(e for _, e in low)} / {len(low)}",
        f"- log-rank p approx：{p_value:.6g}",
        f"- high 3-year survival：{km_survival_at(high, 365.25 * 3):.4f}",
        f"- low 3-year survival：{km_survival_at(low, 365.25 * 3):.4f}",
        f"- high 5-year survival：{km_survival_at(high, 365.25 * 5):.4f}",
        f"- low 5-year survival：{km_survival_at(low, 365.25 * 5):.4f}",
        "",
        "## 分期趋势",
        "",
        "| stage | n | axis mean | axis median |",
        "|---|---:|---:|---:|",
    ]
    for item in stage_rows:
        lines.append(f"| {item['stage_group']} | {item['n']} | {item['mean']} | {item['median']} |")
    lines.extend(
        [
            "",
            "## 结论",
            "",
            "- 这条路线具备热点性和内部数据可行性。",
            "- 若 axis score 在正式 Cox/KM 中继续稳定，下一步最值得补的是 GEO 外部验证，而不是下载 GDC 全量 STAR-counts。",
            "- 如果 axis score 不稳定，则退回到单基因 `MKI67/CA9` + 免疫背景解释，不建议做复杂模型。",
        ]
    )
    (REPORTS / "luad_axis_score_validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote report: {REPORTS / 'luad_axis_score_validation.md'}")


if __name__ == "__main__":
    main()
