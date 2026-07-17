#!/usr/bin/env python3
"""Median-split Kaplan-Meier/log-rank feasibility screen using stdlib only."""

from __future__ import annotations

import csv
import math
from datetime import date
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"
TODAY = date.today().isoformat()

CLINICAL = PROCESSED / "clinical_os_derived.csv"
EXPR = PROCESSED / "luad_tme_immune_gene_expression_subset.csv"


def parse_float(value: object) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def normal_sf_from_chisq1(stat: float) -> float:
    """Survival function for chi-square df=1 via erfc(sqrt(x/2))."""
    if stat < 0:
        return 1.0
    return math.erfc(math.sqrt(stat / 2.0))


def logrank(group_a: list[tuple[float, int]], group_b: list[tuple[float, int]]) -> tuple[float, float]:
    """Return chi-square statistic and approximate p value for two groups."""
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
    return stat, normal_sf_from_chisq1(stat)


def km_survival_at(records: list[tuple[float, int]], day: float) -> float | None:
    if not records:
        return None
    survival = 1.0
    for event_time in sorted({time for time, event in records if event == 1 and time <= day}):
        at_risk = sum(1 for time, _ in records if time >= event_time)
        events = sum(1 for time, event in records if time == event_time and event == 1)
        if at_risk > 0:
            survival *= 1.0 - events / at_risk
    return survival


def load_clinical() -> dict[str, tuple[float, int, str]]:
    rows = {}
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
            rows[sample] = (time, int(event), row.get("stage_group", ""))
    return rows


def load_expression() -> dict[str, dict[str, float]]:
    with EXPR.open(encoding="utf-8", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        samples = header[1:]
        genes = {}
        for row in reader:
            gene = row[0]
            genes[gene] = {
                sample: float(value)
                for sample, value in zip(samples, row[1:])
                if value not in ("", "NA", "nan")
            }
    return genes


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def screen() -> list[dict[str, object]]:
    clinical = load_clinical()
    genes = load_expression()
    rows = []
    for gene, values in sorted(genes.items()):
        usable = [(sample, value, clinical[sample]) for sample, value in values.items() if sample in clinical]
        if len(usable) < 50:
            continue
        cutoff = median([value for _, value, _ in usable])
        high = [(time, event) for _, value, (time, event, _) in usable if value >= cutoff]
        low = [(time, event) for _, value, (time, event, _) in usable if value < cutoff]
        high_events = sum(event for _, event in high)
        low_events = sum(event for _, event in low)
        stat, p_value = logrank(high, low)
        high_s3 = km_survival_at(high, 365.25 * 3)
        low_s3 = km_survival_at(low, 365.25 * 3)
        high_s5 = km_survival_at(high, 365.25 * 5)
        low_s5 = km_survival_at(low, 365.25 * 5)
        rows.append(
            {
                "gene": gene,
                "n": len(usable),
                "median_cutoff": round(cutoff, 4),
                "high_n": len(high),
                "low_n": len(low),
                "high_events": high_events,
                "low_events": low_events,
                "logrank_chisq": round(stat, 4),
                "logrank_p_approx": f"{p_value:.6g}",
                "high_3y_survival": "" if high_s3 is None else round(high_s3, 4),
                "low_3y_survival": "" if low_s3 is None else round(low_s3, 4),
                "high_5y_survival": "" if high_s5 is None else round(high_s5, 4),
                "low_5y_survival": "" if low_s5 is None else round(low_s5, 4),
                "direction": direction(high_s3, low_s3, high_s5, low_s5),
            }
        )
    rows.sort(key=lambda row: float(row["logrank_p_approx"]))
    return rows


def direction(high_s3: float | None, low_s3: float | None, high_s5: float | None, low_s5: float | None) -> str:
    diffs = []
    if high_s3 is not None and low_s3 is not None:
        diffs.append(high_s3 - low_s3)
    if high_s5 is not None and low_s5 is not None:
        diffs.append(high_s5 - low_s5)
    if not diffs:
        return "unknown"
    avg = sum(diffs) / len(diffs)
    if avg > 0.02:
        return "high expression better survival"
    if avg < -0.02:
        return "high expression worse survival"
    return "small survival difference"


def write_report(rows: list[dict[str, object]]) -> None:
    top = rows[:12]
    lines = [
        "# LUAD TME/immune KM Log-rank 初筛",
        "",
        f"- 生成日期：{TODAY}",
        "- 数据：本地 TCGA Xena LUAD 表达矩阵 + 临床 OS 派生表 + 28 个 TME/immune 候选基因。",
        "- 边界：本轮未下载数据、未安装新包；使用标准库实现中位数分组 log-rank 初筛，p 值为近似值，只用于决定下一步是否值得正式建模。",
        "",
        "## 初筛排名",
        "",
        "| 排名 | 基因 | n | log-rank p approx | 方向 | high 3y | low 3y | high 5y | low 5y |",
        "|---:|---|---:|---:|---|---:|---:|---:|---:|",
    ]
    for index, row in enumerate(top, 1):
        lines.append(
            f"| {index} | {row['gene']} | {row['n']} | {row['logrank_p_approx']} | "
            f"{row['direction']} | {row['high_3y_survival']} | {row['low_3y_survival']} | "
            f"{row['high_5y_survival']} | {row['low_5y_survival']} |"
        )
    lines.extend(
        [
            "",
            "## 解释",
            "",
            "- `high` / `low` 按每个基因表达中位数切分。",
            "- `direction` 基于 3 年和 5 年 KM 生存率差的平均方向，不是 hazard ratio。",
            "- 正式论文阶段仍需预先确认统计计划，并用 Cox/KM、校正协变量、FDR 或外部验证来确认。",
            "",
            "## 下一步建议",
            "",
            "1. 不下载数据也可以继续：对 top genes 做正式 KM 图数据、单因素 Cox 近似/或改用 R/Python 统计包。",
            "2. 如果进入外部验证，需要筛 GEO 队列并下载，预计几十 MB。",
            "3. 如果要全量 RNA-seq 原始 STAR-counts，GDC manifest 估计约 2.37 GB，当前仍不建议。",
            "",
            "## 输出文件",
            "",
            "- `data/processed/km_logrank_gene_screen.csv`",
        ]
    )
    (REPORTS / "luad_km_logrank_screen.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    rows = screen()
    write_csv(
        PROCESSED / "km_logrank_gene_screen.csv",
        rows,
        [
            "gene",
            "n",
            "median_cutoff",
            "high_n",
            "low_n",
            "high_events",
            "low_events",
            "logrank_chisq",
            "logrank_p_approx",
            "high_3y_survival",
            "low_3y_survival",
            "high_5y_survival",
            "low_5y_survival",
            "direction",
        ],
    )
    write_report(rows)
    print(f"Wrote report: {REPORTS / 'luad_km_logrank_screen.md'}")


if __name__ == "__main__":
    main()
