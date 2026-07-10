#!/usr/bin/env python3
"""Render the four planned submission figures from the checked source-data CSVs."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
SOURCE = PACKAGE / "source_data"
OUT = PACKAGE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#3B6EA5"
ORANGE = "#C76B32"
TEAL = "#2E8B7A"
GRAY = "#7A7A7A"
LIGHT_GRAY = "#D9D9D9"
BLACK = "#222222"

mpl.rcParams.update(
    {
        "font.family": "Arial",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "axes.labelcolor": BLACK,
        "xtick.color": BLACK,
        "ytick.color": BLACK,
        "legend.frameon": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def rows(filename: str) -> list[dict[str, str]]:
    with (SOURCE / filename).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def value(row: dict[str, str], key: str) -> float:
    return float(row[key])


def save_figure(fig: plt.Figure, basename: str) -> None:
    base = OUT / basename
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white", transparent=False)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", transparent=False)
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight", facecolor="white", transparent=False)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white", transparent=False)
    normalize_rgb_raster(base.with_suffix(".tiff"), 600)
    normalize_rgb_raster(base.with_suffix(".png"), 300)
    plt.close(fig)


def normalize_rgb_raster(path: Path, dpi: int) -> None:
    """Flatten the white canvas so submission rasters use RGB, not RGBA."""
    with Image.open(path) as image:
        if image.mode == "RGB":
            return
        canvas = Image.new("RGB", image.size, "white")
        if "A" in image.getbands():
            canvas.paste(image, mask=image.getchannel("A"))
        else:
            canvas.paste(image)
        save_args = {"dpi": (dpi, dpi)}
        if path.suffix.lower() in {".tif", ".tiff"}:
            save_args["compression"] = "tiff_lzw"
        canvas.save(path, **save_args)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.15, 1.08, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="top")


def km_panel(ax: plt.Axes, filename: str, x_label: str, title: str) -> None:
    data = rows(filename)
    groups = [("axis_group_median=low", "Low axis", GRAY), ("axis_group_median=high", "High axis", ORANGE)]
    for code, label, color in groups:
        subset = [item for item in data if item["strata"] == code]
        if not subset:
            continue
        times = [0.0] + [value(item, "time") for item in subset]
        surv = [1.0] + [value(item, "survival") for item in subset]
        ax.step(times, surv, where="post", lw=1.6, color=color, label=label)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel(x_label)
    ax.set_ylabel("Survival probability")
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold")
    ax.legend(loc="lower left", fontsize=6)


def forest_panel(ax: plt.Axes, entries: list[tuple[str, float, float, float, str]], title: str) -> None:
    y_positions = list(range(len(entries), 0, -1))
    for y, (label, hr, low, high, color) in zip(y_positions, entries):
        ax.errorbar(hr, y, xerr=[[hr - low], [high - hr]], fmt="o", color=color, ms=4, capsize=2, lw=1)
        ax.text(0.92, y, label, ha="right", va="center", fontsize=6)
        ax.text(0.98, y, f"{hr:.2f} ({low:.2f}-{high:.2f})", transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6)
    ax.axvline(1, color=GRAY, lw=0.8, ls="--")
    max_x = max(high for _, _, _, high, _ in entries)
    ax.set_xlim(0.2, max(3.8, max_x * 1.35))
    ax.set_ylim(0.4, len(entries) + 0.6)
    ax.set_yticks([])
    ax.set_xlabel("Hazard ratio")
    ax.set_title(title, loc="left", fontsize=8, fontweight="bold")


def get_cox(filename: str, model: str, term: str) -> dict[str, str]:
    return next(item for item in rows(filename) if item["model"] == model and item["term"] == term)


def figure_1() -> None:
    data = rows("figure1_cohort_hierarchy.csv")
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.text(0.2, 6.65, "Fixed eight-gene axis and evidence hierarchy", fontsize=10, fontweight="bold")
    ax.text(0.2, 6.25, "risk = mean z(CA9, MKI67, STAT1, HIF1A)   minus   immune = mean z(CD4, CD8A, CD8B, TIGIT)", fontsize=7)
    positions = [(0.4, 4.5), (3.7, 4.5), (6.9, 4.5), (3.7, 1.7), (6.9, 1.7)]
    styles = {"discovery": (BLUE, "Discovery OS"), "primary_external": (TEAL, "Primary external OS"), "secondary_external": (ORANGE, "Secondary RFS support"), "sensitivity": (GRAY, "Sensitivity only"), "annotation": (GRAY, "Annotation only")}
    for item, (x, y) in zip(data, positions):
        edge, header = styles[item["role"]]
        box = FancyBboxPatch((x, y), 2.6, 1.1, boxstyle="round,pad=0.03,rounding_size=0.05", facecolor="white", edgecolor=edge, lw=1.5)
        ax.add_patch(box)
        ax.text(x + 0.12, y + 0.82, header, fontsize=6, color=edge, fontweight="bold")
        ax.text(x + 0.12, y + 0.55, item["cohort"], fontsize=7, fontweight="bold")
        detail = "Annotation resource" if item["role"] == "annotation" else f"{item['endpoint']}; n={item['eligible_n']}; {item['events_or_annotation']}"
        ax.text(x + 0.12, y + 0.28, detail, fontsize=5.7)
    for start, end in [((3.0, 5.05), (3.7, 5.05)), ((6.3, 5.05), (6.9, 5.05)), ((2.0, 4.5), (4.8, 2.8)), ((6.3, 4.5), (8.1, 2.8))]:
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "color": GRAY, "lw": 0.9})
    ax.text(0.4, 0.5, "Exact-score cohorts use the same genes, signs, and equal weights. OS and RFS are reported separately.", fontsize=6.5)
    save_figure(fig, "Figure_1_cohort_hierarchy_axis")


def figure_2() -> None:
    stage = rows("figure2_tcga_stage_summary.csv")
    cox = rows("figure2_tcga_cox_results.csv")
    cont = get_cox("figure2_tcga_cox_results.csv", "cox_axis_continuous", "axis_score")
    adj = get_cox("figure2_tcga_cox_results.csv", "cox_axis_adjusted_age_gender_stage", "axis_score")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), gridspec_kw={"width_ratios": [1.0, 1.15, 1.15]})
    ax = axes[0]
    ax.bar([item["stage_group"] for item in stage], [value(item, "mean") for item in stage], color=BLUE, width=0.65)
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.set_ylabel("Mean axis score")
    ax.set_xlabel("Pathological stage")
    ax.set_title("Stage gradient", loc="left", fontsize=8, fontweight="bold")
    panel_label(ax, "a")
    km_panel(axes[1], "figure2_tcga_km_curve.csv", "Overall survival time (days)", "Discovery OS visualization")
    panel_label(axes[1], "b")
    forest_panel(
        axes[2],
        [
            ("Continuous axis", value(cont, "HR"), value(cont, "CI95_low"), value(cont, "CI95_high"), BLUE),
            ("Age/sex/stage adjusted", value(adj, "HR"), value(adj, "CI95_low"), value(adj, "CI95_high"), TEAL),
        ],
        "TCGA-LUAD Cox models",
    )
    panel_label(axes[2], "c")
    fig.tight_layout(w_pad=1.4)
    save_figure(fig, "Figure_2_tcga_discovery")


def figure_3() -> None:
    cont = get_cox("figure3_gse72094_cox_results.csv", "cox_axis_continuous", "axis_score")
    adj = get_cox("figure3_gse72094_cox_results.csv", "cox_axis_adjusted_age_gender_stage", "axis_score")
    risk = get_cox("figure3_gse72094_cox_results.csv", "cox_risk_immune", "risk_score")
    immune = get_cox("figure3_gse72094_cox_results.csv", "cox_risk_immune", "immune_score")
    cindex = rows("figure3_gse72094_cindex.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), gridspec_kw={"width_ratios": [1.1, 1.25, 0.85]})
    km_panel(axes[0], "figure3_gse72094_km_curve.csv", "Overall survival time (days)", "Exact-score OS replication")
    panel_label(axes[0], "a")
    forest_panel(
        axes[1],
        [
            ("Continuous axis", value(cont, "HR"), value(cont, "CI95_low"), value(cont, "CI95_high"), BLUE),
            ("Age/sex/stage adjusted", value(adj, "HR"), value(adj, "CI95_low"), value(adj, "CI95_high"), TEAL),
            ("Risk component", value(risk, "HR"), value(risk, "CI95_low"), value(risk, "CI95_high"), ORANGE),
            ("Immune component", value(immune, "HR"), value(immune, "CI95_low"), value(immune, "CI95_high"), GRAY),
        ],
        "GSE72094 Cox estimates",
    )
    panel_label(axes[1], "b")
    labels = ["Clinical", "Clinical + axis"]
    values = [value(item, "optimism_adjusted_cindex") for item in cindex]
    axes[2].bar(labels, values, color=[GRAY, TEAL], width=0.6)
    axes[2].set_ylim(0.5, 0.75)
    axes[2].set_ylabel("Optimism-adjusted C-index")
    axes[2].set_title("Model-fit increment", loc="left", fontsize=8, fontweight="bold")
    for i, score in enumerate(values):
        axes[2].text(i, score + 0.007, f"{score:.3f}", ha="center", fontsize=6)
    axes[2].text(0.5, 0.52, "LRT p=1.14e-04", ha="center", fontsize=6)
    panel_label(axes[2], "c")
    fig.tight_layout(w_pad=1.4)
    save_figure(fig, "Figure_3_gse72094_external_os")


def figure_4() -> None:
    cont = get_cox("figure4_gse31210_cox_results.csv", "cox_axis_continuous", "axis_score")
    adj = get_cox("figure4_gse31210_cox_results.csv", "cox_axis_adjusted_age_gender_stage", "axis_score")
    strat = get_cox("figure4_gse31210_gender_stratified.csv", "cox_axis_age_stage_stratified_by_gender", "axis_score")
    cindex = rows("figure4_gse31210_cindex.csv")
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.45), gridspec_kw={"width_ratios": [1.1, 1.25, 0.85]})
    km_panel(axes[0], "figure4_gse31210_km_curve.csv", "Relapse-free time (months)", "Exact-score RFS support")
    panel_label(axes[0], "a")
    forest_panel(
        axes[1],
        [
            ("Continuous axis", value(cont, "HR"), value(cont, "CI95_low"), value(cont, "CI95_high"), BLUE),
            ("Age/sex/stage adjusted", value(adj, "HR"), value(adj, "CI95_low"), value(adj, "CI95_high"), ORANGE),
            ("Gender-stratified sensitivity", value(strat, "HR"), value(strat, "CI95_low"), value(strat, "CI95_high"), TEAL),
        ],
        "GSE31210 Cox estimates",
    )
    panel_label(axes[1], "b")
    labels = ["Clinical", "Clinical + axis"]
    values = [value(item, "optimism_adjusted_cindex") for item in cindex]
    axes[2].bar(labels, values, color=[GRAY, TEAL], width=0.6)
    axes[2].set_ylim(0.5, 0.75)
    axes[2].set_ylabel("Optimism-adjusted C-index")
    axes[2].set_title("Exploratory increment", loc="left", fontsize=8, fontweight="bold")
    for i, score in enumerate(values):
        axes[2].text(i, score + 0.007, f"{score:.3f}", ha="center", fontsize=6)
    axes[2].text(0.5, 0.52, "LRT p=0.0458", ha="center", fontsize=6)
    panel_label(axes[2], "c")
    fig.tight_layout(w_pad=1.4)
    save_figure(fig, "Figure_4_gse31210_external_rfs")


def main() -> None:
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    print("Wrote four submission figure bundles")


if __name__ == "__main__":
    main()
