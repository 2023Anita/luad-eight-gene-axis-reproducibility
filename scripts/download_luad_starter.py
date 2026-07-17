#!/usr/bin/env python3
"""Download a small TCGA-LUAD starter pack for topic narrowing."""

from __future__ import annotations

import csv
import gzip
import json
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
REPORTS = ROOT / "reports"

TODAY = date.today().isoformat()

URLS = {
    "xena_luad_hiseqv2": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
    "xena_luad_clinical": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FLUAD_clinicalMatrix",
}

TME_IMMUNE_GENES = [
    "CD274",
    "PDCD1",
    "CTLA4",
    "LAG3",
    "TIGIT",
    "HAVCR2",
    "PDCD1LG2",
    "CD8A",
    "CD8B",
    "CD4",
    "FOXP3",
    "GZMB",
    "PRF1",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "IFNG",
    "STAT1",
    "IDO1",
    "MKI67",
    "VEGFA",
    "HIF1A",
    "CA9",
    "COL1A1",
    "COL1A2",
    "FN1",
    "VIM",
    "EPCAM",
]


def request_json(url: str, timeout: int = 60) -> object:
    req = urllib.request.Request(url, headers={"User-Agent": "codex-luad-starter/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def download(url: str, path: Path, timeout: int = 180) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "codex-luad-starter/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        path.write_bytes(resp.read())


def gdc_luad_star_counts_manifest() -> list[dict[str, object]]:
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-LUAD"]}},
            {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "access", "value": ["open"]}},
        ],
    }
    fields = ",".join(
        [
            "file_id",
            "file_name",
            "file_size",
            "data_type",
            "data_category",
            "experimental_strategy",
            "analysis.workflow_type",
            "cases.submitter_id",
            "cases.samples.sample_type",
        ]
    )
    params = urllib.parse.urlencode(
        {
            "size": 2000,
            "format": "JSON",
            "fields": fields,
            "filters": json.dumps(filters),
        }
    )
    data = request_json(f"https://api.gdc.cancer.gov/files?{params}", timeout=90)
    (RAW / "gdc_luad_open_expression_files.json").write_text(
        json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    rows = []
    for item in data["data"]["hits"]:
        case = (item.get("cases") or [{}])[0]
        sample = ((case.get("samples") or [{}])[0]).get("sample_type", "")
        workflow = (item.get("analysis") or {}).get("workflow_type", "")
        rows.append(
            {
                "file_id": item.get("file_id", ""),
                "file_name": item.get("file_name", ""),
                "file_size": item.get("file_size", 0),
                "workflow_type": workflow,
                "case_submitter_id": case.get("submitter_id", ""),
                "sample_type": sample,
            }
        )
    rows.sort(key=lambda r: (r["sample_type"], r["case_submitter_id"], r["file_name"]))
    write_csv(
        PROCESSED / "gdc_luad_open_expression_manifest.csv",
        rows,
        ["file_id", "file_name", "file_size", "workflow_type", "case_submitter_id", "sample_type"],
    )
    return rows


def cbioportal_luad_studies() -> list[dict[str, object]]:
    url = "https://www.cbioportal.org/api/studies?projection=DETAILED&pageSize=100000&pageNumber=0"
    studies = request_json(url, timeout=90)
    rows = []
    for study in studies:
        if not study.get("publicStudy"):
            continue
        if study.get("cancerTypeId") != "luad":
            continue
        rows.append(
            {
                "study_id": study.get("studyId", ""),
                "name": study.get("name", ""),
                "sample_count": study.get("allSampleCount", 0),
                "sequenced_sample_count": study.get("sequencedSampleCount", 0),
                "rna_seq_sample_count": (study.get("mrnaRnaSeqSampleCount") or 0)
                + (study.get("mrnaRnaSeqV2SampleCount") or 0),
                "cna_sample_count": study.get("cnaSampleCount", 0),
                "pmid": study.get("pmid", ""),
                "citation": study.get("citation", ""),
            }
        )
    rows.sort(key=lambda r: int(r["sample_count"] or 0), reverse=True)
    write_csv(
        PROCESSED / "cbioportal_luad_public_studies.csv",
        rows,
        [
            "study_id",
            "name",
            "sample_count",
            "sequenced_sample_count",
            "rna_seq_sample_count",
            "cna_sample_count",
            "pmid",
            "citation",
        ],
    )
    return rows


def geo_luad_candidates() -> list[dict[str, object]]:
    term = '("lung adenocarcinoma"[All Fields]) AND (expression profiling by high throughput sequencing[DataSet Type] OR "Expression profiling by array"[DataSet Type])'
    params = urllib.parse.urlencode({"db": "gds", "term": term, "retmode": "json", "retmax": 20})
    search = request_json(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}", timeout=60)
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    summary_params = urllib.parse.urlencode({"db": "gds", "id": ",".join(ids), "retmode": "json"})
    summary = request_json(
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{summary_params}", timeout=60
    )
    (RAW / "geo_luad_esummary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    rows = []
    result = summary.get("result", {})
    for uid in result.get("uids", []):
        item = result.get(uid, {})
        rows.append(
            {
                "uid": uid,
                "accession": item.get("accession", ""),
                "title": item.get("title", ""),
                "taxon": item.get("taxon", ""),
                "n_samples": item.get("n_samples", ""),
                "gds_type": item.get("gdsType", ""),
                "entry_type": item.get("entryType", ""),
            }
        )
    write_csv(
        PROCESSED / "geo_luad_candidate_series.csv",
        rows,
        ["uid", "accession", "title", "taxon", "n_samples", "gds_type", "entry_type"],
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def extract_gene_subset(expression_gz: Path) -> dict[str, object]:
    output = PROCESSED / "luad_tme_immune_gene_expression_subset.csv"
    requested = set(TME_IMMUNE_GENES)
    found = []
    samples = 0
    with gzip.open(expression_gz, "rt", encoding="utf-8", newline="") as src, output.open(
        "w", newline="", encoding="utf-8"
    ) as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst)
        header = next(reader)
        samples = max(len(header) - 1, 0)
        writer.writerow(header)
        for row in reader:
            if not row:
                continue
            gene = row[0].split("|", 1)[0]
            if gene in requested:
                writer.writerow(row)
                found.append(gene)
    missing = sorted(requested - set(found))
    summary = {
        "source": str(expression_gz),
        "output": str(output),
        "sample_columns": samples,
        "requested_gene_count": len(TME_IMMUNE_GENES),
        "found_gene_count": len(set(found)),
        "missing_genes": missing,
    }
    (PROCESSED / "luad_tme_immune_gene_subset_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return summary


def count_clinical_rows(path: Path) -> tuple[int, int]:
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = sum(1 for _ in reader)
    return len(header), rows


def write_report(
    manifest_rows: list[dict[str, object]],
    cbio_rows: list[dict[str, object]],
    geo_rows: list[dict[str, object]],
    subset_summary: dict[str, object],
    clinical_shape: tuple[int, int],
) -> None:
    expression_path = RAW / "TCGA.LUAD.HiSeqV2.gz"
    clinical_path = RAW / "TCGA.LUAD.clinicalMatrix.tsv"
    total_gdc_size = sum(int(row["file_size"] or 0) for row in manifest_rows)
    sample_types: dict[str, int] = {}
    for row in manifest_rows:
        sample_types[row["sample_type"]] = sample_types.get(row["sample_type"], 0) + 1
    lines = [
        "# TCGA-LUAD Starter 数据包",
        "",
        f"- 生成日期：{TODAY}",
        "- 方向：肺腺癌 LUAD + 肿瘤微环境/免疫治疗 starter 数据。",
        "- 边界：下载整理后的 TCGA Xena 表达矩阵和临床矩阵；GDC 只保存 open expression manifest，未下载 601 个 STAR-counts 原始 TSV。",
        "",
        "## 本地数据量",
        "",
        f"- 表达矩阵：`{expression_path.name}`，{expression_path.stat().st_size / 1024 / 1024:.2f} MB",
        f"- 临床矩阵：`{clinical_path.name}`，{clinical_path.stat().st_size / 1024:.1f} KB",
        f"- 整个 starter 目录：约 {directory_size_mb(ROOT):.2f} MB",
        "",
        "## 数据概况",
        "",
        f"- TCGA Xena LUAD 表达矩阵样本列数：{subset_summary['sample_columns']}",
        f"- TCGA Xena LUAD 临床矩阵：{clinical_shape[1]} 行，{clinical_shape[0]} 列",
        f"- GDC LUAD open Gene Expression Quantification 文件：{len(manifest_rows)} 个，若全量下载约 {total_gdc_size / 1024 / 1024 / 1024:.2f} GB",
        f"- GDC sample type 分布：{json.dumps(sample_types, ensure_ascii=False)}",
        f"- cBioPortal LUAD public studies：{len(cbio_rows)} 个",
        f"- GEO LUAD 候选 series：{len(geo_rows)} 个初筛条目",
        "",
        "## 已抽取候选基因表达子集",
        "",
        f"- 输出：`data/processed/luad_tme_immune_gene_expression_subset.csv`",
        f"- 已找到基因：{subset_summary['found_gene_count']} / {subset_summary['requested_gene_count']}",
        f"- 缺失基因：{', '.join(subset_summary['missing_genes']) if subset_summary['missing_genes'] else '无'}",
        "",
        "## 下一步建议",
        "",
        "1. 先用表达子集 + 临床矩阵检查样本 ID、OS/DFS 变量、分期和免疫检查点基因表达分布。",
        "2. 再决定是否下载全量 GDC STAR-counts 或只基于 Xena 整理矩阵做差异/预后/免疫相关分析。",
        "3. 研究问题建议收窄为：LUAD 中 TME/immune-checkpoint 相关基因签名与生存、分期、免疫状态的关联，并用 GEO/HPA 做外部支持。",
        "",
        "## 关键文件",
        "",
        "- `data/raw/TCGA.LUAD.HiSeqV2.gz`",
        "- `data/raw/TCGA.LUAD.clinicalMatrix.tsv`",
        "- `data/processed/luad_tme_immune_gene_expression_subset.csv`",
        "- `data/processed/gdc_luad_open_expression_manifest.csv`",
        "- `data/processed/cbioportal_luad_public_studies.csv`",
        "- `data/processed/geo_luad_candidate_series.csv`",
        "",
        "## 来源",
        "",
        "- TCGA Xena Hub: https://xenabrowser.net/datapages/",
        "- GDC API: https://api.gdc.cancer.gov/",
        "- cBioPortal API: https://www.cbioportal.org/api/",
        "- NCBI E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25501/",
    ]
    (REPORTS / "luad_starter_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def directory_size_mb(path: Path) -> float:
    size = sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
    return size / 1024 / 1024


def main() -> None:
    for path in (RAW, PROCESSED, REPORTS):
        path.mkdir(parents=True, exist_ok=True)
    expression = RAW / "TCGA.LUAD.HiSeqV2.gz"
    clinical = RAW / "TCGA.LUAD.clinicalMatrix.tsv"
    download(URLS["xena_luad_hiseqv2"], expression)
    download(URLS["xena_luad_clinical"], clinical)
    manifest_rows = gdc_luad_star_counts_manifest()
    cbio_rows = cbioportal_luad_studies()
    geo_rows = geo_luad_candidates()
    subset_summary = extract_gene_subset(expression)
    clinical_shape = count_clinical_rows(clinical)
    write_report(manifest_rows, cbio_rows, geo_rows, subset_summary, clinical_shape)
    print(f"Wrote report: {REPORTS / 'luad_starter_report.md'}")


if __name__ == "__main__":
    main()
