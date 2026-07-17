#!/usr/bin/env python3
"""Download the two public TCGA Xena files used by the discovery workflow."""

from download_luad_starter import PROCESSED, RAW, REPORTS, URLS, download, extract_gene_subset


def main() -> None:
    for directory in (RAW, PROCESSED, REPORTS):
        directory.mkdir(parents=True, exist_ok=True)

    expression = RAW / "TCGA.LUAD.HiSeqV2.gz"
    clinical = RAW / "TCGA.LUAD.clinicalMatrix.tsv"
    download(URLS["xena_luad_hiseqv2"], expression)
    download(URLS["xena_luad_clinical"], clinical)
    summary = extract_gene_subset(expression)
    print(
        f"Downloaded TCGA-LUAD Xena inputs; recovered "
        f"{summary['found_gene_count']}/{summary['requested_gene_count']} candidate genes"
    )


if __name__ == "__main__":
    main()
