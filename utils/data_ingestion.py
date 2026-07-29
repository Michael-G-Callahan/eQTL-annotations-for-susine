from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from .annotations import (
    download_gtf_if_needed,
    ensure_dir,
    get_gene_info,
    load_gene_annotations_for_gene,
)


def load_vcf(vcf_path: str, *, snp_only: bool = False) -> pd.DataFrame:
    """
    Load VCF file into a DataFrame.
    Handles standard VCF format.
    """
    vcf_path = str(vcf_path)
    if not Path(vcf_path).exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    opener = gzip.open if vcf_path.endswith(".gz") else open

    # Find header line to get column names
    header_line = None
    skip_rows = 0
    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                skip_rows += 1
            elif line.startswith("#CHROM"):
                header_line = line.strip().lstrip("#").split("\t")
                skip_rows += 1
                break
            else:
                break

    if header_line is None:
        header_line = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    df = pd.read_csv(
        vcf_path,
        sep="\t",
        comment="#",
        header=None,
        names=header_line[:8],
        usecols=range(min(8, len(header_line))),
    )

    df = df.rename(
        columns={"CHROM": "chrom", "POS": "pos", "ID": "rsID", "REF": "ref", "ALT": "alt"}
    )

    if not df.empty and not df["chrom"].iloc[0].startswith("chr"):
        df["chrom"] = "chr" + df["chrom"].astype(str)

    df["pos_0based"] = df["pos"] - 1

    if snp_only:
        is_snp = (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)
        n_before = len(df)
        df = df[is_snp].reset_index(drop=True)
        print(f"Filtered {n_before - len(df)} non-SNP variants, {len(df)} SNPs remaining")

    return df
