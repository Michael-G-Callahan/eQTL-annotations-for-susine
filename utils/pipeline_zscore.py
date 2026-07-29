from __future__ import annotations

import pandas as pd
import pyarrow.parquet as pq

from .annotations import download_gtf_if_needed, get_gene_info, load_gene_annotations_for_gene
from .output_layout import LocusOutputPaths
from .paths import ProjectPaths, find_gtex_parquet

SEQ_LEN = 524288
LABEL_LEN = 196608


def run_zscore_export(
    locus_cfg: dict,
    paths: ProjectPaths,
    locus_paths: LocusOutputPaths,
    *,
    write_stage_outputs: bool = True,
    write_prelim: bool = True,
    allow_cache_writes: bool = True,
) -> dict:
    parquet_path = find_gtex_parquet(
        tissue=locus_cfg["gtex_tissue"],
        chrom=locus_cfg["gtex_chrom"],
        anchor=paths.root,
    )
    table = pq.read_table(parquet_path, filters=[("gene_id", "==", locus_cfg["gene_id"])])
    df = table.to_pandas()
    raw_variant_count = len(df)

    df["sample_size"] = df["ma_count"] / (2 * df["af"])
    df = df[(df["af"] >= float(locus_cfg["maf_min"])) & (df["af"] <= float(locus_cfg["maf_max"]))].copy()
    post_af_count = len(df)
    df = df[df["sample_size"] > float(locus_cfg["min_sample_size"])].copy()
    post_sample_size_count = len(df)

    gtf_path = download_gtf_if_needed(
        paths.gtf_cache,
        genome=locus_cfg["reference_genome"],
        allow_download=allow_cache_writes,
    )
    genes_df, exons_df, _ = load_gene_annotations_for_gene(
        gtf_path,
        locus_cfg["gene_name"],
        paths.gtf_shortcuts,
        write_shortcuts=allow_cache_writes,
    )
    gene_info = get_gene_info(locus_cfg["gene_name"], genes_df, exons_df)
    gene_start = gene_info["start"]
    gene_end = gene_info["end"]
    gene_tss = gene_info["tss"]
    gene_len = gene_end - gene_start

    df["snp_pos"] = gene_tss + df["tss_distance"]
    limit_tss_centered = SEQ_LEN // 2
    df["valid_tss_centered"] = df["tss_distance"].abs() <= limit_tss_centered
    output_radius = LABEL_LEN // 2

    def check_snp_centered(row):
        snp_pos = row["snp_pos"]
        window_start = snp_pos - output_radius
        window_end = snp_pos + output_radius
        overlap_start = max(gene_start, window_start)
        overlap_end = min(gene_end, window_end)
        overlap_len = max(0, overlap_end - overlap_start)
        return overlap_len >= (0.5 * gene_len)

    df["valid_snp_centered"] = df.apply(check_snp_centered, axis=1)
    df["valid_any"] = df["valid_tss_centered"] | df["valid_snp_centered"]
    df["z_score"] = df["slope"] / df["slope_se"]

    vcf_data = df["variant_id"].apply(_parse_variant_id).tolist()
    vcf_df = pd.DataFrame(vcf_data, columns=["#CHROM", "POS", "REF", "ALT"])
    vcf_df["ID"] = df["variant_id"].values
    vcf_df["QUAL"] = "."
    vcf_df["FILTER"] = "PASS"
    vcf_df["INFO"] = "."
    vcf_df = vcf_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]]

    result = {
        "locus_id": locus_cfg["locus_id"],
        "gene_name": locus_cfg["gene_name"],
        "gene_id": locus_cfg["gene_id"],
        "gtex_tissue": locus_cfg["gtex_tissue"],
        "gtex_chrom": locus_cfg["gtex_chrom"],
        "raw_gene_variants_loaded": raw_variant_count,
        "post_af_filter": post_af_count,
        "post_sample_size_filter": post_sample_size_count,
        "z_score_variants_exported": len(df),
        "parquet_path": str(parquet_path),
        "_z_scores_df": df,
        "_vcf_variants_df": vcf_df,
    }

    if write_stage_outputs:
        vcf_path = locus_paths.variants_vcf(locus_cfg["gene_name"])
        csv_path = locus_paths.z_score_csv(locus_cfg["gene_name"])
        with vcf_path.open("w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##source=GTEx_Analysis_v10_eQTL\n")
            handle.write(f"##reference={locus_cfg['reference_genome']}\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf_df.to_csv(vcf_path, mode="a", sep="\t", index=False, header=False)
        df.to_csv(csv_path, index=False)
        result["z_score_csv"] = str(csv_path)
        result["variants_vcf"] = str(vcf_path)

    if write_prelim:
        funnel_path = locus_paths.count_funnel_csv(locus_cfg["gene_name"])
        count_funnel = pd.DataFrame(
            [
                {"step": "raw_gene_variants_loaded", "count": raw_variant_count},
                {"step": "post_af_filter", "count": post_af_count},
                {"step": "post_sample_size_filter", "count": post_sample_size_count},
                {"step": "exported_z_score_set", "count": len(df)},
            ]
        )
        count_funnel.to_csv(funnel_path, index=False)
        result["count_funnel_csv"] = str(funnel_path)
    return result


def _parse_variant_id(variant_id: str) -> tuple[str, int, str, str]:
    chrom, pos, ref, alt = variant_id.split("_")[:4]
    return chrom, int(pos), ref, alt
