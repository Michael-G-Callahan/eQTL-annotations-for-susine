from __future__ import annotations

import csv
import urllib.request
from io import StringIO

import numpy as np
import pandas as pd
import pysam

from .output_layout import LocusOutputPaths
from .paths import ProjectPaths


def run_ld_phase1(
    locus_cfg: dict,
    paths: ProjectPaths,
    locus_paths: LocusOutputPaths,
    *,
    z_scores_df: pd.DataFrame | None = None,
    variants_df: pd.DataFrame | None = None,
    write_stage_outputs: bool = True,
    write_prelim: bool = True,
    allow_cache_writes: bool = True,
) -> dict:
    if variants_df is None:
        vcf_path = locus_paths.variants_vcf(locus_cfg["gene_name"])
        if not vcf_path.exists():
            raise FileNotFoundError(f"Missing input VCF: {vcf_path}")
        variants = _load_variants_from_vcf(vcf_path)
    else:
        vcf_path = None
        variants = _load_variants_from_df(variants_df)
    if not variants:
        source_desc = str(vcf_path) if vcf_path is not None else "variants_df"
        raise ValueError(f"No variants found in input source: {source_desc}")

    eur_samples = _load_eur_samples(paths.data, allow_cache_writes=allow_cache_writes)
    G, variants, sample_names, hits = _fetch_1000g_genotypes(variants, eur_samples)

    G_mean = np.nanmean(G, axis=1, keepdims=True)
    G_filled = np.where(np.isnan(G), G_mean, G)
    X = G_filled - G_mean
    ddof = 1 if G.shape[1] > 1 else 0
    stds = X.std(axis=1, ddof=ddof, keepdims=True)
    stds[stds == 0] = np.nan
    X = X / stds
    X = np.nan_to_num(X, nan=0.0)
    R = (X @ X.T) / max(G.shape[1] - 1, 1)
    np.fill_diagonal(R, 1.0)
    R = np.round(R, 3).astype(np.float32)

    variant_labels = [v["id"] for v in variants]
    if z_scores_df is None:
        z_scores = pd.read_csv(locus_paths.z_score_csv(locus_cfg["gene_name"]))
    else:
        z_scores = z_scores_df.copy()
    ld_index = {variant_id: idx for idx, variant_id in enumerate(variant_labels)}
    master_df = z_scores[z_scores["variant_id"].isin(ld_index)].copy()
    master_df["ld_included"] = True
    master_df["ld_matrix_index"] = master_df["variant_id"].map(ld_index)
    master_df = master_df.sort_values("ld_matrix_index").reset_index(drop=True)
    if len(master_df) != len(variant_labels):
        missing_from_z = [variant_id for variant_id in variant_labels if variant_id not in set(master_df["variant_id"])]
        raise ValueError(f"LD matched variants missing from z-score table: {missing_from_z[:10]}")

    variant_map_df = pd.DataFrame(
        {
            "snp_index": np.arange(len(variants), dtype=np.int32),
            "variant_id": [v["id"] for v in variants],
            "chrom": [v["chrom"] for v in variants],
            "pos": [v["pos"] for v in variants],
            "ref": [v["ref"] for v in variants],
            "alt": [v["alt"] for v in variants],
        }
    )
    z_compact_df = (
        master_df[["variant_id", "z_score", "sample_size"]]
        .merge(variant_map_df[["snp_index", "variant_id"]], on="variant_id", how="inner")
        .sort_values("snp_index")
        .reset_index(drop=True)
    )

    tri_i, tri_j = np.triu_indices(len(variant_labels), k=1)
    ld_long_df = pd.DataFrame(
        {
            "snp_index_1": tri_i.astype(np.int32),
            "snp_index_2": tri_j.astype(np.int32),
            "r": R[tri_i, tri_j].astype(np.float32),
        }
    )

    z = master_df["z_score"].to_numpy(dtype=np.float64)
    A = np.abs(R.astype(np.float64, copy=False))
    ut_mask = np.triu(np.ones(A.shape, dtype=bool), k=1)
    ut_abs = A[ut_mask]
    centered_r = R.astype(np.float64, copy=False) - 0.5
    ut_centered = centered_r[ut_mask]
    z2 = np.square(z)
    z4 = np.square(z2)

    metrics = {
        "M1": float(2.0 * np.mean(ut_abs * (1.0 - ut_abs))),
        "M2": float(2.0 * np.mean(np.square(ut_centered))),
        "M4": float(2.0 * np.mean(np.power(ut_centered, 4))),
        "z_count_abs_gt_3": int(np.sum(np.abs(z) > 3.0)),
        "z_eff_signals": float((np.sum(z2) ** 2) / np.sum(z4)) if np.sum(z4) > 0 else np.nan,
    }
    metrics_df = pd.DataFrame(
        [
            {
                "gene_name": locus_cfg["gene_name"],
                "n_variants_z": len(z_scores),
                "n_variants_ld": len(variant_labels),
                "n_variants_master": len(master_df),
                **metrics,
            }
        ]
    )
    result = {
        "locus_id": locus_cfg["locus_id"],
        "ld_input_variants": len(z_scores),
        "ld_matched_variants": len(variant_labels),
        "ld_unmatched_variants": len(z_scores) - len(variant_labels),
        "phase1_master_variants": len(master_df),
        "M1": metrics["M1"],
        "M2": metrics["M2"],
        "M4": metrics["M4"],
        "z_count_abs_gt_3": metrics["z_count_abs_gt_3"],
        "z_eff_signals": metrics["z_eff_signals"],
        "eur_samples_used": len(sample_names),
        "matched_variants_in_1000g": hits,
    }
    if write_stage_outputs:
        order_path = locus_paths.ld_variant_order_tsv(locus_cfg["gene_name"])
        with order_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["index", "id", "chrom", "pos", "ref", "alt"])
            for i, variant in enumerate(variants):
                writer.writerow([i, variant["id"], variant["chrom"], variant["pos"], variant["ref"], variant["alt"]])

        master_path = locus_paths.phase1_master_variants_csv(locus_cfg["gene_name"])
        master_df.to_csv(master_path, index=False)

        variant_map_path = locus_paths.phase1_variant_map_parquet(locus_cfg["gene_name"])
        variant_map_df.to_parquet(variant_map_path, index=False, compression="zstd")

        z_compact_path = locus_paths.phase1_z_scores_parquet(locus_cfg["gene_name"])
        z_compact_df.to_parquet(z_compact_path, index=False, compression="zstd")

        ld_long_path = locus_paths.phase1_ld_long_parquet(locus_cfg["gene_name"])
        ld_long_df.to_parquet(ld_long_path, index=False, compression="zstd")

        result["ld_variant_order_tsv"] = str(order_path)
        result["phase1_master_variants_csv"] = str(master_path)
        result["phase1_variant_map_parquet"] = str(variant_map_path)
        result["phase1_z_scores_parquet"] = str(z_compact_path)
        result["phase1_ld_long_parquet"] = str(ld_long_path)
    if write_prelim:
        metrics_path = locus_paths.dataset_metrics_csv(locus_cfg["gene_name"])
        metrics_df.to_csv(metrics_path, index=False)

        funnel_path = locus_paths.count_funnel_csv(locus_cfg["gene_name"])
        if funnel_path.exists():
            count_funnel = pd.read_csv(funnel_path)
        else:
            count_funnel = pd.DataFrame(columns=["step", "count"])
        ld_steps = pd.DataFrame(
            [
                {"step": "ld_matched", "count": len(variant_labels)},
                {"step": "ld_unmatched", "count": len(z_scores) - len(variant_labels)},
                {"step": "phase1_master_z_intersect_ld", "count": len(master_df)},
            ]
        )
        count_funnel = pd.concat(
            [count_funnel[~count_funnel["step"].isin(ld_steps["step"])], ld_steps],
            ignore_index=True,
        )
        count_funnel.to_csv(funnel_path, index=False)
        result["dataset_metrics_csv"] = str(metrics_path)
        result["count_funnel_csv"] = str(funnel_path)
    return result


def load_existing_phase1_result(locus_cfg: dict, locus_paths: LocusOutputPaths) -> dict:
    funnel_df = pd.read_csv(locus_paths.count_funnel_csv(locus_cfg["gene_name"]))
    metrics_df = pd.read_csv(locus_paths.dataset_metrics_csv(locus_cfg["gene_name"]))
    counts = dict(zip(funnel_df["step"], funnel_df["count"]))
    metrics_row = metrics_df.iloc[0].to_dict()
    return {
        "locus_id": locus_cfg["locus_id"],
        "raw_gene_variants_loaded": int(counts.get("raw_gene_variants_loaded", 0)),
        "post_af_filter": int(counts.get("post_af_filter", 0)),
        "post_sample_size_filter": int(counts.get("post_sample_size_filter", 0)),
        "z_score_variants_exported": int(counts.get("exported_z_score_set", 0)),
        "ld_input_variants": int(counts.get("exported_z_score_set", counts.get("post_sample_size_filter", 0))),
        "ld_matched_variants": int(counts.get("ld_matched", 0)),
        "ld_unmatched_variants": int(counts.get("ld_unmatched", 0)),
        "phase1_master_variants": int(counts.get("phase1_master_z_intersect_ld", 0)),
        "M1": metrics_row.get("M1"),
        "M2": metrics_row.get("M2"),
        "M4": metrics_row.get("M4"),
        "z_count_abs_gt_3": metrics_row.get("z_count_abs_gt_3"),
        "z_eff_signals": metrics_row.get("z_eff_signals"),
        "ld_variant_order_tsv": str(locus_paths.ld_variant_order_tsv(locus_cfg["gene_name"])),
        "phase1_master_variants_csv": str(locus_paths.phase1_master_variants_csv(locus_cfg["gene_name"])),
        "phase1_variant_map_parquet": str(locus_paths.phase1_variant_map_parquet(locus_cfg["gene_name"])),
        "phase1_z_scores_parquet": str(locus_paths.phase1_z_scores_parquet(locus_cfg["gene_name"])),
        "phase1_ld_long_parquet": str(locus_paths.phase1_ld_long_parquet(locus_cfg["gene_name"])),
        "dataset_metrics_csv": str(locus_paths.dataset_metrics_csv(locus_cfg["gene_name"])),
        "count_funnel_csv": str(locus_paths.count_funnel_csv(locus_cfg["gene_name"])),
    }


def _load_variants_from_vcf(vcf_path):
    variants = []
    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, vid, ref, alt = fields[:5]
            variants.append(
                {"chrom": chrom, "pos": int(pos), "id": vid, "ref": ref, "alt": alt.split(",")[0]}
            )
    return variants


def _load_variants_from_df(variants_df: pd.DataFrame) -> list[dict]:
    required_columns = ("#CHROM", "POS", "ID", "REF", "ALT")
    missing_columns = set(required_columns).difference(variants_df.columns)
    if missing_columns:
        raise ValueError(f"Missing columns in variants_df: {sorted(missing_columns)}")
    variants = []
    for row in variants_df[list(required_columns)].to_dict("records"):
        variants.append(
            {
                "chrom": str(row["#CHROM"]),
                "pos": int(row["POS"]),
                "id": row["ID"],
                "ref": row["REF"],
                "alt": str(row["ALT"]).split(",")[0],
            }
        )
    return variants


def _load_eur_samples(data_dir, *, allow_cache_writes: bool = True):
    panel_url = (
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
        "integrated_call_samples_v3.20130502.ALL.panel"
    )
    panel_path = data_dir / "integrated_call_samples_v3.20130502.ALL.panel"
    if panel_path.exists():
        with panel_path.open() as handle:
            handle.readline()
            return [
                fields[0]
                for line in handle
                if (fields := line.strip().split()) and len(fields) >= 3 and fields[2] == "EUR"
            ]

    if allow_cache_writes:
        urllib.request.urlretrieve(panel_url, panel_path)
        with panel_path.open() as handle:
            handle.readline()
            return [
                fields[0]
                for line in handle
                if (fields := line.strip().split()) and len(fields) >= 3 and fields[2] == "EUR"
            ]

    with urllib.request.urlopen(panel_url) as response:
        content = response.read().decode("utf-8")
    handle = StringIO(content)
    handle.readline()
    return [
        fields[0]
        for line in handle
        if (fields := line.strip().split()) and len(fields) >= 3 and fields[2] == "EUR"
    ]


def _fetch_1000g_genotypes(variants, eur_samples):
    chroms = sorted({v["chrom"] for v in variants})
    min_pos = min(v["pos"] for v in variants)
    max_pos = max(v["pos"] for v in variants)
    base_url = (
        "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
        "data_collections/1000_genomes_project/release/"
        "20190312_biallelic_SNV_and_INDEL"
    )
    chrom_no_chr = chroms[0].replace("chr", "")
    vcf_url = (
        f"{base_url}/ALL.chr{chrom_no_chr}."
        "shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
    )

    vcf_in = pysam.VariantFile(vcf_url)
    contigs = set(vcf_in.header.contigs)
    chrom_with_chr = f"chr{chrom_no_chr}"
    if chroms[0] in contigs:
        fetch_chrom = chroms[0]
    elif chrom_with_chr in contigs:
        fetch_chrom = chrom_with_chr
    elif chrom_no_chr in contigs:
        fetch_chrom = chrom_no_chr
    else:
        raise ValueError(f"Chromosome not found in VCF contigs: {chroms[0]}")

    eur_set = set(eur_samples)
    sample_names = [sample for sample in vcf_in.header.samples if sample in eur_set]
    if not sample_names:
        raise ValueError("No EUR samples from the panel were found in the VCF header.")

    key_to_index = {}
    for idx, variant in enumerate(variants):
        chrom_key = variant["chrom"]
        chrom_no_chr = chrom_key.replace("chr", "")
        chrom_with_chr = f"chr{chrom_no_chr}"
        for candidate in {chrom_key, chrom_no_chr, chrom_with_chr}:
            key_to_index[(candidate, variant["pos"], variant["ref"], variant["alt"])] = idx

    G = np.full((len(variants), len(sample_names)), np.nan, dtype=np.float32)
    hits = 0
    for rec in vcf_in.fetch(fetch_chrom, min_pos - 1, max_pos):
        if not rec.alts:
            continue
        key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
        idx = key_to_index.get(key)
        if idx is None:
            continue
        hits += 1
        for j, sample in enumerate(sample_names):
            gt = rec.samples[sample].get("GT")
            if gt is None or None in gt or -1 in gt:
                continue
            G[idx, j] = gt[0] + gt[1]

    found_mask = ~np.isnan(G).all(axis=1)
    if not found_mask.all():
        G = G[found_mask]
        variants = [variant for variant, ok in zip(variants, found_mask) if ok]

    try:
        vcf_in.close()
    except OSError:
        pass

    return G, variants, sample_names, hits
