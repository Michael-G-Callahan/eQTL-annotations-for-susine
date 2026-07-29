from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pandas as pd

from .pipeline_defaults import (
    ANNOTATION_SELECTION_OPTIONAL_COLUMNS,
    ANNOTATION_SELECTION_REQUIRED_COLUMNS,
    DEFAULT_ALPHAGENOME_BATCH_SIZE,
    DEFAULT_ALPHAGENOME_MAX_WORKERS,
    DEFAULT_ALPHAGENOME_RETRY_WAIT_SECONDS,
    DEFAULT_ALPHAGENOME_SEQUENCE_LENGTH,
    DEFAULT_ALPHAGENOME_TARGET_DATA_SOURCE,
    DEFAULT_ENABLED,
    DEFAULT_MAF_MAX,
    DEFAULT_MAF_MIN,
    DEFAULT_MIN_SAMPLE_SIZE,
    DEFAULT_REFERENCE_GENOME,
    MANIFEST_OPTIONAL_COLUMNS,
    MANIFEST_REQUIRED_COLUMNS,
)


def _parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        raise ValueError("Boolean value is missing")
    normalized = str(value).strip().lower()
    if normalized in {"true", "t", "1", "yes", "y"}:
        return True
    if normalized in {"false", "f", "0", "no", "n"}:
        return False
    raise ValueError(f"Could not parse boolean value: {value!r}")


def _ensure_columns(df: pd.DataFrame, required: list[str], optional: list[str]) -> pd.DataFrame:
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = df.copy()
    for col in optional:
        if col not in out.columns:
            out[col] = pd.NA
    return out


def load_loci_manifest(path) -> pd.DataFrame:
    df = pd.read_csv(Path(path))
    df = _ensure_columns(df, MANIFEST_REQUIRED_COLUMNS, MANIFEST_OPTIONAL_COLUMNS)
    if "enabled" not in df.columns:
        df["enabled"] = DEFAULT_ENABLED

    df["enabled"] = df["enabled"].apply(lambda x: DEFAULT_ENABLED if pd.isna(x) else _parse_bool(x))
    df["reference_genome"] = df["reference_genome"].fillna(DEFAULT_REFERENCE_GENOME)
    df["maf_min"] = df["maf_min"].fillna(DEFAULT_MAF_MIN).astype(float)
    df["maf_max"] = df["maf_max"].fillna(DEFAULT_MAF_MAX).astype(float)
    df["min_sample_size"] = df["min_sample_size"].fillna(DEFAULT_MIN_SAMPLE_SIZE).astype(float)
    df["alphagenome_sequence_length"] = df["alphagenome_sequence_length"].fillna(
        DEFAULT_ALPHAGENOME_SEQUENCE_LENGTH
    )
    df["alphagenome_target_data_source"] = df["alphagenome_target_data_source"].fillna(
        DEFAULT_ALPHAGENOME_TARGET_DATA_SOURCE
    )
    df["alphagenome_target_gtex_tissue"] = df["alphagenome_target_gtex_tissue"].fillna(df["gtex_tissue"])
    df["alphagenome_batch_size"] = df["alphagenome_batch_size"].fillna(
        DEFAULT_ALPHAGENOME_BATCH_SIZE
    ).astype(int)
    df["alphagenome_max_workers"] = df["alphagenome_max_workers"].fillna(
        DEFAULT_ALPHAGENOME_MAX_WORKERS
    ).astype(int)
    df["alphagenome_retry_wait_seconds"] = df["alphagenome_retry_wait_seconds"].fillna(
        DEFAULT_ALPHAGENOME_RETRY_WAIT_SECONDS
    ).astype(int)
    return df


def validate_loci_manifest(df: pd.DataFrame) -> None:
    if df["locus_id"].duplicated().any():
        duplicates = sorted(df.loc[df["locus_id"].duplicated(), "locus_id"].unique().tolist())
        raise ValueError(f"Duplicate locus_id values: {duplicates}")

    enabled_df = df[df["enabled"]].copy()
    for col in MANIFEST_REQUIRED_COLUMNS:
        if enabled_df[col].astype(str).str.strip().eq("").any() or enabled_df[col].isna().any():
            raise ValueError(f"Enabled rows must have non-empty values for {col}")

    chrom_ok = enabled_df["gtex_chrom"].astype(str).str.match(r"^chr[0-9XYM]+$")
    if not chrom_ok.all():
        bad = enabled_df.loc[~chrom_ok, ["locus_id", "gtex_chrom"]].to_dict("records")
        raise ValueError(f"Invalid gtex_chrom values: {bad}")

    if (enabled_df["maf_min"] > enabled_df["maf_max"]).any():
        raise ValueError("maf_min cannot exceed maf_max")

    if (enabled_df["min_sample_size"] <= 0).any():
        raise ValueError("min_sample_size must be positive")

    if (enabled_df["alphagenome_batch_size"] <= 0).any():
        raise ValueError("alphagenome_batch_size must be positive")

    if (enabled_df["alphagenome_max_workers"] <= 0).any():
        raise ValueError("alphagenome_max_workers must be positive")

    if (enabled_df["alphagenome_retry_wait_seconds"] < 0).any():
        raise ValueError("alphagenome_retry_wait_seconds must be non-negative")


def load_annotation_selection(path) -> pd.DataFrame:
    df = pd.read_csv(Path(path))
    df = _ensure_columns(
        df, ANNOTATION_SELECTION_REQUIRED_COLUMNS, ANNOTATION_SELECTION_OPTIONAL_COLUMNS
    )
    df["annotate"] = df["annotate"].apply(_parse_bool)
    return df


def validate_annotation_selection(df: pd.DataFrame, loci_df: pd.DataFrame) -> None:
    if df["locus_id"].duplicated().any():
        duplicates = sorted(df.loc[df["locus_id"].duplicated(), "locus_id"].unique().tolist())
        raise ValueError(f"Duplicate locus_id values in annotation selection: {duplicates}")

    manifest_ids = set(loci_df["locus_id"])
    unknown = sorted(set(df["locus_id"]) - manifest_ids)
    if unknown:
        raise ValueError(f"Annotation selection contains unknown locus_id values: {unknown}")


def iter_enabled_loci(df: pd.DataFrame) -> Iterable[dict]:
    for row in df[df["enabled"]].to_dict("records"):
        yield row


def merge_annotation_selection(loci_df: pd.DataFrame, selection_df: pd.DataFrame) -> pd.DataFrame:
    merged = loci_df.merge(selection_df, on="locus_id", how="inner", suffixes=("", "_selection"))
    merged = merged[merged["annotate"]].copy()
    if merged.empty:
        return merged

    if "annotation_gene_name_override" in merged.columns:
        merged["gene_name"] = merged["annotation_gene_name_override"].fillna(merged["gene_name"])
    if "annotation_gene_id_override" in merged.columns:
        merged["gene_id"] = merged["annotation_gene_id_override"].fillna(merged["gene_id"])
    if "annotation_tissue_override" in merged.columns:
        merged["alphagenome_target_gtex_tissue"] = merged["annotation_tissue_override"].fillna(
            merged["alphagenome_target_gtex_tissue"]
        )
    return merged
