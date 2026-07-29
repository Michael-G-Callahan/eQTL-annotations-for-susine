from __future__ import annotations

from typing import Any

import pandas as pd

from .alphagenome import classify_variant_interval_mode


def parse_variant_id(variant_id: str) -> tuple[str, int, str, str]:
    parts = variant_id.split("_")
    if len(parts) < 4:
        raise ValueError(f"Cannot parse variant_id: {variant_id}")
    chrom, pos, ref, alt = parts[:4]
    return chrom, int(pos), ref, alt


def ensure_variant_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    required = {"variant_id", "z_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required variant columns: {sorted(missing)}")

    if not {"chrom", "pos", "ref", "alt"}.issubset(df.columns):
        parsed = df["variant_id"].apply(parse_variant_id)
        parsed_df = pd.DataFrame(parsed.tolist(), columns=["chrom", "pos", "ref", "alt"])
        for col in ["chrom", "pos", "ref", "alt"]:
            if col not in df.columns:
                df[col] = parsed_df[col]

    return df


def chunked(iterable_df: pd.DataFrame, chunk_size: int):
    for start in range(0, len(iterable_df), chunk_size):
        yield iterable_df.iloc[start : start + chunk_size].copy()


def annotate_variant_window_eligibility(
    df: pd.DataFrame,
    gene_info: dict,
    sequence_length_label: str,
) -> pd.DataFrame:
    annotated_df = ensure_variant_columns(df)
    annotations = annotated_df.apply(
        lambda row: classify_variant_interval_mode(row, gene_info, sequence_length_label),
        axis=1,
        result_type="expand",
    )
    return pd.concat([annotated_df.reset_index(drop=True), annotations.reset_index(drop=True)], axis=1)


def filter_variants_to_tss_radius(df: pd.DataFrame, tss: int, radius_bp: int) -> pd.DataFrame:
    filtered_df = ensure_variant_columns(df)
    filtered_df["distance_to_tss_bp"] = (filtered_df["pos"].astype(int) - int(tss)).abs()
    filtered_df["within_tss_radius"] = filtered_df["distance_to_tss_bp"] <= int(radius_bp)
    filtered_df["fits_tss_centered"] = filtered_df["within_tss_radius"]
    return filtered_df[filtered_df["within_tss_radius"]].copy()


def sanitize_for_parquet(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    sanitized = df.copy()
    for col in sanitized.columns:
        series = sanitized[col]
        if series.dtype != "object":
            continue

        non_null = series.dropna()
        if non_null.empty:
            continue

        sample = non_null.iloc[0]
        if isinstance(sample, (str, bytes, int, float, bool)):
            continue

        sanitized[col] = series.map(_coerce_complex_value)

    return sanitized


def _coerce_complex_value(value: Any):
    if pd.isna(value):
        return None
    return str(value)
