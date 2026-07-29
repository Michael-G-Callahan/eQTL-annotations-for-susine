from __future__ import annotations

import time
from typing import Sequence

import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

TSS_CENTERED_MODE = "tss_centered"
MIDPOINT_CENTERED_MODE = "midpoint_centered"
INELIGIBLE_MODE = "ineligible"
INELIGIBLE_EXON_VARIANT_SPAN = "outside_tss_interval_and_gene_variant_span_exceeds_sequence_length"
INELIGIBLE_NO_EXONS = "missing_gene_exons"
INELIGIBLE_UNKNOWN = "unknown"


def resolve_sequence_length_bp(sequence_length_label: str) -> int:
    return int(dna_client.SUPPORTED_SEQUENCE_LENGTHS[f"SEQUENCE_LENGTH_{sequence_length_label}"])


def tss_radius_bp(sequence_length_label: str) -> int:
    return resolve_sequence_length_bp(sequence_length_label) // 2


def interval_bounds_1based(interval: genome.Interval) -> tuple[int, int]:
    return int(interval.start) + 1, int(interval.end)


def variant_allele_span_1based(row) -> dict:
    variant_start_1based = int(row["pos"])
    ref_length = len(str(row["ref"]))
    alt_length = len(str(row["alt"]))
    ref_end_1based = variant_start_1based + ref_length - 1
    alt_end_1based = variant_start_1based + alt_length - 1
    return {
        "variant_start_1based": variant_start_1based,
        "ref_end_1based": ref_end_1based,
        "alt_end_1based": alt_end_1based,
        "variant_min_1based": variant_start_1based,
        "variant_max_1based": max(ref_end_1based, alt_end_1based),
        "ref_length": ref_length,
        "alt_length": alt_length,
    }


def build_centered_interval(
    chromosome: str,
    center_1based: int,
    sequence_length_label: str,
    strand: str,
    *,
    name: str = "",
) -> genome.Interval:
    anchor = genome.Interval(
        chromosome=chromosome,
        start=int(center_1based) - 1,
        end=int(center_1based),
        strand=strand,
        name=name,
    )
    return anchor.resize(resolve_sequence_length_bp(sequence_length_label))


def build_interval_covering_span(
    chromosome: str,
    span_min_1based: int,
    span_max_1based: int,
    sequence_length_label: str,
    strand: str,
    *,
    name: str = "",
) -> genome.Interval:
    sequence_length_bp = resolve_sequence_length_bp(sequence_length_label)
    span_min_1based = int(span_min_1based)
    span_max_1based = int(span_max_1based)
    span_bp = span_max_1based - span_min_1based + 1
    if span_bp > sequence_length_bp:
        raise ValueError(
            f"Span {chromosome}:{span_min_1based}-{span_max_1based} ({span_bp} bp) exceeds sequence length "
            f"{sequence_length_bp} bp."
        )

    padding_left = (sequence_length_bp - span_bp) // 2
    start_1based = span_min_1based - padding_left
    end_1based = start_1based + sequence_length_bp - 1
    return genome.Interval(
        chromosome=chromosome,
        start=start_1based - 1,
        end=end_1based,
        strand=strand,
        name=name,
    )


def build_tss_centered_interval(
    gene_info: dict,
    sequence_length_label: str,
    *,
    name: str = "",
) -> genome.Interval:
    return build_centered_interval(
        gene_info["chrom"],
        int(gene_info["tss"]),
        sequence_length_label,
        gene_info.get("strand", "."),
        name=name,
    )


def gene_exon_span_1based(gene_info: dict) -> tuple[int, int]:
    exons = gene_info.get("exons")
    if exons is None or exons.empty:
        raise ValueError("Target gene has no exon coordinates available for AlphaGenome interval selection.")
    exon_min_1based = int(exons["start"].astype(int).min()) + 1
    exon_max_1based = int(exons["end"].astype(int).max())
    return exon_min_1based, exon_max_1based


def build_midpoint_centered_interval(
    gene_info: dict,
    row,
    sequence_length_label: str,
    *,
    name: str = "",
) -> genome.Interval:
    exon_min_1based, exon_max_1based = gene_exon_span_1based(gene_info)
    allele_span = variant_allele_span_1based(row)
    combined_min_1based = min(exon_min_1based, allele_span["variant_min_1based"])
    combined_max_1based = max(exon_max_1based, allele_span["variant_max_1based"])
    return build_interval_covering_span(
        gene_info["chrom"],
        combined_min_1based,
        combined_max_1based,
        sequence_length_label,
        gene_info.get("strand", "."),
        name=name,
    )


def interval_fully_contains_variant(interval: genome.Interval, row) -> bool:
    interval_start_1based, interval_end_1based = interval_bounds_1based(interval)
    allele_span = variant_allele_span_1based(row)
    return (
        allele_span["variant_min_1based"] >= interval_start_1based
        and allele_span["variant_max_1based"] <= interval_end_1based
    )


def interval_fully_contains_gene_exons(interval: genome.Interval, gene_info: dict) -> bool:
    exon_min_1based, exon_max_1based = gene_exon_span_1based(gene_info)
    interval_start_1based, interval_end_1based = interval_bounds_1based(interval)
    return exon_min_1based >= interval_start_1based and exon_max_1based <= interval_end_1based


def classify_variant_interval_mode(row, gene_info: dict, sequence_length_label: str) -> dict:
    allele_span = variant_allele_span_1based(row)
    distance_to_tss_bp = abs(int(row["pos"]) - int(gene_info["tss"]))
    tss_interval = build_tss_centered_interval(gene_info, sequence_length_label, name=f"{gene_info['gene_name']}_tss")
    tss_interval_start_1based, tss_interval_end_1based = interval_bounds_1based(tss_interval)
    fits_tss_centered = interval_fully_contains_variant(tss_interval, row)

    try:
        gene_exon_min_1based, gene_exon_max_1based = gene_exon_span_1based(gene_info)
    except ValueError:
        return {
            **allele_span,
            "distance_to_tss_bp": distance_to_tss_bp,
            "tss_interval_start_1based": tss_interval_start_1based,
            "tss_interval_end_1based": tss_interval_end_1based,
            "fits_tss_centered": fits_tss_centered,
            "within_tss_radius": fits_tss_centered,
            "midpoint_span_min_1based": pd.NA,
            "midpoint_span_max_1based": pd.NA,
            "midpoint_span_bp": pd.NA,
            "fits_midpoint_centered": False,
            "scoring_mode": INELIGIBLE_MODE,
            "excluded_from_scoring": True,
            "exclusion_reason": INELIGIBLE_NO_EXONS,
            "scoring_interval": None,
            "scoring_interval_start_1based": pd.NA,
            "scoring_interval_end_1based": pd.NA,
        }

    combined_min_1based = min(gene_exon_min_1based, allele_span["variant_min_1based"])
    combined_max_1based = max(gene_exon_max_1based, allele_span["variant_max_1based"])
    combined_span_bp = combined_max_1based - combined_min_1based + 1
    sequence_length_bp = resolve_sequence_length_bp(sequence_length_label)

    if fits_tss_centered:
        return {
            **allele_span,
            "distance_to_tss_bp": distance_to_tss_bp,
            "tss_interval_start_1based": tss_interval_start_1based,
            "tss_interval_end_1based": tss_interval_end_1based,
            "fits_tss_centered": True,
            "within_tss_radius": True,
            "midpoint_span_min_1based": combined_min_1based,
            "midpoint_span_max_1based": combined_max_1based,
            "midpoint_span_bp": combined_span_bp,
            "fits_midpoint_centered": False,
            "scoring_mode": TSS_CENTERED_MODE,
            "excluded_from_scoring": False,
            "exclusion_reason": "",
            "scoring_interval": tss_interval,
            "scoring_interval_start_1based": tss_interval_start_1based,
            "scoring_interval_end_1based": tss_interval_end_1based,
        }

    midpoint_interval = None
    midpoint_interval_start_1based = pd.NA
    midpoint_interval_end_1based = pd.NA
    fits_midpoint_centered = False
    if combined_span_bp <= sequence_length_bp:
        midpoint_interval = build_midpoint_centered_interval(
            gene_info,
            row,
            sequence_length_label,
            name=f"{gene_info['gene_name']}_midpoint",
        )
        midpoint_interval_start_1based, midpoint_interval_end_1based = interval_bounds_1based(midpoint_interval)
        fits_midpoint_centered = interval_fully_contains_variant(
            midpoint_interval, row
        ) and interval_fully_contains_gene_exons(midpoint_interval, gene_info)

    if fits_midpoint_centered and midpoint_interval is not None:
        return {
            **allele_span,
            "distance_to_tss_bp": distance_to_tss_bp,
            "tss_interval_start_1based": tss_interval_start_1based,
            "tss_interval_end_1based": tss_interval_end_1based,
            "fits_tss_centered": False,
            "within_tss_radius": False,
            "midpoint_span_min_1based": combined_min_1based,
            "midpoint_span_max_1based": combined_max_1based,
            "midpoint_span_bp": combined_span_bp,
            "fits_midpoint_centered": True,
            "scoring_mode": MIDPOINT_CENTERED_MODE,
            "excluded_from_scoring": False,
            "exclusion_reason": "",
            "scoring_interval": midpoint_interval,
            "scoring_interval_start_1based": midpoint_interval_start_1based,
            "scoring_interval_end_1based": midpoint_interval_end_1based,
        }

    exclusion_reason = (
        INELIGIBLE_EXON_VARIANT_SPAN if combined_span_bp > sequence_length_bp else INELIGIBLE_UNKNOWN
    )
    return {
        **allele_span,
        "distance_to_tss_bp": distance_to_tss_bp,
        "tss_interval_start_1based": tss_interval_start_1based,
        "tss_interval_end_1based": tss_interval_end_1based,
        "fits_tss_centered": False,
        "within_tss_radius": False,
        "midpoint_span_min_1based": combined_min_1based,
        "midpoint_span_max_1based": combined_max_1based,
        "midpoint_span_bp": combined_span_bp,
        "fits_midpoint_centered": False,
        "scoring_mode": INELIGIBLE_MODE,
        "excluded_from_scoring": True,
        "exclusion_reason": exclusion_reason,
        "scoring_interval": None,
        "scoring_interval_start_1based": midpoint_interval_start_1based,
        "scoring_interval_end_1based": midpoint_interval_end_1based,
    }


def build_variant_inputs(
    batch_df: pd.DataFrame,
    *,
    interval_col: str = "scoring_interval",
):
    variants = []
    intervals = []
    for _, row in batch_df.iterrows():
        interval = row[interval_col]
        if interval is None:
            raise ValueError(f"Missing scoring interval for variant {row['variant_id']}")
        variants.append(
            genome.Variant(
                chromosome=row["chrom"],
                position=int(row["pos"]),
                reference_bases=row["ref"],
                alternate_bases=row["alt"],
            )
        )
        intervals.append(interval)
    return variants, intervals


def score_batch(
    dna_model,
    batch_df: pd.DataFrame,
    *,
    gene_info: dict,
    gene_name: str,
    gene_id: str,
    sequence_length_label: str,
    active_variant_scorers: Sequence,
    max_workers: int = 1,
    progress_bar: bool = False,
    retry_wait_seconds: int = 5,
    interval_col: str = "scoring_interval",
) -> tuple[pd.DataFrame, pd.DataFrame, float, dict]:
    variants, scoring_intervals = build_variant_inputs(batch_df, interval_col=interval_col)
    interval_bounds = [interval_bounds_1based(interval) for interval in scoring_intervals]
    shared_interval = len(set(interval_bounds)) == 1
    grouped_indices: dict[tuple[int, int], list[int]] = {}
    for idx, bounds in enumerate(interval_bounds):
        grouped_indices.setdefault(bounds, []).append(idx)

    execution_info = {
        "requested_max_workers": int(max_workers),
        "used_parallel": False,
        "interval_group_count": len(set(interval_bounds)),
        "largest_interval_group_size": max((len(indices) for indices in grouped_indices.values()), default=0),
        "retry_count": 0,
        "last_retry_reason": "",
    }

    started = time.perf_counter()
    while True:
        try:
            if max_workers > 1 and shared_interval:
                execution_info["used_parallel"] = True
                score_groups = dna_model.score_variants(
                    intervals=scoring_intervals[0],
                    variants=variants,
                    variant_scorers=active_variant_scorers,
                    progress_bar=progress_bar,
                    max_workers=max_workers,
                )
            elif max_workers > 1:
                score_groups = [None] * len(variants)
                for indices in grouped_indices.values():
                    if len(indices) == 1:
                        idx = indices[0]
                        score_groups[idx] = dna_model.score_variant(
                            interval=scoring_intervals[idx],
                            variant=variants[idx],
                            variant_scorers=active_variant_scorers,
                        )
                        continue

                    execution_info["used_parallel"] = True
                    grouped_scores = dna_model.score_variants(
                        intervals=scoring_intervals[indices[0]],
                        variants=[variants[idx] for idx in indices],
                        variant_scorers=active_variant_scorers,
                        progress_bar=progress_bar,
                        max_workers=max_workers,
                    )
                    for idx, grouped_score in zip(indices, grouped_scores):
                        score_groups[idx] = grouped_score
            else:
                score_groups = [
                    dna_model.score_variant(
                        interval=scoring_interval,
                        variant=variant,
                        variant_scorers=active_variant_scorers,
                    )
                    for scoring_interval, variant in zip(scoring_intervals, variants)
                ]
            break
        except Exception as exc:
            error_text = str(exc)
            if "RESOURCE_EXHAUSTED" not in error_text and "Quota exceeded" not in error_text:
                raise

            execution_info["retry_count"] += 1
            execution_info["last_retry_reason"] = error_text
            time.sleep(retry_wait_seconds)
    elapsed_seconds = time.perf_counter() - started

    tidy_frames = []
    variant_records = []
    for (_, row), variant_scores in zip(batch_df.iterrows(), score_groups):
        variant_record = {
            "variant_id": row["variant_id"],
            "chrom": row["chrom"],
            "pos": int(row["pos"]),
            "ref": row["ref"],
            "alt": row["alt"],
            "z_score": row["z_score"],
            "distance_to_tss_bp": int(row["distance_to_tss_bp"]),
            "within_tss_radius": bool(row["within_tss_radius"]),
            "fits_tss_centered": bool(row["fits_tss_centered"]),
            "fits_midpoint_centered": bool(row["fits_midpoint_centered"]),
            "scoring_mode": row["scoring_mode"],
            "excluded_from_scoring": bool(row["excluded_from_scoring"]),
            "exclusion_reason": row["exclusion_reason"],
            "scoring_interval_start_1based": int(row["scoring_interval_start_1based"]),
            "scoring_interval_end_1based": int(row["scoring_interval_end_1based"]),
            "n_scores_total": 0,
            "score_status": "api_error",
            "score_status_reason": "",
        }

        try:
            tidy_df = variant_scorers.tidy_scores([variant_scores], match_gene_strand=True)
            if tidy_df is not None:
                tidy_df["source_variant_id"] = row["variant_id"]
                tidy_df["source_gene_name"] = gene_name
                tidy_df["source_gene_id"] = gene_id
                tidy_df["distance_to_tss_bp"] = int(row["distance_to_tss_bp"])
                tidy_df["source_scoring_mode"] = row["scoring_mode"]
                tidy_df["source_scoring_interval_start_1based"] = int(row["scoring_interval_start_1based"])
                tidy_df["source_scoring_interval_end_1based"] = int(row["scoring_interval_end_1based"])
                tidy_frames.append(tidy_df)
                variant_record["n_scores_total"] = len(tidy_df)
            variant_record["score_status"] = "scored"
        except Exception as exc:
            variant_record["score_status"] = "api_error"
            variant_record["score_status_reason"] = str(exc)
        variant_records.append(variant_record)

    full_scores_df = pd.concat(tidy_frames, ignore_index=True) if tidy_frames else pd.DataFrame()
    scored_variants_df = pd.DataFrame(variant_records)
    return full_scores_df, scored_variants_df, elapsed_seconds, execution_info


def build_run_summary(
    *,
    gene_name: str,
    gene_id: str,
    variant_source_file: str,
    sequence_length: str,
    centering_mode: str,
    gene_info: dict,
    scorer_mode: str,
    batch_size: int,
    scoring_max_workers: int,
    max_variants,
    variants_loaded_count: int,
    variants_df: pd.DataFrame,
    scored_variants_df: pd.DataFrame,
    full_scores_df: pd.DataFrame,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "variant_source_file": variant_source_file,
                "sequence_length": sequence_length,
                "sequence_length_bp": resolve_sequence_length_bp(sequence_length),
                "centering_mode": centering_mode,
                "tss_bp": int(gene_info["tss"]),
                "tss_radius_bp": tss_radius_bp(sequence_length),
                "scorer_mode": scorer_mode,
                "batch_size": batch_size,
                "scoring_max_workers": scoring_max_workers,
                "max_variants": max_variants if max_variants is not None else "",
                "variants_loaded": int(variants_loaded_count),
                "variants_within_tss_radius": int(variants_df["fits_tss_centered"].sum()),
                "variants_attempted": len(scored_variants_df),
                "variants_scored_ok": int((scored_variants_df["score_status"] == "scored").sum())
                if not scored_variants_df.empty
                else 0,
                "variants_api_error": int((scored_variants_df["score_status"] == "api_error").sum())
                if not scored_variants_df.empty
                else 0,
                "full_score_rows": len(full_scores_df),
                "distinct_variants_in_full_scores": int(full_scores_df["source_variant_id"].nunique())
                if not full_scores_df.empty and "source_variant_id" in full_scores_df.columns
                else 0,
                "distinct_gene_ids_returned": int(full_scores_df["gene_id"].nunique())
                if not full_scores_df.empty and "gene_id" in full_scores_df.columns
                else 0,
            }
        ]
    )


def filter_scores_to_target_gene(
    full_scores_df: pd.DataFrame,
    *,
    gene_name: str,
    gene_id: str,
) -> pd.DataFrame:
    if full_scores_df.empty:
        return full_scores_df

    gene_id_base = gene_id.split(".", 1)[0]
    mask = pd.Series(False, index=full_scores_df.index)

    if "gene_id" in full_scores_df.columns:
        gene_ids = full_scores_df["gene_id"].astype(str)
        mask = gene_ids.eq(gene_id) | gene_ids.str.replace(r"\.\d+$", "", regex=True).eq(gene_id_base)

    if "gene_name" in full_scores_df.columns:
        mask = mask | full_scores_df["gene_name"].astype(str).str.upper().eq(gene_name.upper())

    return full_scores_df[mask].copy()


def extract_interval_genes(full_scores_df: pd.DataFrame) -> pd.DataFrame:
    if full_scores_df.empty:
        return pd.DataFrame(columns=["gene_id", "gene_name"])

    cols = [col for col in ["gene_id", "gene_name"] if col in full_scores_df.columns]
    if not cols:
        return pd.DataFrame(columns=["gene_id", "gene_name"])

    genes_df = full_scores_df[cols].drop_duplicates().sort_values(cols).reset_index(drop=True)
    for col in ["gene_id", "gene_name"]:
        if col not in genes_df.columns:
            genes_df[col] = ""
    return genes_df[["gene_id", "gene_name"]]


def filter_scores_to_target_context(
    full_scores_df: pd.DataFrame,
    *,
    gene_name: str,
    gene_id: str,
    data_source: str | None = None,
    gtex_tissue: str | None = None,
) -> pd.DataFrame:
    filtered_df = filter_scores_to_target_gene(
        full_scores_df,
        gene_name=gene_name,
        gene_id=gene_id,
    )
    if filtered_df.empty:
        return filtered_df

    if data_source is not None and "data_source" in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df["data_source"].astype(str).str.lower() == data_source.lower()
        ].copy()

    if gtex_tissue is not None and "gtex_tissue" in filtered_df.columns:
        filtered_df = filtered_df[
            filtered_df["gtex_tissue"].astype(str).str.lower() == gtex_tissue.lower()
        ].copy()

    return filtered_df


def aggregate_mean_raw_score_by_variant(filtered_scores_df: pd.DataFrame) -> pd.DataFrame:
    if filtered_scores_df.empty:
        return pd.DataFrame(
            columns=[
                "source_variant_id",
                "alphagenome_raw_mean",
                "alphagenome_quantile_mean",
                "alphagenome_track_count",
            ]
        )

    required_columns = {"source_variant_id", "raw_score", "quantile_score"}
    missing_columns = required_columns.difference(filtered_scores_df.columns)
    if missing_columns:
        raise ValueError(f"Filtered scores must include columns: {sorted(required_columns)}")

    aggregated_df = (
        filtered_scores_df.groupby("source_variant_id", as_index=False)
        .agg(
            alphagenome_raw_mean=("raw_score", "mean"),
            alphagenome_quantile_mean=("quantile_score", "mean"),
            alphagenome_track_count=("raw_score", "size"),
        )
        .sort_values("source_variant_id")
        .reset_index(drop=True)
    )
    return aggregated_df
