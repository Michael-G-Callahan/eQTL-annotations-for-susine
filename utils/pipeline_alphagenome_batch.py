from __future__ import annotations

import gc
import shutil
import time
from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pyarrow.parquet as pq
from alphagenome.models import dna_client, variant_scorers

from .alphagenome import (
    INELIGIBLE_MODE,
    MIDPOINT_CENTERED_MODE,
    TSS_CENTERED_MODE,
    build_interval_covering_span,
    extract_interval_genes,
    filter_scores_to_target_context,
    gene_exon_span_1based,
    interval_bounds_1based,
    resolve_sequence_length_bp,
    score_batch,
    tss_radius_bp,
)
from .annotations import download_gtf_if_needed, get_gene_info, load_gene_annotations_for_gene
from .output_layout import LocusOutputPaths
from .paths import ProjectPaths
from .pipeline_defaults import DEFAULT_ALPHAGENOME_BATCH_SIZE
from .variant_processing import annotate_variant_window_eligibility, ensure_variant_columns, sanitize_for_parquet


def estimate_alphagenome_workload(
    locus_cfg: dict,
    paths: ProjectPaths,
    locus_paths: LocusOutputPaths,
    *,
    max_variants: int | None = None,
) -> dict:
    raw_variants_df, gene_info = _load_annotation_inputs(locus_cfg, paths, locus_paths)
    eligibility_df = _build_eligibility_df(
        raw_variants_df,
        gene_info,
        locus_cfg["alphagenome_sequence_length"],
        batch_size=int(locus_cfg.get("alphagenome_batch_size", DEFAULT_ALPHAGENOME_BATCH_SIZE)),
        max_variants=max_variants,
    )
    counts = _summarize_eligibility(eligibility_df)
    return {
        "locus_id": locus_cfg["locus_id"],
        "gene_name": locus_cfg["gene_name"],
        "variants_loaded": len(raw_variants_df),
        **counts,
        "tss_radius_bp": tss_radius_bp(locus_cfg["alphagenome_sequence_length"]),
        "sequence_length_bp": resolve_sequence_length_bp(locus_cfg["alphagenome_sequence_length"]),
    }


def run_alphagenome_annotation(
    locus_cfg: dict,
    paths: ProjectPaths,
    locus_paths: LocusOutputPaths,
    api_key: str,
    *,
    max_variants: int | None = None,
    write_batch_checkpoints: bool = True,
    verbose: bool = False,
    progress_callback: Callable[[dict], None] | None = None,
) -> dict:
    raw_variants_df, gene_info = _load_annotation_inputs(locus_cfg, paths, locus_paths)
    if gene_info["exons"].empty:
        raise ValueError(f"Target gene {locus_cfg['gene_name']} has no exon coordinates; cannot build midpoint windows.")

    variants_loaded_count = len(raw_variants_df)
    sequence_length = locus_cfg["alphagenome_sequence_length"]
    sequence_length_bp = resolve_sequence_length_bp(sequence_length)
    eligibility_df = _build_eligibility_df(
        raw_variants_df,
        gene_info,
        sequence_length,
        batch_size=int(locus_cfg.get("alphagenome_batch_size", DEFAULT_ALPHAGENOME_BATCH_SIZE)),
        max_variants=max_variants,
    )
    eligibility_path = locus_paths.alphagenome_variant_window_eligibility_csv(locus_cfg["gene_name"])
    _write_eligibility_csv(eligibility_df, eligibility_path)
    counts = _summarize_eligibility(eligibility_df)

    eligible_variants_df = eligibility_df[~eligibility_df["excluded_from_scoring"]].copy()
    if eligible_variants_df.empty:
        reason_counts = (
            eligibility_df["exclusion_reason"].fillna("unknown").value_counts(dropna=False).to_dict()
            if not eligibility_df.empty
            else {}
        )
        raise ValueError(
            "No variants remain after AlphaGenome interval eligibility filtering. "
            f"Reasons: {reason_counts}"
        )

    if verbose:
        print(f"Gene: {gene_info['gene_name']}")
        print(f"  Location: {gene_info['chrom']}:{gene_info['start']}-{gene_info['end']}")
        print(f"  Strand: {gene_info['strand']}")
        print(f"  TSS: {gene_info['tss']}")
        print(
            "  Variants loaded: "
            f"{variants_loaded_count}; tss_centered={counts['variants_tss_centered']}; "
            f"midpoint_centered={counts['variants_midpoint_centered']}; "
            f"ineligible={counts['variants_ineligible']}"
        )
        if max_variants is not None:
            print(f"  MAX_VARIANTS applied: {int(max_variants)}")

    dna_model = dna_client.create(api_key)
    target_output = dna_client.OutputType.RNA_SEQ
    active_variant_scorers = [variant_scorers.GeneMaskLFCScorer(requested_output=target_output)]

    preflight_variant_df = _select_preflight_variants(eligible_variants_df)
    preflight_scoring_mode = preflight_variant_df.iloc[0]["scoring_mode"]
    preflight_full_scores_df, _, _, _ = score_batch(
        dna_model,
        preflight_variant_df,
        gene_info=gene_info,
        gene_name=locus_cfg["gene_name"],
        gene_id=locus_cfg["gene_id"],
        sequence_length_label=sequence_length,
        active_variant_scorers=active_variant_scorers,
        max_workers=1,
        progress_bar=False,
        retry_wait_seconds=int(locus_cfg["alphagenome_retry_wait_seconds"]),
    )
    interval_genes_df = extract_interval_genes(preflight_full_scores_df)
    interval_genes_path = locus_paths.interval_genes_csv(locus_cfg["gene_name"])
    interval_genes_df.to_csv(interval_genes_path, index=False)
    if verbose:
        print(f"Wrote interval gene manifest: {interval_genes_path}")
        print(
            f"Genes returned in {sequence_length} {preflight_scoring_mode} preflight interval: "
            f"{len(interval_genes_df)}"
        )

    gene_id_base = locus_cfg["gene_id"].split(".", 1)[0]
    target_present = (
        interval_genes_df["gene_name"].astype(str).str.upper().eq(locus_cfg["gene_name"].upper()).any()
        or interval_genes_df["gene_id"].astype(str).eq(locus_cfg["gene_id"]).any()
        or interval_genes_df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True).eq(gene_id_base).any()
    )
    if not target_present:
        raise ValueError(
            f"Target gene {locus_cfg['gene_name']} ({locus_cfg['gene_id']}) was not found in interval genes."
        )

    preflight_filtered_scores_df = filter_scores_to_target_context(
        preflight_full_scores_df,
        gene_name=locus_cfg["gene_name"],
        gene_id=locus_cfg["gene_id"],
        data_source=locus_cfg["alphagenome_target_data_source"],
        gtex_tissue=locus_cfg["alphagenome_target_gtex_tissue"],
    )
    if verbose:
        print(
            f"Rows for {locus_cfg['gene_name']} after "
            f"{locus_cfg['alphagenome_target_data_source']}/{locus_cfg['alphagenome_target_gtex_tissue']} "
            f"filters on one preflight variant: {len(preflight_filtered_scores_df)}"
        )
    if preflight_filtered_scores_df.empty:
        raise ValueError(
            f"No rows remained for {locus_cfg['gene_name']} after filtering to "
            f"data_source={locus_cfg['alphagenome_target_data_source']!r} and "
            f"gtex_tissue={locus_cfg['alphagenome_target_gtex_tissue']!r}."
        )

    shard_dir = locus_paths.alphagenome_batch_shard_dir(locus_cfg["gene_name"])
    if write_batch_checkpoints:
        if shard_dir.exists():
            shutil.rmtree(shard_dir)
        shard_dir.mkdir(parents=True, exist_ok=True)

    aggregate_state = {}
    variant_records = []
    filtered_row_count = 0
    written_shard_count = 0
    filtered_batches: list[pd.DataFrame] = []

    if verbose:
        print(
            f"Scoring {len(eligible_variants_df)} variants "
            f"(tss_centered={counts['variants_tss_centered']}, "
            f"midpoint_centered={counts['variants_midpoint_centered']}) "
            f"with sequence length {sequence_length}"
        )

    scoring_batches = list(
        _chunk_variants_for_scoring(eligible_variants_df, int(locus_cfg["alphagenome_batch_size"]))
    )
    total_batches = len(scoring_batches)
    total_variants_to_score = len(eligible_variants_df)
    variants_scored_so_far = 0
    locus_started = time.perf_counter()

    for batch_idx, batch_df in enumerate(scoring_batches, start=1):
        if verbose:
            batch_modes = batch_df["scoring_mode"].value_counts().to_dict()
            print(
                f"[{locus_cfg['gene_name']}] Starting batch {batch_idx}/{total_batches}: "
                f"{len(batch_df)} variants; modes={batch_modes}",
                flush=True,
            )
        batch_full_scores_df, batch_scored_variants_df, batch_elapsed_seconds, batch_exec_info = score_batch(
            dna_model,
            batch_df,
            gene_info=gene_info,
            gene_name=locus_cfg["gene_name"],
            gene_id=locus_cfg["gene_id"],
            sequence_length_label=sequence_length,
            active_variant_scorers=active_variant_scorers,
            max_workers=int(locus_cfg["alphagenome_max_workers"]),
            progress_bar=False,
            retry_wait_seconds=int(locus_cfg["alphagenome_retry_wait_seconds"]),
        )
        variant_records.extend(batch_scored_variants_df.to_dict("records"))
        filtered_batch_scores_df = filter_scores_to_target_context(
            batch_full_scores_df,
            gene_name=locus_cfg["gene_name"],
            gene_id=locus_cfg["gene_id"],
            data_source=locus_cfg["alphagenome_target_data_source"],
            gtex_tissue=locus_cfg["alphagenome_target_gtex_tissue"],
        )
        filtered_row_count += len(filtered_batch_scores_df)
        _update_aggregate_state(aggregate_state, filtered_batch_scores_df)
        write_started = pd.Timestamp.utcnow()
        shard_path = None
        if write_batch_checkpoints:
            shard_path = _write_filtered_batch_shard(filtered_batch_scores_df, shard_dir, batch_idx)
            if shard_path is not None:
                written_shard_count += 1
        else:
            filtered_batches.append(filtered_batch_scores_df.copy())
        write_elapsed_seconds = (pd.Timestamp.utcnow() - write_started).total_seconds()
        if verbose:
            variants_scored_so_far += len(batch_df)
            locus_elapsed_seconds = time.perf_counter() - locus_started
            progress = _format_locus_progress(
                variants_done=variants_scored_so_far,
                variants_total=total_variants_to_score,
                elapsed_seconds=locus_elapsed_seconds,
            )
            print(
                f"[{locus_cfg['gene_name']}] Finished batch {batch_idx}/{total_batches}; "
                f"{progress}; API scoring={batch_elapsed_seconds:.2f}s "
                f"({batch_elapsed_seconds / len(batch_df):.3f}s/variant); "
                f"parallel_used={batch_exec_info['used_parallel']}; "
                f"interval_groups={batch_exec_info['interval_group_count']}; "
                f"largest_group={batch_exec_info['largest_interval_group_size']}; "
                f"retries={batch_exec_info['retry_count']}",
                flush=True,
            )
            if batch_exec_info["last_retry_reason"]:
                print(
                    f"[{locus_cfg['gene_name']}] RATE-LIMIT/QUOTA retry observed; "
                    f"waited {int(locus_cfg['alphagenome_retry_wait_seconds'])}s and retried. "
                    f"Reason: {_shorten_message(batch_exec_info['last_retry_reason'])}",
                    flush=True,
                )
            print(
                f"[{locus_cfg['gene_name']}] Batch write={write_elapsed_seconds:.2f}s; "
                f"filtered_rows={len(filtered_batch_scores_df)}; "
                f"shard_written={shard_path is not None if write_batch_checkpoints else False}",
                flush=True,
            )
        else:
            variants_scored_so_far += len(batch_df)
        if progress_callback is not None:
            progress_callback(
                {
                    "event": "batch_complete",
                    "locus_id": locus_cfg["locus_id"],
                    "gene_name": locus_cfg["gene_name"],
                    "batch_idx": batch_idx,
                    "total_batches": total_batches,
                    "batch_n_variants": len(batch_df),
                    "variants_done": variants_scored_so_far,
                    "variants_total": total_variants_to_score,
                    "locus_elapsed_seconds": time.perf_counter() - locus_started,
                    "batch_elapsed_seconds": batch_elapsed_seconds,
                    "retry_count": int(batch_exec_info["retry_count"]),
                    "last_retry_reason": batch_exec_info["last_retry_reason"],
                }
            )
        del batch_full_scores_df, batch_scored_variants_df, filtered_batch_scores_df
        gc.collect()

    filtered_scores_path = locus_paths.alphagenome_filtered_scores_parquet(locus_cfg["gene_name"])
    if write_batch_checkpoints:
        final_filtered_row_count = _finalize_filtered_scores(shard_dir, filtered_scores_path)
    else:
        if filtered_batches:
            filtered_scores_df = pd.concat(filtered_batches, ignore_index=True)
        else:
            filtered_scores_df = pd.DataFrame()
        sanitize_for_parquet(filtered_scores_df).to_parquet(filtered_scores_path, index=False)
        final_filtered_row_count = len(filtered_scores_df)
    aggregated_scores_df = _finalize_aggregated_scores(aggregate_state)
    aggregated_scores_path = locus_paths.alphagenome_variant_scores_csv(locus_cfg["gene_name"])
    aggregated_scores_df.to_csv(aggregated_scores_path, index=False)

    histogram_path = locus_paths.alphagenome_histogram_png(locus_cfg["gene_name"])
    trimmed_histogram_path = locus_paths.alphagenome_trimmed_histogram_png(locus_cfg["gene_name"])
    _write_histograms(aggregated_scores_df, locus_cfg["gene_name"], histogram_path, trimmed_histogram_path)

    scored_variants_df = pd.DataFrame(variant_records)
    if verbose:
        print("Scoring complete.")
        print(f"Variant records: {len(scored_variants_df)}")
        print(f"Filtered rows accumulated across batches: {filtered_row_count}")
        print(f"Batch shard files written: {written_shard_count if write_batch_checkpoints else 0}")
        print(f"Wrote: {filtered_scores_path}")
        print(f"Wrote: {aggregated_scores_path}")
        print(f"Wrote: {eligibility_path}")
    return {
        "locus_id": locus_cfg["locus_id"],
        "variants_loaded": variants_loaded_count,
        "variants_within_tss_radius": counts["variants_tss_centered"],
        "variants_tss_centered": counts["variants_tss_centered"],
        "variants_midpoint_centered": counts["variants_midpoint_centered"],
        "variants_ineligible": counts["variants_ineligible"],
        "variants_to_score": counts["variants_to_score"],
        "sequence_length_bp": sequence_length_bp,
        "preflight_scoring_mode": preflight_scoring_mode,
        "variants_attempted": len(scored_variants_df),
        "variants_scored_ok": int((scored_variants_df["score_status"] == "scored").sum())
        if not scored_variants_df.empty
        else 0,
        "variants_api_error": int((scored_variants_df["score_status"] == "api_error").sum())
        if not scored_variants_df.empty
        else 0,
        "filtered_score_rows": final_filtered_row_count,
        "aggregated_variants": len(aggregated_scores_df),
        "filtered_scores_parquet": str(filtered_scores_path),
        "variant_scores_csv": str(aggregated_scores_path),
        "histogram_png": str(histogram_path),
        "trimmed_histogram_png": str(trimmed_histogram_path),
        "interval_genes_csv": str(interval_genes_path),
        "variant_window_eligibility_csv": str(eligibility_path),
        "write_batch_checkpoints": bool(write_batch_checkpoints),
        "batch_shard_count": written_shard_count if write_batch_checkpoints else 0,
    }


def _format_seconds(seconds: float | int | None) -> str:
    if seconds is None or not pd.notna(seconds) or seconds < 0:
        return "unknown"
    seconds = int(round(float(seconds)))
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h{minutes:02d}m{secs:02d}s"
    if minutes:
        return f"{minutes}m{secs:02d}s"
    return f"{secs}s"


def _format_locus_progress(*, variants_done: int, variants_total: int, elapsed_seconds: float) -> str:
    pct = 100.0 * variants_done / variants_total if variants_total else 100.0
    rate = variants_done / elapsed_seconds if elapsed_seconds > 0 else 0.0
    eta_seconds = (variants_total - variants_done) / rate if rate > 0 else None
    return (
        f"{variants_done}/{variants_total} variants ({pct:.1f}%); "
        f"elapsed={_format_seconds(elapsed_seconds)}; "
        f"locus_eta={_format_seconds(eta_seconds)}"
    )


def _shorten_message(message: str, max_chars: int = 220) -> str:
    message = " ".join(str(message).split())
    if len(message) <= max_chars:
        return message
    return message[: max_chars - 3] + "..."


def load_existing_annotation_result(locus_cfg: dict, locus_paths: LocusOutputPaths) -> dict:
    scores_df = pd.read_csv(locus_paths.alphagenome_variant_scores_csv(locus_cfg["gene_name"]))
    filtered_rows = pq.read_table(
        locus_paths.alphagenome_filtered_scores_parquet(locus_cfg["gene_name"])
    ).num_rows
    counts = {
        "variants_tss_centered": pd.NA,
        "variants_midpoint_centered": pd.NA,
        "variants_ineligible": pd.NA,
        "variants_to_score": pd.NA,
        "sequence_length_bp": pd.NA,
        "preflight_scoring_mode": pd.NA,
    }
    eligibility_path = locus_paths.alphagenome_variant_window_eligibility_csv(locus_cfg["gene_name"])
    if eligibility_path.exists():
        eligibility_df = pd.read_csv(eligibility_path)
        counts.update(_summarize_eligibility(eligibility_df))
    return {
        "locus_id": locus_cfg["locus_id"],
        "variants_loaded": pd.NA,
        "variants_within_tss_radius": counts["variants_tss_centered"],
        **counts,
        "variants_attempted": pd.NA,
        "variants_scored_ok": pd.NA,
        "variants_api_error": pd.NA,
        "filtered_score_rows": filtered_rows,
        "aggregated_variants": len(scores_df),
        "filtered_scores_parquet": str(locus_paths.alphagenome_filtered_scores_parquet(locus_cfg["gene_name"])),
        "variant_scores_csv": str(locus_paths.alphagenome_variant_scores_csv(locus_cfg["gene_name"])),
        "histogram_png": str(locus_paths.alphagenome_histogram_png(locus_cfg["gene_name"])),
        "trimmed_histogram_png": str(
            locus_paths.alphagenome_trimmed_histogram_png(locus_cfg["gene_name"])
        ),
        "interval_genes_csv": str(locus_paths.interval_genes_csv(locus_cfg["gene_name"])),
        "variant_window_eligibility_csv": str(eligibility_path),
    }


def _load_annotation_inputs(
    locus_cfg: dict,
    paths: ProjectPaths,
    locus_paths: LocusOutputPaths,
) -> tuple[pd.DataFrame, dict]:
    variant_source_file = locus_paths.phase1_master_variants_csv(locus_cfg["gene_name"])
    if not variant_source_file.exists():
        raise FileNotFoundError(
            f"Missing phase-1 master variants for locus {locus_cfg['locus_id']} "
            f"({locus_cfg['gene_name']}) at {variant_source_file}. "
            "This usually means the locus was only run in screening mode. "
            "AlphaGenome annotation requires phase-1 processing outputs under output/ld."
        )
    gtf_path = download_gtf_if_needed(paths.gtf_cache, genome=locus_cfg["reference_genome"])
    genes_df, exons_df, _ = load_gene_annotations_for_gene(
        gtf_path, locus_cfg["gene_name"], paths.gtf_shortcuts
    )
    gene_info = get_gene_info(locus_cfg["gene_name"], genes_df, exons_df)

    raw_variants_df = pd.read_csv(variant_source_file)
    raw_variants_df = ensure_variant_columns(raw_variants_df)
    raw_variants_df = raw_variants_df[["variant_id", "chrom", "pos", "ref", "alt", "z_score"]].copy()
    return raw_variants_df, gene_info


def _build_eligibility_df(
    raw_variants_df: pd.DataFrame,
    gene_info: dict,
    sequence_length: str,
    *,
    batch_size: int,
    max_variants: int | None = None,
) -> pd.DataFrame:
    working_df = raw_variants_df.copy()
    if max_variants is not None:
        working_df = working_df.head(int(max_variants)).copy()
    eligibility_df = annotate_variant_window_eligibility(working_df, gene_info, sequence_length)
    eligibility_df = _coalesce_midpoint_intervals(
        eligibility_df,
        gene_info,
        sequence_length,
        batch_size=batch_size,
    )
    return _sort_variants_for_scoring(eligibility_df)


def _summarize_eligibility(eligibility_df: pd.DataFrame) -> dict:
    if eligibility_df.empty:
        return {
            "variants_tss_centered": 0,
            "variants_midpoint_centered": 0,
            "variants_ineligible": 0,
            "variants_to_score": 0,
        }
    mode_counts = eligibility_df["scoring_mode"].value_counts().to_dict()
    return {
        "variants_tss_centered": int(mode_counts.get(TSS_CENTERED_MODE, 0)),
        "variants_midpoint_centered": int(mode_counts.get(MIDPOINT_CENTERED_MODE, 0)),
        "variants_ineligible": int(mode_counts.get(INELIGIBLE_MODE, 0)),
        "variants_to_score": int(mode_counts.get(TSS_CENTERED_MODE, 0) + mode_counts.get(MIDPOINT_CENTERED_MODE, 0)),
    }


def _select_preflight_variants(eligible_variants_df: pd.DataFrame) -> pd.DataFrame:
    for mode in (TSS_CENTERED_MODE, MIDPOINT_CENTERED_MODE):
        mode_df = eligible_variants_df[eligible_variants_df["scoring_mode"] == mode]
        if not mode_df.empty:
            return mode_df.head(1).copy()
    raise ValueError("No eligible variants were available for AlphaGenome preflight.")


def _sort_variants_for_scoring(eligibility_df: pd.DataFrame) -> pd.DataFrame:
    if eligibility_df.empty:
        return eligibility_df

    mode_rank = {
        TSS_CENTERED_MODE: 0,
        MIDPOINT_CENTERED_MODE: 1,
        INELIGIBLE_MODE: 2,
    }
    sorted_df = eligibility_df.copy()
    sorted_df["_mode_rank"] = sorted_df["scoring_mode"].map(mode_rank).fillna(99).astype(int)
    sorted_df = sorted_df.sort_values(
        by=[
            "_mode_rank",
            "scoring_interval_start_1based",
            "scoring_interval_end_1based",
            "variant_min_1based",
            "variant_id",
        ],
        kind="stable",
        na_position="last",
    ).drop(columns="_mode_rank")
    return sorted_df.reset_index(drop=True)


def _coalesce_midpoint_intervals(
    eligibility_df: pd.DataFrame,
    gene_info: dict,
    sequence_length: str,
    *,
    batch_size: int,
) -> pd.DataFrame:
    if eligibility_df.empty:
        return eligibility_df

    midpoint_mask = eligibility_df["scoring_mode"] == MIDPOINT_CENTERED_MODE
    midpoint_df = eligibility_df[midpoint_mask].copy()
    if midpoint_df.empty:
        return eligibility_df

    sequence_length_bp = resolve_sequence_length_bp(sequence_length)
    exon_min_1based, exon_max_1based = gene_exon_span_1based(gene_info)
    midpoint_df = midpoint_df.sort_values(["variant_min_1based", "variant_max_1based", "variant_id"], kind="stable")

    bucket_assignments: list[tuple[list[int], object, int, int]] = []
    current_indices: list[int] = []
    current_span_min: int | None = None
    current_span_max: int | None = None

    for row in midpoint_df.itertuples():
        proposed_span_min = row.variant_min_1based if current_span_min is None else min(current_span_min, row.variant_min_1based)
        proposed_span_max = row.variant_max_1based if current_span_max is None else max(current_span_max, row.variant_max_1based)
        required_min = min(exon_min_1based, proposed_span_min)
        required_max = max(exon_max_1based, proposed_span_max)
        required_span_bp = required_max - required_min + 1

        if current_indices and (len(current_indices) >= batch_size or required_span_bp > sequence_length_bp):
            shared_interval = build_interval_covering_span(
                gene_info["chrom"],
                min(exon_min_1based, current_span_min),
                max(exon_max_1based, current_span_max),
                sequence_length,
                gene_info.get("strand", "."),
                name=f"{gene_info['gene_name']}_midpoint_coalesced",
            )
            interval_start, interval_end = interval_bounds_1based(shared_interval)
            bucket_assignments.append((current_indices.copy(), shared_interval, interval_start, interval_end))
            current_indices = []
            current_span_min = None
            current_span_max = None
            proposed_span_min = row.variant_min_1based
            proposed_span_max = row.variant_max_1based

        current_indices.append(row.Index)
        current_span_min = proposed_span_min
        current_span_max = proposed_span_max

    if current_indices:
        shared_interval = build_interval_covering_span(
            gene_info["chrom"],
            min(exon_min_1based, current_span_min),
            max(exon_max_1based, current_span_max),
            sequence_length,
            gene_info.get("strand", "."),
            name=f"{gene_info['gene_name']}_midpoint_coalesced",
        )
        interval_start, interval_end = interval_bounds_1based(shared_interval)
        bucket_assignments.append((current_indices.copy(), shared_interval, interval_start, interval_end))

    updated_df = eligibility_df.copy()
    for indices, shared_interval, interval_start, interval_end in bucket_assignments:
        updated_df.loc[indices, "scoring_interval"] = shared_interval
        updated_df.loc[indices, "scoring_interval_start_1based"] = interval_start
        updated_df.loc[indices, "scoring_interval_end_1based"] = interval_end
    return updated_df


def _write_eligibility_csv(eligibility_df: pd.DataFrame, output_path: Path) -> None:
    export_df = eligibility_df.copy()
    if "scoring_interval" in export_df.columns:
        export_df = export_df.drop(columns=["scoring_interval"])
    export_df.to_csv(output_path, index=False)


def _chunked(df: pd.DataFrame, chunk_size: int):
    for start in range(0, len(df), chunk_size):
        yield df.iloc[start : start + chunk_size].copy()


def _chunk_variants_for_scoring(df: pd.DataFrame, chunk_size: int):
    if df.empty:
        return

    group_columns = [
        "scoring_mode",
        "scoring_interval_start_1based",
        "scoring_interval_end_1based",
    ]
    sorted_df = df.reset_index(drop=True).copy()

    group_start = 0
    groups: list[pd.DataFrame] = []
    for row_idx in range(1, len(sorted_df) + 1):
        if row_idx == len(sorted_df):
            groups.append(sorted_df.iloc[group_start:row_idx].copy())
            break
        current_key = tuple(sorted_df.loc[row_idx - 1, group_columns].tolist())
        next_key = tuple(sorted_df.loc[row_idx, group_columns].tolist())
        if current_key != next_key:
            groups.append(sorted_df.iloc[group_start:row_idx].copy())
            group_start = row_idx

    pending_groups: list[pd.DataFrame] = []
    pending_size = 0

    def flush_pending():
        nonlocal pending_groups, pending_size
        if not pending_groups:
            return None
        batch_df = pd.concat(pending_groups, ignore_index=True)
        pending_groups = []
        pending_size = 0
        return batch_df

    for group_df in groups:
        if len(group_df) > chunk_size:
            flushed = flush_pending()
            if flushed is not None:
                yield flushed
            for split_df in _chunked(group_df, chunk_size):
                yield split_df.reset_index(drop=True)
            continue

        if pending_size + len(group_df) > chunk_size:
            flushed = flush_pending()
            if flushed is not None:
                yield flushed

        pending_groups.append(group_df)
        pending_size += len(group_df)

    flushed = flush_pending()
    if flushed is not None:
        yield flushed


def _write_filtered_batch_shard(filtered_scores_df: pd.DataFrame, shard_dir: Path, batch_idx: int) -> Path | None:
    if filtered_scores_df.empty:
        return None
    shard_path = shard_dir / f"part-{batch_idx:04d}.parquet"
    sanitize_for_parquet(filtered_scores_df).to_parquet(shard_path, index=False)
    return shard_path


def _update_aggregate_state(aggregate_state: dict[str, dict[str, float]], filtered_scores_df: pd.DataFrame) -> None:
    if filtered_scores_df.empty:
        return
    batch_stats = (
        filtered_scores_df.groupby("source_variant_id", as_index=False)
        .agg(
            raw_score_sum=("raw_score", "sum"),
            quantile_score_sum=("quantile_score", "sum"),
            raw_score_count=("raw_score", "size"),
        )
    )
    for row in batch_stats.itertuples(index=False):
        state = aggregate_state.setdefault(
            row.source_variant_id,
            {"raw_score_sum": 0.0, "quantile_score_sum": 0.0, "raw_score_count": 0},
        )
        state["raw_score_sum"] += float(row.raw_score_sum)
        state["quantile_score_sum"] += float(row.quantile_score_sum)
        state["raw_score_count"] += int(row.raw_score_count)


def _finalize_aggregated_scores(aggregate_state: dict[str, dict[str, float]]) -> pd.DataFrame:
    records = []
    for variant_id, state in aggregate_state.items():
        count = int(state["raw_score_count"])
        mean_score = float(state["raw_score_sum"]) / count if count else float("nan")
        mean_quantile_score = float(state["quantile_score_sum"]) / count if count else float("nan")
        records.append(
            {
                "source_variant_id": variant_id,
                "alphagenome_raw_mean": mean_score,
                "alphagenome_quantile_mean": mean_quantile_score,
                "alphagenome_track_count": count,
            }
        )
    if not records:
        return pd.DataFrame(
            columns=[
                "source_variant_id",
                "alphagenome_raw_mean",
                "alphagenome_quantile_mean",
                "alphagenome_track_count",
            ]
        )
    return pd.DataFrame(records).sort_values("source_variant_id").reset_index(drop=True)


def _finalize_filtered_scores(shard_dir: Path, output_path: Path) -> int:
    shard_paths = sorted(shard_dir.glob("part-*.parquet"))
    if not shard_paths:
        pd.DataFrame().to_parquet(output_path, index=False)
        return 0
    writer = None
    total_rows = 0
    try:
        for shard_path in shard_paths:
            table = pq.read_table(shard_path)
            total_rows += table.num_rows
            if writer is None:
                writer = pq.ParquetWriter(output_path, table.schema, compression="zstd")
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    return total_rows


def _write_histograms(
    aggregated_scores_df: pd.DataFrame,
    gene_name: str,
    histogram_path: Path,
    trimmed_histogram_path: Path,
) -> None:
    if aggregated_scores_df.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(aggregated_scores_df["alphagenome_raw_mean"], bins=40, color="#2f6c8f", edgecolor="white")
    ax.set_title(f"{gene_name} AlphaGenome Mean Raw Score Distribution")
    ax.set_xlabel("Mean raw score")
    ax.set_ylabel("Variant count")
    fig.tight_layout()
    fig.savefig(histogram_path, dpi=150)
    plt.close(fig)

    trimmed_scores_df = aggregated_scores_df[
        ~aggregated_scores_df["alphagenome_raw_mean"].between(-0.002, 0.002, inclusive="both")
    ].copy()
    fig, ax = plt.subplots(figsize=(8, 5))
    if trimmed_scores_df.empty:
        ax.text(0.5, 0.5, "No variants with |score| > 0.002", ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.hist(trimmed_scores_df["alphagenome_raw_mean"], bins=40, color="#b85c38", edgecolor="white")
        ax.set_xlabel("Mean raw score")
        ax.set_ylabel("Variant count")
    ax.set_title(f"{gene_name} AlphaGenome Mean Raw Score Distribution (|score| > 0.002)")
    fig.tight_layout()
    fig.savefig(trimmed_histogram_path, dpi=150)
    plt.close(fig)
