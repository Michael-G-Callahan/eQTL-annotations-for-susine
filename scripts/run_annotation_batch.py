from __future__ import annotations

import argparse
import os
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils.locus_manifest import (
    load_annotation_selection,
    load_loci_manifest,
    merge_annotation_selection,
    validate_annotation_selection,
    validate_loci_manifest,
)
from utils.output_layout import get_locus_output_paths
from utils.paths import add_project_root_to_sys_path, configure_runtime_env, find_project_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run AlphaGenome annotation for selected manifest loci.")
    parser.add_argument("--manifest", required=True, help="Path to loci_manifest.csv")
    parser.add_argument("--selection", required=True, help="Path to annotation_selection.csv")
    parser.add_argument(
        "--summary-path",
        default=None,
        help=(
            "Optional output CSV for the batch summary. Defaults to "
            "output/annotation/alphagenome/annotation_batch_summary.csv"
        ),
    )
    parser.add_argument("--api-key-env-var", default="ALPHAGENOME_API_KEY", help="Environment variable holding the AlphaGenome API key")
    parser.add_argument("--locus-id", action="append", dest="locus_ids", help="Optional locus_id filter")
    parser.add_argument("--force", action="store_true", help="Rerun loci even if outputs already exist")
    parser.add_argument("--fail-fast", action="store_true", help="Stop on the first locus failure")
    parser.add_argument("--quiet", action="store_true", help="Suppress per-locus and per-batch progress logging")
    return parser.parse_args()


def main() -> int:
    project_root = add_project_root_to_sys_path(find_project_root())
    paths = configure_runtime_env(project_root)
    args = parse_args()
    api_key = os.environ.get(args.api_key_env_var)
    if not api_key:
        raise RuntimeError(f"Missing AlphaGenome API key in environment variable {args.api_key_env_var}")

    from utils.pipeline_alphagenome_batch import load_existing_annotation_result, run_alphagenome_annotation

    manifest_df = load_loci_manifest(args.manifest)
    validate_loci_manifest(manifest_df)
    selection_df = load_annotation_selection(args.selection)
    validate_annotation_selection(selection_df, manifest_df)

    loci_df = merge_annotation_selection(manifest_df, selection_df)
    if args.locus_ids:
        loci_df = loci_df[loci_df["locus_id"].isin(set(args.locus_ids))].copy()

    summary_rows = []
    summary_path = (
        Path(args.summary_path)
        if args.summary_path
        else paths.output_annotation_alphagenome / "annotation_batch_summary.csv"
    )
    if not summary_path.is_absolute():
        summary_path = project_root / summary_path
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    total_loci = len(loci_df)
    batch_started = perf_counter()
    progress_enabled = not args.quiet
    if progress_enabled:
        print(
            f"AlphaGenome batch starting: {total_loci} locus/loci; "
            f"summary_path={summary_path}",
            flush=True,
        )
    for locus_idx, locus_cfg in enumerate(loci_df.to_dict("records"), start=1):
        row = {
            "locus_id": locus_cfg["locus_id"],
            "gene_name": locus_cfg["gene_name"],
            "gene_id": locus_cfg["gene_id"],
            "gtex_tissue": locus_cfg["gtex_tissue"],
            "gtex_chrom": locus_cfg["gtex_chrom"],
            "status": "failed",
            "error_message": "",
            "stage": "annotation",
        }
        started_at = datetime.now(timezone.utc)
        row["started_at"] = started_at.isoformat()
        locus_timer = perf_counter()
        if progress_enabled:
            print(
                f"[locus {locus_idx}/{total_loci}] Starting {locus_cfg['locus_id']} "
                f"({locus_cfg['gene_name']})",
                flush=True,
            )
        progress_callback = _build_progress_callback(
            locus_idx=locus_idx,
            total_loci=total_loci,
            batch_started=batch_started,
        ) if progress_enabled else None
        try:
            locus_paths = get_locus_output_paths(paths, locus_cfg["locus_id"])
            if not locus_paths.phase1_alphagenome_ready(locus_cfg["gene_name"]):
                raise FileNotFoundError(
                    f"AlphaGenome inputs are incomplete for locus {locus_cfg['locus_id']} ({locus_cfg['gene_name']})"
                )
            if locus_paths.annotation_complete(locus_cfg["gene_name"]) and not args.force:
                row.update(load_existing_annotation_result(locus_cfg, locus_paths))
                row["status"] = "skipped_existing"
            else:
                row.update(
                    run_alphagenome_annotation(
                        locus_cfg,
                        paths,
                        locus_paths,
                        api_key,
                        verbose=progress_enabled,
                        progress_callback=progress_callback,
                    )
                )
                row["status"] = "completed"
        except Exception as exc:  # noqa: BLE001
            row["error_message"] = f"{type(exc).__name__}: {exc}"
            row["traceback"] = traceback.format_exc()
            if args.fail_fast:
                summary_rows.append(_finalize_row(row, started_at))
                pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
                raise
        summary_rows.append(_finalize_row(row, started_at))
        pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
        if progress_enabled:
            completed_loci = len(summary_rows)
            elapsed = perf_counter() - batch_started
            per_locus = elapsed / completed_loci if completed_loci else 0.0
            remaining_loci = max(total_loci - completed_loci, 0)
            print(
                f"[locus {locus_idx}/{total_loci}] {row['status']} "
                f"{locus_cfg['locus_id']} ({locus_cfg['gene_name']}); "
                f"locus_elapsed={_format_seconds(perf_counter() - locus_timer)}; "
                f"overall_elapsed={_format_seconds(elapsed)}; "
                f"overall_eta={_format_seconds(per_locus * remaining_loci)}; "
                f"summary_rows={completed_loci}",
                flush=True,
            )
    return 0


def _finalize_row(row: dict, started_at: datetime) -> dict:
    finished_at = datetime.now(timezone.utc)
    row["finished_at"] = finished_at.isoformat()
    row["elapsed_seconds"] = round((finished_at - started_at).total_seconds(), 3)
    return row


def _build_progress_callback(*, locus_idx: int, total_loci: int, batch_started: float):
    def _callback(event: dict) -> None:
        if event.get("event") != "batch_complete":
            return
        locus_fraction = (
            float(event["variants_done"]) / float(event["variants_total"])
            if event.get("variants_total")
            else 1.0
        )
        completed_equiv_loci = (locus_idx - 1) + locus_fraction
        elapsed = perf_counter() - batch_started
        overall_eta = None
        if completed_equiv_loci > 0:
            overall_eta = elapsed * (total_loci - completed_equiv_loci) / completed_equiv_loci
        msg = (
            f"[overall {completed_equiv_loci:.2f}/{total_loci} loci] "
            f"overall_elapsed={_format_seconds(elapsed)}; "
            f"overall_eta={_format_seconds(overall_eta)}"
        )
        if event.get("retry_count", 0) > 0:
            msg += (
                f"; batch_retries={event['retry_count']}; "
                "rate_limit_or_quota_retry=yes"
            )
        print(msg, flush=True)

    return _callback


def _format_seconds(seconds) -> str:
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


if __name__ == "__main__":
    raise SystemExit(main())
