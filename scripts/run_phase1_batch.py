from __future__ import annotations

import argparse
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils.locus_manifest import iter_enabled_loci, load_loci_manifest, validate_loci_manifest
from utils.output_layout import get_locus_output_paths
from utils.paths import add_project_root_to_sys_path, configure_runtime_env, find_project_root
from utils.pipeline_ld import load_existing_phase1_result, run_ld_phase1
from utils.pipeline_zscore import run_zscore_export


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 1 z-score + LD pipeline for manifest loci.")
    parser.add_argument("--manifest", required=True, help="Path to loci_manifest.csv")
    parser.add_argument("--locus-id", action="append", dest="locus_ids", help="Optional locus_id filter")
    parser.add_argument("--force", action="store_true", help="Rerun loci even if outputs already exist")
    parser.add_argument("--fail-fast", action="store_true", help="Stop on the first locus failure")
    parser.add_argument("--max-loci", type=int, default=None, help="Limit number of loci processed")
    parser.add_argument(
        "--run-mode",
        choices=("screening", "processing"),
        default="processing",
        help="Use 'screening' to emit only aggregate prelim CSVs, or 'processing' to save all phase-1 outputs.",
    )
    return parser.parse_args()


def main() -> int:
    project_root = add_project_root_to_sys_path(find_project_root())
    args = parse_args()
    dir_keys = ("output", "output_prelim") if args.run_mode == "screening" else None
    paths = configure_runtime_env(project_root, dir_keys=dir_keys)

    manifest_df = load_loci_manifest(args.manifest)
    validate_loci_manifest(manifest_df)

    if args.locus_ids:
        manifest_df = manifest_df[manifest_df["locus_id"].isin(set(args.locus_ids))].copy()

    loci = list(iter_enabled_loci(manifest_df))
    if args.max_loci is not None:
        loci = loci[: args.max_loci]

    summary_rows = []
    phase1_summary_path = paths.output_prelim / "phase1_batch_summary.csv"
    metrics_summary_path = paths.output_prelim / "phase1_dataset_metrics_all_loci.csv"

    for locus_cfg in loci:
        row = {
            "locus_id": locus_cfg["locus_id"],
            "gene_name": locus_cfg["gene_name"],
            "gene_id": locus_cfg["gene_id"],
            "gtex_tissue": locus_cfg["gtex_tissue"],
            "gtex_chrom": locus_cfg["gtex_chrom"],
            "status": "failed",
            "error_message": "",
            "stage": "",
        }
        started_at = datetime.now(timezone.utc)
        row["started_at"] = started_at.isoformat()
        try:
            locus_paths = get_locus_output_paths(
                paths,
                locus_cfg["locus_id"],
                ensure_dirs=args.run_mode == "processing",
            )
            if locus_paths.phase1_complete(locus_cfg["gene_name"]) and not args.force:
                row.update(load_existing_phase1_result(locus_cfg, locus_paths))
                row["status"] = "skipped_existing"
                row["stage"] = "phase1"
            else:
                row["stage"] = "zscore"
                zscore_result = run_zscore_export(
                    locus_cfg,
                    paths,
                    locus_paths,
                    write_stage_outputs=args.run_mode == "processing",
                    write_prelim=args.run_mode == "processing",
                    allow_cache_writes=args.run_mode == "processing",
                )
                z_scores_df = zscore_result.pop("_z_scores_df")
                variants_df = zscore_result.pop("_vcf_variants_df")
                row.update(zscore_result)
                row["stage"] = "ld"
                row.update(
                    run_ld_phase1(
                        locus_cfg,
                        paths,
                        locus_paths,
                        z_scores_df=z_scores_df,
                        variants_df=variants_df,
                        write_stage_outputs=args.run_mode == "processing",
                        write_prelim=args.run_mode == "processing",
                        allow_cache_writes=args.run_mode == "processing",
                    )
                )
                row["status"] = "completed"
                row["stage"] = "phase1"
        except Exception as exc:  # noqa: BLE001
            row["error_message"] = f"{type(exc).__name__}: {exc}"
            row["traceback"] = traceback.format_exc()
            if args.fail_fast:
                summary_rows.append(_finalize_row(row, started_at))
                _write_phase1_summaries(
                    summary_rows,
                    phase1_summary_path,
                    metrics_summary_path,
                    write_status_summary=args.run_mode == "processing",
                )
                raise
        summary_rows.append(_finalize_row(row, started_at))
        _write_phase1_summaries(
            summary_rows,
            phase1_summary_path,
            metrics_summary_path,
            write_status_summary=args.run_mode == "processing",
        )

    return 0


def _finalize_row(row: dict, started_at: datetime) -> dict:
    finished_at = datetime.now(timezone.utc)
    row["finished_at"] = finished_at.isoformat()
    row["elapsed_seconds"] = round((finished_at - started_at).total_seconds(), 3)
    return row


def _write_phase1_summaries(
    rows: list[dict],
    summary_path: Path,
    metrics_path: Path,
    *,
    write_status_summary: bool,
) -> None:
    summary_df = pd.DataFrame(rows)
    if write_status_summary:
        summary_df.to_csv(summary_path, index=False)
    metric_cols = [
        "locus_id",
        "gene_name",
        "gene_id",
        "gtex_tissue",
        "gtex_chrom",
        "raw_gene_variants_loaded",
        "post_af_filter",
        "post_sample_size_filter",
        "z_score_variants_exported",
        "ld_matched_variants",
        "ld_unmatched_variants",
        "phase1_master_variants",
        "M1",
        "M2",
        "M4",
        "z_count_abs_gt_3",
        "z_eff_signals",
        "status",
        "error_message",
    ]
    metrics_df = summary_df[[col for col in metric_cols if col in summary_df.columns]].copy()
    metrics_df.to_csv(metrics_path, index=False)


if __name__ == "__main__":
    raise SystemExit(main())
