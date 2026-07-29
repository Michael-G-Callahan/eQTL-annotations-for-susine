from __future__ import annotations

import argparse
import csv
import fcntl
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils.locus_manifest import iter_enabled_loci, load_loci_manifest, validate_loci_manifest
from utils.output_layout import get_locus_output_paths
from utils.paths import add_project_root_to_sys_path, configure_runtime_env, find_project_root
from utils.pipeline_ld import run_ld_phase1
from utils.pipeline_zscore import run_zscore_export

METRIC_COLUMNS = [
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
    "started_at",
    "finished_at",
    "elapsed_seconds",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Phase 1 z-score + LD for one manifest chunk by task index.")
    parser.add_argument("--manifest", required=True, help="Path to loci_manifest.csv")
    parser.add_argument(
        "--task-id",
        type=int,
        required=True,
        help="1-based task index into the enabled manifest chunks (matches Slurm array task IDs).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1,
        help="Number of enabled loci to process per task.",
    )
    parser.add_argument("--force", action="store_true", help="Rerun the locus even if outputs already exist")
    parser.add_argument(
        "--run-mode",
        choices=("screening", "processing"),
        default="screening",
        help="Use 'screening' to emit only the aggregate prelim CSV, or 'processing' to save all phase-1 outputs.",
    )
    return parser.parse_args()


def main() -> int:
    project_root = add_project_root_to_sys_path(find_project_root())
    args = parse_args()
    dir_keys = ("output", "output_prelim") if args.run_mode == "screening" else None
    paths = configure_runtime_env(project_root, dir_keys=dir_keys)

    manifest_df = load_loci_manifest(args.manifest)
    validate_loci_manifest(manifest_df)
    loci = list(iter_enabled_loci(manifest_df))
    if args.chunk_size < 1:
        raise ValueError(f"chunk-size must be >= 1, got {args.chunk_size}")

    chunk_count = (len(loci) + args.chunk_size - 1) // args.chunk_size
    if args.task_id < 1 or args.task_id > chunk_count:
        raise IndexError(
            f"task-id {args.task_id} is out of range for {chunk_count} enabled chunks in {args.manifest}"
        )

    batch_metrics_path = paths.output_prelim / f"{Path(args.manifest).stem}_phase1_dataset_metrics.csv"
    chunk_start = (args.task_id - 1) * args.chunk_size
    chunk_loci = loci[chunk_start : chunk_start + args.chunk_size]

    if not chunk_loci:
        raise RuntimeError(f"No loci resolved for task-id {args.task_id} with chunk-size {args.chunk_size}")

    summaries = []
    overall_ok = True
    ensure_locus_dirs = args.run_mode == "processing"

    for locus_cfg in chunk_loci:
        locus_started_at = datetime.now(timezone.utc)
        locus_paths = get_locus_output_paths(paths, locus_cfg["locus_id"], ensure_dirs=ensure_locus_dirs)
        summary = {
            "task_id": args.task_id,
            "chunk_size": args.chunk_size,
            "locus_id": locus_cfg["locus_id"],
            "gene_name": locus_cfg["gene_name"],
            "gene_id": locus_cfg["gene_id"],
            "gtex_tissue": locus_cfg["gtex_tissue"],
            "gtex_chrom": locus_cfg["gtex_chrom"],
            "status": "failed",
            "stage": "",
            "error_message": "",
            "started_at": locus_started_at.isoformat(),
        }

        try:
            existing_row = _load_summary_row(batch_metrics_path, locus_cfg["locus_id"])
            if existing_row and not args.force:
                summary.update(existing_row)
                summary["status"] = "skipped_existing"
                summary["stage"] = "phase1"
            else:
                summary["stage"] = "zscore"
                zscore_result = run_zscore_export(
                    locus_cfg,
                    paths,
                    locus_paths,
                    write_prelim=False,
                    write_stage_outputs=args.run_mode == "processing",
                    allow_cache_writes=args.run_mode == "processing",
                )
                z_scores_df = zscore_result.pop("_z_scores_df")
                variants_df = zscore_result.pop("_vcf_variants_df")
                summary.update(zscore_result)
                summary["stage"] = "ld"
                summary.update(
                    run_ld_phase1(
                        locus_cfg,
                        paths,
                        locus_paths,
                        z_scores_df=z_scores_df,
                        variants_df=variants_df,
                        write_stage_outputs=args.run_mode == "processing",
                        write_prelim=False,
                        allow_cache_writes=args.run_mode == "processing",
                    )
                )
                summary["status"] = "completed"
                summary["stage"] = "phase1"
        except Exception as exc:  # noqa: BLE001
            summary["error_message"] = f"{type(exc).__name__}: {exc}"
            summary["traceback"] = traceback.format_exc()
            overall_ok = False
        finally:
            finished_at = datetime.now(timezone.utc)
            summary["finished_at"] = finished_at.isoformat()
            summary["elapsed_seconds"] = round((finished_at - locus_started_at).total_seconds(), 3)

        _append_summary_row(batch_metrics_path, summary)
        summaries.append(summary)

    print(json.dumps(summaries, indent=2))
    return 0 if overall_ok else 1


def _load_summary_row(summary_path: Path, locus_id: str) -> dict | None:
    if not summary_path.exists():
        return None
    try:
        with summary_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                if row.get("locus_id") == locus_id:
                    return row
            return None
    except Exception:
        return None


def _append_summary_row(summary_path: Path, summary: dict) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    row = {column: summary.get(column, "") for column in METRIC_COLUMNS}

    with summary_path.open("a+", newline="") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        existing_rows = []
        if handle.tell() == 0:
            pass
        reader = csv.DictReader(handle)
        if reader.fieldnames:
            existing_rows = list(reader)
        existing_rows = [existing for existing in existing_rows if existing.get("locus_id") != row["locus_id"]]
        existing_rows.append(row)
        existing_rows.sort(key=lambda existing: existing.get("locus_id", ""))

        handle.seek(0)
        handle.truncate(0)
        writer = csv.DictWriter(handle, fieldnames=METRIC_COLUMNS)
        writer.writeheader()
        writer.writerows(existing_rows)
        handle.flush()
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


if __name__ == "__main__":
    raise SystemExit(main())
