from __future__ import annotations

import argparse
import os
import sys
import traceback
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from alphagenome.models import dna_client, variant_scorers

from utils.alphagenome import TSS_CENTERED_MODE, build_variant_inputs
from utils.locus_manifest import (
    load_annotation_selection,
    load_loci_manifest,
    merge_annotation_selection,
    validate_annotation_selection,
    validate_loci_manifest,
)
from utils.output_layout import get_locus_output_paths
from utils.paths import add_project_root_to_sys_path, configure_runtime_env, find_project_root
from utils.pipeline_alphagenome_batch import _build_eligibility_df, _load_annotation_inputs


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="One-off AlphaGenome quota test: submit a single raw score_variants batch and print the raw result or exception."
    )
    parser.add_argument(
        "--manifest",
        default="config/loci_manifest_sample_100_per_chrom.csv",
        help="Path to loci manifest.",
    )
    parser.add_argument(
        "--selection",
        default="output/prelim/phase1_metrics_screening_review/representative_gene_sample_n10_annotation_selection.csv",
        help="Path to annotation selection CSV.",
    )
    parser.add_argument("--locus-id", default=None, help="Optional locus_id to test.")
    parser.add_argument("--batch-size", type=int, default=200, help="Number of variants to submit in one API batch.")
    parser.add_argument("--max-workers", type=int, default=32, help="max_workers passed to score_variants.")
    parser.add_argument(
        "--api-key-env-var",
        default="ALPHAGENOME_API_KEY",
        help="Environment variable holding the AlphaGenome API key.",
    )
    return parser.parse_args()


def main() -> int:
    project_root = add_project_root_to_sys_path(find_project_root())
    paths = configure_runtime_env(project_root)
    args = parse_args()
    api_key = os.environ.get(args.api_key_env_var)
    if not api_key:
        raise RuntimeError(f"Missing AlphaGenome API key in environment variable {args.api_key_env_var}")

    manifest_df = load_loci_manifest(args.manifest)
    validate_loci_manifest(manifest_df)
    selection_df = load_annotation_selection(args.selection)
    validate_annotation_selection(selection_df, manifest_df)
    loci_df = merge_annotation_selection(manifest_df, selection_df)
    if loci_df.empty:
        raise ValueError("No loci remain after applying the annotation selection file.")

    locus_cfg = _select_locus(args, loci_df.to_dict("records"), paths)
    locus_paths = get_locus_output_paths(paths, locus_cfg["locus_id"], ensure_dirs=False)
    raw_variants_df, gene_info = _load_annotation_inputs(locus_cfg, paths, locus_paths)
    eligibility_df = _build_eligibility_df(
        raw_variants_df,
        gene_info,
        locus_cfg["alphagenome_sequence_length"],
        batch_size=int(locus_cfg["alphagenome_batch_size"]),
    )
    tss_df = eligibility_df[eligibility_df["scoring_mode"] == TSS_CENTERED_MODE].copy()
    if len(tss_df) < args.batch_size:
        raise ValueError(
            f"Locus {locus_cfg['locus_id']} has only {len(tss_df)} TSS-centered variants; "
            f"cannot build a {args.batch_size}-variant shared-interval batch."
        )

    batch_df = tss_df.head(int(args.batch_size)).copy().reset_index(drop=True)
    variants, intervals = build_variant_inputs(batch_df)
    shared_interval = intervals[0]
    scorer = variant_scorers.GeneMaskLFCScorer(requested_output=dna_client.OutputType.RNA_SEQ)
    dna_model = dna_client.create(api_key)

    print(f"Manifest: {args.manifest}")
    print(f"Selection: {args.selection}")
    print(f"Locus: {locus_cfg['locus_id']} ({locus_cfg['gene_name']})")
    print(f"Submitting one raw score_variants batch with batch_size={args.batch_size}, max_workers={args.max_workers}")
    print(f"Scoring mode: {TSS_CENTERED_MODE}")
    print(f"Interval: {shared_interval}")
    print(f"First variant: {batch_df.iloc[0]['variant_id']}")
    print(f"Last variant: {batch_df.iloc[-1]['variant_id']}")
    print("Calling AlphaGenome API now...")

    result = None
    try:
        result = dna_model.score_variants(
            intervals=shared_interval,
            variants=variants,
            variant_scorers=[scorer],
            progress_bar=False,
            max_workers=int(args.max_workers),
        )
        print("API call returned successfully.")
        print(f"Python type: {type(result)}")
        try:
            print(f"len(result): {len(result)}")
        except Exception:
            print("len(result): <not available>")
        print("Result summary:")
        _print_result_summary(result)
    except Exception as exc:  # noqa: BLE001
        print("API call raised an exception.")
        print(f"Exception type: {type(exc).__name__}")
        print(f"Exception str: {exc}")
        for attr in ("code", "details", "debug_error_string"):
            if hasattr(exc, attr):
                try:
                    value = getattr(exc, attr)
                    value = value() if callable(value) else value
                    print(f"{attr}: {value}")
                except Exception as attr_exc:  # noqa: BLE001
                    print(f"{attr}: <error while reading attribute: {attr_exc}>")
        print("traceback:")
        print(traceback.format_exc())
        print(f"result variable after exception is None: {result is None}")
        return 1

    return 0


def _select_locus(args: argparse.Namespace, loci: list[dict], paths) -> dict:
    if args.locus_id:
        matches = [locus for locus in loci if locus["locus_id"] == args.locus_id]
        if not matches:
            raise ValueError(f"Requested locus_id {args.locus_id!r} was not present in the selected loci inventory.")
        return matches[0]

    for locus_cfg in loci:
        locus_paths = get_locus_output_paths(paths, locus_cfg["locus_id"], ensure_dirs=False)
        try:
            raw_variants_df, gene_info = _load_annotation_inputs(locus_cfg, paths, locus_paths)
            eligibility_df = _build_eligibility_df(
                raw_variants_df,
                gene_info,
                locus_cfg["alphagenome_sequence_length"],
                batch_size=int(locus_cfg["alphagenome_batch_size"]),
            )
        except Exception:
            continue
        if int((eligibility_df["scoring_mode"] == TSS_CENTERED_MODE).sum()) >= int(args.batch_size):
            return locus_cfg

    raise ValueError(
        f"Could not find any selected locus with at least {args.batch_size} TSS-centered variants for a raw shared-interval batch test."
    )

def _print_result_summary(result: Any) -> None:
    if not result:
        print("  empty result")
        return

    first_item = result[0]
    print(f"  outer container type: {type(result).__name__}")
    print(f"  item type: {type(first_item).__name__}")

    if hasattr(first_item, "n_obs") and hasattr(first_item, "n_vars"):
        print(f"  first item shape: n_obs={first_item.n_obs}, n_vars={first_item.n_vars}")
        print(f"  first item obs columns: {list(first_item.obs.columns)}")
        print(f"  first item var columns: {list(first_item.var.columns)}")
        print(f"  first item layers: {list(first_item.layers.keys())}")
        print(f"  first item uns keys: {list(first_item.uns.keys())}")

        try:
            gene_names = first_item.obs["gene_name"].dropna().astype(str).unique().tolist()
            print(f"  first item gene_name sample: {gene_names[:10]}")
        except Exception:
            pass

        try:
            gtex_tissues = first_item.var["gtex_tissue"].dropna().astype(str).unique().tolist()
            print(f"  first item GTEx tissue sample: {gtex_tissues[:10]}")
        except Exception:
            pass
        return

    print(f"  first item repr: {repr(first_item)}")


if __name__ == "__main__":
    raise SystemExit(main())
