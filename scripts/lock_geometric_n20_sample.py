from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.append(str(PROJECT_ROOT))

from utils.paths import find_project_root


DEFAULT_SAMPLE_REL = (
    "output/baseline_susie_screen/collected/"
    "baseline_susie_pip_k_eff_geometric_sample_n20.csv"
)
DEFAULT_MANIFEST_REL = "config/loci_manifest_sample_100_per_chrom.csv"
DEFAULT_OUT_REL = "output/baseline_susie_screen/selected_geometric_n20"
DEFAULT_CONFIG_SELECTION_REL = "config/annotation_selection_geometric_n20.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Lock the predicted-vs-actual geometric N20 sample for real-data AlphaGenome annotation."
    )
    parser.add_argument("--sample", default=DEFAULT_SAMPLE_REL, help="Geometric sample CSV from the 0.4 workbook")
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST_REL, help="Source loci manifest CSV")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_REL, help="Directory for locked sample artifacts")
    parser.add_argument(
        "--config-selection",
        default=DEFAULT_CONFIG_SELECTION_REL,
        help="Optional config-level copy of the annotation selection CSV",
    )
    parser.add_argument("--sample-n", type=int, default=20, help="Expected number of locked loci")
    return parser.parse_args()


def resolve_path(root: Path, path_like: str) -> Path:
    path = Path(path_like)
    return path if path.is_absolute() else root / path


def main() -> int:
    root = find_project_root(PROJECT_ROOT)
    args = parse_args()
    sample_path = resolve_path(root, args.sample)
    manifest_path = resolve_path(root, args.manifest)
    out_dir = resolve_path(root, args.out_dir)
    config_selection_path = resolve_path(root, args.config_selection) if args.config_selection else None

    if not sample_path.exists():
        raise FileNotFoundError(f"Missing geometric sample CSV: {sample_path}")
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing loci manifest CSV: {manifest_path}")

    sample = pd.read_csv(sample_path)
    manifest = pd.read_csv(manifest_path)

    required_sample = {"geometric_sample_order", "locus_id", "gene_name"}
    missing_sample = required_sample.difference(sample.columns)
    if missing_sample:
        raise ValueError(f"Sample is missing required columns: {sorted(missing_sample)}")
    required_manifest = {"locus_id", "gene_name"}
    missing_manifest = required_manifest.difference(manifest.columns)
    if missing_manifest:
        raise ValueError(f"Manifest is missing required columns: {sorted(missing_manifest)}")
    if sample["locus_id"].duplicated().any():
        dupes = sorted(sample.loc[sample["locus_id"].duplicated(), "locus_id"].unique())
        raise ValueError(f"Sample contains duplicate locus_id values: {dupes}")
    if len(sample) != args.sample_n:
        raise ValueError(f"Expected {args.sample_n} loci, found {len(sample)} in {sample_path}")

    sample = sample.sort_values("geometric_sample_order").reset_index(drop=True)
    expected_order = list(range(1, args.sample_n + 1))
    observed_order = sample["geometric_sample_order"].astype(int).tolist()
    if observed_order != expected_order:
        raise ValueError(
            "geometric_sample_order must be contiguous 1.."
            f"{args.sample_n}; observed {observed_order}"
        )

    manifest_ids = set(manifest["locus_id"].astype(str))
    unknown = sorted(set(sample["locus_id"].astype(str)) - manifest_ids)
    if unknown:
        raise ValueError(f"Locked sample has locus_id values absent from manifest: {unknown}")

    locked_with_manifest = sample.merge(manifest, on="locus_id", how="left", suffixes=("_sample", ""))
    gene_sample_col = "gene_name_sample"
    if gene_sample_col in locked_with_manifest.columns:
        mismatched_gene = locked_with_manifest[
            locked_with_manifest[gene_sample_col].astype(str) != locked_with_manifest["gene_name"].astype(str)
        ]
        if not mismatched_gene.empty:
            raise ValueError(
                "Sample gene_name does not match manifest for: "
                f"{mismatched_gene[['locus_id', gene_sample_col, 'gene_name']].to_dict('records')}"
            )

    annotation_selection = pd.DataFrame(
        {
            "locus_id": sample["locus_id"].astype(str),
            "annotate": True,
            "notes": "baseline_susie_predicted_actual_geometric_n20",
            "priority": sample["geometric_sample_order"].astype(int),
            "annotation_gene_name_override": pd.NA,
            "annotation_gene_id_override": pd.NA,
            "annotation_tissue_override": pd.NA,
        }
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    locked_loci_path = out_dir / "locked_loci.csv"
    locked_manifest_path = out_dir / "locked_loci_with_manifest.csv"
    out_selection_path = out_dir / "annotation_selection_geometric_n20.csv"

    sample.to_csv(locked_loci_path, index=False)
    locked_with_manifest.to_csv(locked_manifest_path, index=False)
    annotation_selection.to_csv(out_selection_path, index=False)

    if config_selection_path is not None:
        config_selection_path.parent.mkdir(parents=True, exist_ok=True)
        annotation_selection.to_csv(config_selection_path, index=False)

    print(f"Wrote locked loci: {locked_loci_path}")
    print(f"Wrote locked manifest join: {locked_manifest_path}")
    print(f"Wrote annotation selection: {out_selection_path}")
    if config_selection_path is not None:
        print(f"Wrote config annotation selection: {config_selection_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
