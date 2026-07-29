from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from statistics import NormalDist

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
from utils.paths import configure_runtime_env, find_project_root


BOUNDARY_EPS = 1e-4
WINSOR_LIMIT = 2.5
DEFAULT_MANIFEST_REL = "config/loci_manifest_sample_100_per_chrom.csv"
DEFAULT_SELECTION_REL = "config/annotation_selection_geometric_n20.csv"
DEFAULT_OUTPUT_REL = "output/susine_mu0/geometric_n20"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare standardized AlphaGenome annotation vectors for the geometric N20 real-data sample."
    )
    parser.add_argument("--manifest", default=DEFAULT_MANIFEST_REL, help="Source loci manifest CSV")
    parser.add_argument("--selection", default=DEFAULT_SELECTION_REL, help="Geometric N20 annotation selection CSV")
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_REL, help="Output directory for mu0 artifacts")
    parser.add_argument("--boundary-eps", type=float, default=BOUNDARY_EPS, help="Signed-quantile boundary clip epsilon")
    parser.add_argument("--winsor-limit", type=float, default=WINSOR_LIMIT, help="Inverse-normal winsorization limit")
    parser.add_argument("--fail-fast", action="store_true", help="Stop on first locus failure")
    return parser.parse_args()


def resolve_path(root: Path, path_like: str) -> Path:
    path = Path(path_like)
    return path if path.is_absolute() else root / path


def rms(values: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce")
    x = x[~x.isna()]
    if x.empty:
        return math.nan
    return float(math.sqrt((x**2).mean()))


def q95_abs(values: pd.Series) -> float:
    x = pd.to_numeric(values, errors="coerce").abs()
    x = x[~x.isna()]
    if x.empty:
        return math.nan
    return float(x.quantile(0.95))


def transform_annotation(q: pd.Series, *, boundary_eps: float, winsor_limit: float) -> pd.DataFrame:
    normal = NormalDist()
    q_num = pd.to_numeric(q, errors="coerce").fillna(0.0)
    q_star = q_num.clip(lower=-1.0 + boundary_eps, upper=1.0 - boundary_eps)
    p = (q_star + 1.0) / 2.0
    a_raw = p.map(normal.inv_cdf)
    a_clip = a_raw.clip(lower=-winsor_limit, upper=winsor_limit)
    a_rms = rms(a_clip)
    if math.isfinite(a_rms) and a_rms > 0:
        a = a_clip / a_rms
    else:
        a = pd.Series([0.0] * len(a_clip), index=a_clip.index)
        a_rms = 0.0
    return pd.DataFrame(
        {
            "q_star": q_star,
            "a_raw": a_raw,
            "a_clip": a_clip,
            "annotation_a": a,
            "a_rms": a_rms,
        }
    )


def finite_min(x: float, y: float) -> float:
    if math.isfinite(x) and math.isfinite(y):
        return min(x, y)
    return math.nan


def safe_divide(num: pd.Series, den: pd.Series) -> pd.Series:
    num = pd.to_numeric(num, errors="coerce")
    den = pd.to_numeric(den, errors="coerce")
    out = num / den
    out = pd.to_numeric(out, errors="coerce")
    finite = out.map(lambda value: math.isfinite(value) if pd.notna(value) else False)
    return out.mask(~finite, math.nan)


def prepare_one_locus(locus_cfg: dict, paths, out_dir: Path, *, boundary_eps: float, winsor_limit: float) -> tuple[pd.DataFrame, dict]:
    locus_paths = get_locus_output_paths(paths, locus_cfg["locus_id"], ensure_dirs=False)
    gene_name = locus_cfg["gene_name"]

    master_path = locus_paths.phase1_master_variants_csv(gene_name)
    order_path = locus_paths.ld_variant_order_tsv(gene_name)
    score_path = locus_paths.alphagenome_variant_scores_csv(gene_name)
    for label, path in {
        "phase1_master_variants": master_path,
        "ld_variant_order": order_path,
        "alphagenome_variant_scores": score_path,
    }.items():
        if not path.exists():
            raise FileNotFoundError(f"Missing {label} for {locus_cfg['locus_id']}: {path}")

    master = pd.read_csv(master_path)
    if "ld_included" in master.columns:
        master = master[master["ld_included"].astype(bool)].copy()
    order = pd.read_csv(order_path, sep="\t").sort_values("index").reset_index(drop=True)
    scores = pd.read_csv(score_path)

    required_master = {"variant_id", "z_score", "sample_size"}
    missing_master = required_master.difference(master.columns)
    if missing_master:
        raise ValueError(f"{master_path} is missing columns: {sorted(missing_master)}")
    required_order = {"index", "id"}
    missing_order = required_order.difference(order.columns)
    if missing_order:
        raise ValueError(f"{order_path} is missing columns: {sorted(missing_order)}")
    required_scores = {"source_variant_id", "alphagenome_quantile_mean"}
    missing_scores = required_scores.difference(scores.columns)
    if missing_scores:
        raise ValueError(f"{score_path} is missing columns: {sorted(missing_scores)}")

    expected_index = list(range(len(order)))
    observed_index = order["index"].astype(int).tolist()
    if observed_index != expected_index:
        raise ValueError(f"LD order index is not contiguous 0..p-1 for {locus_cfg['locus_id']}")
    if master["variant_id"].duplicated().any():
        raise ValueError(f"Duplicate variant_id values in {master_path}")
    if scores["source_variant_id"].duplicated().any():
        raise ValueError(f"Duplicate source_variant_id values in {score_path}")

    aligned = (
        order.merge(master, left_on="id", right_on="variant_id", how="left")
        .merge(scores, left_on="id", right_on="source_variant_id", how="left")
        .sort_values("index")
        .reset_index(drop=True)
    )
    if aligned["z_score"].isna().any():
        missing = aligned.loc[aligned["z_score"].isna(), "id"].head(10).tolist()
        raise ValueError(f"Missing z_score rows for {locus_cfg['locus_id']}: {missing}")

    aligned["annotation_missing"] = aligned["alphagenome_quantile_mean"].isna()
    aligned["alphagenome_quantile_mean_filled"] = aligned["alphagenome_quantile_mean"].fillna(0.0)
    aligned = pd.concat(
        [
            aligned,
            transform_annotation(
                aligned["alphagenome_quantile_mean_filled"],
                boundary_eps=boundary_eps,
                winsor_limit=winsor_limit,
            ),
        ],
        axis=1,
    )

    n = pd.to_numeric(aligned["sample_size"], errors="coerce")
    z = pd.to_numeric(aligned["z_score"], errors="coerce")
    adj = (n - 1.0) / ((z**2) + n - 2.0)
    aligned["adj"] = adj
    aligned["beta_hat_std"] = z * (adj / (n - 1.0)).map(math.sqrt)
    aligned["shat2_std"] = 1.0 / (n - 1.0)

    if {"slope", "slope_se"}.issubset(aligned.columns):
        aligned["beta_hat_slope"] = safe_divide(pd.to_numeric(aligned["slope"], errors="coerce"), adj.map(math.sqrt))
    else:
        aligned["beta_hat_slope"] = math.nan

    if {"slope", "slope_se", "af"}.issubset(aligned.columns):
        af = pd.to_numeric(aligned["af"], errors="coerce")
        var_x = 2.0 * af * (1.0 - af)
        r2 = (z**2) / ((z**2) + n - 2.0)
        aligned["var_y_hat_from_slope"] = safe_divide((pd.to_numeric(aligned["slope"], errors="coerce") ** 2) * var_x, r2)
        aligned["var_y_hat_from_se"] = safe_divide(
            (pd.to_numeric(aligned["slope_se"], errors="coerce") ** 2) * (n - 1.0) * var_x,
            1.0 - r2,
        )
    else:
        aligned["var_y_hat_from_slope"] = math.nan
        aligned["var_y_hat_from_se"] = math.nan

    c_rms_l = rms(aligned["beta_hat_std"])
    q95_abs_beta_hat_std = q95_abs(aligned["beta_hat_std"])
    max_abs_a = float(pd.to_numeric(aligned["annotation_a"], errors="coerce").abs().max())
    c_cap_l = q95_abs_beta_hat_std / max_abs_a if math.isfinite(max_abs_a) and max_abs_a > 0 else 0.0
    baseline_c_l = finite_min(c_rms_l, c_cap_l)
    aligned["baseline_c_l"] = baseline_c_l
    aligned["mu0"] = baseline_c_l * aligned["annotation_a"]

    export = pd.DataFrame(
        {
            "locus_id": locus_cfg["locus_id"],
            "gene_name": gene_name,
            "variant_id": aligned["id"].astype(str),
            "ld_matrix_index": aligned["index"].astype(int),
            "chrom": aligned.get("chrom_x", aligned.get("chrom", pd.Series([pd.NA] * len(aligned)))),
            "pos": aligned.get("pos_x", aligned.get("pos", pd.Series([pd.NA] * len(aligned)))),
            "ref": aligned.get("ref_x", aligned.get("ref", pd.Series([pd.NA] * len(aligned)))),
            "alt": aligned.get("alt_x", aligned.get("alt", pd.Series([pd.NA] * len(aligned)))),
            "z_score": pd.to_numeric(aligned["z_score"], errors="coerce"),
            "sample_size": pd.to_numeric(aligned["sample_size"], errors="coerce"),
            "slope": pd.to_numeric(aligned["slope"], errors="coerce") if "slope" in aligned.columns else math.nan,
            "slope_se": pd.to_numeric(aligned["slope_se"], errors="coerce") if "slope_se" in aligned.columns else math.nan,
            "af": pd.to_numeric(aligned["af"], errors="coerce") if "af" in aligned.columns else math.nan,
            "alphagenome_quantile_mean": pd.to_numeric(aligned["alphagenome_quantile_mean"], errors="coerce"),
            "alphagenome_quantile_mean_filled": aligned["alphagenome_quantile_mean_filled"],
            "annotation_missing": aligned["annotation_missing"].astype(bool),
            "q_star": aligned["q_star"],
            "a_raw": aligned["a_raw"],
            "a_clip": aligned["a_clip"],
            "annotation_a": aligned["annotation_a"],
            "adj": aligned["adj"],
            "beta_hat_std": aligned["beta_hat_std"],
            "shat2_std": aligned["shat2_std"],
            "beta_hat_slope": aligned["beta_hat_slope"],
            "var_y_hat_from_slope": aligned["var_y_hat_from_slope"],
            "var_y_hat_from_se": aligned["var_y_hat_from_se"],
            "baseline_c_l": aligned["baseline_c_l"],
            "mu0": aligned["mu0"],
        }
    )

    per_locus_dir = out_dir / "per_locus_annotations"
    per_locus_dir.mkdir(parents=True, exist_ok=True)
    per_locus_path = per_locus_dir / f"{gene_name}_mu0_variant_annotations.csv"
    export.to_csv(per_locus_path, index=False)

    summary = {
        "locus_id": locus_cfg["locus_id"],
        "gene_name": gene_name,
        "n_variants": len(export),
        "n_annotation_missing": int(export["annotation_missing"].sum()),
        "missing_annotation_rate": float(export["annotation_missing"].mean()),
        "q95_abs_beta_hat_std": q95_abs_beta_hat_std,
        "max_abs_a": max_abs_a,
        "a_rms_pre_normalization": float(aligned["a_rms"].iloc[0]) if len(aligned) else math.nan,
        "c_rms_l": c_rms_l,
        "c_cap_l": c_cap_l,
        "baseline_c_l": baseline_c_l,
        "capped": bool(math.isfinite(c_rms_l) and math.isfinite(c_cap_l) and c_cap_l < c_rms_l),
        "mean_abs_mu0": float(export["mu0"].abs().mean()),
        "max_abs_mu0": float(export["mu0"].abs().max()),
        "cor_annotation_z": float(export["annotation_a"].corr(export["z_score"])),
        "cor_annotation_beta_hat_std": float(export["annotation_a"].corr(export["beta_hat_std"])),
        "per_locus_annotation_path": str(per_locus_path),
    }
    return export, summary


def main() -> int:
    root = find_project_root(PROJECT_ROOT)
    paths = configure_runtime_env(root)
    args = parse_args()
    manifest_path = resolve_path(root, args.manifest)
    selection_path = resolve_path(root, args.selection)
    out_dir = resolve_path(root, args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_loci_manifest(manifest_path)
    validate_loci_manifest(manifest)
    selection = load_annotation_selection(selection_path)
    validate_annotation_selection(selection, manifest)
    loci = merge_annotation_selection(manifest, selection)
    if "priority" in loci.columns:
        loci = loci.sort_values(["priority", "locus_id"], na_position="last")
    else:
        loci = loci.sort_values("locus_id")
    if len(loci) != 20:
        raise ValueError(f"Expected 20 selected loci, found {len(loci)} in {selection_path}")

    variant_frames: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    error_rows: list[dict] = []
    for locus_cfg in loci.to_dict("records"):
        try:
            variants, summary = prepare_one_locus(
                locus_cfg,
                paths,
                out_dir,
                boundary_eps=args.boundary_eps,
                winsor_limit=args.winsor_limit,
            )
            variant_frames.append(variants)
            summary_rows.append(summary)
            print(f"Prepared {locus_cfg['locus_id']} ({len(variants)} variants)")
        except Exception as exc:  # noqa: BLE001
            row = {
                "locus_id": locus_cfg.get("locus_id"),
                "gene_name": locus_cfg.get("gene_name"),
                "status": "failed",
                "error_message": f"{type(exc).__name__}: {exc}",
            }
            error_rows.append(row)
            print(f"FAILED {row['locus_id']}: {row['error_message']}")
            if args.fail_fast:
                raise

    if error_rows:
        pd.DataFrame(error_rows).to_csv(out_dir / "mu0_locus_errors.csv", index=False)
    if not variant_frames:
        raise RuntimeError("No loci were prepared successfully.")

    variant_table = pd.concat(variant_frames, ignore_index=True)
    summary_table = pd.DataFrame(summary_rows)
    variant_table.to_parquet(out_dir / "mu0_variant_table.parquet", index=False)
    variant_table.to_csv(out_dir / "mu0_variant_annotations_all.csv", index=False)
    summary_table.to_csv(out_dir / "mu0_locus_summary.csv", index=False)

    print(f"Wrote locus summary: {out_dir / 'mu0_locus_summary.csv'}")
    print(f"Wrote variant parquet: {out_dir / 'mu0_variant_table.parquet'}")
    print(f"Wrote per-locus annotations: {out_dir / 'per_locus_annotations'}")
    if error_rows:
        raise RuntimeError(f"{len(error_rows)} locus/loci failed; see {out_dir / 'mu0_locus_errors.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
