from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    data: Path
    output: Path
    output_z_score: Path
    output_ld: Path
    output_annotation: Path
    output_annotation_alphagenome: Path
    output_annotation_legacy_borzoi: Path
    output_prelim: Path
    vignettes: Path
    utils: Path
    wandb: Path
    genomes: Path
    gtex: Path
    gtf_cache: Path
    gtf_shortcuts: Path
    fasta_shortcuts: Path


def _candidate_starts(anchor: Optional[Path] = None) -> list[Path]:
    candidates: list[Path] = []

    if anchor is not None:
        candidates.append(Path(anchor).resolve())

    if "__file__" in globals():
        candidates.append(Path(__file__).resolve())

    candidates.append(Path.cwd().resolve())
    return candidates


def find_project_root(anchor: Optional[Path] = None) -> Path:
    """
    Locate the repository root by walking upward until expected project
    directories are present.
    """
    for start in _candidate_starts(anchor):
        current = start if start.is_dir() else start.parent
        for candidate in [current, *current.parents]:
            if (
                (candidate / "utils").is_dir()
                and (candidate / "vignettes").is_dir()
                and (candidate / "README.md").exists()
            ):
                return candidate

    raise FileNotFoundError("Could not locate eQTL_annotations_for_susine project root")


def get_project_paths(anchor: Optional[Path] = None) -> ProjectPaths:
    root = find_project_root(anchor)
    return ProjectPaths(
        root=root,
        data=root / "data",
        output=root / "output",
        output_z_score=root / "output" / "z_score",
        output_ld=root / "output" / "ld",
        output_annotation=root / "output" / "annotation",
        output_annotation_alphagenome=root / "output" / "annotation" / "alphagenome",
        output_annotation_legacy_borzoi=root / "output" / "annotation" / "legacy_borzoi",
        output_prelim=root / "output" / "prelim",
        vignettes=root / "vignettes",
        utils=root / "utils",
        wandb=root / "wandb",
        genomes=root / "data" / "genomes",
        gtex=root / "data" / "gtex",
        gtf_cache=root / "data" / "gtf_cache",
        gtf_shortcuts=root / "data" / "gtf_shortcuts",
        fasta_shortcuts=root / "data" / "fasta_shortcuts",
    )


def ensure_project_dirs(
    anchor: Optional[Path] = None,
    *,
    dir_keys: Optional[tuple[str, ...]] = None,
) -> ProjectPaths:
    paths = get_project_paths(anchor)
    path_map = {
        "data": paths.data,
        "output": paths.output,
        "output_z_score": paths.output_z_score,
        "output_ld": paths.output_ld,
        "output_annotation": paths.output_annotation,
        "output_annotation_alphagenome": paths.output_annotation_alphagenome,
        "output_annotation_legacy_borzoi": paths.output_annotation_legacy_borzoi,
        "output_prelim": paths.output_prelim,
        "wandb": paths.wandb,
        "genomes": paths.genomes,
        "gtex": paths.gtex,
        "gtf_cache": paths.gtf_cache,
        "gtf_shortcuts": paths.gtf_shortcuts,
        "fasta_shortcuts": paths.fasta_shortcuts,
    }
    selected_keys = tuple(path_map) if dir_keys is None else dir_keys
    for key in selected_keys:
        path = path_map[key]
        path.mkdir(parents=True, exist_ok=True)
    return paths


def add_project_root_to_sys_path(anchor: Optional[Path] = None) -> Path:
    root = find_project_root(anchor)
    root_str = str(root)
    if root_str not in sys.path:
        sys.path.append(root_str)
    return root


def configure_runtime_env(
    anchor: Optional[Path] = None,
    *,
    dir_keys: Optional[tuple[str, ...]] = None,
) -> ProjectPaths:
    paths = ensure_project_dirs(anchor, dir_keys=dir_keys)
    os.environ.setdefault("GENOMEPY_CACHE_DIR", str(paths.genomes))
    os.environ.setdefault("WANDB_CACHE_DIR", str(paths.wandb))
    return paths


def resolve_repo_path(path_like: os.PathLike[str] | str, anchor: Optional[Path] = None) -> Path:
    path = Path(path_like).expanduser()
    if path.is_absolute():
        return path
    return find_project_root(anchor) / path


def find_gtex_parquet(
    tissue: str = "Lung",
    chrom: Optional[str] = None,
    anchor: Optional[Path] = None,
) -> Path:
    """
    Find a GTEx parquet file under data/gtex by tissue and optional chromosome.
    """
    paths = get_project_paths(anchor)
    matches = []
    tissue_lower = tissue.lower()
    chrom_lower = chrom.lower() if chrom else None

    for candidate in sorted(paths.gtex.glob("*.parquet")):
        name_lower = candidate.name.lower()
        if tissue_lower not in name_lower:
            continue

        if chrom_lower:
            chrom_match = re.search(r"\.(chr[^.]+)\.parquet$", name_lower)
            if chrom_match is None:
                continue
            candidate_chrom = chrom_match.group(1)
            if candidate_chrom != chrom_lower:
                continue
        matches.append(candidate)

    if len(matches) == 1:
        return matches[0]
    if not matches:
        chrom_msg = f" and chromosome {chrom!r}" if chrom else ""
        raise FileNotFoundError(
            f"No GTEx parquet matched tissue {tissue!r}{chrom_msg} in {paths.gtex}"
        )
    raise FileExistsError(
        f"Multiple GTEx parquet files matched tissue {tissue!r}: {[str(m.name) for m in matches]}"
    )
