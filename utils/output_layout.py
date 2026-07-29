from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .paths import ProjectPaths


@dataclass(frozen=True)
class LocusOutputPaths:
    locus_id: str
    z_score_dir: Path
    ld_dir: Path
    prelim_dir: Path
    annotation_alphagenome_dir: Path
    annotation_legacy_borzoi_dir: Path

    def ensure_dirs(self, dir_keys: tuple[str, ...] | None = None) -> "LocusOutputPaths":
        path_map = {
            "z_score": self.z_score_dir,
            "ld": self.ld_dir,
            "prelim": self.prelim_dir,
            "annotation_alphagenome": self.annotation_alphagenome_dir,
            "annotation_legacy_borzoi": self.annotation_legacy_borzoi_dir,
        }
        selected_keys = tuple(path_map) if dir_keys is None else dir_keys
        for key in selected_keys:
            path_map[key].mkdir(parents=True, exist_ok=True)
        return self

    def z_score_csv(self, gene_name: str) -> Path:
        return self.z_score_dir / f"{gene_name}_GTEx_z_scores.csv"

    def variants_vcf(self, gene_name: str) -> Path:
        return self.z_score_dir / f"{gene_name}_variants.vcf"

    def count_funnel_csv(self, gene_name: str) -> Path:
        return self.prelim_dir / f"{gene_name}_phase1_count_funnel.csv"

    def dataset_metrics_csv(self, gene_name: str) -> Path:
        return self.prelim_dir / f"{gene_name}_phase1_dataset_metrics.csv"

    def interval_genes_csv(self, gene_name: str) -> Path:
        return self.prelim_dir / f"{gene_name}_alphagenome_interval_genes.csv"

    def alphagenome_variant_window_eligibility_csv(self, gene_name: str) -> Path:
        return self.prelim_dir / f"{gene_name}_alphagenome_variant_window_eligibility.csv"

    def ld_variant_order_tsv(self, gene_name: str) -> Path:
        return self.ld_dir / f"{gene_name}_LD_variant_order.tsv"

    def phase1_master_variants_csv(self, gene_name: str) -> Path:
        return self.ld_dir / f"{gene_name}_phase1_master_variants.csv"

    def phase1_variant_map_parquet(self, gene_name: str) -> Path:
        return self.ld_dir / f"{gene_name}_phase1_variant_map.parquet"

    def phase1_z_scores_parquet(self, gene_name: str) -> Path:
        return self.ld_dir / f"{gene_name}_phase1_z_scores.parquet"

    def phase1_ld_long_parquet(self, gene_name: str) -> Path:
        return self.ld_dir / f"{gene_name}_phase1_LD_R_long.parquet"

    def alphagenome_filtered_scores_parquet(self, gene_name: str) -> Path:
        return self.annotation_alphagenome_dir / f"{gene_name}_alphagenome_filtered_scores.parquet"

    def alphagenome_variant_scores_csv(self, gene_name: str) -> Path:
        return self.annotation_alphagenome_dir / f"{gene_name}_alphagenome_variant_scores.csv"

    def alphagenome_histogram_png(self, gene_name: str) -> Path:
        return self.annotation_alphagenome_dir / f"{gene_name}_alphagenome_variant_scores_histogram.png"

    def alphagenome_trimmed_histogram_png(self, gene_name: str) -> Path:
        return self.annotation_alphagenome_dir / f"{gene_name}_alphagenome_variant_scores_histogram_trimmed.png"

    def alphagenome_batch_shard_dir(self, gene_name: str) -> Path:
        return self.annotation_alphagenome_dir / f"{gene_name}_alphagenome_filtered_scores_batches"

    def legacy_borzoi_variant_effects_csv(self, gene_name: str) -> Path:
        return self.annotation_legacy_borzoi_dir / f"{gene_name}_variant_effects.csv"

    def legacy_borzoi_centering_png(self, gene_name: str) -> Path:
        return self.annotation_legacy_borzoi_dir / f"{gene_name}_centering_comparison.png"

    def legacy_borzoi_hist_png(self, gene_name: str) -> Path:
        return self.annotation_legacy_borzoi_dir / f"{gene_name}_expression_delta_histograms.png"

    def phase1_complete(self, gene_name: str) -> bool:
        return all(
            path.exists()
            for path in [
                self.z_score_csv(gene_name),
                self.variants_vcf(gene_name),
                self.ld_variant_order_tsv(gene_name),
                self.phase1_master_variants_csv(gene_name),
                self.phase1_variant_map_parquet(gene_name),
                self.phase1_z_scores_parquet(gene_name),
                self.phase1_ld_long_parquet(gene_name),
                self.count_funnel_csv(gene_name),
                self.dataset_metrics_csv(gene_name),
            ]
        )

    def phase1_alphagenome_ready(self, gene_name: str) -> bool:
        return all(
            path.exists()
            for path in [
                self.phase1_master_variants_csv(gene_name),
            ]
        )

    def annotation_complete(self, gene_name: str) -> bool:
        return all(
            path.exists()
            for path in [
                self.alphagenome_filtered_scores_parquet(gene_name),
                self.alphagenome_variant_scores_csv(gene_name),
                self.alphagenome_histogram_png(gene_name),
            ]
        )


def get_locus_output_paths(
    paths: ProjectPaths,
    locus_id: str,
    *,
    ensure_dirs: bool = True,
    dir_keys: tuple[str, ...] | None = None,
) -> LocusOutputPaths:
    result = LocusOutputPaths(
        locus_id=locus_id,
        z_score_dir=paths.output_z_score / locus_id,
        ld_dir=paths.output_ld / locus_id,
        prelim_dir=paths.output_prelim / locus_id,
        annotation_alphagenome_dir=paths.output_annotation_alphagenome / locus_id,
        annotation_legacy_borzoi_dir=paths.output_annotation_legacy_borzoi / locus_id,
    )
    return result.ensure_dirs(dir_keys=dir_keys) if ensure_dirs else result
