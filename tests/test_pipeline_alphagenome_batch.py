from __future__ import annotations

import pandas as pd

from utils.pipeline_alphagenome_batch import (
    _build_eligibility_df,
    _chunk_variants_for_scoring,
    estimate_alphagenome_workload,
)


def test_estimate_alphagenome_workload_splits_tss_midpoint_and_ineligible(monkeypatch):
    raw_variants_df = pd.DataFrame(
        [
            {"variant_id": "chr1_1000000_A_G_b38", "chrom": "chr1", "pos": 1_000_000, "ref": "A", "alt": "G", "z_score": 1.0},
            {"variant_id": "chr1_1530000_A_G_b38", "chrom": "chr1", "pos": 1_530_000, "ref": "A", "alt": "G", "z_score": 1.0},
            {"variant_id": "chr1_2100000_A_G_b38", "chrom": "chr1", "pos": 2_100_000, "ref": "A", "alt": "G", "z_score": 1.0},
        ]
    )
    gene_info = {
        "gene_name": "GENE1",
        "chrom": "chr1",
        "start": 900_000,
        "end": 1_200_000,
        "strand": "+",
        "tss": 1_000_000,
        "exons": pd.DataFrame(
            [
                {"chrom": "chr1", "start": 949_999, "end": 950_050},
                {"chrom": "chr1", "start": 1_049_999, "end": 1_050_050},
            ]
        ),
    }

    monkeypatch.setattr(
        "utils.pipeline_alphagenome_batch._load_annotation_inputs",
        lambda locus_cfg, paths, locus_paths: (raw_variants_df, gene_info),
    )

    workload = estimate_alphagenome_workload(
        {"locus_id": "locus1", "gene_name": "GENE1", "alphagenome_sequence_length": "1MB"},
        paths=None,
        locus_paths=None,
    )
    assert workload["variants_loaded"] == 3
    assert workload["variants_tss_centered"] == 1
    assert workload["variants_midpoint_centered"] == 1
    assert workload["variants_ineligible"] == 1
    assert workload["variants_to_score"] == 2


def test_build_eligibility_df_coalesces_nearby_midpoint_variants():
    raw_variants_df = pd.DataFrame(
        [
            {"variant_id": "chr1_1530000_A_G_b38", "chrom": "chr1", "pos": 1_530_000, "ref": "A", "alt": "G", "z_score": 1.0},
            {"variant_id": "chr1_1535000_A_G_b38", "chrom": "chr1", "pos": 1_535_000, "ref": "A", "alt": "G", "z_score": 1.0},
        ]
    )
    gene_info = {
        "gene_name": "GENE1",
        "chrom": "chr1",
        "start": 900_000,
        "end": 1_200_000,
        "strand": "+",
        "tss": 1_000_000,
        "exons": pd.DataFrame(
            [
                {"chrom": "chr1", "start": 949_999, "end": 950_050},
                {"chrom": "chr1", "start": 1_049_999, "end": 1_050_050},
            ]
        ),
    }

    eligibility_df = _build_eligibility_df(
        raw_variants_df,
        gene_info,
        "1MB",
        batch_size=100,
    )
    assert eligibility_df["scoring_mode"].eq("midpoint_centered").all()
    assert eligibility_df["scoring_interval_start_1based"].nunique() == 1
    assert eligibility_df["scoring_interval_end_1based"].nunique() == 1


def test_finalize_aggregated_scores_includes_quantile_mean():
    from utils.pipeline_alphagenome_batch import _finalize_aggregated_scores

    aggregate_state = {
        "var1": {"raw_score_sum": 3.0, "quantile_score_sum": 0.5, "raw_score_count": 2},
    }
    df = _finalize_aggregated_scores(aggregate_state)
    assert list(df.columns) == [
        "source_variant_id",
        "alphagenome_raw_mean",
        "alphagenome_quantile_mean",
        "alphagenome_track_count",
    ]
    assert df.loc[0, "alphagenome_raw_mean"] == 1.5
    assert df.loc[0, "alphagenome_quantile_mean"] == 0.25


def test_chunk_variants_for_scoring_preserves_interval_group_boundaries():
    df = pd.DataFrame(
        [
            {
                "variant_id": f"g1_{idx}",
                "scoring_mode": "midpoint_centered",
                "scoring_interval_start_1based": 1000,
                "scoring_interval_end_1based": 2000,
            }
            for idx in range(78)
        ]
        + [
            {
                "variant_id": f"g2_{idx}",
                "scoring_mode": "midpoint_centered",
                "scoring_interval_start_1based": 3000,
                "scoring_interval_end_1based": 4000,
            }
            for idx in range(22)
        ]
        + [
            {
                "variant_id": f"g3_{idx}",
                "scoring_mode": "midpoint_centered",
                "scoring_interval_start_1based": 5000,
                "scoring_interval_end_1based": 6000,
            }
            for idx in range(78)
        ]
    )

    batches = list(_chunk_variants_for_scoring(df, 100))
    assert [len(batch) for batch in batches] == [100, 78]
    first_batch_keys = {
        tuple(row)
        for row in batches[0][["scoring_interval_start_1based", "scoring_interval_end_1based"]].drop_duplicates().to_numpy()
    }
    second_batch_keys = {
        tuple(row)
        for row in batches[1][["scoring_interval_start_1based", "scoring_interval_end_1based"]].drop_duplicates().to_numpy()
    }
    assert first_batch_keys == {(1000, 2000), (3000, 4000)}
    assert second_batch_keys == {(5000, 6000)}
