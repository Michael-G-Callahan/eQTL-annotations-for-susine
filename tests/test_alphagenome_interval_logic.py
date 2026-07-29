from __future__ import annotations

import pandas as pd

from utils.alphagenome import (
    INELIGIBLE_MODE,
    MIDPOINT_CENTERED_MODE,
    TSS_CENTERED_MODE,
    build_tss_centered_interval,
    classify_variant_interval_mode,
    interval_bounds_1based,
)
from utils.variant_processing import annotate_variant_window_eligibility


def _gene_info(*, strand: str = "+") -> dict:
    return {
        "gene_name": "GENE1",
        "chrom": "chr1",
        "start": 900_000,
        "end": 1_200_000,
        "strand": strand,
        "tss": 1_000_000,
        "exons": pd.DataFrame(
            [
                {"chrom": "chr1", "start": 949_999, "end": 950_050},
                {"chrom": "chr1", "start": 1_049_999, "end": 1_050_050},
            ]
        ),
    }


def _row(pos: int, ref: str = "A", alt: str = "G") -> pd.Series:
    return pd.Series(
        {
            "variant_id": f"chr1_{pos}_{ref}_{alt}_b38",
            "chrom": "chr1",
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "z_score": 1.0,
        }
    )


def test_snv_fully_inside_tss_interval_is_tss_centered():
    gene_info = _gene_info()
    result = classify_variant_interval_mode(_row(1_000_000), gene_info, "1MB")
    assert result["scoring_mode"] == TSS_CENTERED_MODE
    assert result["fits_tss_centered"] is True
    assert result["excluded_from_scoring"] is False


def test_alt_allele_spilling_past_edge_is_not_tss_centered():
    gene_info = _gene_info()
    interval = build_tss_centered_interval(gene_info, "1MB")
    _, interval_end = interval_bounds_1based(interval)
    result = classify_variant_interval_mode(_row(interval_end, ref="A", alt="GG"), gene_info, "1MB")
    assert result["fits_tss_centered"] is False
    assert result["scoring_mode"] != TSS_CENTERED_MODE


def test_variant_can_fall_back_to_midpoint_centered():
    gene_info = _gene_info()
    result = classify_variant_interval_mode(_row(1_530_000), gene_info, "1MB")
    assert result["fits_tss_centered"] is False
    assert result["fits_midpoint_centered"] is True
    assert result["scoring_mode"] == MIDPOINT_CENTERED_MODE


def test_variant_exceeding_gene_variant_span_is_ineligible():
    gene_info = _gene_info()
    result = classify_variant_interval_mode(_row(2_100_000), gene_info, "1MB")
    assert result["fits_tss_centered"] is False
    assert result["fits_midpoint_centered"] is False
    assert result["scoring_mode"] == INELIGIBLE_MODE
    assert result["excluded_from_scoring"] is True


def test_negative_strand_gene_uses_same_containment_logic():
    gene_info = _gene_info(strand="-")
    result = classify_variant_interval_mode(_row(1_530_000), gene_info, "1MB")
    assert result["scoring_mode"] == MIDPOINT_CENTERED_MODE


def test_annotate_variant_window_eligibility_adds_expected_columns():
    gene_info = _gene_info()
    df = pd.DataFrame(
        [
            _row(1_000_000).to_dict(),
            _row(1_530_000).to_dict(),
            _row(2_100_000).to_dict(),
        ]
    )
    annotated = annotate_variant_window_eligibility(df, gene_info, "1MB")
    for column in [
        "ref_length",
        "alt_length",
        "variant_min_1based",
        "variant_max_1based",
        "fits_tss_centered",
        "fits_midpoint_centered",
        "scoring_mode",
        "excluded_from_scoring",
        "exclusion_reason",
        "scoring_interval_start_1based",
        "scoring_interval_end_1based",
        "midpoint_span_bp",
    ]:
        assert column in annotated.columns
