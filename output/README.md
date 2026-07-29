# Output Guide

This directory contains artifacts produced by the notebook and script pipeline in `eQTL_annotations_for_susine`. Most files here are generated outputs, not hand-edited source files.

## Pipeline Map

The repo currently has three connected stages.

1. Manifest creation and curation
   - `config/loci_manifest*.csv` defines which loci are eligible to run.
   - `vignettes/0_1_sample_gtex_loci.ipynb` creates sampled manifests such as `config/loci_manifest_sample_100_per_chrom.csv`.

2. Phase 1: GTEx z-scores plus 1000G LD
   - Prototype notebooks:
     - `vignettes/1_get_z_scores.ipynb`
     - `vignettes/2_cbx8_1000g_ld.ipynb`
   - Batch entry points:
     - `vignettes/1_2_batch_controller.ipynb`
     - `scripts/run_phase1_task.py` for one Slurm-array task / manifest chunk
     - `scripts/run_phase1_batch.py` for a local or serial batch run
   - Phase 1 produces either:
     - screening-only aggregate metrics under `output/prelim/`, or
     - full per-locus z-score, LD, and prelim outputs under `output/z_score/`, `output/ld/`, and `output/prelim/`.

3. Screening review and AlphaGenome annotation
   - `vignettes/0_2_phase1_dataset_metrics_screening.ipynb` reads the manifest-level Phase 1 metrics CSV and writes review plots plus selected-locus CSVs under `output/prelim/phase1_metrics_screening_review/`.
   - `vignettes/3_Get_annotations_alphagenome_selected_loci.ipynb` consumes:
     - a loci manifest from `config/`
     - an annotation selection CSV from `output/prelim/phase1_metrics_screening_review/`
   - `scripts/run_annotation_batch.py` is the script equivalent for selected-locus AlphaGenome runs.
   - `vignettes/3_Get_annotations_alphagenome.ipynb` is the earlier single-locus/manual AlphaGenome workflow.
   - `vignettes/3_Get_annotations_legacy_borzoi.ipynb` is historical reference output from the earlier Borzoi-based approach.
   - `vignettes/4_AlphaGenome_output_walkthrough.ipynb` inspects AlphaGenome outputs after scoring completes.

## Directory Layout

### `output/prelim/`

This is the hub for lightweight summaries and screening artifacts.

- Manifest-level screening output:
  - `<manifest_stem>_phase1_dataset_metrics.csv`
  - Example: `loci_manifest_sample_100_per_chrom_phase1_dataset_metrics.csv`
- Per-locus prelim output when Phase 1 runs in `processing` mode:
  - `<locus_id>/<GENE>_phase1_count_funnel.csv`
  - `<locus_id>/<GENE>_phase1_dataset_metrics.csv`
  - `<locus_id>/<GENE>_alphagenome_interval_genes.csv`
  - `<locus_id>/<GENE>_alphagenome_variant_window_eligibility.csv`
- Screening review output from notebook 6:
  - `phase1_metrics_screening_review/*.png`
  - `phase1_metrics_screening_review/*representative*.csv`
  - `phase1_metrics_screening_review/*annotation_selection.csv`

Current examples in this repo:

- `representative_gene_sample_n10_manifest_phase1_dataset_metrics.csv`
- `phase1_metrics_screening_review/representative_gene_sample_n10_manifest.csv`
- `phase1_metrics_screening_review/representative_gene_sample_n10_annotation_selection.csv`

### `output/z_score/`

Phase 1 z-score export artifacts.

- New layout from the shared runners:
  - `<locus_id>/<GENE>_GTEx_z_scores.csv`
  - `<locus_id>/<GENE>_variants.vcf`
- Legacy single-locus notebook outputs may still appear at the top level:
  - `CBX8_GTEx_z_scores.csv`
  - `CBX8_variants.vcf`

These files come from `utils.pipeline_zscore.run_zscore_export()` and are the direct inputs to the LD stage.

### `output/ld/`

Phase 1 LD-derived outputs built from the z-score export and 1000 Genomes EUR genotypes.

- New layout from the shared runners:
  - `<locus_id>/<GENE>_LD_variant_order.tsv`
  - `<locus_id>/<GENE>_phase1_master_variants.csv`
  - `<locus_id>/<GENE>_phase1_variant_map.parquet`
  - `<locus_id>/<GENE>_phase1_z_scores.parquet`
  - `<locus_id>/<GENE>_phase1_LD_R_long.parquet`
- Legacy notebook outputs may still appear at the top level for CBX8.

The key bridge file for downstream annotation is:

- `<GENE>_phase1_master_variants.csv`

That file is the Phase 1 `z ∩ LD` variant set used by AlphaGenome scoring.

### `output/annotation/alphagenome/`

AlphaGenome outputs for selected loci.

Scoring now mixes two interval strategies inside the same locus when needed:

- `tss_centered`: preferred first pass when both alleles fit fully inside the canonical 1 Mb TSS-centered receptive field
- `midpoint_centered`: fallback when the minimal span covering all target-gene exons plus both alleles still fits inside 1 Mb

Variants that fit neither interval are excluded before API submission and recorded in the per-locus eligibility CSV under `output/prelim/`.

- Current manual / older notebook outputs may be flat files at the folder root:
  - `CBX8_alphagenome_filtered_scores.parquet`
  - `CBX8_alphagenome_variant_scores.csv`
  - histogram PNGs
- The new shared output layout expects per-locus subdirectories:
  - `<locus_id>/<GENE>_alphagenome_filtered_scores.parquet`
  - `<locus_id>/<GENE>_alphagenome_variant_scores.csv`
  - `<locus_id>/<GENE>_alphagenome_variant_scores_histogram.png`
  - `<locus_id>/<GENE>_alphagenome_variant_scores_histogram_trimmed.png`
  - `<locus_id>/<GENE>_alphagenome_filtered_scores_batches/part-*.parquet`
- Batch summary files may also be written at the folder root:
  - `annotation_batch_summary.csv`
  - `*_annotation_selection_summary.csv`

Inputs:

- `output/ld/<locus_id>/<GENE>_phase1_master_variants.csv`
- a manifest row from `config/loci_manifest*.csv`
- a selection row from an annotation selection CSV

### `output/annotation/legacy_borzoi/`

Historical reference artifacts from the earlier Borzoi workflow.

- `<GENE>_variant_effects.csv`
- `<GENE>_centering_comparison.png`
- `<GENE>_expression_delta_histograms.png`

These are not the current production path for selected-locus runs.

### `output/slurm_scripts/`

Generated helper files from the batch-controller notebook.

- `phase1_*.slurm` job scripts
- downloaded `.tbi` indexes used by the LD stage

Treat this folder as generated orchestration output, not source code.

## Screening vs Processing

The same Phase 1 code can run in two modes.

- `screening`
  - used to score many loci cheaply
  - writes the manifest-level aggregate metrics CSV in `output/prelim/`
  - does not retain per-locus z-score and LD artifacts
- `processing`
  - writes the full per-locus z-score, LD, and prelim outputs
  - required before AlphaGenome annotation can run from the script-based pipeline

If notebook `3_Get_annotations_alphagenome_selected_loci.ipynb` is running against a selection emitted by notebook `0_2`, the dependency chain is:

`config/loci_manifest*.csv`
-> `scripts/run_phase1_task.py` or `scripts/run_phase1_batch.py`
-> `output/prelim/<manifest_stem>_phase1_dataset_metrics.csv`
-> `vignettes/0_2_phase1_dataset_metrics_screening.ipynb`
-> `output/prelim/phase1_metrics_screening_review/*annotation_selection.csv`
-> `vignettes/3_Get_annotations_alphagenome_selected_loci.ipynb` or `scripts/run_annotation_batch.py`
-> `output/annotation/alphagenome/`

## Practical Notes

- The repo is mid-transition from older single-locus flat outputs to newer per-locus subdirectories keyed by `locus_id`. Both layouts are present in the current worktree.
- Files in `output/` should generally be regenerated from notebooks or scripts rather than edited manually.
- The schema for `config/loci_manifest*.csv` and `annotation_selection.csv` is documented in `config/manifest_schema.md`.
