# eQTL_annotations_for_susine

This repo prepares GTEx Lung loci, 1000 Genomes LD, and AlphaGenome-derived signed
annotations for annotation-informed **SuSiNE** RSS analyses. The production path is the
**20-locus geometric-N20 panel** selected from the Phase 1 screening and baseline SuSiE
workflow. It builds the `mu_0 = c * a` handoff objects that the downstream `test_susine`
real-data pipeline ingests; it does **not** run the SuSiNE model grid itself.

> **Canonical reproduction guide:** for the paper's real-data inputs, follow
> [`real_data_geometric_n20_workflow.md`](real_data_geometric_n20_workflow.md). It is the
> authoritative provenance record and command sequence (candidate sampling → z-scores → LD
> → baseline SuSiE screen → geometric-N20 selection → AlphaGenome scoring → `mu_0` prep).
> The output inventory is documented in [output/README.md](output/README.md).

## Setup

```bash
pip install -r requirements.txt
```

**Prerequisites:**

- **`ALPHAGENOME_API_KEY`** environment variable — the AlphaGenome client (Google DeepMind;
  PyPI `alphagenome`, source `google-deepmind/alphagenome`) reads it at runtime. No key is
  stored in this repo.
- **Raw GTEx v10 Lung allpairs parquet** under `data/gtex/` — not shipped (large, third-party;
  obtain from the GTEx Portal / dbGaP). `data/` is git-ignored.
- **Linux/HPC assumed for LD generation** — `pysam` (VCF/tabix reads of 1000G genotypes) is
  POSIX-oriented; the LD stage is not expected to run on native Windows.
- **Downstream fit backend:** the vanilla-SuSiE screen and RSS validation R workbooks use
  upstream [`stephenslab/susieR`](https://github.com/stephenslab/susieR), pinned to tag
  `0.15.53` for the real-data analyses:
  `remotes::install_github("stephenslab/susieR", ref = "0.15.53")`.

## Repo Role

The repo is responsible for:

- sampling and screening GTEx loci,
- exporting Phase 1 GTEx summary-stat inputs,
- building LD-aligned per-locus variant sets,
- running a baseline SuSiE screen and selecting the geometric-N20 panel,
- scoring selected loci with AlphaGenome,
- transforming AlphaGenome outputs into annotation vectors `a`,
- calibrating the annotation scale for `mu_0 = c * a`,
- exporting ready-to-fit RSS handoff objects for later SuSiNE fitting.

The full annotation-informed SuSiNE model grid runs downstream in the SuSiNE repos.

## Selected Panel — geometric-N20

The paper uses the **20-locus geometric-N20 panel**, selected by a predicted-vs-actual
maximin (spread-out) sample over the baseline SuSiE posterior-diffuseness plot. The
selection is defined by `config/annotation_selection_geometric_n20.csv` and the locked-loci
artifacts under `output/baseline_susie_screen/selected_geometric_n20/`, produced by
`scripts/lock_geometric_n20_sample.py`. See workflow doc §4 for the algorithm.

> **Superseded — earlier representative-n10 selection.** An earlier version of this repo used
> a 10-locus representative panel (`aoc3`, `cdc42se2`, `gal3st1`, `gpr89b`, `hdac5`, `mei1`,
> `pus7l`, `spring1`, `tnpo1`, `tspan1`). It is retained only as selection history; the paper's
> real-data study uses the geometric-N20 panel above.

## Workflow

### 0.1 Sample candidate loci

`vignettes/0_1_sample_gtex_loci.ipynb` — scans local GTEx Lung allpairs parquet, samples
protein-coding genes by chromosome, writes a candidate manifest
(`config/loci_manifest_sample_100_per_chrom.csv`).

### 0.2 Review Phase 1 screening metrics

`vignettes/0_2_phase1_dataset_metrics_screening.ipynb` — loads aggregate screening metrics
for the sampled loci and reviews the metric distributions.

### 0.3 / 0.4 Baseline SuSiE screen and geometric-N20 selection (R)

`vignettes/0_3_qc_and_baseline_susie_screen.Rmd` runs QC and a vanilla `susieR::susie_rss()`
screen on each eligible locus; `vignettes/0_4_explore_baseline_susie_metrics.Rmd` fits the
predicted-vs-actual model and selects the geometric-N20 sample. The lock step is
`scripts/lock_geometric_n20_sample.py`.

### 1. Export GTEx summary-stat inputs

`vignettes/1_get_z_scores.ipynb` — reads GTEx allpairs parquet; exports `slope`, `slope_se`,
`z_score`, `af`, `sample_size`, and variant lists for the LD stage.

### 1.2 Batch controller for Phase 1

`vignettes/1_2_batch_controller.ipynb` — builds SLURM job-array scripts for the shared Phase 1
runners (z-score export, LD construction, metric generation).

### 2. LD construction

The LD stage builds the canonical Phase 1 variant set and LD artifacts per locus from 1000G
EUR genotypes. Key outputs:

- `output/ld/<locus_id>/<GENE>_phase1_master_variants.csv`
- `output/ld/<locus_id>/<GENE>_LD_variant_order.tsv`
- `output/ld/<locus_id>/<GENE>_phase1_LD_R_long.parquet`

`vignettes/2_cbx8_1000g_ld.ipynb` is the readable single-locus (CBX8) prototype.
`phase1_master_variants.csv` is the bridge into AlphaGenome scoring and `mu_0` prep.

### 3. AlphaGenome scoring for selected loci

`vignettes/3_Get_annotations_alphagenome_selected_loci.ipynb` (script equivalent:
`scripts/run_annotation_batch.py`) scores the selected Phase 1 loci and writes per-locus
filtered parquet under `output/annotation/alphagenome/<locus_id>/`. Scoring uses a two-pass
1 Mb interval policy (`tss_centered`, then `midpoint_centered` fallback); ineligible variants
are excluded before API submission and logged under `output/prelim/`.

`vignettes/3_Get_annotations_alphagenome.ipynb` is the earlier single-locus/manual flow.
`vignettes/3_Get_annotations_legacy_borzoi.ipynb` is the **legacy Borzoi path — superseded,
retained for provenance** (see `annotation_scaling_decision_2026-04-09.md`).

### 4. Review AlphaGenome outputs

`vignettes/4_AlphaGenome_output_walkthrough.ipynb` — scored-variant counts, TSS- vs
midpoint-centered usage, and raw/quantile score distributions.

### 5. Build annotation-derived `mu_0` inputs

`vignettes/5_build_annotation_mu0_selected_loci.ipynb` (script equivalent for the N20 panel:
`scripts/prepare_geometric_n20_mu0.py`) aligns the selected loci across Phase 1 summary stats,
LD order, and AlphaGenome outputs; converts signed `quantile_score` into per-locus annotation
templates `a`; computes standardized-path calibration quantities; estimates per-locus and
pooled annotation scales `c`; and exports `mu_0 = c * a` grids and ready-to-fit RSS handoff
objects under `output/susine_mu0/<selection_name>/` (`geometric_n20` for the paper).

Key artifacts: `mu0_locus_summary.csv`, `mu0_variant_table.parquet`, `mu0_scale_grid.csv`,
`mu0_ready_rss_inputs.pkl`, plus diagnostic plots.

### Legacy vanilla RSS baseline (R)

`vignettes/5_run_susine_rss_selected_loci.Rmd` runs plain `susieR::susie_rss()` on the
selected loci (reference baseline, no `mu_0`); `vignettes/5_validate_susie_rss.Rmd` validates
z/LD/annotation alignment; `vignettes/6_run_susine_baseline_selected_loci.Rmd` is the baseline
run workbook.

## How this feeds `test_susine`

The Step 5 exports under `output/susine_mu0/geometric_n20/` are the cross-repo handoff:

- `test_susine::sync_real_data_inputs()` pulls the prepared per-locus annotation files into
  `test_susine`.
- `test_susine:::load_real_data_locus_bundle()` validates that the annotation file, LD order,
  master-variant table, and LD matrix all share one variant order, then exposes
  `z`, `R`, `n_sample`, `a`, `baseline_c_l`, and `variant_map`.
- The real-data task runner passes `mu_0 = c_value * bundle$a` to `susine::susine_rss()`.

See workflow doc §7 for the full handoff contract.

## Important Interfaces

- canonical SNP order comes from `output/ld/<locus_id>/<GENE>_LD_variant_order.tsv`
- the per-locus annotation join key is `source_variant_id`
- final `a` and all exported `mu_0` vectors are in LD order
- downstream SuSiNE RSS consumers should use the Step 5 exports rather than recomputing transforms
- the Python Step 5 workbook serializes the ready-to-fit handoff object as pickle, not `.rds`

## Repo Layout Notes

- `vignettes/` mixes Python prep (`.ipynb`) and R validation/screen (`.Rmd`); the `.Rproj` is
  kept for the R workbooks.
- Shared path helpers live in `utils/paths.py`; the **Python core** resolves paths from the
  project root (via `find_project_root()`) rather than machine-specific absolute paths. The
  `.Rmd` workbooks and SLURM templates contain example HPC paths/emails to edit for your
  environment.
- Outputs are written under git-ignored `data/`, `output/`, and `wandb/`; only
  `output/README.md` is tracked. Config CSVs under `config/` (the study definition) are tracked.
- Provenance docs to read: [`real_data_geometric_n20_workflow.md`](real_data_geometric_n20_workflow.md)
  (canonical), `annotation_scaling_decision_2026-04-09.md` (Borzoi → AlphaGenome + scaling
  rationale), and `config/manifest_schema.md` (manifest field schema).

## License

MIT — see [LICENSE](LICENSE).
