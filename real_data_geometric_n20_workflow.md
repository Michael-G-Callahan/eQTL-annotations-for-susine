# Geometric N20 Real-Data Workflow

This document records the data provenance and preprocessing choices for the
20-locus real-data workflow used to generate AlphaGenome annotations and
annotation-informed SuSiNE RSS inputs. It is a reproducibility and methods
handoff note, not a polished manuscript Methods section.

## 1. Candidate Locus Source

Candidate loci were drawn from local GTEx v10 Lung allpairs parquet files under
`data/gtex/`. The locus sampling workflow is defined in
`vignettes/0_1_sample_gtex_loci.ipynb`.

The candidate manifest was built by scanning the available chromosome-specific
GTEx Lung allpairs files, mapping GTEx gene IDs to GENCODE gene metadata, and
sampling protein-coding genes by chromosome. GENCODE annotations were read from
the cached GTF file under `data/gtf_cache/`; when missing, the code downloads
GENCODE v44 for hg38 from EBI and stores shortcuts under `data/gtf_shortcuts/`.

The main manifest used for the baseline screen and geometric N20 selection is:

```text
config/loci_manifest_sample_100_per_chrom.csv
```

Important manifest fields are:

- `locus_id`: stable per-locus identifier used for output directories.
- `gene_name`, `gene_id`: GTEx/GENCODE target gene identity.
- `gtex_tissue`, `gtex_chrom`: source tissue and chromosome.
- `reference_genome`: default `hg38`.
- `maf_min`, `maf_max`: default allele-frequency bounds `0.01` and `0.99`.
- `min_sample_size`: default effective sample-size threshold `50`.
- `alphagenome_sequence_length`: default `1MB`.
- `alphagenome_target_data_source`: default `gtex`.
- `alphagenome_target_gtex_tissue`: default inherited from `gtex_tissue`.
- `alphagenome_batch_size`, `alphagenome_max_workers`,
  `alphagenome_retry_wait_seconds`: default AlphaGenome execution controls.

The schema and defaults are documented in `config/manifest_schema.md` and
implemented in `utils/locus_manifest.py` and `utils/pipeline_defaults.py`.

## 2. GTEx Z-Score Preprocessing

The z-score export logic is implemented in `utils/pipeline_zscore.py` and was
first prototyped in `vignettes/1_get_z_scores.ipynb`.

For each manifest locus, the workflow reads the GTEx Lung allpairs parquet rows
for the target `gene_id`. The source GTEx fields include `slope`, `slope_se`,
`af`, and `ma_count`.

The effective sample size is computed as:

```text
sample_size = ma_count / (2 * af)
```

The marginal z-score is computed as:

```text
z_score = slope / slope_se
```

Variants are then filtered by:

```text
maf_min <= af <= maf_max
sample_size > min_sample_size
```

For the geometric N20 workflow, the default manifest thresholds are:

```text
0.01 <= af <= 0.99
sample_size > 50
```

The z-score stage writes:

```text
output/z_score/<locus_id>/<GENE>_GTEx_z_scores.csv
output/z_score/<locus_id>/<GENE>_variants.vcf
output/prelim/<locus_id>/<GENE>_phase1_count_funnel.csv
```

The exported z-score CSV retains the summary-statistic inputs needed later:
`variant_id`, `slope`, `slope_se`, `z_score`, `af`, and `sample_size`.

## 3. LD Reference Panel

LD is computed with the 1000 Genomes Project GRCh38 reference panel. The LD
logic is implemented in `utils/pipeline_ld.py` and was first prototyped in
`vignettes/2_cbx8_1000g_ld.ipynb`.

The workflow uses:

- the 1000 Genomes integrated sample panel:
  `integrated_call_samples_v3.20130502.ALL.panel`
- EUR samples from that panel
- chromosome-specific 1000 Genomes GRCh38 biallelic SNV/INDEL VCFs from EBI:

```text
https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/
1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/
ALL.chr<chrom>.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
```

For each locus, the variant list from the GTEx z-score export is matched to the
1000 Genomes VCF by chromosome, position, reference allele, and alternate
allele. Dosage genotypes are extracted for EUR samples, centered, scaled, and
used to compute the Pearson correlation LD matrix:

```text
R = cor(G_EUR)
```

The LD matrix is rounded to three decimals and stored compactly as upper-triangle
long-form parquet. The canonical variant order for all downstream work is the LD
order:

```text
output/ld/<locus_id>/<GENE>_LD_variant_order.tsv
```

The LD stage also writes:

```text
output/ld/<locus_id>/<GENE>_phase1_master_variants.csv
output/ld/<locus_id>/<GENE>_phase1_variant_map.parquet
output/ld/<locus_id>/<GENE>_phase1_z_scores.parquet
output/ld/<locus_id>/<GENE>_phase1_LD_R_long.parquet
output/prelim/<locus_id>/<GENE>_phase1_dataset_metrics.csv
```

`phase1_master_variants.csv` is the bridge between the GTEx z-score inputs,
AlphaGenome scoring, baseline SuSiE screening, and final SuSiNE handoff.

## 4. Baseline SuSiE Screen And Geometric Sampling

The baseline screen is controlled by:

```text
vignettes/0_3_qc_and_baseline_susie_screen.Rmd
scripts/run_baseline_susie_screen_task.R
```

Each QC-passing locus is fit with vanilla SuSiE RSS:

```r
susieR::susie_rss(
  z = inputs$z,
  R = inputs$R,
  n = inputs$n_rss,
  L = 10,
  estimate_residual_variance = TRUE,
  check_prior = FALSE
)
```

The screen writes compact per-locus metrics only. It does not save full SuSiE
fit objects. The key baseline metrics include:

- `max_pip`
- `pip_sum`
- `pip_entropy`
- `pip_k_eff`
- PIP threshold counts
- z-score concentration metrics
- per-effect alpha concentration and credible-set purity summaries

The exploration workbook is:

```text
vignettes/0_4_explore_baseline_susie_metrics.Rmd
```

That workbook joins baseline SuSiE metrics to Phase 1 dataset metrics and fits a
five-fold cross-validated PCA-plus-linear-regression predictor of:

```text
log1p(pip_k_eff)
```

from dataset-level metrics such as `M1`, `M2`, `M4`, z-score effective signal
count, z-score count above 3, and variant-count/matching summaries.

The official geometric N20 sample is the first predicted-vs-actual maximin
sample. It is selected on the two-dimensional scatterplot:

```text
x = log1p(predicted pip_k_eff from dataset metrics)
y = log1p(actual baseline SuSiE pip_k_eff)
```

The selection algorithm robustly scales these two coordinates, starts from the
forced anchor locus if present, and then repeatedly adds the locus whose minimum
Euclidean distance to the already-selected set is largest. In plain terms, it
chooses 20 loci that spread out over the observed predicted-vs-actual posterior
diffuseness plot, without imposing hard region quotas.

The locked sample artifacts are written by the workbook and can also be
regenerated with `scripts/lock_geometric_n20_sample.py`:

```text
output/baseline_susie_screen/selected_geometric_n20/locked_loci.csv
output/baseline_susie_screen/selected_geometric_n20/locked_loci_with_manifest.csv
output/baseline_susie_screen/selected_geometric_n20/annotation_selection_geometric_n20.csv
config/annotation_selection_geometric_n20.csv
```

The residual-geometry sample is retained only as a diagnostic comparison. The
official sample for this workflow is the predicted-vs-actual geometric maximin
sample.

## 5. AlphaGenome Annotation Generation

AlphaGenome annotation is run for the locked 20 loci using:

```text
scripts/run_annotation_batch.py
utils/pipeline_alphagenome_batch.py
utils/alphagenome.py
```

The run uses the main manifest and the geometric N20 annotation selection:

```bash
python scripts/run_annotation_batch.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --summary-path output/annotation/alphagenome/geometric_n20_annotation_batch_summary.csv
```

The AlphaGenome API key is read from `ALPHAGENOME_API_KEY` unless
`--api-key-env-var` is supplied.

For each selected locus, the annotation runner reads:

```text
output/ld/<locus_id>/<GENE>_phase1_master_variants.csv
```

and scores variants against the target gene in the GTEx/Lung RNA-seq context.
The default sequence length is `1MB`.

The interval policy is:

- `tss_centered`: use a canonical TSS-centered interval when both alleles fit
  inside the 1 Mb receptive field.
- `midpoint_centered`: use a fallback interval centered around the variant/gene
  span when TSS-centered scoring is not possible but the target exon span plus
  both alleles still fits.
- ineligible variants are excluded before API submission and recorded in the
  eligibility CSV.

The runner performs a preflight request to confirm that the target gene and
target GTEx/Lung rows are present in the returned AlphaGenome output. It then
scores eligible variants in batches, filters returned rows to the target gene,
data source, and tissue, and writes both detailed filtered scores and compact
per-variant aggregates.

Important AlphaGenome outputs are:

```text
output/annotation/alphagenome/<locus_id>/<GENE>_alphagenome_filtered_scores.parquet
output/annotation/alphagenome/<locus_id>/<GENE>_alphagenome_variant_scores.csv
output/prelim/<locus_id>/<GENE>_alphagenome_variant_window_eligibility.csv
output/prelim/<locus_id>/<GENE>_alphagenome_interval_genes.csv
```

The compact variant score file contains:

```text
source_variant_id
alphagenome_raw_mean
alphagenome_quantile_mean
alphagenome_track_count
```

The geometric N20 `mu_0` prep uses `alphagenome_quantile_mean` as the signed
annotation source.

## 6. Annotation Standardization And Scaling

The geometric N20 annotation prep is implemented in:

```text
scripts/prepare_geometric_n20_mu0.py
```

It writes:

```text
output/susine_mu0/geometric_n20/mu0_locus_summary.csv
output/susine_mu0/geometric_n20/mu0_variant_table.parquet
output/susine_mu0/geometric_n20/mu0_variant_annotations_all.csv
output/susine_mu0/geometric_n20/per_locus_annotations/<GENE>_mu0_variant_annotations.csv
```

For each locus, the prep script aligns AlphaGenome scores to the LD order. The
join key is:

```text
LD order id == AlphaGenome source_variant_id
```

Missing AlphaGenome aggregate scores are treated as neutral annotation evidence:

```text
alphagenome_quantile_mean_filled = 0
annotation_missing = TRUE
```

Let `q_j` be the filled signed AlphaGenome quantile score for variant `j`.
The transform is:

```text
q*_j = clip(q_j, -1 + 1e-4, 1 - 1e-4)
a_raw_j = qnorm((q*_j + 1) / 2)
a_clip_j = clip(a_raw_j, -2.5, 2.5)
a_j = a_clip_j / sqrt(mean(a_clip^2))
```

The annotation is RMS-normalized within locus and is not mean-centered, because
zero has a semantic meaning: no directional AlphaGenome evidence.

The standardized RSS calibration target is computed from the z-only path. For
variant `j`:

```text
adj_j = (n_j - 1) / (z_j^2 + n_j - 2)
beta_hat_std_j = z_j * sqrt(adj_j / (n_j - 1))
shat2_std_j = 1 / (n_j - 1)
```

The per-locus baseline annotation scale is:

```text
c_rms_l = sqrt(mean(beta_hat_std^2))
c_cap_l = quantile(abs(beta_hat_std), 0.95) / max(abs(a))
baseline_c_l = min(c_rms_l, c_cap_l)
mu0_j = baseline_c_l * a_j
```

The safety cap prevents the RMS scale from assigning a prior mean whose maximum
absolute value is too large relative to the observed marginal standardized
effect scale.

## 7. Model Handoff To SuSiE And SuSiNE

The canonical downstream inputs are:

```text
z: GTEx marginal z-score vector in LD order
R: 1000 Genomes EUR LD correlation matrix in the same order
n: rounded or median effective sample size
a: RMS-normalized AlphaGenome annotation vector in the same order
mu_0: c * a
```

Vanilla SuSiE receives only the summary-statistic inputs:

```r
susieR::susie_rss(
  z = z,
  R = R,
  n = n,
  L = 10
)
```

Annotation-informed SuSiNE RSS receives the same `z`, `R`, and `n`, plus the
annotation-derived prior mean:

```r
susine::susine_rss(
  z = z,
  R = R,
  n = n,
  L = 10,
  mu_0 = c_value * a,
  sigma_0_2 = sigma_0_2,
  prior_update_method = "none"
)
```

The `test_susine` real-data pipeline loads the per-locus annotation CSV through
`load_real_data_locus_bundle()`. It validates that the annotation file, LD order,
master variant table, and LD matrix all use the same variant order. The bundle
then exposes:

```text
bundle$z
bundle$R
bundle$n_sample
bundle$a
bundle$baseline_c_l
bundle$variant_map
```

The real-data task runner passes:

```r
mu_0 = as.numeric(run_row$c_value) * bundle$a
```

to SuSiNE RSS. Vanilla SuSiE anchor runs use the same `z`, `R`, and `n` but no
annotation vector.

For the geometric N20 handoff, sync the prepared annotation files into
`test_susine` with:

```r
test_susine::sync_real_data_inputs(
  # example -- replace with your clone path
  source_repo_root = "/path/to/eQTL_annotations_for_susine",
  dest_root = here::here("data", "real_case_studies", "geometric_n20_loci"),
  source_mu0_name = "geometric_n20"
)
```

Then build and run the real-data job config against the resulting
`locus_manifest.csv`.

## 8. Reproducibility Commands

Lock the official sample after running the baseline exploration workbook:

```bash
python scripts/lock_geometric_n20_sample.py
```

Run AlphaGenome scoring:

```bash
python scripts/run_annotation_batch.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --summary-path output/annotation/alphagenome/geometric_n20_annotation_batch_summary.csv
```

Prepare standardized annotation and `mu_0` files:

```bash
python scripts/prepare_geometric_n20_mu0.py \
  --manifest config/loci_manifest_sample_100_per_chrom.csv \
  --selection config/annotation_selection_geometric_n20.csv \
  --output-dir output/susine_mu0/geometric_n20
```

Expected validation checks:

- `locked_loci.csv` has exactly 20 unique loci and `geometric_sample_order`
  equals `1..20`.
- all locked loci appear in `config/loci_manifest_sample_100_per_chrom.csv`.
- the AlphaGenome batch summary reports each selected locus as `completed` or
  `skipped_existing`.
- each per-locus annotation file has `variant_id` exactly matching LD order.
- `test_susine:::load_real_data_locus_bundle()` can load at least one synced
  geometric N20 locus without schema changes.
