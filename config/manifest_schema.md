# Manifest Schema

## `loci_manifest.csv`

Required columns:

- `locus_id`: unique stable identifier used for per-locus output subfolders
- `gene_name`
- `gene_id`
- `gtex_tissue`
- `gtex_chrom`

Optional columns and defaults:

- `reference_genome`: default `hg38`
- `maf_min`: default `0.01`
- `maf_max`: default `0.99`
- `min_sample_size`: default `50`
- `alphagenome_sequence_length`: default `1MB`
- `alphagenome_target_data_source`: default `gtex`
- `alphagenome_target_gtex_tissue`: default value from `gtex_tissue`
- `alphagenome_batch_size`: default `100`
- `alphagenome_max_workers`: default `8`
- `alphagenome_retry_wait_seconds`: default `5`
- `enabled`: default `TRUE`

Validation rules:

- `locus_id` must be unique
- required columns must be non-empty for enabled rows
- `gtex_chrom` must match `^chr[0-9XYM]+$`
- `enabled` must parse as boolean
- numeric override fields must parse as numeric
- `maf_min <= maf_max`
- positive batch size / worker count

## `annotation_selection.csv`

Required columns:

- `locus_id`
- `annotate`

Optional columns:

- `notes`
- `priority`
- `annotation_gene_name_override`
- `annotation_gene_id_override`
- `annotation_tissue_override`

Behavior:

- only rows with `annotate == TRUE` are passed to the annotation batch runner
- `locus_id` values must already exist in `loci_manifest.csv`
- optional overrides supersede the main manifest only for annotation
