from __future__ import annotations

DEFAULT_REFERENCE_GENOME = "hg38"
DEFAULT_MAF_MIN = 0.01
DEFAULT_MAF_MAX = 0.99
DEFAULT_MIN_SAMPLE_SIZE = 50

DEFAULT_ALPHAGENOME_SEQUENCE_LENGTH = "1MB"
DEFAULT_ALPHAGENOME_TARGET_DATA_SOURCE = "gtex"
DEFAULT_ALPHAGENOME_BATCH_SIZE = 100
DEFAULT_ALPHAGENOME_MAX_WORKERS = 8
DEFAULT_ALPHAGENOME_RETRY_WAIT_SECONDS = 5

DEFAULT_ENABLED = True

MANIFEST_REQUIRED_COLUMNS = [
    "locus_id",
    "gene_name",
    "gene_id",
    "gtex_tissue",
    "gtex_chrom",
]

MANIFEST_OPTIONAL_COLUMNS = [
    "reference_genome",
    "maf_min",
    "maf_max",
    "min_sample_size",
    "alphagenome_sequence_length",
    "alphagenome_target_data_source",
    "alphagenome_target_gtex_tissue",
    "alphagenome_batch_size",
    "alphagenome_max_workers",
    "alphagenome_retry_wait_seconds",
    "enabled",
]

ANNOTATION_SELECTION_REQUIRED_COLUMNS = [
    "locus_id",
    "annotate",
]

ANNOTATION_SELECTION_OPTIONAL_COLUMNS = [
    "notes",
    "priority",
    "annotation_gene_name_override",
    "annotation_gene_id_override",
    "annotation_tissue_override",
]
