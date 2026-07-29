from __future__ import annotations

import sys
import types
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
project_root_str = str(PROJECT_ROOT)
if project_root_str not in sys.path:
    sys.path.insert(0, project_root_str)


if "alphagenome" not in sys.modules:
    alphagenome_module = types.ModuleType("alphagenome")
    data_module = types.ModuleType("alphagenome.data")
    genome_module = types.ModuleType("alphagenome.data.genome")
    models_module = types.ModuleType("alphagenome.models")
    dna_client_module = types.ModuleType("alphagenome.models.dna_client")
    variant_scorers_module = types.ModuleType("alphagenome.models.variant_scorers")

    class Interval:
        def __init__(self, chromosome, start, end, strand=".", name=""):
            self.chromosome = chromosome
            self.start = int(start)
            self.end = int(end)
            self.strand = strand
            self.name = name

        def resize(self, length: int):
            center = (self.start + self.end) / 2.0
            new_start = int(round(center - (length / 2.0)))
            new_end = new_start + int(length)
            return Interval(self.chromosome, new_start, new_end, self.strand, self.name)

    class Variant:
        def __init__(self, chromosome, position, reference_bases, alternate_bases):
            self.chromosome = chromosome
            self.position = position
            self.reference_bases = reference_bases
            self.alternate_bases = alternate_bases

    dna_client_module.SUPPORTED_SEQUENCE_LENGTHS = {"SEQUENCE_LENGTH_1MB": 1_048_576}

    class OutputType:
        RNA_SEQ = "RNA_SEQ"

    dna_client_module.OutputType = OutputType
    dna_client_module.create = lambda api_key: None
    variant_scorers_module.GeneMaskLFCScorer = lambda requested_output=None: ("GeneMaskLFCScorer", requested_output)
    variant_scorers_module.tidy_scores = lambda scores, match_gene_strand=True: None

    genome_module.Interval = Interval
    genome_module.Variant = Variant
    data_module.genome = genome_module
    models_module.dna_client = dna_client_module
    models_module.variant_scorers = variant_scorers_module
    alphagenome_module.data = data_module
    alphagenome_module.models = models_module

    sys.modules["alphagenome"] = alphagenome_module
    sys.modules["alphagenome.data"] = data_module
    sys.modules["alphagenome.data.genome"] = genome_module
    sys.modules["alphagenome.models"] = models_module
    sys.modules["alphagenome.models.dna_client"] = dna_client_module
    sys.modules["alphagenome.models.variant_scorers"] = variant_scorers_module


try:
    import pyarrow  # noqa: F401  -- prefer the real library when it is installed
except ImportError:
    pyarrow_module = types.ModuleType("pyarrow")
    pyarrow_module.__version__ = "0.0.0"
    parquet_module = types.ModuleType("pyarrow.parquet")

    class _DummyArray:  # pandas 2.x accesses pyarrow.Array during DataFrame construction
        pass

    pyarrow_module.Array = _DummyArray

    class _DummyTable:
        def __init__(self, num_rows=0, schema=None):
            self.num_rows = num_rows
            self.schema = schema

    parquet_module.read_table = lambda *args, **kwargs: _DummyTable()

    class _DummyWriter:
        def __init__(self, *args, **kwargs):
            pass

        def write_table(self, table):
            return None

        def close(self):
            return None

    parquet_module.ParquetWriter = _DummyWriter
    pyarrow_module.parquet = parquet_module
    sys.modules["pyarrow"] = pyarrow_module
    sys.modules["pyarrow.parquet"] = parquet_module
