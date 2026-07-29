#!/bin/bash
set -euo pipefail

# Run this on a node with internet access and 'tabix' installed.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"
OUTPUT_DIR="${PROJECT_ROOT}/output"

mkdir -p "${OUTPUT_DIR}"

tabix -h \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  17:79530699-80073804 \
  > "${OUTPUT_DIR}/1kg_chr17_79530699_80073804.vcf"
