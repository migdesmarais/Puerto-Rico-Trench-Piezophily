#!/usr/bin/env bash
set -euo pipefail

# CoverM quantification of reads recruited to Puerto Rico Trench MAGs

THREADS=12
OUT_DIR="mag_recruitment"

coverm genome \
    --bam-files "${OUT_DIR}"/bam/*.bam \
    --genome-fasta-directory "${OUT_DIR}/renamed_MAGs" \
    --genome-fasta-extension fa \
    --methods covered_bases covered_fraction mean rpkm relative_abundance \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 75 \
    --output-file "${OUT_DIR}/mag_coverage_summary.tsv" \
    --threads "${THREADS}"

echo "CoverM quantification complete."
