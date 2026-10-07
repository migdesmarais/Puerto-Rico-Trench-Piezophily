#!/usr/bin/env bash
set -euo pipefail

# Pairwise average amino acid identity among JS1 genomes
# CompareM v0.1.2

PROTEIN_DIR="predicted_proteins"
OUTDIR="comparem_aai"

mkdir -p "${OUTDIR}"

comparem aai_wf \
    "${PROTEIN_DIR}" \
    "${OUTDIR}" \
    --cpus 24
