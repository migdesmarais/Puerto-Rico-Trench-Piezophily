#!/usr/bin/env bash
set -euo pipefail

# JS1 16S rRNA gene sequence clustering
#
# JS1-affiliated sequence variants are clustered at 99%
# nucleotide identity using VSEARCH.

INPUT="JS1_zotus.fa"
OUTDIR="vsearch_99"

mkdir -p "${OUTDIR}"

vsearch \
    --cluster_fast "${INPUT}" \
    --id 0.99 \
    --centroids "${OUTDIR}/JS1_centroids_99.fa" \
    --uc "${OUTDIR}/JS1_clusters_99.uc" \
    --threads 8

echo "JS1 99%-identity clustering complete."
