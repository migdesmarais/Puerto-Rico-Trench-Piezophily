#!/usr/bin/env bash
set -euo pipefail

# Phylogenomic reconstruction of JS1/Atribacterota genomes.
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"

THREADS=24
GENOME_DIR="genomes_filtered"
OUTDIR="phylogeny"

mkdir -p "${OUTDIR}"

# Generate list of genomes retained for phylogenomic analysis
find "$(realpath "${GENOME_DIR}")" \
    -maxdepth 1 \
    -type f \
    -name "*.fna" \
    | sort \
    > "${OUTDIR}/genomes.list"

# Extract and align conserved bacterial single-copy genes
# using GToTree and the bacterial HMM set.
GToTree \
    -f "${OUTDIR}/genomes.list" \
    -H Bacteria \
    -B \
    -o "${OUTDIR}/gtotree" \
    -t "${THREADS}"

# Infer maximum-likelihood phylogeny using ModelFinder
# and 1,000 ultrafast bootstrap replicates.
iqtree \
    -s "${OUTDIR}/gtotree/Aligned_SCGs.faa" \
    -m MFP \
    -B 1000 \
    -T "${THREADS}" \
    --prefix "${OUTDIR}/JS1_iqtree"

echo "Phylogenomic reconstruction complete."
