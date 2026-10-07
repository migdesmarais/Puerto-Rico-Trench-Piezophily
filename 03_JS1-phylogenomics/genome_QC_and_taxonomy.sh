#!/usr/bin/env bash
set -euo pipefail

# Genome quality assessment and taxonomic classification of
# JS1/Atribacterota genomes.
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"

THREADS=24
GENOME_DIR="genomes"
OUTDIR="qc_taxonomy"

mkdir -p "${OUTDIR}/checkm"
mkdir -p "${OUTDIR}/gtdbtk"

# ---------------------------------------------------------
# CheckM genome quality assessment
# ---------------------------------------------------------

checkm lineage_wf \
    -x fna \
    -t "${THREADS}" \
    "${GENOME_DIR}" \
    "${OUTDIR}/checkm"

checkm qa \
    "${OUTDIR}/checkm/lineage.ms" \
    "${OUTDIR}/checkm" \
    -o 2 \
    > "${OUTDIR}/checkm/qa.tsv"

# ---------------------------------------------------------
# GTDB-Tk taxonomic classification
# ---------------------------------------------------------

gtdbtk classify_wf \
    --genome_dir "${GENOME_DIR}" \
    --out_dir "${OUTDIR}/gtdbtk" \
    --extension fna \
    --cpus "${THREADS}"

echo "Genome quality assessment and taxonomy complete."
