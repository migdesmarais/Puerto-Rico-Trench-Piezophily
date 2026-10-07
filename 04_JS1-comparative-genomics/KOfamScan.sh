#!/usr/bin/env bash
set -euo pipefail

# Functional annotation of JS1/Atribacterota genomes with KOfamScan
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"

THREADS=24
GENOME_DIR="../03_JS1-phylogenomics/genomes_filtered"
PROTEIN_DIR="proteins"
OUTDIR="kofamscan"

mkdir -p "${PROTEIN_DIR}"
mkdir -p "${OUTDIR}"

# ---------------------------------------------------------
# Predict protein-coding sequences with Prodigal
# ---------------------------------------------------------

for GENOME in "${GENOME_DIR}"/*.fna; do

    ID=$(basename "${GENOME}" .fna)

    echo "Predicting proteins: ${ID}"

    prodigal \
        -i "${GENOME}" \
        -a "${PROTEIN_DIR}/${ID}.faa" \
        -d "${PROTEIN_DIR}/${ID}.fna" \
        -p single \
        -q

done

# ---------------------------------------------------------
# Annotate predicted proteins with KOfamScan
# ---------------------------------------------------------

for PROTEINS in "${PROTEIN_DIR}"/*.faa; do

    ID=$(basename "${PROTEINS}" .faa)

    echo "Running KOfamScan: ${ID}"

    exec_annotation \
        -o "${OUTDIR}/${ID}_kofam.tsv" \
        --cpu "${THREADS}" \
        --format mapper \
        "${PROTEINS}"

done

echo "KOfamScan annotation complete."
