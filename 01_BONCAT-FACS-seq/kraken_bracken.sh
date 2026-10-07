#!/usr/bin/env bash
set -euo pipefail

# Kraken2/Bracken taxonomic classification of BONCAT-FACS-seq reads
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"
#
# Kraken2 database: GTDB release 226
# Kraken2 confidence threshold: 0.05
# Bracken read length: 150 bp
# Bracken minimum abundance threshold: 10 reads

# -----------------------------
# User settings
# -----------------------------

DB="path/to/kraken_gtdb_r226_128g"
READS="trimmed_reads"
OUT="kraken_bracken"
THREADS=24
CONF=0.05
READLEN=150

mkdir -p "${OUT}/kraken" "${OUT}/bracken" "${OUT}/logs"

# -----------------------------
# Check for paired reads
# -----------------------------

shopt -s nullglob
R1S=( "${READS}"/*_paired_R1.fastq.gz )
shopt -u nullglob

if [[ ${#R1S[@]} -eq 0 ]]; then
    echo "ERROR: No *_paired_R1.fastq.gz files found in ${READS}"
    exit 1
fi

echo "Found ${#R1S[@]} samples."

# -----------------------------
# Build Bracken database
# for 150-bp reads if necessary
# -----------------------------

KMER=$(kraken2-inspect --db "${DB}" 2>/dev/null |
       awk -F: '/k-mer length/{gsub(/ /,""); print $2}')

KMER=${KMER:-35}

if [[ ! -f "${DB}/database${READLEN}mers.kmer_distrib" ]]; then
    echo "Building Bracken database for ${READLEN}-bp reads..."
    bracken-build \
        -d "${DB}" \
        -t "${THREADS}" \
        -k "${KMER}" \
        -l "${READLEN}"
fi

# -----------------------------
# Kraken2 classification
# -----------------------------

for R1 in "${R1S[@]}"; do

    R2="${R1/_paired_R1/_paired_R2}"
    SAMPLE=$(basename "${R1}" _paired_R1.fastq.gz)

    if [[ ! -f "${R2}" ]]; then
        echo "WARNING: Missing R2 for ${SAMPLE}; skipping."
        continue
    fi

    echo "Running Kraken2: ${SAMPLE}"

    kraken2 \
        --db "${DB}" \
        --threads "${THREADS}" \
        --paired \
        --gzip-compressed \
        --confidence "${CONF}" \
        --report "${OUT}/kraken/${SAMPLE}.kreport" \
        --output "${OUT}/kraken/${SAMPLE}.kraken" \
        "${R1}" "${R2}" \
        2> "${OUT}/logs/${SAMPLE}.kraken2.stderr"

done

# -----------------------------
# Bracken genus-level estimates
# -----------------------------

for KREPORT in "${OUT}"/kraken/*.kreport; do

    SAMPLE=$(basename "${KREPORT}" .kreport)

    echo "Running Bracken (genus): ${SAMPLE}"

    bracken \
        -d "${DB}" \
        -i "${KREPORT}" \
        -o "${OUT}/bracken/${SAMPLE}.bracken_genus.txt" \
        -w "${OUT}/bracken/${SAMPLE}.bracken_genus.report" \
        -r "${READLEN}" \
        -l G \
        -t 10

done

# -----------------------------
# Bracken phylum-level estimates
# -----------------------------

for KREPORT in "${OUT}"/kraken/*.kreport; do

    SAMPLE=$(basename "${KREPORT}" .kreport)

    echo "Running Bracken (phylum): ${SAMPLE}"

    bracken \
        -d "${DB}" \
        -i "${KREPORT}" \
        -o "${OUT}/bracken/${SAMPLE}.bracken_phylum.txt" \
        -w "${OUT}/bracken/${SAMPLE}.bracken_phylum.report" \
        -r "${READLEN}" \
        -l P \
        -t 10

done

# -----------------------------
# Generate combined summary tables
# -----------------------------

GENUS_SUMMARY="${OUT}/bracken/all_genus_summary.tsv"
PHYLUM_SUMMARY="${OUT}/bracken/all_phylum_summary.tsv"

printf "sample\ttaxon\test_reads\tfraction\n" > "${GENUS_SUMMARY}"
printf "sample\ttaxon\test_reads\tfraction\n" > "${PHYLUM_SUMMARY}"

for FILE in "${OUT}"/bracken/*.bracken_genus.txt; do
    SAMPLE=$(basename "${FILE}" .bracken_genus.txt)

    awk -v S="${SAMPLE}" \
        'NR > 1 {print S "\t" $1 "\t" $6 "\t" $7}' \
        "${FILE}" >> "${GENUS_SUMMARY}"
done

for FILE in "${OUT}"/bracken/*.bracken_phylum.txt; do
    SAMPLE=$(basename "${FILE}" .bracken_phylum.txt)

    awk -v S="${SAMPLE}" \
        'NR > 1 {print S "\t" $1 "\t" $6 "\t" $7}' \
        "${FILE}" >> "${PHYLUM_SUMMARY}"
done

echo "Kraken2/Bracken analysis complete."
echo "Genus summary:  ${GENUS_SUMMARY}"
echo "Phylum summary: ${PHYLUM_SUMMARY}"
