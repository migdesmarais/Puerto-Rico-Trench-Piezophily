#!/usr/bin/env bash
set -euo pipefail

# Read quality control for BONCAT-FACS-seq libraries
# Desmarais et al., "Piezophilic activity dominates the Puerto Rico
# Trench deep subseafloor biosphere"

THREADS=12
RAW_DIR="path/to/raw_reads"
TRIM_DIR="trimmed_reads"
QC_DIR="quality_control"

mkdir -p "${TRIM_DIR}" "${QC_DIR}/raw" "${QC_DIR}/trimmed"

# Initial read quality assessment
fastqc --threads "${THREADS}" \
    --outdir "${QC_DIR}/raw" \
    "${RAW_DIR}"/*.fastq.gz

# Paired-end trimming
for R1 in "${RAW_DIR}"/*_R1_001.fastq.gz; do

    SAMPLE=$(basename "${R1}" _R1_001.fastq.gz)
    R2="${RAW_DIR}/${SAMPLE}_R2_001.fastq.gz"

    trimmomatic PE -phred33 -threads "${THREADS}" \
        "${R1}" "${R2}" \
        "${TRIM_DIR}/${SAMPLE}_paired_R1.fastq.gz" \
        "${TRIM_DIR}/${SAMPLE}_unpaired_R1.fastq.gz" \
        "${TRIM_DIR}/${SAMPLE}_paired_R2.fastq.gz" \
        "${TRIM_DIR}/${SAMPLE}_unpaired_R2.fastq.gz" \
        ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36
done

# Post-trimming QC of paired reads
fastqc --threads "${THREADS}" \
    --outdir "${QC_DIR}/trimmed" \
    "${TRIM_DIR}"/*_paired_R1.fastq.gz \
    "${TRIM_DIR}"/*_paired_R2.fastq.gz

# Aggregate QC reports
multiqc "${QC_DIR}" -o "${QC_DIR}/multiqc"
