#!/usr/bin/env bash
set -euo pipefail

# Recruitment of BONCAT-FACS-seq reads to Puerto Rico Trench MAGs
#
# Desmarais et al.
# "Piezophilic activity dominates the Puerto Rico Trench
# deep subseafloor biosphere"

THREADS=12
MAPQ=20

READS_DIR="../01_BONCAT-FACS-seq/trimmed_reads"
MAG_DIR="path/to/MAGs"
OUT_DIR="mag_recruitment"

REF="${OUT_DIR}/all_MAGs_unique.fa"
IDX="${OUT_DIR}/all_MAGs_index"

mkdir -p "${OUT_DIR}"/{renamed_MAGs,bam,logs,counts}

# ---------------------------------------------------------
# Prepare MAG reference
# ---------------------------------------------------------

for MAG in "${MAG_DIR}"/*.fa; do
    ID=$(basename "${MAG}" .fa)

    awk -v prefix="${ID}" \
        '/^>/{print ">" prefix "_" substr($0,2)} !/^>/' \
        "${MAG}" > "${OUT_DIR}/renamed_MAGs/${ID}.fa"
done

cat "${OUT_DIR}"/renamed_MAGs/*.fa > "${REF}"

# ---------------------------------------------------------
# Build Bowtie2 index
# ---------------------------------------------------------

bowtie2-build "${REF}" "${IDX}"
samtools faidx "${REF}"

# ---------------------------------------------------------
# Map paired-end reads
# ---------------------------------------------------------

for R1 in "${READS_DIR}"/*_paired_R1.fastq.gz; do

    R2="${R1/_paired_R1/_paired_R2}"
    SAMPLE=$(basename "${R1}" _paired_R1.fastq.gz)

    echo "Mapping ${SAMPLE}"

    bowtie2 \
        --very-sensitive-local \
        -p "${THREADS}" \
        -k 1 \
        -X 2000 \
        -x "${IDX}" \
        -1 "${R1}" \
        -2 "${R2}" \
        2> "${OUT_DIR}/logs/${SAMPLE}_bowtie2.log" \
    | samtools view \
        -h \
        -b \
        -q "${MAPQ}" \
        -F 0x904 \
    | samtools sort \
        -@ "${THREADS}" \
        -o "${OUT_DIR}/bam/${SAMPLE}.q20.primary.bam"

    samtools index \
        "${OUT_DIR}/bam/${SAMPLE}.q20.primary.bam"

    # Add MD/NM tags for identity calculations
    samtools calmd \
        -bAr \
        "${OUT_DIR}/bam/${SAMPLE}.q20.primary.bam" \
        "${REF}" \
        > "${OUT_DIR}/bam/${SAMPLE}.tmp.bam"

    mv \
        "${OUT_DIR}/bam/${SAMPLE}.tmp.bam" \
        "${OUT_DIR}/bam/${SAMPLE}.q20.primary.bam"

    samtools index \
        "${OUT_DIR}/bam/${SAMPLE}.q20.primary.bam"

    samtools idxstats \
        "${OUT_DIR}/bam/${SAMPLE}.q20.primary.bam" \
        > "${OUT_DIR}/counts/${SAMPLE}_idxstats.tsv"

done

echo "MAG read recruitment complete."
