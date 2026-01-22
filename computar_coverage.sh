#!/bin/bash
set -euo pipefail

# --- Verifica parámetros ---
if [ "$#" -lt 4 ]; then
    echo "Uso: $0 <BAM> <KIT_BED> <FASTA_REF> <OUTFILE>"
    exit 1
fi

BAM=$1
KIT_BED=$2
FASTA=$3
OUTFILE=$4

# --- Archivos intermedios ---
FAI="${FASTA}.fai"
KITSORTED="kit_sorted.bed"
ONTARGET_BAM="ontarget.bam"
SORTBAM="ontarget_sorted.bam"
SORTBAM_BAI="${SORTBAM}.bai"
GENOMESIZES="sizes.genome"

THREADS=16
TMP_PREFIX="temp_sort_$$"

########################################
# 1. Index del FASTA
########################################
if [ ! -f "$FAI" ]; then
    echo "Indexando FASTA..."
    samtools faidx "$FASTA"
else
    echo "FASTA ya indexado, se omite."
fi

########################################
# 2. Ordenar BED del kit
########################################
if [ ! -f "$KITSORTED" ]; then
    echo "Ordenando BED del kit..."
    sort -k1,1 -k2,2n -k3,3n "$KIT_BED" > "$KITSORTED"
else
    echo "BED del kit ya ordenado, se omite."
fi

########################################
# 3. Recortar BAM al kit (on-target)
########################################
if [ ! -f "$ONTARGET_BAM" ]; then
    echo "Recortando BAM al kit de exoma..."
    samtools view -b -h -L "$KITSORTED" "$BAM" > "$ONTARGET_BAM"
else
    echo "BAM on-target ya existe, se omite."
fi

########################################
# 4. Ordenar BAM on-target
########################################
if [ ! -f "$SORTBAM" ]; then
    echo "Ordenando BAM on-target..."
    samtools sort -@ "$THREADS" -T "$TMP_PREFIX" -o "$SORTBAM" "$ONTARGET_BAM"
else
    echo "BAM on-target ya ordenado, se omite."
fi

########################################
# 5. Indexar BAM
########################################
if [ ! -f "$SORTBAM_BAI" ]; then
    echo "Indexando BAM on-target..."
    samtools index "$SORTBAM"
else
    echo "Índice BAM ya existe, se omite."
fi

########################################
# 6. Genome sizes (orden real del BAM)
########################################
if [ ! -f "$GENOMESIZES" ]; then
    echo "Generando genome sizes desde BAM..."
    samtools idxstats "$SORTBAM" | cut -f1,2 > "$GENOMESIZES"
else
    echo "Genome sizes ya existe, se omite."
fi

########################################
# 7. Cobertura
########################################
if [ ! -f "$OUTFILE" ]; then
    echo "Calculando cobertura con bedtools..."
    bedtools coverage \
        -a "$KITSORTED" \
        -b "$SORTBAM" \
        -g "$GENOMESIZES" \
        -sorted \
        -hist > "$OUTFILE"
else
    echo "Archivo de cobertura ya existe, se omite."
fi

echo "¡Pipeline finalizado correctamente!"
