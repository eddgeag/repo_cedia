#!/bin/bash

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
CHROMSIZES="chrom.sizes"
BEDSPLIT="bed_split.bed"
OUTBAM="out.bam"
SORTBAM="myfile_sorted.bam"
KITSORTED="kit_sorted.bed"
GENOMESIZES="sizes.genome"

# 1. Index del fasta si no existe
if [ ! -f "$FAI" ]; then
    echo "Indexando $FASTA..."
    samtools faidx "$FASTA"
fi

# 2. chrom.sizes
cut -f1,2 "$FAI" > "$CHROMSIZES"

# 3. bed_split.bed
awk '{printf("%s\t0\t%s\n",$1,$2);}' "$CHROMSIZES" > "$BEDSPLIT"

# 4. out.bam (solo si NO existe)
if [ ! -f "$OUTBAM" ]; then
    echo "Extrayendo regiones del BAM..."
    samtools view -b -h -L "$BEDSPLIT" "$BAM" > "$OUTBAM"
else
    echo "Ya existe $OUTBAM, saltando extracción."
fi

# 5. myfile_sorted.bam (solo si NO existe)
if [ ! -f "$SORTBAM" ]; then
    echo "Ordenando BAM..."
    samtools sort -@ 16 -T temp_sort -o "$SORTBAM" "$OUTBAM"
else
    echo "Ya existe $SORTBAM, saltando ordenamiento."
fi

# 6. BED ordenado
sort -k1,1 -k2,2n -k3,3n "$KIT_BED" > "$KITSORTED"

# 7. Copiar sizes.genome (mismo orden que BAM)
cp "$CHROMSIZES" "$GENOMESIZES"

# 8. Cobertura
echo "Calculando cobertura con bedtools..."
bedtools coverage -a "$KITSORTED" -b "$SORTBAM" -g "$GENOMESIZES" -sorted -hist > "$OUTFILE"

echo "¡Listo! Resultado guardado en $OUTFILE"

