#!/bin/bash
set -euo pipefail

# =================================================
# Uso (NO CAMBIA)
# =================================================
if [ "$#" -lt 4 ]; then
    echo "Uso: $0 <BAM> <KIT_BED> <FASTA_REF> <OUTFILE>"
    exit 1
fi

BAM=$1
KIT_BED=$2
FASTA=$3
OUTFILE=$4   # <-- sigue siendo el .hist

# =================================================
# Archivos derivados (internos)
# =================================================
FAI="${FASTA}.fai"

KITSORTED="kit_sorted.bed"
ONTARGET_BAM="ontarget.bam"
SORTBAM="ontarget_sorted.bam"
SORTBAM_BAI="${SORTBAM}.bai"
GENOMESIZES="sizes.genome"

# Derivar prefijo desde OUTFILE (sin tocar R)
OUT_PREFIX="${OUTFILE%.hist.txt}"
DEPTH_STATS_OUT="${OUT_PREFIX}.depth.stats.txt"

THREADS=16
TMP_PREFIX="temp_sort_$$"

echo "=== Pipeline de cobertura on-target ==="

########################################
# 1. Index FASTA
########################################
if [ ! -f "$FAI" ]; then
    echo "[1/8] Indexando FASTA..."
    samtools faidx "$FASTA"
else
    echo "[1/8] FASTA ya indexado."
fi

########################################
# 2. Ordenar BED
########################################
if [ ! -f "$KITSORTED" ]; then
    echo "[2/8] Ordenando BED del kit..."
    sort -k1,1 -k2,2n -k3,3n "$KIT_BED" > "$KITSORTED"
else
    echo "[2/8] BED del kit ya ordenado."
fi

########################################
# 3. Recortar BAM al kit
########################################
if [ ! -f "$ONTARGET_BAM" ]; then
    echo "[3/8] Recortando BAM al kit..."
    samtools view -b -h -L "$KITSORTED" "$BAM" > "$ONTARGET_BAM"
else
    echo "[3/8] BAM on-target ya existe."
fi

########################################
# 4. Ordenar BAM
########################################
if [ ! -f "$SORTBAM" ]; then
    echo "[4/8] Ordenando BAM on-target..."
    samtools sort -@ "$THREADS" -T "$TMP_PREFIX" -o "$SORTBAM" "$ONTARGET_BAM"
else
    echo "[4/8] BAM on-target ya ordenado."
fi

########################################
# 5. Indexar BAM
########################################
if [ ! -f "$SORTBAM_BAI" ]; then
    echo "[5/8] Indexando BAM..."
    samtools index "$SORTBAM"
else
    echo "[5/8] Índice BAM ya existe."
fi

########################################
# 6. Genome sizes (para bedtools)
########################################
if [ ! -f "$GENOMESIZES" ]; then
    echo "[6/8] Generando genome sizes..."
    samtools idxstats "$SORTBAM" | cut -f1,2 > "$GENOMESIZES"
else
    echo "[6/8] Genome sizes ya existe."
fi

########################################
# 7. Histograma de cobertura (QC)
########################################
if [ ! -f "$OUTFILE" ]; then
    echo "[7/8] Calculando histograma de cobertura (bedtools)..."
    bedtools coverage \
        -a "$KITSORTED" \
        -b "$SORTBAM" \
        -g "$GENOMESIZES" \
        -sorted \
        -hist > "$OUTFILE"
else
    echo "[7/8] Histograma ya existe."
fi

########################################
# 8. Cobertura REAL estilo proveedor
#    (samtools depth, sin ceros)
########################################
if [ ! -f "$DEPTH_STATS_OUT" ]; then
    echo "[8/8] Calculando cobertura real (samtools depth)..."

    samtools depth \
        -b "$KITSORTED" \
        "$SORTBAM" \
    | awk '
        {
          sum += $3;
          n++;
          if ($3 >= 20) n20++;
          if ($3 >= 30) n30++;
        }
        END {
          printf("Mean_depth\t%.2f\n", (n>0 ? sum/n : 0));
          printf("Pct_ge_20X\t%.2f\n", (n>0 ? 100*n20/n : 0));
          printf("Pct_ge_30X\t%.2f\n", (n>0 ? 100*n30/n : 0));
          printf("Bases_used\t%d\n", n);
        }
    ' > "$DEPTH_STATS_OUT"

else
    echo "[8/8] Métricas depth ya existen."
fi

echo "=== Pipeline de cobertura finalizado ==="
echo "Histograma QC:      $OUTFILE"
echo "Cobertura real:     $DEPTH_STATS_OUT"
