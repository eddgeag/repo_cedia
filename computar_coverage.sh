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
OUTFILE=$4   # viene con path .../metrics/XXX_coverage.hist.txt

# =================================================
# Derivar directorio y prefijo CORRECTOS
# =================================================
OUT_DIR="$(dirname "$OUTFILE")"
mkdir -p "$OUT_DIR"

BASENAME="$(basename "$OUTFILE")"
OUT_PREFIX="${OUT_DIR}/${BASENAME%.hist.txt}"

# =================================================
# Archivos de salida (TODOS en metrics/)
# =================================================
DEPTH_TXT="${OUT_PREFIX}.depth.txt"
DEPTH_STATS="${OUT_PREFIX}.depth.stats.txt"
DEPTH_HIST="${OUT_PREFIX}.depth.hist.tsv"

THREADS=16

echo "=== Cobertura on-target (samtools depth) ==="
echo "Output dir: $OUT_DIR"

########################################
# 1. Cobertura por base (streaming)
########################################
if [ ! -f "$DEPTH_TXT" ]; then
    echo "[1/3] Calculando profundidad por base..."
    samtools depth \
        -b "$KIT_BED" \
        "$BAM" > "$DEPTH_TXT"
else
    echo "[1/3] Depth por base ya existe."
fi

########################################
# 2. Métricas clínicas
########################################
if [ ! -f "$DEPTH_STATS" ]; then
    echo "[2/3] Calculando métricas clínicas..."
    awk '
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
    ' "$DEPTH_TXT" > "$DEPTH_STATS"
else
    echo "[2/3] Métricas clínicas ya existen."
fi

########################################
# 3. Histograma ligero (para gráficas)
########################################
if [ ! -f "$DEPTH_HIST" ]; then
    echo "[3/3] Generando histograma ligero..."
    awk '{ h[$3]++ }
         END { for (d in h) printf("%d\t%d\n", d, h[d]) }' \
      "$DEPTH_TXT" | sort -n > "$DEPTH_HIST"
else
    echo "[3/3] Histograma ya existe."
fi

echo "=== Cobertura finalizada correctamente ==="
echo "Depth por base:   $DEPTH_TXT"
echo "Stats clínicas:   $DEPTH_STATS"
echo "Histograma:      $DEPTH_HIST"
