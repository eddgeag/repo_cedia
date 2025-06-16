#!/bin/bash
# Uso: bash script.sh <BAM> <KIT_BED> <FASTA_REF>
# Ejemplo: bash script.sh sample.bam xgen-exome.bed ~/datos_exomas/datos_gatk/hg38/hg38.fasta

# Control de errores y argumentos
if [ "$#" -lt 4 ]; then
    echo "Uso: $0 <BAM> <KIT_BED> <FASTA_REF> <OUTFILE>"
    exit 1
fi

BAM=$1
KIT_BED=$2
FASTA=$3
OUTFILE=$4

# Chequeo de existencia de archivos
for f in "$BAM" "$KIT_BED" "$FASTA"; do
    if [ ! -f "$f" ]; then
        echo "Archivo no encontrado: $f"
        exit 1
    fi
done

# 1. Obtener chrom.sizes y bed_split.bed
echo "Generando chrom.sizes..."
faidx "$FASTA" -i chromsizes >  chrom.sizes
echo "Generando bed_split.bed para cada cromosoma..."

awk '/^chr[0-9XY]*\t/ {printf("%s\t0\t%s\n",$1,$2);}' "${FASTA}.fai" > bed_split.bed
# 2. Extraer regiones del BAM (por cromosoma)
echo "Extrayendo regiones del BAM..."
samtools view -b -h -L bed_split.bed $BAM > out.bam
# 3. Preparar archivo genome de tamaños (para bedtools)
cat sizes.genome|sort -V > sizes.genome.sort

echo "sortear el bam a 16 gb" 
samtools sort out.bam -@ 16 -o myfile_sorted.bam

# 5. Ordenar el archivo BED del kit
echo "Ordenando el archivo BED del kit..."

cat "$KIT_BED"| sort -k1,1 -k2,2n -k3,3n > kit_sorted.bed

echo "Calculando cobertura con bedtools..."

bedtools coverage -a kit_sorted.bed  -b myfile_sorted.bam -g sizes.genome.sort -sorted -hist > "$OUTFILE"

echo "¡Listo! Resultado guardado en $OUTFILE"
