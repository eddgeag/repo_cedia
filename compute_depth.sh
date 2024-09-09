#!/bin/bash
faidx ~/datos_exomas/datos_gatk/hg38/hg38.fasta -i chromsizes >  chrom.sizes
awk '/^chr[0-9XY]*\t/ {printf("%s\t0\t%s\n",$1,$2);}' ~/datos_exomas/datos_gatk/hg38/hg38.fasta.fai > bed_split.bed
samtools view -b -h -L bed_split.bed $1 > out.bam
awk '/^chr[0-9XY]*\t/ {printf("%s\t%s\n",$1,$2);}' chrom.sizes > sizes.genome
cat sizes.genome|sort -V > sizes.genome.sort
echo "sortear el bam a 16 gb" 
samtools sort out.bam -@ 32 -o myfile_sorted.bam
echo "sorteando xgen"
cat xgen-exome.bed| sort -k1,1 -k2,2n -k3,3n > xgen-exome_sorted.bed
echo "computando cobertura"
bedtools coverage -a xgen-exome_sorted.bed  -b myfile_sorted.bam -g sizes.genome -sorted -hist 
