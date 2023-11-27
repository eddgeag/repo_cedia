#!/bin/bash
faidx ~/datos_exomas/datos_gatk/referencia/ensembl/referencia_filtered.fasta -i chromsizes >  chrom.sizes
awk '/^chr[0-9XY]*\t/ {printf("%s\t0\t%s\n",$1,$2);}' ~/datos_exomas/datos_gatk/referencia/ensembl/referencia_filtered.fasta.fai > bed_split.bed
samtools view -L bed_split.bed -o out.bam $1
cat chrom.sizes|sort -V > sizes.genome.sort
bedtools coverage -a ~/datos_exomas/coverage/hg38/xgen-exome.bed  -b out.bam -g sizes.genome.sort -sorted -hist 
