#!/usr/bin/bash

# Ruta al archivo gtf
gtf=./Homo_sapiens.GRCh38.114.gtf

# Crear un directorio de salida
mkdir -p ./htseq_counts

# Ejecutar htseq para todos los archivos
for bam in /Users/tu_usuario/Documents/mi_proyecto/mapped/*_sorted.bam; do

	sample=$(basename "$bam" _sorted.bam)
	htseq-count -t exon -i gene_id --stranded=no -f bam -r pos -s no "$bam" "$gtf" > ./htseq_counts/"${sample}_counts.tsv"

done
