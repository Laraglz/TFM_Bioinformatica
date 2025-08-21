#!/usr/bin/bash
#
# Ruta al ref_genome
ref_genome=/Users/tu_usuario/Documents/mi_proyecto

# Ruta a los archivos
archivos=/Users/tu_usuario/Documents/mi_proyecto

# Crear un directorio de salida
mkdir -p ../mapped

# Iterar sobre todos los archivos pares
for archivo1 in $archivos/*_001_val_1.fq.gz; do
	archivo2="${archivo1/_001_val_1.fq.gz/_002_val_2.fq.gz}"	
        if [ -f "$archivo2" ]; then

		# Ejecutar HISAT2
		nombre_base=$(basename "${archivo1%.part_001_val_1.fq.gz}")
		salida=../mapped/"${nombre_base}.sam"
		salida_bam=../mapped/"${nombre_base}.bam"
		hisat2 -k1 -x $ref_genome -1 $archivo1 -2 $archivo2 -S $salida
		
		# Convertir el archivo SAM a BAM
		samtools view -bS $salida > $salida_bam

		# Ordenar el archivo BAM
		samtools sort $salida_bam -o ../mapped/${nombre_base}_sorted.bam

		# Indexar el archivo BAM
		samtools index ../mapped/${nombre_base}_sorted.bam
		
		# Eliminar el archivo SAM
		rm $salida
	fi
done	
