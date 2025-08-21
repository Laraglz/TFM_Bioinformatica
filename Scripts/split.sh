#!usr/bin/bash

for file in *.fastq.gz; do
	seqkit split2 -p 2 "$file"
done
