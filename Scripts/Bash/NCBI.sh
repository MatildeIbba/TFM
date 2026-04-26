#!/bin/bash
# This script was wrote by Matilde Ibba
# NCBI AMRFinderPlus on original .fasta files

input_dir=$1

for f in "$input_dir"/*.fasta
do
	base=$(basename "$f" .fasta)
	abricate --db ncbi "$f" > "${base}_NCBI.tab"
done
