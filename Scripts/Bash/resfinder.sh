#!/bin/bash
# This script was wrote by Matilde Ibba
# ResFinder on original .fasta files

input_dir=$1

for f in "$input_dir"/*.fasta
do
	base=$(basename "$f" .fasta)
	abricate --db resfinder "$f" > "${base}_resfinder.tab"
done
