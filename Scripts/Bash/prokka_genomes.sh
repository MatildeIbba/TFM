#!/bin/bash
# Valencia, 13/02/2026
# This script was wrote by Matilde Ibba
# Simple script for analyzing genome of C.Acnes with Prokka

for genome in *.fasta
do
 base=$(basename $genome .fasta)
	prokka $genome --outdir prokka_$base --prefix $base --genus Cutibacterium --species acnes --usegenus --cpus 2
done
