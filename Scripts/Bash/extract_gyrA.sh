#!/bin/bash
# This script was wrote by Matilde Ibba
# 1. Automatic extraction of gyrA from BLASTn
# 2. Extraction of the sequences taking into account the orientation

input_fasta=$1
blast_file=$2

# 1. Generate unique coordinates (best hit for contig)

awk '{
  contig=$2;
  start=($9<$10?$9:$10);
  end=($9>$10?$9:$10);
  strand=($9<$10?"+":"-");
  print contig"\t"start"\t"end"\t"strand"\t"$4
}' $blast_file | sort -k1,1 -k5,5nr | awk '!seen[$1]++' > gyrA_coords.tsv

# 2. Extraction of all the sequences

> gyrA_all.fna
while read contig start end strand length; do
  
  echo ">gyrA_${i}_${contig}_${strand}_${end}" >> gyrA_all.fna

  if [ "$strand" = "+" ]; then
      samtools faidx $input_fasta ${contig}:${start}-${end} | tail -n +2 >> gyrA_all.fna
  else
      samtools faidx $input_fasta ${contig}:${start}-${end} | tail -n +2 | tr "ACGTacgt" "TGCAtgca" | rev >> gyrA_all.fna
  fi

  i=$((i+1))

done < gyrA_coords.tsv
