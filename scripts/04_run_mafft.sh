#!/bin/bash
# 04_run_mafft.sh
# Align each candidate's group fasta (candidate + top 5 Arabidopsis hits) with MAFFT (L-INS-i)
cd ../dataset/msa/groups

mkdir -p ../aligned
for f in *.fasta; do
    name=$(basename "$f")
    mafft --localpair --maxiterate 1000 "$f" > "../aligned/${name%.fasta}_aligned.fasta"
done