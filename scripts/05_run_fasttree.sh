#!/bin/bash
# 05_run_fasttree.sh
# Build a phylogenetic tree (Newick format) from the aligned fasta files
cd ../dataset/msa/aligned

mkdir -p ../../trees
for f in *.fasta; do
    name=$(basename "$f")
    fasttree "$f" > "../../trees/${name%_aligned.fasta}.nwk"
done