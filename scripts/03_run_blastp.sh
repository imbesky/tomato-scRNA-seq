#!/bin/bash
# 03_run_blastp.sh
# Search tomato candidate gene sequences against the Arabidopsis BLAST database
cd ../dataset/blast

blastp \
  -query tomato_candidates.fasta \
  -db araport_db \
  -out blast_results.tsv \
  -outfmt 6 \
  -evalue 1e-5 \
  -max_target_seqs 5