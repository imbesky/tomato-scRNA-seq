#!/bin/bash
# 02_build_blast_db.sh
# Build a BLAST database from the Arabidopsis protein sequences
cd ../dataset/blast

makeblastdb -in GCF_000001735.4_TAIR10.1_protein.faa -dbtype prot -out araport_db -parse_seqids