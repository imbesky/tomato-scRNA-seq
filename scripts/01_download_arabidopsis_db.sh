#!/bin/bash
# 01_download_arabidopsis_db.sh
# Download the Arabidopsis reference protein sequences (NCBI TAIR10.1)
mkdir -p ../dataset/blast
cd ../dataset/blast

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz"
gunzip GCF_000001735.4_TAIR10.1_protein.faa.gz