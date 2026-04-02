#!/bin/bash
# Download reference files required by this pipeline
# Run this script once before running main.nf
# Creates: reference/ directory with genome, GTF, and STAR index

set -euo pipefail
mkdir -p reference

echo "Downloading GRCh38 genome..."
wget -P reference/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
gunzip reference/GRCh38.primary_assembly.genome.fa.gz

echo "Downloading GENCODE v44 annotation..."
wget -P reference/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip reference/gencode.v44.annotation.gtf.gz

echo "Done. Now build STAR index:"
echo "  sbatch build_star_index.slurm"
