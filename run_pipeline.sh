#!/bin/bash
#SBATCH --job-name=rnaseq_pipeline
#SBATCH --output=logs/pipeline_%j.log
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --account=r00270
#SBATCH --partition=general
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=saipamid@iu.edu 

# Exit immediately on any error
set -euo pipefail

# Create log directory if it doesn't exist
mkdir -p logs

# ---------------------------------------------
# Load conda and activate Nextflow environment
# ---------------------------------------------
module load conda
conda activate nextflow_env

# Load Apptainer for containers
module load apptainer

# ---------------------------------------------
# Project directory (confirmed path)
# ---------------------------------------------
cd /N/slate/bkota/slate_project

echo "========================================"
echo "Starting RNA-seq pipeline"
echo "Date: $(date)"
echo "Project dir: $(pwd)"
echo "Nextflow version: $(nextflow -version)"
echo "========================================"

# Run Nextflow pipeline
nextflow run main.nf \
    -resume \
    -with-report logs/report.html \
    -with-trace logs/trace.txt \
    -with-timeline logs/timeline.html

echo "========================================"
echo "Pipeline finished"
echo "Date: $(date)"
echo "========================================"
