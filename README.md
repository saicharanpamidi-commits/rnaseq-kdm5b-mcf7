# RNA-seq Analysis: KDM5B Knockdown in MCF-7 Breast Cancer Cells

[![Nextflow](https://img.shields.io/badge/nextflow-25.10.4-brightgreen)](https://nextflow.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview
End-to-end bulk RNA-seq pipeline analyzing transcriptome-wide effects of KDM5B
knockdown in MCF-7 human breast cancer cells using Nextflow DSL2 on HPC.

## Biological Question
What genes and pathways are dysregulated when KDM5B (histone H3K4me3 demethylase)
is knocked down in MCF-7 breast cancer cells?

## Key Findings
- **5,906 differentially expressed genes** (34% of expressed transcriptome)
- **Warburg metabolic shift**: OXPHOS suppressed (NES -2.33, p=3.9e-14), glycolysis activated
- **TFF1** top upregulated gene (padj ~10^-300) — canonical estrogen receptor target
- **ER timing dysregulation**: early ER response suppressed, late ER targets de-repressed
- **NF-kB inflammatory activation** with paradoxical apoptosis suppression
- Automated strandedness detection confirmed **reverse-stranded library** (-s 2)

## Pipeline Architecture
```
FASTQ
  ↓
FastQC (raw QC)           [FastQC 0.12.1]
  ↓
FASTP (trimming)          [FASTP 0.23.4]
  ↓
STAR (alignment)          [STAR 2.7.11a → GRCh38]
  ↓
samtools (BAM indexing)   [samtools 1.17]
  ↓
RSeQC (strandedness)      [RSeQC 5.0.4 — auto-detects -s flag]
  ↓
featureCounts             [Subread 2.0.6 — GENCODE v44]
  ↓
MultiQC                   [MultiQC 1.27.1]
  ↓
DESeq2 + Pathway Analysis [GO | KEGG | Reactome | GSEA Hallmarks]
```

## Dataset
| Sample | Condition | Description |
|--------|-----------|-------------|
| SRR31992838 | shSCR | Scramble control replicate 1 |
| SRR31992839 | shSCR | Scramble control replicate 2 |
| SRR31992840 | shSCR | Scramble control replicate 3 |
| SRR31992841 | shKDM5B | KDM5B knockdown replicate 1 |
| SRR31992842 | shKDM5B | KDM5B knockdown replicate 2 |
| SRR31992843 | shKDM5B | KDM5B knockdown replicate 3 |

- **Cell line**: MCF-7 (ER+ luminal breast cancer)
- **Platform**: Illumina NovaSeq 6000, paired-end
- **Reference**: GRCh38 / GENCODE v44

## How to Run

### 1. Download references
```bash
bash scripts/download_references.sh
```

### 2. Build STAR index
```bash
# Edit build_star_index.slurm with your paths, then:
sbatch build_star_index.slurm
```

### 3. Configure pipeline
Edit `nextflow.config` — update account and partition for your HPC.

### 4. Submit pipeline
```bash
mkdir -p logs
sbatch run_pipeline.sh
```

### 5. Run DESeq2 analysis
```bash
Rscript scripts/deseq2_full_analysis.R
```

## Requirements
```bash
# Conda environment (Nextflow)
conda env create -f environment.yml
conda activate nextflow_env

# R packages
BiocManager::install(c("DESeq2","apeglm","org.Hs.eg.db",
                       "clusterProfiler","ReactomePA","enrichplot"))
install.packages(c("ggplot2","ggrepel","pheatmap","msigdbr","dplyr","tibble"))
```

## Repository Structure
```
.
├── README.md
├── environment.yml              # Conda environment
├── main.nf                      # Nextflow DSL2 pipeline (9 processes)
├── nextflow.config              # SLURM + resource configuration
├── run_pipeline.sh              # HPC sbatch submission script
├── scripts/
│   ├── deseq2_full_analysis.R   # DESeq2 + pathway analysis
│   ├── download_references.sh   # Reference genome download script
│   ├── samplesheet.csv          # Sample metadata
│   ├── run_deseq2.R
│   └── check_counts.R
└── results/
    ├── counts/
    │   ├── gene_counts.txt      # featureCounts output
    │   └── gene_counts.txt.summary
    ├── deseq2/
    │   ├── deseq2_results_annotated.csv
    │   ├── pca_plot.png
    │   ├── volcano_plot.png
    │   ├── ma_plot.png
    │   ├── heatmap_top50.png
    │   └── pathway/
    │       ├── gsea_hallmarks_results.csv
    │       ├── go_bp_results.csv
    │       ├── kegg_results.csv
    │       └── reactome_results.csv
    ├── multiqc/
    │   └── multiqc_report.html
    └── rseqc/
        └── strandedness/
            └── strandedness_consensus.txt
```

## HPC Notes
- Scheduler: SLURM
- Container runtime: Apptainer 1.3.6
- Several biocontainers lack `ps` (procps) — resolved by using HPC modules
  for STAR, samtools, FastQC, and featureCounts

## Author
**Sai Charan Pamidi**
Bioinformatics / Computational Biology

## License
MIT
