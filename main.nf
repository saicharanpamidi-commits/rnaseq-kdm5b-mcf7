// =============================================================================
// RNA-seq Pipeline  Production Grade
// Processes: FASTQC, FASTP, STAR, INDEX_BAM, INFER_STRANDEDNESS,
//            READ_DISTRIBUTION, PARSE_STRANDEDNESS, FEATURECOUNTS, MULTIQC
// HPC: SLURM + Apptainer + HPC modules
// =============================================================================


// -----------------------------------------------------------------------------
// PROCESS 1: FASTQC  Raw read quality control
// Run on raw FASTQs before trimming to capture baseline quality
// Uses container (FastQC biocontainer includes ps)
// If ps error occurs: run `module spider fastqc` and switch to module approach
// -----------------------------------------------------------------------------
process FASTQC {

    tag "$sample_id"

    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    path "*.zip",                         emit: zip

    script:
    """
    source /etc/profile.d/modules.sh
    module load fastqc/0.12.1

    fastqc ${reads[0]} ${reads[1]} --threads ${task.cpus} --outdir .
    """
}

// -----------------------------------------------------------------------------
// PROCESS 2: FASTP  Adapter trimming and QC
// Produces trimmed reads + HTML/JSON QC reports for MultiQC
// -----------------------------------------------------------------------------
process FASTP {

    tag "$sample_id"

    container "quay.io/biocontainers/fastp:0.23.4--h125f33a_5"

    publishDir "results/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: trimmed
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_1.trimmed.fastq.gz \
        -O ${sample_id}_2.trimmed.fastq.gz \
        -h ${sample_id}.fastp.html \
        -j ${sample_id}.fastp.json \
        -w ${task.cpus}
    """
}


// -----------------------------------------------------------------------------
// PROCESS 3: STAR_ALIGN  Splice-aware alignment to reference genome
// Uses HPC module (biocontainer lacks ps on this cluster)
// -----------------------------------------------------------------------------
process STAR_ALIGN {

    tag "$sample_id"

    publishDir "results/star", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    source /etc/profile.d/modules.sh
    module load star/2.7.11a

    STAR \
        --genomeDir ${params.genome_index} \
        --readFilesIn ${read1} ${read2} \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${sample_id}.

    mv ${sample_id}.Aligned.sortedByCoord.out.bam ${sample_id}.sorted.bam
    """
}


// -----------------------------------------------------------------------------
// PROCESS 4: INDEX_BAM  Index sorted BAM files with samtools
// Emits BAM + BAI together so downstream processes stage both files
// Uses HPC module (samtools biocontainer lacks ps on this cluster)
// -----------------------------------------------------------------------------
process INDEX_BAM {

    tag "$sample_id"

    publishDir "results/star", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path(bam), path("${bam}.bai")

    script:
    """
    source /etc/profile.d/modules.sh
    module load gnu/12.2.0
    module load samtools/1.17

    samtools index $bam
    """
}


// -----------------------------------------------------------------------------
// PROCESS 5: INFER_STRANDEDNESS  Detect library strandedness per sample
// Uses RSeQC infer_experiment.py + BED12 gene model
// Output parsed to numeric strandedness: 0=unstranded, 1=forward, 2=reverse
// Uses conda rseqc_env (installed via pip to avoid libdeflate conflict)
// -----------------------------------------------------------------------------
process INFER_STRANDEDNESS {

    tag "$sample_id"

    publishDir "results/rseqc/strandedness", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed

    output:
    path "${sample_id}.infer_experiment.txt", emit: txt

    script:
    """
    source /etc/profile.d/modules.sh
    module load conda
    conda activate ${params.rseqc_env}

    infer_experiment.py \
        -r ${bed} \
        -i ${bam} \
        > ${sample_id}.infer_experiment.txt 2>&1
    """
}


// -----------------------------------------------------------------------------
// PROCESS 6: READ_DISTRIBUTION  Quantify where reads map in the genome
// Reports fraction of reads in CDS, UTR, intron, intergenic regions
// Output captured by MultiQC automatically
// -----------------------------------------------------------------------------
process READ_DISTRIBUTION {

    tag "$sample_id"

    publishDir "results/rseqc/read_distribution", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed

    output:
    path "${sample_id}.read_distribution.txt", emit: txt

    script:
    """
    source /etc/profile.d/modules.sh
    module load conda
    conda activate ${params.rseqc_env}

    read_distribution.py \
        -r ${bed} \
        -i ${bam} \
        > ${sample_id}.read_distribution.txt 2>&1
    """
}


// -----------------------------------------------------------------------------
// PROCESS 7: PARSE_STRANDEDNESS  Determine consensus strandedness from all samples
// Reads all infer_experiment outputs, votes across samples, emits single integer
// 0 = unstranded, 1 = forward stranded, 2 = reverse stranded
// Feeds directly into featureCounts -s flag
// -----------------------------------------------------------------------------
process PARSE_STRANDEDNESS {
    input:
    path txt_files

    output:
    path "strandedness_consensus.txt", emit: full_report
    env strand_int,                    emit: strand_int

    script:
    """
    python3 << 'PYEOF'
import os
files = "${txt_files}".split()
votes = []
for f in files:
    fwd = rev = 0.0
    with open(f) as fh:
        for line in fh:
            if '1++,1--,2+-,2-+' in line:
                try:
                    fwd = float(line.strip().split()[-1])
                except:
                    pass
            if '1+-,1-+,2++,2--' in line:
                try:
                    rev = float(line.strip().split()[-1])
                except:
                    pass
    if rev > 0.6:
        votes.append(2)
    elif fwd > 0.6:
        votes.append(1)
    else:
        votes.append(0)
from collections import Counter
consensus = Counter(votes).most_common(1)[0][0] if votes else 0
label = {0: "unstranded", 1: "forward_stranded", 2: "reverse_stranded"}
print(f"Strandedness votes per sample: {votes}")
print(f"Consensus: {consensus} ({label.get(consensus, 'unknown')})")
with open("strandedness_consensus.txt", "w") as out:
    out.write(f"{consensus}\\n")
    out.write(f"Label: {label.get(consensus, 'unknown')}\\n")
    out.write(f"Votes: {votes}\\n")
PYEOF
    strand_int=\$(head -1 strandedness_consensus.txt | tr -d '[:space:]')
    """
}


// -----------------------------------------------------------------------------
// PROCESS 8: FEATURECOUNTS  Gene-level read quantification
// Strandedness is passed dynamically from PARSE_STRANDEDNESS
// No hardcoded -s flag  fully automated
// Uses HPC module (subread biocontainer lacks ps on this cluster)
// -----------------------------------------------------------------------------
process FEATURECOUNTS {

    tag "featureCounts"

    publishDir "results/featurecounts", mode: 'copy'

    input:
    path bam_files
    path bai_files
    path gtf
    val strandedness

    output:
    path "gene_counts.txt",         emit: counts
    path "gene_counts.txt.summary", emit: summary

    script:
    """
    source /etc/profile.d/modules.sh
    module load subread/2.0.6

    echo "Using strandedness: ${strandedness}"

    featureCounts \
        -T ${task.cpus} \
        -p -B -C \
        -s ${strandedness} \
        -t exon \
        -g gene_id \
        -a ${gtf} \
        -o gene_counts.txt \
        ${bam_files}
    """
}


// -----------------------------------------------------------------------------
// PROCESS 9: MULTIQC  Aggregate all QC reports into one HTML report
// Collects: FastQC + FASTP + STAR logs + RSeQC outputs + featureCounts summary
// -----------------------------------------------------------------------------
process MULTIQC {

    tag "multiqc"

    container "quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0"

    publishDir "results/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . -o . --force
    """
}


// =============================================================================
// WORKFLOW
// =============================================================================

workflow {

    //  Inputs 
    fastq_ch = Channel.fromFilePairs(params.reads)
    bed_ch   = file(params.bed)

    //  Step 1: Raw QC 
    FASTQC(fastq_ch)

    //  Step 2: Trimming 
    FASTP(fastq_ch)

    //  Step 3: Alignment 
    star_out = STAR_ALIGN(FASTP.out.trimmed)

    //  Step 4: Index BAM 
    indexed = INDEX_BAM(star_out)

    //  Step 5: RSeQC per-sample QC 
    INFER_STRANDEDNESS(indexed, bed_ch)
    READ_DISTRIBUTION(indexed, bed_ch)

    //  Step 6: Determine consensus strandedness 
    // Collect all infer_experiment outputs  parse  emit single integer
    strand_val = PARSE_STRANDEDNESS(
        INFER_STRANDEDNESS.out.txt.collect()
    ).strand_int.map { it.trim() }

    //  Step 7: Collect BAMs and BAIs for featureCounts 
    bam_files = indexed.map { sample_id, bam, bai -> bam }.collect()
    bai_files = indexed.map { sample_id, bam, bai -> bai }.collect()

    //  Step 8: Gene quantification with auto-detected strandedness 
    FEATURECOUNTS(
        bam_files,
        bai_files,
        file(params.gtf),
        strand_val
    )

    //  Step 9: MultiQC  aggregate all QC 
    qc_ch = Channel.empty()
        .mix(FASTQC.out.zip.collect())
        .mix(FASTP.out.json.collect())
        .mix(INFER_STRANDEDNESS.out.txt.collect())
        .mix(READ_DISTRIBUTION.out.txt.collect())
        .collect()

    MULTIQC(qc_ch)
}
