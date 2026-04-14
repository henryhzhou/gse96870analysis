#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ============================================================
// main.nf — GSE96870 RNA-seq Pipeline
// ============================================================

// 打印运行信息
log.info """
    GSE96870 RNA-seq Pipeline
    =========================
    fastq_dir  : ${params.fastq_dir}
    results_dir: ${params.results_dir}
    star_index : ${params.star_index}
    """

// ============================================================
// PROCESS 1: FastQC
// ============================================================
process FASTQC {
    tag "${sample_id}"
    publishDir "${params.results_dir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{zip,html}"

    script:
    """
    fastqc -t ${task.cpus} -o . ${reads}
    """
}

// ============================================================
// PROCESS 2: Trim Galore
// ============================================================
process TRIM {
    tag "${sample_id}"
    publishDir "${params.results_dir}/trim", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_val_{1,2}.fq.gz"), emit: trimmed_reads
    path "*_trimming_report.txt"

    script:
    """
    trim_galore --paired \
        --cores ${task.cpus} \
        --quality 20 \
        --length 20 \
        ${reads[0]} ${reads[1]}
    """
}

// ============================================================
// PROCESS 3: STAR Index（只需要建一次）
// ============================================================
process STAR_INDEX {
    publishDir "${params.refs_dir}", mode: 'copy'

    output:
    path "star_index", emit: star_index

    script:
    """
    mkdir -p star_index
    STAR --runMode genomeGenerate \
        --runThreadN ${task.cpus} \
        --genomeDir star_index \
        --genomeFastaFiles ${params.genome_fa} \
        --sjdbGTFfile ${params.genome_gtf} \
        --genomeSAindexNbases 14
    """
}

// ============================================================
// PROCESS 4: STAR Alignment
// ============================================================
process STAR_ALIGN {
    tag "${sample_id}"
    publishDir "${params.results_dir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path star_index

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.sortedByCoord.out.bam"), emit: bam
    path "${sample_id}.Log.final.out"

    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${star_index} \
        --readFilesIn ${reads[0]} ${reads[1]} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM \
        --outFileNamePrefix ${sample_id}. \
        --quantMode GeneCounts
    """
}

// ============================================================
// PROCESS 5: featureCounts
// ============================================================
process FEATURE_COUNTS {
    publishDir "${params.results_dir}/counts", mode: 'copy'

    input:
    path bam_files

    output:
    path "counts.txt"
    path "counts.txt.summary"

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -p \
        -a ${params.genome_gtf} \
        -o counts.txt \
        ${bam_files}
    """
}

// ============================================================
// PROCESS 6: MultiQC
// ============================================================
process MULTIQC {
    publishDir "${params.results_dir}/multiqc", mode: 'copy'

    input:
    path qc_files

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc . -o .
    """
}

// ============================================================
// WORKFLOW — 把所有process串联起来
// ============================================================
workflow {

    // 读取样本：自动配对 _1 和 _2 文件
    reads_ch = Channel
        .fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz", checkIfExists: true)

    // 步骤1: FastQC（并行 × 45）
    FASTQC(reads_ch)

    // 步骤2: Trim（并行 × 45）
    TRIM(reads_ch)

    // 步骤3: 建STAR index（只跑一次）
    STAR_INDEX()

    // 步骤4: Alignment（并行 × 45）
    STAR_ALIGN(TRIM.out.trimmed_reads, STAR_INDEX.out.star_index)

    // 步骤5: featureCounts（收集所有bam，跑一次）
    all_bams = STAR_ALIGN.out.bam
        .map { sample_id, bam -> bam }
        .collect()
    FEATURE_COUNTS(all_bams)

    // 步骤6: MultiQC（收集所有QC文件，跑一次）
    all_qc = FASTQC.out
        .mix(TRIM.out[1])
        .collect()
    MULTIQC(all_qc)
}
