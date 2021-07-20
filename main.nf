#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Define help message

def helpMessage() {
    log.info"""
    RedC-nextflow is a pipeline to map RNA-DNA interactions in paired-end sequencing data.

    Usage:
    The typical command for launching the pipeline:
      nextflow run main.nf -profile test,conda,debug

    All the parameters should be listed in the profile file.
    """.stripIndent()
}

// Show help message
if (params.get('help', 'false').toBoolean()) {
    helpMessage()
    exit 0
}

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
genome = params.get('genome', [:])
if (params.genome.assembly) { assembly = params.genome.assembly } else { exit 1, 'Genome assembly is not specified!' }
if (params.genome.genes_gtf) { genes_gtf = params.genome.genes_gtf } else { exit 1, 'Genome RNA annotation GTF is not specified!' }

// Check optional parameters
def chunksize = params.get('chunksize', 0)
def dedup_crop_length = params.get('dedup_crop_length', 50)

include { INPUT_CHECK_DOWNLOAD } from './subworkflows/local/input_check' addParams( options: [:] )
include { INPUT_SPLIT          } from './subworkflows/local/input_split' addParams( options: [chunksize: chunksize] )

include { GENOME_PREPARE } from './modules/local/genome_prepare/main'  addParams( options: [args: [genome: [chromsizes: genome.get("chromsizes", ""), index_prefix: genome.get("index_prefix", ""), fasta: genome.get("fasta", "")], auto_download_genome: genome.get("auto_download_genome", true)] ] )
include { GENOME_RESTRICT } from './modules/local/genome_restrict/main' addParams( options: [:])
include { GENOME_PREPARE_RNA_ANNOTATIONS } from './modules/local/genome_prepare_rna_annotations/main' addParams( options: [args: [genes_gtf: genes_gtf, rna_annotation_suffix: genome.get('rna_annotation_suffix', '')]])

include { DEDUP } from './modules/local/dedup/main' addParams( options: [ args: [dedup_crop_length: dedup_crop_length]] )
include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { TABLE_CONVERT } from './modules/local/table_convert/main' addParams( options: [:] )
include { TRIM_TRIMMOMATIC} from './modules/local/trimmomatic/main' addParams( options: [args: [params_trimmomatic: params.get('params_trimmomatic', '')]])
include { HISAT2_ALIGN as HISAT2_ALIGN_RNA } from './modules/local/hisat2/main' addParams( options: [args:'--known-splicesite-infile', suffix:'.rna'] )

workflow REDC {

    /* Prepare input */
    // Check input FASTQ files
    INPUT_CHECK_DOWNLOAD (
            ch_input
        )
        .map {
            meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }
        .set { ch_fastq }

    // fastqc quality of the sequencing
    FASTQC (
        ch_fastq
    )

    // deduplication of input sequences
    DEDUP (
        ch_fastq
    )

    // split input into chunks
    if (chunksize) {
        INPUT_SPLIT(
                ch_fastq
        ).set { ch_fastq_chunks }
    } else {
        ch_fastq_chunks = ch_fastq
    }

    /* Prepare genome and annotations */
    // Prepare genome (index, chromosome sizes and fasta)
    GENOME_PREPARE (
            assembly
    )
    hisat2_index = GENOME_PREPARE.out.genome_index
    genome_fasta = GENOME_PREPARE.out.genome_fasta

    // Restrict the genome
    if (params.get('check_restriction', false)) {
        list_renzymes = Channel.fromList(params.protocol.renzymes)
        GENOME_RESTRICT (
            assembly,
            list_renzymes,
            genome_fasta
        )
    }

    // Get the splice sites
    GENOME_PREPARE_RNA_ANNOTATIONS (
        assembly
    )
    hisat2_splicesites = GENOME_PREPARE_RNA_ANNOTATIONS.out.genome_splicesites


    /* Start of the reads processing */

    TABLE_CONVERT (
        ch_fastq_chunks
    )

    TRIM_TRIMMOMATIC (
        ch_fastq_chunks
    )

    //OLIGOS_ALIGNMENT (
    //    ch_fastq_chunks
    //)

    //OLIGOS_CHECK (
    //    ch_fastq_chunks
    //)

    //COMPLEMENTARY_CHECK (
    //    ch_fastq_chunks
    //)

    //SUBSTRINGS_GET (
    //    ch_fastq_chunks
    //)

    //EXTEND (
    //    ch_parsed.dna
    //)

    //HISAT2_ALIGN_RNA (
    //        ch_fastq,
    //        hisat2_index,
    //        hisat2_splicesites
    //)

    //HISAT2_ALIGN_RNA.out.bam.view()

}

workflow {

    REDC ( )

}
workflow.onComplete {
    log.info "Done!"
}