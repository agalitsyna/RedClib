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

// Check required parameters
if (params.input) { InputSamplesheet = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
Genome = params.get('genome', [:])
if (params.genome.assembly) { Assembly = params.genome.assembly } else { exit 1, 'Genome assembly is not specified!' }
if (params.genome.genes_gtf) { GenesGtf = params.genome.genes_gtf } else { exit 1, 'Genome RNA annotation GTF is not specified!' }

// Check optional parameters
def chunksize = params.get('chunksize', 0)
def dedup_crop = params.get('dedup_crop', 50)
def trimmomatic_params = params.get('trimmomatic_params', '')
def check_restriction = params.get('check_restriction', false)
def RenzymesPreloaded = Genome.get('restricted', {})

// Include modules and subworkflows
include { INPUT_CHECK_DOWNLOAD } from './subworkflows/local/input_check' addParams( options: [:] )
include { INPUT_SPLIT          } from './subworkflows/local/input_split' addParams( options: [chunksize: chunksize] )

include { GENOME_PREPARE } from './modules/local/genome_prepare/main'  addParams( options: [args: [
                                                            genome: [chromsizes: Genome.get("chromsizes", ""),
                                                            index_prefix: Genome.get("index_prefix", ""),
                                                            fasta: Genome.get("fasta", "")],
                                                            auto_download_genome: Genome.get("auto_download_genome", true)
                                                            ] ] )
include { GENOME_RESTRICT } from './modules/local/genome_restrict/main' addParams( options: [:])
include { GENOME_PREPARE_RNA_ANNOTATIONS } from './modules/local/genome_prepare_rna_annotations/main' addParams( options: [args: [
                                                            genes_gtf: GenesGtf,
                                                            rna_annotation_suffix: Genome.get('rna_annotation_suffix', '')
                                                            ]])

include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { TABLE_CONVERT } from './modules/local/table_convert/main' addParams( options: [:] )

include { TRIMMOMATIC as FASTQ_CROP } from './modules/local/trimmomatic/main' addParams( options: [args: 'CROP:'+dedup_crop, suffix:'.crop', args2: [gzip: false]])
include { FASTUNIQ as TABLE_FASTUNIQ} from './modules/local/fastuniq/main' addParams( options: [:] )

include { TRIMMOMATIC as FASTQ_TRIM } from './modules/local/trimmomatic/main' addParams( options: [args: trimmomatic_params, suffix:'.trim', args2: [gzip: false]])
include { HISAT2_ALIGN as HISAT2_ALIGN_RNA } from './modules/local/hisat2/main' addParams( options: [args:'--known-splicesite-infile', suffix:'.rna'] )

include { BIN_ENCODE as BIN_OLIGOS } from './modules/local/bin_encode/main' addParams( options: [args: [mode: 'fasta']] ) // Bin input oligos from fasta file
include { BIN_ENCODE as BIN_FASTQ }  from './modules/local/bin_encode/main' addParams( options: [args: [mode: 'fastq']] ) // Bin input fastq

include { OLIGOS_ALIGN } from './modules/local/align_encoded/main' addParams( options: [args: [:]] ) // Align oligos

// Define workflow
workflow REDC {

    /* Prepare input */
    // Check input FASTQ files
    Fastq = INPUT_CHECK_DOWNLOAD(InputSamplesheet)
        .map {meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }

    // fastqc quality of the sequencing
    FASTQC(Fastq)

    // deduplication of input sequences
    //FastqCropped = Fastq | FASTQ_CROP
    //Nodup = FastqCropped.fastq | TABLE_FASTUNIQ
    FASTQ_CROP(Fastq)
    Nodup = FASTQ_CROP.out.fastq | TABLE_FASTUNIQ

    // split input into chunks if chunksize was passed
    FastqChunks = chunksize ? INPUT_SPLIT(Fastq) : Fastq

    /* Prepare genome and annotations */
    // Prepare genome (index, chromosome sizes and fasta)
    GENOME_PREPARE(Assembly)
    Hisat2Index = GENOME_PREPARE.out.genome_index
    GenomeFasta = GENOME_PREPARE.out.genome_fasta

    // Restrict the genome or load pre-computed restriction sites
    if (check_restriction) {
        Renzymes = Channel
            .fromList(params.protocol.renzymes)
            .branch { it ->
                    forRestriction : !RenzymesPreloaded.containsKey(it) // Restrict only missing renzymes
                        return [renzyme: it, assembly: Assembly]
                    loaded : RenzymesPreloaded.containsKey(it) // Load files if restriction is pre-computed
                        return [[renzyme: it, assembly: Assembly], file(RenzymesPreloaded[it])]
                    }
        GENOME_RESTRICT(
            Renzymes.forRestriction,
            GenomeFasta
        )
        Restricted = GENOME_RESTRICT.out.genome_restricted.mix( Renzymes.loaded ) // Concatenate two outputs
    }

    // Get the splice sites by hisat2 script
    GENOME_PREPARE_RNA_ANNOTATIONS(Assembly)
    SpliceSites = GENOME_PREPARE_RNA_ANNOTATIONS.out.genome_splicesites


    /* Start of the reads processing */
    TableChunks = TABLE_CONVERT(FastqChunks)
    TrimmedChunks = FASTQ_TRIM(FastqChunks)


    /* Check for the presence of oligos with a custom Rabin-Karp-based algorithm */
    Oligos = Channel
               .fromList(params.oligos.keySet())
               .map { it -> [it, params.oligos[it]] } // construct meta for each oligo
               .map { oligo_name, meta ->
                       meta.single_end = true
                       meta.id = oligo_name
                   [meta, file(meta.file)]
               }

    IndexOligos = BIN_OLIGOS(Oligos)
    IndexFastqs = BIN_FASTQ(FastqChunks)

    MappedOligos = IndexOligos.bin
        .flatMap {
            meta, bin -> // Groovy list comprehesion to get one channel emission per requested side of oligo:
            meta.sides.collect{ it -> meta.side=it; [meta, bin] }
        }
        .combine(IndexFastqs.bin) \
        | OLIGOS_ALIGN
    MappedOligos.aligned.view()

    /* Check presence of GA dunucleotide in the bridge: */
    // OLIGOS_CHECK_GA(TableChunks, MappedOligos)

    /* Check complementary RNA ends: */
    // OLIGOS_CHECK_COMPLEMENTARY(TableChunks, IndexFastqs, MappedOligos) // should be a subworkflow


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