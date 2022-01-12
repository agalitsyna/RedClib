#!/usr/bin/env nextflow

// TODOs:
// 1. Check the specifications for all the modules. Restructure irrelevant as simple scripts.
// 2. Process low/high check.
// 3. Add checks for the read length.
// 4. verify intermediate file names
// 5. do not store the intermediary files
// 6. Add hdf5 backend instead of pyarrow/parquet

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
include { TABLE_CONVERT as FASTQ_TABLE_CONVERT } from './modules/local/table_convert/main' addParams( options: [:] )

include { TRIMMOMATIC as FASTQ_CROP } from './modules/local/trimmomatic/main' addParams( options: [args: 'CROP:'+dedup_crop, suffix:'.crop', args2: [gzip: false]])
include { FASTUNIQ as TABLE_FASTUNIQ} from './modules/local/fastuniq/main' addParams( options: [:] )

include { TRIMMOMATIC as FASTQ_TRIM } from './modules/local/trimmomatic/main' addParams( options: [args: trimmomatic_params, suffix:'.trim', args2: [gzip: false]])
include { TABLE_TRIM } from './modules/local/table_trimmomatic' addParams( options: [:] )

include { OLIGOS_MAP } from './subworkflows/local/oligos_map' addParams( options: [:] )

include { TSV_MERGE as TABLE_MERGE } from './modules/local/tsv_merge/main' addParams( options: [args: [:], suffix: '.table'] )

include { PARQUET_CONVERT as TABLE_CONVERT } from './modules/local/parquet_convert/main' addParams( options: [args: [:]] )
//include { PARQUET_MERGE as TABLE_MERGE } from './modules/local/parquet_merge/main' addParams( options: [args: [suffixes:['']]] )

include { PARQUET_EVALUATE as TABLE_EVALUATE_FRAGMENTS } from './modules/local/parquet_evaluate'  addParams( options: [args: [format:'int'], suffix:'.fragments'] )
include { PARQUET2FASTQ as FRAGMENTS_TO_FASTQ } from './modules/local/parquet2fastq'  addParams( options: [args: [:], suffix:'.fragments'] )
def dna_extension = params.fragments.dna.get("extension_suffix", "")
include { FASTQ_EXTEND as FASTQ_EXTEND_DNA } from './modules/local/extend_fastq/main' addParams( options: [args: [suffix: dna_extension], args2: [gzip: true]] )

include { HISAT2_ALIGN as HISAT2_ALIGN_RNA1 } from './modules/local/hisat2/main' addParams( options: [args:'-k 10 --no-softclip --dta-cufflinks --known-splicesite-infile', suffix:'.rna1'] )
include { HISAT2_ALIGN as HISAT2_ALIGN_RNA2 } from './modules/local/hisat2/main' addParams( options: [args:'-k 10 --no-softclip --dta-cufflinks --known-splicesite-infile', suffix:'.rna2'] )
include { HISAT2_ALIGN as HISAT2_ALIGN_DNA } from './modules/local/hisat2/main' addParams( options: [args:'-k 10 --no-softclip --no-spliced-alignment', suffix:'.dna'] )

include { BAM2BED as BAM2BED_RNA1 } from './modules/local/bam2bed/main' addParams( options: [filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2', filter2: 'samtools view -h -d NH:1 -', suffix: '.rna1'])
include { BAM2BED as BAM2BED_RNA2 } from './modules/local/bam2bed/main' addParams( options: [filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2', filter2: 'samtools view -h -d NH:1 -', suffix: '.rna2'])
include { BAM2BED as BAM2BED_DNA  } from './modules/local/bam2bed/main' addParams( options: [filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2', filter2: 'samtools view -h -d NH:1 -', suffix: '.dna'])


include { COOLER_MAKE  } from './modules/local/cooler_make/main' addParams( options: [assembly: Assembly, resolution: params.get("cooler_resolution", 1000000)])


// Define workflow
workflow REDC {

    /* Prepare input */
    // Check input FASTQ files
    Fastq = INPUT_CHECK_DOWNLOAD(InputSamplesheet)
        .map {meta, fastq ->
                meta.id = meta.id.split('_')[0..-2].join('_')
                [ meta, fastq ] }

//    // fastqc quality of the sequencing
//    FASTQC(Fastq)

    // deduplication of input sequences
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
    TableChunks = FASTQ_TABLE_CONVERT(FastqChunks).table

    /* Trim reads by quality */
    TrimmedChunks = FASTQ_TRIM(FastqChunks).fastq
    TrimTable = TABLE_TRIM( join_2_channels(TableChunks, TrimmedChunks, 'id') ).table

    /* Map and check oligos */
    Hits = OLIGOS_MAP(FastqChunks, TableChunks)

    /* One of the important tables with merged statistics for each read: */
    TablesGrouped = TableChunks.map{ it -> [ [it[0]['id']], it ] }
                             .combine( TrimTable.map{ it -> [ [it[0]['id']], it ] }, by: 0 )
                             .combine( Hits.HitsOligos.map{ it -> [ [it[0]['id']], it ] }, by: 0 )
                             .combine( Hits.HitsSmallOligos.map{ it -> [ [it[0]['id']], it ] }, by: 0 )
                             .combine( Hits.HitsComplementary.map{ it -> [ [it[0]['id']], it ] }, by: 0 )
                             .multiMap{id, ch1, ch2, ch3, ch4, ch5 ->
                                dataset: [ch1[0], [ch1[1], ch2[1], ch3[1], ch4[1], ch5[1]]]
                                suffixes: ['', '', '', '', '']
                           }
    ResultingTable = TABLE_MERGE( TablesGrouped ).output

    ResultingParquetTable = TABLE_CONVERT( ResultingTable ).parquet

    /* Read metadata about the fragments */
    def FragmentsColumns = [:] // todo: add necessary check here
    for (fragment_name in params.fragments.keySet()){ // select all fragments
        fragment_params = params.fragments[fragment_name] // take the parameters for each fragment
        for( col in fragment_params.new_columns.keySet() ) { // add new columns to the list
            FragmentsColumns[col] = fragment_params.new_columns[col]
        }
    }
    FragmentsFilters = Channel.from(FragmentsColumns)

    ResultingParquetTableWithFragments = TABLE_EVALUATE_FRAGMENTS (ResultingParquetTable
                                                                    .combine( FragmentsFilters )
                                                                    .multiMap{meta, id, filt ->
                                                                              pq: [meta, id]
                                                                              filters: filt }
                                                                  ).parquet

    FragmentsSelection = ResultingParquetTableWithFragments
                        .combine( Channel.fromList(params.fragments.collect { fragment, fragment_params -> [fragment: fragment, params:fragment_params] }) )
                        .map { meta_pq, pq1, pq2, fragment ->
                                        bin_reads: update_meta( [meta_pq, [pq1, pq2]],
                                                                [fragment: fragment.fragment,
                                                                 side: fragment.params.side,
                                                                 single_end: true,
                                                                 selection_criteria: fragment.params.selection_criteria] )
                                   }

    FragmentsFastq = FRAGMENTS_TO_FASTQ( FragmentsSelection ).output

    /* Define channels for Fastq files with fragments: */
    FastqRNA1 = FragmentsFastq.filter{ it[0].fragment == 'rna1' }
    FastqRNA2 = FragmentsFastq.filter{ it[0].fragment == 'rna2' }
    if (dna_extension) {
        FastqDNA = FASTQ_EXTEND_DNA (
            FragmentsFastq.filter{ it[0].fragment == 'dna' }
        ).fastq
    } else {
        FastqDNA = FragmentsFastq.filter{ it[0].fragment == 'dna' }
    }

    HISAT2_ALIGN_RNA1 (
            FastqRNA1,
            Hisat2Index,
            SpliceSites
    )

    HISAT2_ALIGN_RNA2 (
            FastqRNA2,
            Hisat2Index,
            SpliceSites
    )

    HISAT2_ALIGN_DNA (
            FastqDNA,
            Hisat2Index,
            SpliceSites
    )

    HISAT2_ALIGN_RNA1.out.bam.view()
    BAM2BED_RNA1 ( HISAT2_ALIGN_RNA1.out.bam)
    BAM2BED_RNA2 ( HISAT2_ALIGN_RNA2.out.bam)
    BAM2BED_DNA  ( HISAT2_ALIGN_DNA.out.bam)

    // TODO: make custom table as output
//    ChromSizes = GENOME_PREPARE.out.chromsizes
//    COOLER_MAKE (
//        FinalTable,
//        ChromSizes
//    )

}

workflow {

    REDC ( )

}
workflow.onComplete {
    log.info "Done!"
}

def join_2_channels(channel_l, channel_r, k){
    // Take two channels and a meta key and return two channels combined by that key:
    channel_l
        .map{ it -> [it[0][k], it] }
        .combine(channel_r.map{ it -> [it[0][k], it] }, by: 0)
        .multiMap { key, left, right ->
                    left: left
                    right: right
         }
}

def update_meta( it, hashMap ) {
    def meta = [:]
    keys = it[0].keySet() as String[]
    for( key in keys ) {
        meta[key] = it[0][key]
    }

    if (!meta.single_end) {
        file2 = it[1][1]
    }

    keys_new = hashMap.keySet() as String[]
    for( key in keys_new ) {
        meta[key] = hashMap[key]
    }

   array = [ meta, it[1] ]
    return array
}