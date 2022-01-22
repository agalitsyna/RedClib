#!/usr/bin/env nextflow

// TODOs:
// 1. Check the specifications for all the modules. Restructure irrelevant as simple scripts.
// 2. Process low/high check.
// 3. Add checks for the read length.
// 4. Verify intermediate file names
// 5. Do not store the intermediary files
// 6. Add backend choice hdf5/parquet/tsv/csv
// 7. Remove non-deterministic steps to make resume work:
// https://www.nextflow.io/blog/2019/troubleshooting-nextflow-resume.html
// 8. Check cardinality of REDC:HISAT2_ALIGN_RNA1 and REDC:FASTQ_EXTEND_DNA

nextflow.enable.dsl = 2

// Define help message

def helpMessage() {
    log.info"""
    RedC-nextflow is a pipeline to map RNA-DNA interactions in paired-end sequencing data.

    Usage:
    The typical command for launching the pipeline:
      nextflow run main.nf -profile test,conda,debug -params-file params-redc.yml

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
include { RNADNATOOLS_GENOME_RECSITES as GENOME_RESTRICT } from './modules/rnadnatools/genome_recsites/main' addParams( options: [:])
include { GENOME_PREPARE_RNA_ANNOTATIONS } from './modules/local/genome_prepare_rna_annotations/main' addParams( options: [args: [
                                                            genes_gtf: GenesGtf,
                                                            rna_annotation_suffix: Genome.get('rna_annotation_suffix', '')
                                                            ]])

//include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { FASTQ2TSV as TABLE_FASTQ2TSV } from './modules/local/fastq2table/main' addParams( options: [:] )

include { DEDUP_FASTUNIQ as TABLE_DEDUP } from './subworkflows/local/dedup_fastuniq' addParams( options: [:] )
include { TRIMTABLE_CREATE as TABLE_TRIM } from './subworkflows/local/trimtable_create' addParams( options: [:] )

include { OLIGOS_MAP } from './subworkflows/local/oligos_map' addParams( options: [:] )

include { TSV_MERGE as TABLE_MERGE } from './modules/local/tsv_merge/main' addParams( options: [args: [:], suffix: '.table'] )

include { RNADNATOOLS_TABLE_CONVERT as TABLE_CONVERT } from './modules/rnadnatools/table_convert/main' addParams( options: [args: [
                                                             input_format: 'tsv',
                                                             output_format: 'parquet'
                                                             ]])
//include { RNADNATOOLS_TABLE_MERGE as TABLE_MERGE } from './modules/rnadnatools/table_merge/main' addParams( options: [args: [suffixes:['']]] )

include { RNADNATOOLS_TABLE_EVALUATE as TABLE_EVALUATE_FRAGMENTS } from './modules/rnadnatools/table_evaluate/main'  addParams( options: [args: [
                                                             format:'int',
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.fragments'] )
include { RNADNATOOLS_SEGMENT_EXTRACT_FASTQ as FRAGMENTS_TO_FASTQ } from './modules/rnadnatools/segment_extract_fastq/main'  addParams( options: [args: [
                                                             input_format: 'parquet'
                                                             ], suffix:'.fragments'] )

def dna_extension = params.fragments.dna.get("extension_suffix", "")
include { FASTQ_EXTEND as FASTQ_EXTEND_DNA } from './modules/local/extend_fastq/main' addParams( options: [args: [suffix: dna_extension], args2: [gzip: true]] )

include { HISAT2_ALIGN as HISAT2_ALIGN_RNA1 } from './modules/local/hisat2/main' addParams( options: [args:'-k 10 --no-softclip --dta-cufflinks --known-splicesite-infile', suffix:'.rna1'] )
include { HISAT2_ALIGN as HISAT2_ALIGN_RNA2 } from './modules/local/hisat2/main' addParams( options: [args:'-k 10 --no-softclip --dta-cufflinks --known-splicesite-infile', suffix:'.rna2'] )
include { HISAT2_ALIGN as HISAT2_ALIGN_DNA } from './modules/local/hisat2/main' addParams( options: [args:'-k 10 --no-softclip --no-spliced-alignment', suffix:'.dna'] )

include { BAM2BED as BAM2BED_RNA1 } from './modules/local/bam2bed/main' addParams( options: [filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2', filter2: 'samtools view -h -d NH:1 -', suffix: '.rna1'])
include { BAM2BED as BAM2BED_RNA2 } from './modules/local/bam2bed/main' addParams( options: [filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2', filter2: 'samtools view -h -d NH:1 -', suffix: '.rna2'])
include { BAM2BED as BAM2BED_DNA  } from './modules/local/bam2bed/main' addParams( options: [filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2', filter2: 'samtools view -h -d NH:1 -', suffix: '.dna'])

include { RNADNATOOLS_SEGMENT_GETCLOSEST as BED_ANNOTATE_RESTRICTION_RNA1 } from './modules/rnadnatools/segment_getclosest/main'  addParams( options: [args: [:], suffix:'.rna1'] )
include { RNADNATOOLS_SEGMENT_GETCLOSEST as BED_ANNOTATE_RESTRICTION_RNA2 } from './modules/rnadnatools/segment_getclosest/main'  addParams( options: [args: [:], suffix:'.rna2'] )

include { RNADNATOOLS_TABLE_EVALUATE as TABLE_EVALUATE_FILTERS } from './modules/rnadnatools/table_evaluate/main'  addParams( options: [args: [
                                                             format:'bool',
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.filters'] )

//include { COOLER_MAKE  } from './modules/local/cooler_make/main' addParams( options: [assembly: Assembly, resolution: params.get("cooler_resolution", 1000000)])


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

//    // deduplication of input sequences
//    TABLE_DEDUP(Fastq)

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
    TableChunks = TABLE_FASTQ2TSV(FastqChunks).table

    /* Trim reads by quality */
    TrimTable = TABLE_TRIM(FastqChunks, TableChunks).trimtable

//    TrimmedChunks = FASTQ_TRIM(FastqChunks).fastq
//    TrimTable = TABLE_TRIM( join_2_channels(TableChunks, TrimmedChunks, 'id') ).table

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

    ResultingParquetTable = TABLE_CONVERT( ResultingTable ).table

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
                                                                  ).output

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
    if (dna_extension) { // todo: check cardinality of the channel
        FastqDNA = FASTQ_EXTEND_DNA (
            FragmentsFastq.filter{ it[0].fragment == 'dna' }
        ).fastq
    } else {
        FastqDNA = FragmentsFastq.filter{ it[0].fragment == 'dna' }
    }

    HISAT2_ALIGN_RNA1 (  // todo: check cardinality of the channel
            FastqRNA1,
            Hisat2Index,
            SpliceSites
    )

    HISAT2_ALIGN_RNA2 ( // todo: check cardinality of the channel
            FastqRNA2,
            Hisat2Index,
            SpliceSites
    )

    HISAT2_ALIGN_DNA ( // todo: check cardinality of the channel
            FastqDNA,
            Hisat2Index,
            SpliceSites
    )

    BAM2BED_RNA1 ( HISAT2_ALIGN_RNA1.out.bam)
    BAM2BED_RNA2 ( HISAT2_ALIGN_RNA2.out.bam)
    BAM2BED_DNA  ( HISAT2_ALIGN_DNA.out.bam)

    RestrictedPlus = Restricted.map{ update_meta(it, [renz_strand: "+"] ) }
    RestrictedMinus = Restricted.map{ update_meta(it, [renz_strand: "-"] ) }
    RestrictedExtended = RestrictedPlus.mix( RestrictedMinus )

    ParquetTableRestrictionRNA1 = BED_ANNOTATE_RESTRICTION_RNA1 ( join_2_channels(BAM2BED_RNA1.out.bed, ResultingParquetTable, 'id'), RestrictedExtended ).output
    ParquetTableRestrictionRNA2 = BED_ANNOTATE_RESTRICTION_RNA2 ( join_2_channels(BAM2BED_RNA2.out.bed, ResultingParquetTable, 'id'), RestrictedExtended ).output


    def FiltersColumns = [:] // todo: add necessary check here
    for (filter_name in params.filters.keySet()){
        filter_params = params.filters[filter_name]
        if (filter_params instanceof List) {
            filter_params = "(" + filter_params.join(") | (") + ")"
        }
        FiltersColumns[filter_name] = filter_params
    }
    Filters = Channel.from(FiltersColumns)


    ResultingParquetTable.view()
    ResultingParquetTableWithFragments.view()
    ParquetTableRestrictionRNA1.view()
    ParquetTableRestrictionRNA2.view()


//    TABLE_EVALUATE_FILTERS( ResultingParquetTableWithFragments
//                                .combine( Filters )
//                                .multiMap{meta, id, filt ->
//                                          pq: [meta, id]
//                                          filters: filt }
//                              )

    // TODO: add custom tables as output
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
