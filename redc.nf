#!/usr/bin/env nextflow

// TODOs:
// 1. Check the specifications for all the modules. Restructure irrelevant as simple scripts.
// 2. Process low/high check.
// 3. Add checks for the read length.
// 4. Verify intermediate file names
// 5. Do not store the intermediary files
// 6. Add backend choice hdf5/parquet/tsv/csv

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
if (params.input) { InputSamplesheet = file(params.input) }
    else { exit 1, 'Input samplesheet not specified!' }

Genome = params.get('genome', [:])
if (params.genome.assembly) { Assembly = params.genome.assembly }
    else { exit 1, 'Genome assembly is not specified!' }
if (params.genome.genes_gtf) { GenesGtf = params.genome.genes_gtf }
    else { exit 1, 'Genome RNA annotation GTF is not specified!' }

// Check optional parameters
def chunksize = params.get('chunksize', 0)
def check_restriction = params.get('check_restriction', false)
def RenzymesPreloaded = Genome.get('restricted', {})

// Include modules and subworkflows
include { INPUT_CHECK_DOWNLOAD } from './subworkflows/local/input_check' addParams( options: [:] )
include { INPUT_SPLIT          } from './subworkflows/local/input_split' addParams( options: [chunksize: chunksize] )

include { GENOME_PREPARE } from './modules/local/genome_prepare/main'  addParams( options: [
                                                             args: [
                                                            genome: [chromsizes: Genome.get("chromsizes", ""),
                                                            index_prefix: Genome.get("index_prefix", ""),
                                                            fasta: Genome.get("fasta", "")],
                                                            auto_download_genome: Genome.get("auto_download_genome", true)
                                                            ] ] )
include { RNADNATOOLS_GENOME_RECSITES as GENOME_RESTRICT } from './modules/rnadnatools/genome_recsites/main' addParams( options: [:])
include { GENOME_PREPARE_RNA_ANNOTATIONS } from './modules/local/genome_prepare_rna_annotations/main' addParams( options: [
                                                             args: [
                                                            genes_gtf: GenesGtf,
                                                            rna_annotation_suffix: Genome.get('rna_annotation_suffix', '')
                                                            ]])

//include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { FASTQ2TSV as TABLE_FASTQ2TSV } from './modules/local/fastq2table/main' addParams( options: [:] )

include { DEDUP_FASTUNIQ as TABLE_DEDUP } from './subworkflows/local/dedup_fastuniq' addParams( options: [:] )
include { TRIMTABLE_CREATE as TABLE_TRIM } from './subworkflows/local/trimtable_create' addParams( options: [:] )

include { OLIGOS_MAP } from './subworkflows/local/oligos_map' addParams( options: [:] )

include { TSV_MERGE as TABLE_MERGE } from './modules/local/tsv_merge/main' addParams( options: [
                                                             args: [:],
                                                             suffix: '.table'] )

include { RNADNATOOLS_TABLE_CONVERT as TABLE_CONVERT } from './modules/rnadnatools/table_convert/main' addParams( options: [
                                                             args: [
                                                             input_format: 'tsv',
                                                             output_format: 'parquet'
                                                             ]])
//include { RNADNATOOLS_TABLE_MERGE as TABLE_MERGE } from './modules/rnadnatools/table_merge/main' addParams( options: [args: [suffixes:['']]] )

include { RNADNATOOLS_TABLE_EVALUATE as TABLE_EVALUATE_FRAGMENTS } from './modules/rnadnatools/table_evaluate/main'  addParams( options: [
                                                             args: [
                                                             format:'int',
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.fragments'] )
include { RNADNATOOLS_SEGMENT_EXTRACT_FASTQ as FRAGMENTS_TO_FASTQ } from './modules/rnadnatools/segment_extract_fastq/main'  addParams( options: [
                                                             args: [
                                                             input_format: 'parquet'
                                                             ], suffix:'.fragments'] )

def dna_extension = params.fragments.dna.get("extension_suffix", "")
include { FASTQ_EXTEND as FASTQ_EXTEND_DNA } from './modules/local/extend_fastq/main' addParams( options: [
                                                             args: [suffix: dna_extension],
                                                             args2: [gzip: true]] )

include { HISAT2_ALIGN as HISAT2_ALIGN_RNA1 } from './modules/local/hisat2/main' addParams( options: [
                                                             args:'-k 10 --no-softclip --dta-cufflinks --known-splicesite-infile',
                                                             suffix:'.rna1'] )
include { HISAT2_ALIGN as HISAT2_ALIGN_RNA2 } from './modules/local/hisat2/main' addParams( options: [
                                                             args:'-k 10 --no-softclip --dta-cufflinks --known-splicesite-infile',
                                                             suffix:'.rna2'] )
include { HISAT2_ALIGN as HISAT2_ALIGN_DNA } from './modules/local/hisat2/main' addParams( options: [
                                                             args:'-k 10 --no-softclip --no-spliced-alignment',
                                                             suffix:'.dna'] )

include { BAM2BED as BAM2BED_RNA1 } from './modules/local/bam2bed/main' addParams( options: [
                                                             filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2',
                                                             filter2: 'samtools view -h -d NH:1 -',
                                                             suffix: '.rna1'])
include { BAM2BED as BAM2BED_RNA2 } from './modules/local/bam2bed/main' addParams( options: [
                                                             filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2',
                                                             filter2: 'samtools view -h -d NH:1 -',
                                                             suffix: '.rna2'])
include { BAM2BED as BAM2BED_DNA  } from './modules/local/bam2bed/main' addParams( options: [
                                                             filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2',
                                                             filter2: 'samtools view -h -d NH:1 -',
                                                             suffix: '.dna'])

include { RNADNATOOLS_SEGMENT_GETCLOSEST as BED_ANNOTATE_RESTRICTION_RNA1 } from './modules/rnadnatools/segment_getclosest/main'  addParams( options: [
                                                             args: [:],
                                                             suffix:'.rna1'] )
include { RNADNATOOLS_SEGMENT_GETCLOSEST as BED_ANNOTATE_RESTRICTION_RNA2 } from './modules/rnadnatools/segment_getclosest/main'  addParams( options: [
                                                             args: [:],
                                                             suffix:'.rna2'] )

include { RNADNATOOLS_TABLE_EVALUATE as TABLE_EVALUATE_FILTERS } from './modules/rnadnatools/table_evaluate/main'  addParams( options: [
                                                             args: [
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
    if (chunksize) {
        FastqChunks = INPUT_SPLIT(Fastq)
    }
    else {
        FastqChunks = Fastq
    }

//test-download1_T1,sra:SRR10010324:start=0&end=10000,,125,false
//test-download2_T2,sra:SRR10010324:start=0&end=1000,,125,false

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
        RestrictedPlus = Restricted.map{ update_meta(it, [renz_strand: "+"] ) }
        RestrictedMinus = Restricted.map{ update_meta(it, [renz_strand: "-"] ) }
        RestrictedExtended = RestrictedPlus.mix( RestrictedMinus ).map{ it -> [it] }
    }

    // Get the splice sites by hisat2 script
    GENOME_PREPARE_RNA_ANNOTATIONS(Assembly)
    SpliceSites = GENOME_PREPARE_RNA_ANNOTATIONS.out.genome_splicesites

    /* Start of the reads processing */
    TableChunks = TABLE_FASTQ2TSV(FastqChunks).table

    /* Trim reads by quality */
    TrimTable = TABLE_TRIM(FastqChunks, TableChunks).trimtable

    /* Map and check oligos */
    Hits = OLIGOS_MAP(FastqChunks, TableChunks)

    /* Table with read, trimming and oligos mapping data: */
    // Create synchronized multichannel with tables (each table emitted separately):
    TablesGrouped_mult = TableChunks.combine( TrimTable, by:0 ).map{it -> [it[0], it[1..-1]]}
        .cross(Hits.HitsOligos) { it -> it[0].id }.map { it -> [it[0][0], [*it[0][1], it[1][1]] ]}
        .cross(Hits.HitsSmallOligos) { it -> it[0].id }.map { it -> [it[0][0], [*it[0][1], it[1][1]] ]}
        .cross(Hits.HitsComplementary) { it -> it[0].id }.map { it -> [it[0][0], [*it[0][1], it[1][1]] ]}
         .multiMap{it ->
                dataset: it // Channel structure: [meta, [paths]]
                suffixes: ['', '', '', '', ''] // No suffixes
            }
    MergedTable = TABLE_MERGE( TablesGrouped_mult ).table
    Table = TABLE_CONVERT( MergedTable ).table

    /* Table with fragments data.
      Information about the RNA/DNA fragments will be stored in columns of table
      Columns are evaluated based on input expressions */
    def FragmentColumns = [:] // todo: add necessary check here
    for (fragment_name in params.fragments.keySet()){ // select all fragments
        fragment_params = params.fragments[fragment_name] // take the parameters for each fragment
        for( col in fragment_params.new_columns.keySet() ) { // add new columns to the list
            FragmentColumns[col] = fragment_params.new_columns[col]
        }
    }
    FragmentFilters = Channel.from(FragmentColumns)
    // Create multichannel with fragment filters:
    TableFragments_mult = Table.combine( FragmentFilters )
                      .multiMap { meta, id, filt ->
                          pq: [meta, id]
                          filters: filt
                          }
    TableFragments = TABLE_EVALUATE_FRAGMENTS( TableFragments_mult ).table

    /* Extract fragments FASTQ */
    def FragmentData = params.fragments.collect {
            fragment, fragment_params -> [fragment: fragment, params:fragment_params]
            }
    FragmentFilterData = Channel.fromList(FragmentData)
    FragmentsSelected = Table
                        .combine(TableFragments, by:0)
                        .combine( FragmentFilterData )
                        .map { meta_pq, pq1, pq2, fragment ->
                                        bin_reads: update_meta( [meta_pq, [pq1, pq2]],
                                                                [fragment: fragment.fragment,
                                                                 side: fragment.params.side,
                                                                 single_end: true,
                                                                 selection_criteria: fragment.params.selection_criteria] )
                                   }

    FragmentsFastq = FRAGMENTS_TO_FASTQ( FragmentsSelected ).fastq

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

    /* Map the reads */
    BamRNA1 = HISAT2_ALIGN_RNA1 (
            FastqRNA1,
            Hisat2Index,
            SpliceSites
    ).bam

    BamRNA2 = HISAT2_ALIGN_RNA2 (
            FastqRNA2,
            Hisat2Index,
            SpliceSites
    ).bam

    BamDNA = HISAT2_ALIGN_DNA (
            FastqDNA,
            Hisat2Index,
            SpliceSites
    ).bam

    /* Convert to BED */
    BedRNA1 = BAM2BED_RNA1 ( BamRNA1 ).bed
    BedRNA2 = BAM2BED_RNA2 ( BamRNA2 ).bed
    BedDNA  = BAM2BED_DNA  ( BamDNA  ).bed

    /* Annotate restriction recognition sites */
    TableRestrictionRNA1_mult = BedRNA1.cross( Table ) { it -> it[0].id }
                .combine( RestrictedExtended )
                .multiMap{ it ->
                        ch1: it[0]
                        ch2: it[1]
                        ch3: it[2] }
    TableRestrictionRNA1 = BED_ANNOTATE_RESTRICTION_RNA1( TableRestrictionRNA1_mult ).table

    TableRestrictionRNA2_mult = BedRNA2.cross( Table ) { it -> it[0].id }
                .combine( RestrictedExtended )
                .multiMap{ it ->
                            ch1: it[0]
                            ch2: it[1]
                            ch3: it[2] }
    TableRestrictionRNA2 = BED_ANNOTATE_RESTRICTION_RNA2( TableRestrictionRNA2_mult ).table

    // Combined table of restriction:
    // Channel strcuture: [meta, [tables]]
    TableRestriction = TableRestrictionRNA1.groupTuple()
        .cross( TableRestrictionRNA2.groupTuple() ){ it -> it[0].id }
        .map{t1, t2 -> [t1[0], [*t1[1], *t2[1]]]}

    // Collect final filters:
    def FiltersColumns = [:] // todo: add necessary check here
    for (filter_name in params.filters.keySet()){
        def filter_params = params.filters[filter_name]
        if (filter_params instanceof List) {
            filter_params = "(" + filter_params.join(") | (") + ")"
        }
        FiltersColumns[filter_name] = filter_params
    }
    Filters = Channel.from(FiltersColumns)

    TablesCombined = Table.combine(TableFragments, by:0)
         .cross(TableRestriction){ it -> it[0].id }
         .map{t1, t2 -> [t1[0], [*t1[1..-1], *t2[1]]]}

//    // Evaluate filters on collected tables:
//    // Multichannel with collected tables and filters:
//    TableCombined_mult = TablesCombined
//                            .combine( Filters.map{ it -> [it] } ).view()
//                            .multiMap{ it ->
//                                  pq: [it[0], it[1]]
//                                  filters: it[2]
//                            }
//    TableFilters = TABLE_EVALUATE_FILTERS( TableCombined_mult ).table
//    TableFilters.view()

//    // TODO: add custom tables as output
////    ChromSizes = GENOME_PREPARE.out.chromsizes
////    COOLER_MAKE (
////        FinalTable,
////        ChromSizes
////    )

}

workflow {

    REDC ( )

}
workflow.onComplete {
    log.info "Done!"
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
