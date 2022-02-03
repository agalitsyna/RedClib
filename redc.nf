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
def protocol = params.get('protocol', [:])
def chunksize = protocol.get('chunksize', 1000000000000)
def check_restriction = protocol.get('check_restriction', false)
def RenzymesPreloaded = Genome.get('restricted', [:])

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

include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
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

include { RNADNATOOLS_TABLE_ALIGN as TABLE_BED_RNA1 } from './modules/rnadnatools/table_align/main'  addParams( options: [
                                                            args: [
                                                             input_format: 'tsv',
                                                             ref_format: 'parquet',
                                                             output_format: 'parquet',
                                                             params: '--no-input-header --fill-values .,-1,-1,.,0,.,. \
                                                             --new-colnames chrom_rna1,start_rna1,end_rna1,readID,mapq_rna1,strand_rna1,cigar_rna1 \
                                                             --key-column 3 --ref-colname readID'
                                                            ], suffix:'.bed-rna1'] )
include { RNADNATOOLS_TABLE_ALIGN as TABLE_BED_RNA2 } from './modules/rnadnatools/table_align/main'  addParams( options: [
                                                            args: [
                                                             input_format: 'tsv',
                                                             ref_format: 'parquet',
                                                             output_format: 'parquet',
                                                             params: '--no-input-header --fill-values .,-1,-1,.,0,.,. \
                                                             --new-colnames chrom_rna2,start_rna2,end_rna2,readID,mapq_rna2,strand_rna2,cigar_rna2 \
                                                             --key-column 3 --ref-colname readID'
                                                            ], suffix:'.bed-rna2'] )
include { RNADNATOOLS_TABLE_ALIGN as TABLE_BED_DNA } from './modules/rnadnatools/table_align/main'  addParams( options: [
                                                            args: [
                                                             input_format: 'tsv',
                                                             ref_format: 'parquet',
                                                             output_format: 'parquet',
                                                             params: '--no-input-header --fill-values .,-1,-1,.,0,.,. \
                                                             --new-colnames chrom_dna,start_dna,end_dna,readID,mapq_dna,strand_dna,cigar_dna \
                                                             --key-column 3 --ref-colname readID'
                                                            ], suffix:'.bed-dna'] )

include { RNADNATOOLS_TABLE_ALIGN as TABLE_DEDUP_ALIGN } from './modules/rnadnatools/table_align/main'  addParams( options: [
                                                            args: [
                                                             input_format: 'tsv',
                                                             ref_format: 'parquet',
                                                             output_format: 'parquet',
                                                             params: '--no-input-header --key-column 0 --ref-colname readID \
                                                             --fill-values ".",0 --drop-key --new-colnames readID,isUnique'
                                                            ], suffix:'.fastuniq'] )

include { RNADNATOOLS_TABLE_STACK as RESULTS_STACK } from './modules/rnadnatools/table_stack/main' addParams( options: [
                                                             args: [
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ]])

//include { COOLER_MAKE  } from './modules/local/cooler_make/main' addParams( options: [assembly: Assembly, resolution: params.get("cooler_resolution", 1000000)])


// Define workflow
workflow REDC {

    /* Prepare input */
    // Check input FASTQ files
    Fastq = INPUT_CHECK_DOWNLOAD( InputSamplesheet )

    /* Check quality with FASTQC */
    FASTQC( Fastq )

    /* Split into chunks */
    if (chunksize) {
        FastqChunks = INPUT_SPLIT( Fastq ) /* Subworkflow that takes fastq stream and outputs chunks */
    }
    else {
        FastqChunks = Fastq
    }

    /* Prepare genome and annotations */

    // Index, chromosome sizes and fasta:
    GenomeComplete = GENOME_PREPARE( Assembly )
    Hisat2Index = GenomeComplete.genome_index
    GenomeFasta = GenomeComplete.genome_fasta

    // Restrict the genome or load pre-computed restriction sites
    if (check_restriction) {
        Renzymes = Channel.fromList(params.protocol.renzymes)
                            .branch { it ->
                                    forRestriction : !RenzymesPreloaded.containsKey(it) // Restrict only missing renzymes
                                        return [renzyme: it, assembly: Assembly]
                                    loaded : RenzymesPreloaded.containsKey(it) // Load files if restriction is pre-computed
                                        return [[renzyme: it, assembly: Assembly], file(RenzymesPreloaded[it])]
                                    }

        Restricted = GENOME_RESTRICT(Renzymes.forRestriction, GenomeFasta).genome_restricted.mix( Renzymes.loaded )
        // Create separate channels for future searchs:
        RestrictedPlus = Restricted.map{ update_meta(it, [renz_strand: "+"] ) }
        RestrictedMinus = Restricted.map{ update_meta(it, [renz_strand: "-"] ) }
        RestrictedBoth = Restricted.map{ update_meta(it, [renz_strand: "b"] ) }
        RestrictedChannel = RestrictedPlus.mix(RestrictedMinus).mix(RestrictedBoth).map{ it -> [it] }
    }

    // Get the splice sites by hisat2 script:
    SpliceSites = GENOME_PREPARE_RNA_ANNOTATIONS( Assembly ).genome_splicesites

    /* Convert input reads to text table */
    TableChunks = TABLE_FASTQ2TSV( FastqChunks ).table

    /* Trim reads by quality */
    TrimTable = TABLE_TRIM( FastqChunks, TableChunks ).trimtable

    /* Map and check oligos */
    Hits = OLIGOS_MAP( FastqChunks, TableChunks ) /* Subworkflow for oligos mapping */

    /* Tables with read, trimming and mapped oligo: */
    // Create synchronized multichannel with tables:
    InputTablesGrouped = TableChunks.mix(TrimTable)
                                    .mix(Hits.HitsOligos)
                                    .mix(Hits.HitsShortOligos)
                                    .mix(Hits.HitsComplementary)
                                    .map{ removeKeys(it, ['oligo', 'side', 'idx'] ) }
                                    .groupTuple(by:0)
                                    .combine(Channel.from(['']))
                                    .multiMap{it ->
                                        table: [it[0], it[1]]
                                        suffixes: it[2]
                                    }

    MergedTable = TABLE_MERGE( InputTablesGrouped ).table
    Table = TABLE_CONVERT( MergedTable ).table

    /* Deduplicate input sequences */
    FastuniqOut = TABLE_DEDUP( Fastq )
    TableDedup = TABLE_DEDUP_ALIGN( Table.combine(FastuniqOut)
                .filter{ it[0].original_id==it[2].id }
                .multiMap{ it ->
                    reference: [ it[0], it[3] ]
                    table: [ it[0], it[1] ]
                } ).table

    /* Create table with information about fragments. */
    def FragmentColumns = [:]
    for (fragment_name in params.fragments.keySet()){ // select all fragments
        fragment_params = params.fragments[fragment_name] // take the parameters for each fragment
        for( col in fragment_params.new_columns.keySet() ) { // add new columns to the list
            FragmentColumns[col] = fragment_params.new_columns[col]
        }
    }
    FragmentFilters = Channel.from(FragmentColumns)
    // Create multichannel by combining with filters:
    TableFragments_mult = Table.combine( FragmentFilters )
                        .multiMap { meta, id, filt ->
                            table: [meta, id]
                            filters: filt
                        }
    TableFragments = TABLE_EVALUATE_FRAGMENTS( TableFragments_mult ).table

    /* Extract FASTQs with fragments */
    def FragmentData = params.fragments.collect {
            fragment, fragment_params -> [fragment: fragment, params:fragment_params]
            }
    FragmentFilterData = Channel.fromList(FragmentData)
    FragmentsSelected = Table.combine( TableFragments, by:0 )
                        .combine( FragmentFilterData )
                        .map { meta, table1, table2, fragment -> update_meta([meta, [table1, table2]],
                                        [fragment: fragment.fragment, side: fragment.params.side,
                                        single_end: true, selection_criteria: fragment.params.selection_criteria] ) }

    FragmentsFastqNoExtended = FRAGMENTS_TO_FASTQ( FragmentsSelected ).fastq.map{ update_meta( it, [extended:false] ) }

    /* Extend fragments, if needed: */

    FragmentFilterData.filter{it.params.get('extension_suffix', '')}.view()


//    /* Define channels for Fastq files with fragments: */
//    FastqRNA1 = FragmentsFastq.filter{ it[0].fragment == 'rna1' }
//    FastqRNA2 = FragmentsFastq.filter{ it[0].fragment == 'rna2' }
//    if (dna_extension) {
//        FastqDNA = FASTQ_EXTEND_DNA( FragmentsFastq.filter{ it[0].fragment == 'dna' } ).fastq
//    } else {
//        FastqDNA = FragmentsFastq.filter{ it[0].fragment == 'dna' }
//    }
//
//    /* Map the fragments */
//    BamRNA1 = HISAT2_ALIGN_RNA1 (
//            FastqRNA1,
//            Hisat2Index,
//            SpliceSites
//    ).bam
//
//    BamRNA2 = HISAT2_ALIGN_RNA2 (
//            FastqRNA2,
//            Hisat2Index,
//            SpliceSites
//    ).bam
//
//    BamDNA = HISAT2_ALIGN_DNA (
//            FastqDNA,
//            Hisat2Index,
//            SpliceSites
//    ).bam
//
//    /* Convert to BED */
//    BedRNA1 = BAM2BED_RNA1( BamRNA1 ).bed
//    BedRNA2 = BAM2BED_RNA2( BamRNA2 ).bed
//    BedDNA  = BAM2BED_DNA( BamDNA ).bed
//
//    TableBedRNA1 = TABLE_BED_RNA1( BedRNA1.cross(Table){ it -> it[0].id }.multiMap{ it ->
//        input: [ it[0], it[1] ]
//        reference: [ it[2], it[3] ]
//         } ).table
//    TableBedRNA2 = TABLE_BED_RNA2( BedRNA2.cross(Table){ it -> it[0].id }.multiMap{ it ->
//        input: [ it[0], it[1] ]
//        reference: [ it[2], it[3] ]
//         } ).table
//    TableBedDNA = TABLE_BED_DNA( BedDNA.cross(Table){ it -> it[0].id }.view().multiMap{ it ->
//        input: [ it[0], it[1] ]
//        reference: [ it[2], it[3] ]
//         } ).table
//
//    TableBed = TableBedRNA1.mix( TableBedRNA2 ).mix( TableBedDNA )
//         .map{ removeKeys(it, ['fragment', 'side', 'selection_criteria', 'singe-end'] )}
//         .groupTuple(by:0)
//
//    /* Annotate restriction recognition sites */
//    TableRestrictionRNA1_mult = BedRNA1.cross(Table){ it -> it[0].id }
//                .combine( RestrictedExtended )
//                .multiMap{ it ->
//                        bed: it[0]
//                        table_input: it[1]
//                        rsites: it[2] }
//    TableRestrictionRNA1 = BED_ANNOTATE_RESTRICTION_RNA1( TableRestrictionRNA1_mult ).table
//
//    TableRestrictionRNA2_mult = BedRNA2.cross( Table ) { it -> it[0].id }
//                .combine( RestrictedExtended )
//                .multiMap{ it ->
//                        bed: it[0]
//                        table_input: it[1]
//                        rsites: it[2] }
//    TableRestrictionRNA2 = BED_ANNOTATE_RESTRICTION_RNA2( TableRestrictionRNA2_mult ).table
//
//    /* Collect final tables */
//    TablesCombined = Table
//        .mix(TableBedRNA1)
//        .mix(TableBedRNA2)
//        .mix(TableBedDNA)
//        .mix(TableRestrictionRNA1)
//        .mix(TableRestrictionRNA2)
//        .mix(TableDedup)
//        .mix(TableFragments)
//        .map{ removeKeys(it, ['fragment', 'side', 'selection_criteria', 'single_end'] ) }.groupTuple(by:0)
//
//    /*  Collect final filters: */
//    def FiltersColumns = [:] // todo: add necessary check here
//    for (filter_name in params.filters.keySet()){
//        def filter_params = params.filters[filter_name]
//        if (filter_params instanceof List) {
//            filter_params = "(" + filter_params.join(") | (") + ")"
//        }
//        FiltersColumns[filter_name] = filter_params
//    }
//
//    /* Evaluate filters on collected tables: */
//    TableFinalChunked = TABLE_EVALUATE_FILTERS( TablesCombined, Channel.from(FiltersColumns) ).table
//
//    /* Merge final tables between replicates and chunks: */
//    TableFinalChunked.view() //.map{}.groupTuple(by:0, sort  { it1, it2 -> return it1[0].id==it2[0].id  })
//
//    // TableFinal = RESULTS_STACK( TableFinalChunked ).table
//
////    if (params.get('output', [:]).get('make_final_table', false)){
////        TablesCollection = params.output.tables // collect // filter
////        RESULTS_DUMP_TABLE(TableFilters, FinalColumns)
////
////        if (params.get('output', [:]).get('make_cooler', false)){
////            ChromSizes = GenomeComplete.chromsizes
////            COOLER_MAKE(TableFilters, ChromSizes)
////        }
////    }



}

workflow {

    REDC ( )

}
workflow.onComplete {
    log.info "Done!"
}

def update_meta( it, hashMap ) {
    def meta = [:]
    def keys = it[0].keySet() as String[]
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
def removeKeys( it, ks ) {
    def meta = [:]
    def keys = it[0].keySet() as String[]
    for( key in keys ) {
        if (!(key in ks)) {
            meta[key] = it[0][key]
        }
    }

    def array = [ meta, it[1] ]
    return array
}