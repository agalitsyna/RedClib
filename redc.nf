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
if (params.getOrDefault('help', 'false').toBoolean()) {
    helpMessage()
    exit 0
}

// Check required parameters
if (params.input) { InputSamplesheet = params.input }
    else { exit 1, 'Input samplesheet not specified!' }

Genome = params.getOrDefault('genome', [:])
if (params.genome.assembly) { Assembly = params.genome.assembly }
    else { exit 1, 'Genome assembly is not specified!' }
if (params.genome.genes_gtf) { GenesGtf = params.genome.genes_gtf }
    else { exit 1, 'Genome RNA annotation GTF is not specified!' }

// Check optional parameters
def protocol = params.getOrDefault('protocol', [:])
def chunksize = protocol.getOrDefault('chunksize', 100000000000)
def check_restriction = protocol.getOrDefault('check_restriction', false)
def RenzymesPreloaded = Genome.getOrDefault('restricted', [:])

// Include modules and subworkflows
include { INPUT_CHECK_DOWNLOAD } from './subworkflows/local/input_check' addParams( options: [:] )
include { INPUT_SPLIT          } from './subworkflows/local/input_split' addParams( options: [chunksize: chunksize] )

include { GENOME_PREPARE } from './modules/local/genome_prepare/main'  addParams( options: [
                             args: [
                                genome: [chromsizes: Genome.getOrDefault("chromsizes", ""),
                                index_prefix: Genome.getOrDefault("index_prefix", ""),
                                fasta: Genome.getOrDefault("fasta", "")],
                                auto_download_genome: Genome.getOrDefault("auto_download_genome", true)
                            ] ] )
include { RNADNATOOLS_GENOME_RECSITES as GENOME_RESTRICT } from './modules/rnadnatools/genome_recsites/main' addParams( options: [:])
include { GENOME_PREPARE_RNA_ANNOTATIONS } from './modules/local/genome_prepare_rna_annotations/main' addParams( options: [
                             args: [
                                genes_gtf: GenesGtf,
                                rna_annotation_suffix: Genome.getOrDefault('rna_annotation_suffix', '')
                            ]])

include { FASTQC } from './modules/nf-core/software/fastqc/main' addParams( options: [:] )
include { FASTQ2TSV as TABLE_FASTQ2TSV } from './modules/local/fastq2table/main' addParams( options: [:] )
include { DEDUP_FASTUNIQ as TABLE_DEDUP } from './subworkflows/local/dedup_fastuniq' addParams( options: [:] )
include { TRIMTABLE_CREATE as TABLE_TRIM } from './subworkflows/local/trimtable_create' addParams( options: [:] )

include { RNADNATOOLS_TABLE_ALIGN as TABLE_DEDUP_ALIGN } from './modules/rnadnatools/table_align/main'  addParams( options: [
                                                            args: [
                                                             input_format: 'tsv',
                                                             ref_format: 'parquet',
                                                             output_format: 'parquet',
                                                             chunksize: 10000000,
                                                             chunksize_writer: 100000,
                                                             params: '--no-input-header --key-column 0 --ref-colname readID \
                                                             --fill-values ".",0 --drop-key --new-colnames readID,isUnique'
                                                            ], suffix:'.fastuniq'] )

include { OLIGOS_MAP } from './subworkflows/local/oligos_map' addParams( options: [:] )

include { TSV_MERGE as TABLE_MERGE } from './modules/local/tsv_merge/main' addParams( options: [
                                                             args: [:],
                                                             suffix: '.table'] )

include { RNADNATOOLS_TABLE_CONVERT as TABLE_CONVERT } from './modules/rnadnatools/table_convert/main' addParams( options: [
                                                             args: [
                                                             input_format: 'tsv',
                                                             output_format: 'parquet'
                                                             ]])

include { RNADNATOOLS_TABLE_EVALUATE as TABLE_EVALUATE_FRAGMENTS } from './modules/rnadnatools/table_evaluate/main'  addParams( options: [
                                                             args: [
                                                             format:'int',
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.fragment-columns'] )

include { RNADNATOOLS_SEGMENT_EXTRACT_FASTQ as FRAGMENTS_TO_FASTQ } from './modules/rnadnatools/segment_extract_fastq/main'  addParams( options: [
                                                             args: [
                                                             input_format: 'parquet'
                                                             ], suffix:'.fragments'] )

include { FASTQ_EXTEND } from './modules/local/fastq_extend/main' addParams( options: [
                                                             args: [
                                                             gzip: true
                                                             ], suffix:'.ext'] )

include { HISAT2_ALIGN } from './modules/local/hisat2/main' addParams( options: [args: [:]] )
include { BAM2BED } from './modules/local/bam2bed/main' addParams( options: [args: [extraTags: ['NH', 'XM'], samFilter: '']])

include { RNADNATOOLS_SEGMENT_GETCLOSEST as TABLE_ANNOTATE_RESTRICTION } from './modules/rnadnatools/segment_getclosest/main'  addParams( options: [args: [
                                                             input_format: 'parquet',
                                                             reference_format: 'tsv',
                                                             output_format: 'parquet',
                                                             params: '--key-columns chrom,start,end --ref-columns 0,1,5 --no-ref-header',
                                                             chunksize: 1000000
                                                             ]] )

include { RNADNATOOLS_TABLE_EVALUATE as TABLE_EVALUATE_FILTERS } from './modules/rnadnatools/table_evaluate/main'  addParams( options: [
                                                             args: [
                                                             format:'bool',
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.filters'] )

include { RNADNATOOLS_TABLE_ALIGN as TABLE_BED } from './modules/rnadnatools/table_align/main'  addParams( options: [
                                                            args: [
                                                             input_format: 'tsv',
                                                             ref_format: 'parquet',
                                                             output_format: 'parquet',
                                                             chunksize: 10000000,
                                                             chunksize_writer: 100000,
                                                             params: '--no-input-header --fill-values .,-1,-1,.,0,.,.,-1,-1 \
                                                             --new-colnames chrom,start,end,readID,mapq,strand,cigar,nMultiMap,nSub \
                                                             --key-column 3 --ref-colname readID'
                                                            ]] )

include { RNADNATOOLS_TABLE_MERGE as TABLE_FRAGMENTS_MERGE } from './modules/rnadnatools/table_merge/main' addParams( options: [
                                                             args: [
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.fragments'] )

include { RNADNATOOLS_TABLE_MERGE as TABLE_RESTRICTION_MERGE } from './modules/rnadnatools/table_merge/main' addParams( options: [
                                                             args: [
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ], suffix:'.restriction'] )

include { RNADNATOOLS_TABLE_DUMP as RESULTS_DUMP } from './modules/rnadnatools/table_dump/main' addParams( options: [
                                                             args: [
                                                             input_format: 'parquet',
                                                             output_format: 'parquet'
                                                             ]])

include { RNADNATOOLS_TABLE_STACK as RESULTS_STACK } from './modules/rnadnatools/table_stack/main' addParams( options: [
                                                             args: [
                                                             input_format: 'parquet',
                                                             output_format: 'tsv'
                                                             ]])

include { RNADNATOOLS_TABLE_STATS as RESULTS_STATS } from './modules/rnadnatools/table_stats/main' addParams( options: [
                                                             args: [
                                                             input_format: 'tsv'
                                                             ]])

// Extra keys of the metadata that might be dynamically added in the pipeline:
def extraKeys = ['fragment', 'side', 'selection_criteria', 'single_end', 'mapping_args', 'restriction', 'extended', 'ext_suffix', 'ext_prefix']
def oligoKeys = ['oligo', 'side', 'idx']
def dedupKeys = ['single_end']

include { COOLER_MAKE  } from './modules/local/cooler_make/main' addParams( options: [args:[assembly: Assembly]] )


// Define workflow
workflow REDC {

    /* Prepare input */
    // Check input FASTQ files
    Fastq = INPUT_CHECK_DOWNLOAD( file(InputSamplesheet) )

//    /* Check quality with FASTQC */
//    FASTQC( Fastq )

    /* Split into chunks */
    if (chunksize) {
        /* Run subworkflow that takes fastq stream and outputs chunks: */
        FastqChunks = INPUT_SPLIT( Fastq )
    }
    else {
        FastqChunks = Fastq
    }

    /* Prepare genome and annotations */

    // Retrieve index, chromosome sizes and fasta:
    AssemblyInput = Channel.from([[Assembly, file(InputSamplesheet)]])
    GenomeComplete = GENOME_PREPARE( AssemblyInput )
    Hisat2Index = GenomeComplete.genome_index
    GenomeFasta = GenomeComplete.genome_fasta
    ChromSizes = GenomeComplete.genome_chromsizes

    // Restrict the genome or load pre-computed restriction sites
    if (check_restriction) {
        Renzymes = Channel.fromList(params.protocol.renzymes)
                             // Restrict only missing renzymes, create two channel branches:
                             // forRestriction and loaded (restriction is pre-computed)
                                .branch { it ->
                                    forRestriction : !RenzymesPreloaded.containsKey(it)
                                        return [renzyme: it, assembly: Assembly]
                                    loaded : RenzymesPreloaded.containsKey(it)
                                        return [[renzyme: it, assembly: Assembly], file(RenzymesPreloaded[it])]
                                }

        Restricted = GENOME_RESTRICT(Renzymes.forRestriction, GenomeFasta).genome_restricted.mix( Renzymes.loaded )
        // Create separate channels for +, - and +/- (both) strand orientations of enzymes recognition:
        RestrictedPlus = Restricted.map{ update_meta(it, [renz_strand: "+"] ) }
        RestrictedMinus = Restricted.map{ update_meta(it, [renz_strand: "-"] ) }
        RestrictedBoth = Restricted.map{ update_meta(it, [renz_strand: "b"] ) }
        // Mix all channels with annotated restriction recognition sites to a single channel:
        RestrictedChannel = RestrictedPlus.mix(RestrictedMinus).mix(RestrictedBoth)
    }

    // Get genomic splice sites by hisat2 script:
    SpliceSites = GENOME_PREPARE_RNA_ANNOTATIONS( Assembly ).genome_splicesites

    /* Convert input reads to text table */
    TableChunks = TABLE_FASTQ2TSV( FastqChunks ).table

    /* Trim reads by quality */
    TrimTable = TABLE_TRIM( FastqChunks, TableChunks ).trimtable

    /* Map and check oligos */
    Hits = OLIGOS_MAP( FastqChunks, TableChunks ) /* Subworkflow for oligos mapping */

    /* Tables with read, trimming and mapped oligo: */
    // Create synchronized multichannel with tables:
    TablesGroupedInput = TableChunks.mix(TrimTable) // Mix all oligos-related channels into a single channel:
                                    .mix(Hits.HitsOligos)
                                    .mix(Hits.HitsShortOligos)
                                    .mix(Hits.HitsComplementary)
                               // Remove oligo parameters from metadata:
                                    .map{ removeKeys(it, oligoKeys ) }
                               // Group by metadata, which should be identical for the same sample/chunk
                               // Note that group size depends on the number of oligos and short oligos:
                                    .groupTuple(by:0, sort:true)
                               // Combine with channel with empty suffixes
                                    .combine(Channel.from(['']))
                               // Create restructured multi-channel:
                                    .multiMap{it ->
                                        table: [it[0], it[1]]
                                        suffixes: it[2]
                                    }
    // Merge all the tables in a single file for each sample/chunk:
    MergedTable = TABLE_MERGE( TablesGroupedInput ).table

    // Convert the table file format for better storage/operations:
    Table = TABLE_CONVERT( MergedTable ).table

    /* Deduplicate input sequences */

    // Run subworkflow that trims first basepairs of reads and runs fastuniq on them:
    FastuniqOut = TABLE_DEDUP( Fastq )

    // Align deduplicated reads with read table
    // (we want to guarantee the same order of reads in each table):
    DedupAlignInput = Table.combine(FastuniqOut).filter{ it[0].original_id==it[2].id }
                                .multiMap{ it ->
                                    reference: [ it[0], it[3] ]
                                    table: [ it[0], it[1] ]
                                }
    TableDedup = TABLE_DEDUP_ALIGN( DedupAlignInput ).table.map{ removeKeys(it, dedupKeys) }

    /* Collect fragment columns data. */
    // params.fragments has the list of new columns with expressions that will be evaluated for each fragment
    // Take params.fragments and collect columns:
    def FragmentCollection = [:]
    for (fragment_name in params.fragments.keySet()){ // select all fragments
        fragment_params = params.fragments[fragment_name] // take the parameters for each fragment
        for( col in fragment_params.new_columns.keySet() ) { // add new columns to the list
            FragmentCollection[col] = fragment_params.new_columns[col]
        }
    }

    // Create multichannel by combining table with collection of columns:
    TableFragmentsColumnsInput = Table.combine( Channel.from(FragmentCollection) )
                                .multiMap { meta, id, col ->
                                    table: [meta, id]
                                    filters: col
                                }

    // Evaluate new columns for fragments:
    TableFragmentsColumns = TABLE_EVALUATE_FRAGMENTS( TableFragmentsColumnsInput ).table

    /* Collect fragments data: */
    def FragmentData = params.fragments.collect {
            fragment, fragment_params -> [fragment: fragment, params:fragment_params]
            }
    FragmentsSelected = Table.combine( TableFragmentsColumns, by:0 ).combine( Channel.fromList(FragmentData) )
                        // We now have a combined channel with structure: meta, table1, table2, fragment.
                        // Re-structure this channel to [meta, files] and store fragment data in meta:
                            .map { meta, table1, table2, fragment -> update_meta(
                                    [meta, [table1, table2]],
                                        [fragment: fragment.fragment,
                                        side: fragment.params.side,
                                        single_end: true,
                                        selection_criteria: fragment.params.selection_criteria,
                                        ext_suffix: fragment.params.getOrDefault('extension_suffix', ''),
                                        ext_prefix: fragment.params.getOrDefault('extension_prefix', ''),
                                        mapping_args: fragment.params.mapping_args,
                                        restriction: fragment.params.getOrDefault('annotate_restriction', [['', '']])]
                                    ) }
                        // Tag it with the frag name (allows to recognize separate instances and have unique file names):
                            .map { add_tag(it, it[0].fragment) }

    /* Extract FASTQs with fragments */
    FragmentsFastqPreExtend = FRAGMENTS_TO_FASTQ( FragmentsSelected ).fastq // Convert to fastq
                        // Split channels that will be extended with nucleotide suffix and prefix:
                            .branch { it ->
                                noExtend : (!(it[0].ext_suffix || it[0].ext_prefix))
                                    return update_meta( it, [extended:false] )
                                forExtend : (it[0].ext_suffix || it[0].ext_prefix)
                                    return update_meta( it, [extended:true] )
                                }
    // Extend the fragments for which extension was requested:
    FragmentsFastqExtended = FASTQ_EXTEND( FragmentsFastqPreExtend.forExtend ).fastq
    // Mix channels with and without extension into a single channel:
    FragmentsFastq = FragmentsFastqPreExtend.noExtend.mix(FragmentsFastqExtended)

    /* Map the fragments */
    // Combine fragments with splicesites and index to make the operation deterministic (enables -resume)
    Hisat2Input = FragmentsFastq.combine(SpliceSites).combine(Hisat2Index)
                    // Channel structure: meta, fastq, splice_sites, ...index_files...
                        .multiMap{it ->
                             FragmentsFastq: [it[0], it[1]]
                             Hisat2Index: it[3..-1]
                             SpliceSites: it[2]
                        }
    // Align with hisat2:
    Bam = HISAT2_ALIGN ( Hisat2Input ).bam

    /* Convert to BED */
    Bed = BAM2BED( Bam ).bed

    /* Convert BED to table */
    // Bed file might have some reads missing or unordered.
    // Convert BED to table where the read order is guaranteed:
    TableBedInput = Bed.combine(Table)
                    // Combine input bed with input table by tagless ids:
                        .filter{ it -> tagless(it[0].id) == it[2].id }
                    // Channel structure: meta_bed, bed, meta_ref, ref
                        .multiMap{ it ->
                            bed: [it[0], it[1]]
                            table_ref: [it[2], it[3]]
                        }
    // Note that you can get an error "No columns to parse" if the reads were not mapped at all:
    TableBedFragments = TABLE_BED( TableBedInput ).table

    TableBedMergeInput = TableBedFragments
                    // Remove rfrag tag from id and restructure the channel:
                        .map{ meta, table -> remove_tag([meta, [table, meta.fragment]]) }
                    // Remove any extra parameters of meta:
                        .map{ removeKeys(it, extraKeys) }
                    // groupTuple with sorting to guarantee deterministic output,
                    // note that size depends on the number of fragments:
                        .groupTuple(by:0, sort:{it->it[1]})
                        .map{meta, it -> [meta, it.transpose().collect()]}
                    // Create multi-channel:
                        .multiMap{ meta, it ->
                            table_bed: [meta, it[0]]
                            suffixes: it[1]
                        }

    TableBed = TABLE_FRAGMENTS_MERGE( TableBedMergeInput ).table


    /* Annotate restriction */
    // Restriction annotation takes BED file, reads table and renzyme file as input
    // and for each fragment reports the closest restriction recognition sites (on the left and on the right).
    // First, create formatted input for the TABLE_ANNOTATE_RESTRICTION process:
    TableBedExpanded = TableBedFragments
                    // Filter out non-restricted fragments:
                        .filter{ it[0].getOrDefault('restriction', '') }
                    // Emit for each restriction enzyme and orientation:
                        .map{it -> [it[0], it[1], it[0].restriction] }.transpose(by:2)
                    // Update meta params:
                        .map{ meta, file, restriction -> update_meta( [meta, file], [restriction: restriction] ) }
                    // Update id to include restriction enzyme:
                        .map { add_tag(it, it[0].restriction.join()) }

    TableRestrictionInput = TableBedExpanded.combine( RestrictedChannel )
                    // Channel: meta_bed, bed, meta_restr, restr
                    // Take renzyme annotation for the same renzyme and strand orientation:
                        .filter{ it -> it[0].restriction == [it[2].renzyme, it[2].renz_strand] }
                    // Multi-channel:
                        .multiMap{ it ->
                            bed: [it[0], it[1]]
                            restr: [it[2], it[3]]
                        }
    // Annotate bed files by the sites of restriction recognition:
    TableRestrictionFragments = TABLE_ANNOTATE_RESTRICTION( TableRestrictionInput ).table

    // Merge restriction tables into a single one:
    TableRestrictionMergeInput = TableRestrictionFragments
                    // Remove restriction tag:
                        .map{ remove_tag(it, 2) }
                    // Remove any extra parameters of meta:
                        .map{ removeKeys(it, extraKeys ) }
                    // Group by meta. Note that size depends on the number of fragments:
                        .groupTuple(by:0, sort:true)
                    // Combine with empty suffixes(required input for TABLE_RESTRICTION_MERGE).
                    // We don't have to add specific suffixes here, because fragments/renzymes
                    // are stored in the columns of individual tables:
                        .combine(Channel.from(''))
                    // Multi-channel:
                        .multiMap{ meta, files, suffixes ->
                            table_bed: [meta, files]
                            suffixes: suffixes
                        }

    TableRestriction = TABLE_RESTRICTION_MERGE( TableRestrictionMergeInput ).table

    /* Evaluation of final filters */
    /* Collect final tables */
    TablesCollected = Table.mix(TableBed)
                        .mix(TableRestriction)
                        .mix(TableDedup)
                        .mix(TableFragmentsColumns)
                    // Remove any extra parameters of meta:
                        .map{ removeKeys(it, extraKeys) }
                    // All channels have the same keys in meta, we can group by it.
                    // Channel: meta, [table, table_bed, table_dedup, table_fragments]
                    // (tables might be in some other order, sorting depends on full file names)
                        .groupTuple(by:0, sort:true)

    /*  Collect final filters: */
    def FilterColumns = [:]
    for (filter_name in params.filters.keySet()){
        def filter_params = params.filters[filter_name]
        if (filter_params instanceof List) {
            filter_params = "(" + filter_params.join(") | (") + ")"
        }
        FilterColumns[filter_name] = filter_params
    }

    /* Evaluate filters on collected tables: */
    TableEvaluateInput = TablesCollected.combine(Channel.from(FilterColumns))
                          .multiMap{ it ->
                                tables: [it[0], it[1]]
                                filters: it[2]
                          }
    TableFinalChunked = TABLE_EVALUATE_FILTERS( TableEvaluateInput ).table

    /* Output tables */
    TableAllChunked = TablesCollected.combine(TableFinalChunked, by:0)
                        .map{it -> [it[0], it[1]+[it[2]]]}

    if (params.output['tables']){

        // Channel: [table_name, filter, [header_columns]]
        FilesCollection = Channel.fromList(
                params.output.tables.collect{ k, v -> [k, v['filter'], v['header'].split(" ")]}
            )

        DumpInput = TableAllChunked.combine(FilesCollection)
                    // Channel: [meta, files, table_name, filter, [header_columns]]
                    // Add tag with table name:
                        .map{ add_tag(it, it[2]) }
                    // Add meta keys. Table_name is for identification of table type.
                    // Columns will guarantee the order in STACK:
                        .map{ update_meta(it, [table_name: it[2], columns:it[4]]) }
                    // Multi-channel:
                        .multiMap{ it ->
                            table: [it[0], it[1]]
                            filter: it[3]
                            columns: it[4]
                        }
        TableDumpedChunked = RESULTS_DUMP( DumpInput ).table
        TableChunked = TableDumpedChunked.mix( TableFinalChunked.map{ update_meta(it, [table_name:'filters']) } )
    }
    else {
        TableChunked = TableFinalChunked.map{ update_meta(it, [table_name:'filters']) }
    }

    /* Merge final tables between replicates: */
    TableStackChunksInput = TableChunked.map{ it -> [[sample:it[0].original_id, table_name:it[0].table_name], [it[0], it[1]]] }
                 // Group chunks and sort by chunk number:
                     .groupTuple(by:0, sort:{it -> it[0].chunk})
                 // Channel: [meta, [..tables..]] (tables were sorted by chunk and table type)
                     .map{id, it -> [it[0][0], it.transpose().collect()[1]]}
                 // Replace id with original_id and remove chunk from meta:
                    .map{ update_meta(it, [id: it[0].original_id]) }
                    .map{ removeKeys(it, ['chunk']) }

    TableFinal = RESULTS_STACK( TableStackChunksInput ).table

    if (params.output['stats']){
        /* Write table with stats */

        // Channel: [stats_name, [header_columns]]
        StatsList = Channel.fromList(
                params.output.stats.collect{ k, v -> [k, v] }
        )

        StatsInput = TableFinal.combine(StatsList)
                    // Channel: [meta, file, stats_name, [cols]]
                    // Filter appropriate tables:
                        .filter{ it[0].table_name==it[3].table_name }
                    // Add tag with stats name:
                        .map{ add_tag(it, it[2]) }
                    // Update meta with the name of stats:
                        .map{ update_meta(it, [stats_name:it[2]]) }
                    // Multi-channel:
                        .multiMap{ it ->
                            table: [it[0], it[1]]
                            columns: it[3].filters
                        }
        TableStats = RESULTS_STATS( StatsInput ).table

    }

    if (params.output['tables'] && params.output['cooler']){
        /* Write coolers */

        // Channel: [cool_name, params]
        CoolerList = Channel.fromList(
                params.output.cooler.collect{ k, v -> [k, v] }
        )

        CoolerInput = TableFinal.combine(CoolerList)
                    // Channel: [meta, file, cool_name, params]
                    // Select only requested dumped table:
                        .filter{ it[0].table_name == it[3].table_name }
                    // Add tag with stats name:
                        .map{ add_tag(it, it[2]) }
                    // Update meta with parameters:
                        .map{ update_meta(it, it[3]) }
                        .map{ update_meta(it, params) }
                    // Combine with chromosome sizes:
                        .combine(ChromSizes)
                    // Channel: [meta, file, cool_name, params, chromsizes]
                    // Multi-channel:
                        .multiMap{ it ->
                            table: [it[0], it[1]]
                            chromsizes: it[4]
                        }
        Coolers = COOLER_MAKE( CoolerInput )

    }

}

workflow {

    REDC ( )

}
workflow.onComplete {
    log.info "Done!"
}

def update_meta( it, hashMap ) {
/* Update the meta hashMap. Takes channel of structure [meta, ..data..]. */
    def meta = [:]
    def keys = it[0].keySet() as String[]
    for( def key in keys ) {
        meta[key] = it[0][key]
    }

    def keys_new = hashMap.keySet() as String[]
    for( def key in keys_new ) {
        meta[key] = hashMap[key]
    }

    def array = [ meta, *it[1..-1] ]
    return array
}

def add_tag( it, tag, id_key='id' ) {
/* Add the tag to id (after dot). Takes channel of structure [meta, ..data..].  */
    def meta = [:]
    def keys = it[0].keySet() as String[]
    for( def key in keys ) {
        meta[key] = it[0][key]
    }
    meta[id_key] = meta[id_key]+'.'+tag
    def array = [ meta, *it[1..-1] ]
    return array
}

def tagless( s, n=1 ) {
/* Remove n tags from string s  */
    def res = s
    for( def i in 1..n ) {
        res = res.take(res.lastIndexOf('.'))
    }
    return res
}

def remove_tag( it, n=1, id_key='id' ) {
/* Remove n tags from id (after dot). Takes channel of structure [meta, ..data..]. */
    def meta = [:]
    def keys = it[0].keySet() as String[]
    for( def key in keys ) {
        meta[key] = it[0][key]
    }
    for( def i in 1..n ) {
        meta[id_key] = meta[id_key].take(meta[id_key].lastIndexOf('.')) // same as: meta[id_key] = tagless(meta[id_key], n)
    }
    def array = [ meta, *it[1..-1] ]
    return array
}

def removeKeys( it, ks ) {
/* Remove keys from hashMap. Take channel of structure [meta, ..data..].  */
    def meta = [:]
    def keys = it[0].keySet() as String[]
    for( def key in keys ) {
        if (!(key in ks)) {
            meta[key] = it[0][key]
        }
    }

    def array = [ meta, *it[1..-1] ]
    return array
}