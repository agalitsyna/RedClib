//
// Map oligos and check strict substrings in these oligos for RedC.
//

params.options = [:]

include { RKLIB_SEQ2HASH as BIN_OLIGOS } from '../../modules/rklib/rk_seq2hash/main' addParams( options: [args: [mode: 'fasta']] ) // Bin input oligos from fasta file
include { RKLIB_SEQ2HASH as BIN_FASTQ }  from '../../modules/rklib/rk_seq2hash/main' addParams( options: [args: [mode: 'fastq']] ) // Bin input fastq

include { RKLIB_QUERYSEARCH as OLIGOS_ALIGN} from '../../modules/rklib/rk_querysearch/main' addParams( options: [args: [:]] ) // Align oligos
include { RKLIB_CHECK_COMPLEMENTARY as OLIGOS_CHECK_COMPLEMENTARY } from '../../modules/rklib/rk_check_complementary/main' addParams( options: [args: [:]] ) // Align complementary fragments of RNA

include { RNADNATOOLS_READ_CHECK_NUCLEOTIDES as OLIGOS_CHECK_SHORT } from '../../modules/rnadnatools/read_check_nucleotides/main' addParams( options: [args: [:]] ) // Check oligos presence

include { TSV_MERGE as TABLE_MERGE_OLIGOS } from '../../modules/local/tsv_merge/main' addParams( options: [args: [:], suffix: '.oligos'] )

workflow OLIGOS_MAP {
    take:
        FastqChunks // Channel with loaded FastqChunks
        TableChunks // Channel with loaded TableChunks

    main:

        if (!params.containsKey('oligos')){
            HitsOligos = Channel.empty()
        } else {

            /* I. Collect oligos from the params file */
            Oligos = Channel
                       .from(params.oligos.keySet().withIndex())
                       .map { it, idx -> [it, idx, params.oligos[it]] } // construct meta for each oligo
                       .map { oligo_name, idx, meta ->
                               meta.single_end = true
                               meta.id = oligo_name
                               meta.idx = idx
                           [meta, file(meta.file)]
                       }
            def nOligos = params.oligos.keySet().size()

            /* II. Check for the presence of oligos with a Rabin-Karp (rklib) */
            /* Step 1: Index oligos and input reads */
            IndexedOligos = BIN_OLIGOS( Oligos ).bin
            IndexedFastqs = BIN_FASTQ( FastqChunks ).bin

            IndexedLibrary = IndexedFastqs.combine(IndexedOligos)
                                            .multiMap { it ->
                                                bin_reads: update_meta( [it[0], it[1]],
                                                                        [oligo: it[2].id, side: it[2].side, idx: it[2].idx] )
                                                bin_oligos: [it[2], it[3]]
                                            }
            /* Step 2: Check oligos presence */
            HitsOligosStream = OLIGOS_ALIGN( IndexedLibrary ).aligned

            /* Step 3: Format the tables and channel: */
            HitsOligosGrouped = HitsOligosStream
                                .map{ it -> [ it[0].id, it ] }
                                .groupTuple(by: 0, sort: { a, b -> a[0].idx <=> b[0].idx }, size: nOligos)
                                .map{ it -> it.collect()[1].collect{ item -> [item[0], file(item[1]), item[0].oligo+'_R'+item[0].side] }.transpose() }
                                .multiMap{meta, files, suffixes ->
                                    dataset: [meta[0],  files]
                                    suffixes: suffixes
                               }

            HitsOligos = TABLE_MERGE_OLIGOS( HitsOligosGrouped ).table

        }

        if (!params.containsKey('short_oligos')){
            HitsShortOligos = Channel.empty()
        } else {

            /* III. Collect short oligos and check their presence with rnadnatools */
            /* Step 1: Collect */
            ShortOligos = Channel
                           .from( params.short_oligos.keySet().withIndex() )
                           .map { it, idx -> [it, idx, params.short_oligos[it]] } // meta for each short_oligo
                           .map { short_oligo_name, idx, meta ->
                                    meta.id = short_oligo_name
                                    meta.idx = idx
                                    meta.single_end = true
                                [meta]
                           }

            /* Step 2: Check presence of short oligos with rnadnatools: */
            InputShortOligos = TableChunks.combine(ShortOligos).combine(HitsOligosStream) // Select matching tables
                           .filter{ meta, table, meta_short_oligo, meta_ref, ref ->
                                        (meta.id==meta_ref.id &&
                                        meta_ref.oligo==meta_short_oligo.reference_oligo &&
                                        meta_ref.side==meta_short_oligo.side) }
                           .multiMap{ meta, table, meta_short_oligo, meta_ref, ref ->
                                       input: [meta, table]
                                       meta_short_oligo: meta_short_oligo
                                       reference: [meta_ref, ref]
                           }

            HitsShortOligos = OLIGOS_CHECK_SHORT( InputShortOligos ).hits
        }

        if (!params.containsKey('complementary')){
            HitsComplementary = Channel.empty()
        } else {
            /* Step 1: Collect information about left and right side of complementary regions */
            ComplementaryInput = Channel.from(params.complementary)

            InputComplementary = TableChunks
                .combine(IndexedFastqs)
                .combine(HitsOligosStream)
                .combine(HitsOligosStream)
                .combine(ComplementaryInput)
                .filter{meta, table, meta_reads, bin_reads,
                            meta_oligos_left, aligned_left,
                            meta_oligos_right, aligned_right, meta_compl ->
                        (meta.id==meta_reads.id &&
                        meta.id==meta_oligos_left.id &&
                        meta.id==meta_oligos_right.id &&
                        meta_oligos_left.oligo==meta_compl.left_reference_oligo &&
                        meta_oligos_right.oligo==meta_compl.right_reference_oligo)
                }
                .multiMap{meta, table, meta_reads, bin_reads,
                            meta_oligos_left, aligned_left,
                            meta_oligos_right, aligned_right, meta_compl ->
                    input_table: [meta, table]
                    input_reads: [meta_reads, bin_reads]
                    input_oligos_left: [meta_oligos_left, aligned_left]
                    input_oligos_right: [meta_oligos_right, aligned_right]
                    input_meta_compl: meta_compl
                    }

            /* Step 2: Check complementary RNA ends: */
            HitsComplementary = OLIGOS_CHECK_COMPLEMENTARY( InputComplementary ).hits
        }

    emit:
    HitsOligos // channel: [ val(meta), [ output_table ] ] or empty
    HitsShortOligos // channel: [ val(meta), [ output_table ] ] or empty
    HitsComplementary // channel: [ val(meta), [ output_table ] ] or empty
}

def update_meta( it, hashMap ) {
    def meta = [:]
    keys = it[0].keySet() as String[]
    for( key in keys ) {
        meta[key] = it[0][key]
    }
    def file1 = it[1][0]
    def file2 = ""

    if (!meta.single_end) {
        file2 = it[1][1]
    }

    keys_new = hashMap.keySet() as String[]
    for( key in keys_new ) {
        meta[key] = hashMap[key]
    }

   if (meta.single_end) {
        array = [ meta, [ file(file1) ] ]
    } else {
        array = [ meta, [ file(file1), file(file2) ] ]
    }
    return array
}