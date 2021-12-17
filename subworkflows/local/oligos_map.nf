//
// Map oligos and check strict substrings in these oligos for RedC.
//

params.options = [:]

include { BIN_ENCODE as BIN_OLIGOS } from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fasta']] ) // Bin input oligos from fasta file
include { BIN_ENCODE as BIN_FASTQ }  from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fastq']] ) // Bin input fastq

include { OLIGOS_ALIGN } from '../../modules/local/align_encoded/main' addParams( options: [args: [:]] ) // Align oligos
include { OLIGOS_CHECK as OLIGOS_CHECK_GA } from '../../modules/local/oligos_check/main' addParams( options: [args: [oligo: 'GA', position: 35, orientation: 'F']] ) // Check oligos presence
include { OLIGOS_CHECK_COMPLEMENTARY } from '../../modules/local/oligos_complementary/main' addParams( options: [args: [rna_complementary_length: 14]] ) // Align comnplementary fragments of RNA

include { PARQUET_CONVERT as TABLE_CONVERT } from '../../modules/local/parquet_convert/main' addParams( options: [args: [:]] )
//include { PARQUET_MERGE as TABLE_MERGE } from '../../modules/local/parquet_merge/main' addParams( options: [args: [suffixes:['']]] )
include { TSV_MERGE as TABLE_MERGE } from '../../modules/local/tsv_merge/main' addParams( options: [args: [suffixes:['']]] )

workflow OLIGOS_MAP {
    take:
        FastqChunks // Channel with loaded FastqChunks
        TableChunks // Channel with loaded TableChunks

    main:

        /* Check for the presence of oligos with a custom Rabin-Karp-based algorithm */
        // Load oligos:
        Oligos = Channel
                   .from(params.oligos.keySet().withIndex())
                   .map { it, idx -> [it, idx, params.oligos[it]] } // construct meta for each oligo
                   .map { oligo_name, idx, meta ->
                           meta.single_end = true
                           meta.id = oligo_name
                           meta.idx = idx
                       [meta, file(meta.file)]
                   }

        // Index oligos:
        IndexedOligos = BIN_OLIGOS(Oligos)
        IndexedFastqs = BIN_FASTQ(FastqChunks)

        IndexedLibrary = IndexedFastqs.bin.combine(IndexedOligos.bin).multiMap { it ->
                                                                bin_reads: update_meta( [it[0], it[1]], [oligo: it[2].id, side: it[2].side, idx: it[2].idx] )
                                                                bin_oligos: [it[2], it[3]]
                                                           } // Two branches: bin_reads and bin_oligos
        /* Check oligos presence */
        HitsOligosStream = OLIGOS_ALIGN(IndexedLibrary).aligned //| TABLE_CONVERT

        /* Check presence of GA dinucleotide in the bridge: */
        HitsSmallOligos = OLIGOS_CHECK_GA(
                            join_2_channels(TableChunks, HitsOligosStream.filter { it[0].oligo=='bridge_forward' && it[0].side==1 }, 'id')
                            ).hits

        /* Check complementary RNA ends: */
        HitsComplementary = OLIGOS_CHECK_COMPLEMENTARY( join_4_channels(
                                TableChunks,
                                IndexedLibrary.bin_reads,
                                HitsOligosStream.filter { it[0].oligo=='bridge_forward' && it[0].side==1 },
                                HitsOligosStream.filter { it[0].oligo=='ggg' && it[0].side==2 },
                                'id')
                            ).hits

        HitsOligosGrouped = HitsOligosStream
                            .map{ it -> [ [it[0]['id']], it ] }.groupTuple(by: 0, sort: { a, b -> a[0].idx <=> b[0].idx })
                            .map{ it -> it.collect()[1].collect{ item -> [item[0], file(item[1]), item[0].oligo+'_R'+item[0].side] }.transpose() } //.transpose().collect()
                            .multiMap{meta, files, suffixes ->
                                dataset: [meta[0],  files]
                                suffixes: suffixes
                           }

        HitsOligos = TABLE_MERGE( HitsOligosGrouped ).output

    emit:
    HitsOligos // channel: [ val(meta), [ output_table ] ]
    HitsSmallOligos // channel: [ val(meta), [ output_table ] ]
    HitsComplementary // channel: [ val(meta), [ output_table ] ]
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

def join_4_channels(channel_1, channel_2, channel_3, channel_4, k){
    // Take 4 channels and a meta key and return 4 channels combined by that key:
    channel_1
        .map{ it -> [it[0][k], it] }
        .combine(channel_2.map{ it -> [it[0][k], it] }, by: 0)
        .combine(channel_3.map{ it -> [it[0][k], it] }, by: 0)
        .combine(channel_4.map{ it -> [it[0][k], it] }, by: 0)
        .multiMap { key, ch1, ch2, ch3, ch4 ->
                    ch1: ch1
                    ch2: ch2
                    ch3: ch3
                    ch4: ch4
         }
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