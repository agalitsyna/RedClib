//
// Map oligos and check strict substrings in these oligos for RedC.
//

params.options = [:]

include { BIN_ENCODE as BIN_OLIGOS } from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fasta']] ) // Bin input oligos from fasta file
include { BIN_ENCODE as BIN_FASTQ }  from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fastq']] ) // Bin input fastq

include { OLIGOS_ALIGN } from '../../modules/local/align_encoded/main' addParams( options: [args: [:]] ) // Align oligos
include { OLIGOS_CHECK as OLIGOS_CHECK_GA } from '../../modules/local/oligos_check/main' addParams( options: [args: [oligo: 'GA', position: 35, orientation: 'F']] ) // Check oligos presence
include { OLIGOS_CHECK_COMPLEMENTARY } from '../../modules/local/oligos_complementary/main' addParams( options: [args: [rna_complementary_length: 14]] ) // Align comnplementary fragments of RNA

workflow OLIGOS_MAP {
    take:
        FastqChunks // Channel with loaded FastqChunks
        TableChunks // Channel with loaded TableChunks

    main:

        /* Check for the presence of oligos with a custom Rabin-Karp-based algorithm */
        // Load oligos:
        Oligos = Channel
                   .fromList(params.oligos.keySet())
                   .map { it -> [it, params.oligos[it]] } // construct meta for each oligo
                   .map { oligo_name, meta ->
                           meta.single_end = true
                           meta.id = oligo_name
                       [meta, file(meta.file)]
                   }

        // Index oligos:
        IndexedOligos = BIN_OLIGOS(Oligos)
        IndexedFastqs = BIN_FASTQ(FastqChunks)

        IndexedLibrary = IndexedFastqs.bin.combine(IndexedOligos.bin).multiMap { it ->
                                                                bin_reads: update_meta( [it[0], it[1]], [oligo: it[2].id, side: it[2].side] )
                                                                bin_oligos: [it[2], it[3]]
                                                           } // Two branches: bin_reads and bin_oligos
        /* Check oligos presence */
        aligned_oligos = OLIGOS_ALIGN(IndexedLibrary).aligned

        /* Check presence of GA dinucleotide in the bridge: */
        hits_small_oligos = OLIGOS_CHECK_GA(
                            join_2_channels(TableChunks, aligned_oligos.filter { it[0].oligo=='bridge_forward' && it[0].side==1 }, 'id')
                            ).hits
        hits_small_oligos

        /* Check complementary RNA ends: */
        hits_complementary = OLIGOS_CHECK_COMPLEMENTARY( join_4_channels(
                                TableChunks,
                                IndexedLibrary.bin_reads,
                                aligned_oligos.filter { it[0].oligo=='bridge_forward' && it[0].side==1 },
                                aligned_oligos.filter { it[0].oligo=='ggg' && it[0].side==2 },
                                'id')
                            ).hits

    emit:
    aligned_oligos // channel: [ val(meta), [ aligned ] ]
    hits_small_oligos // channel: [ val(meta), [ hits ] ]
    hits_complementary // channel: [ val(meta), [ hits ] ]

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