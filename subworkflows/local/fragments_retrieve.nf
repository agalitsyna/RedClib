//
// Get RNA/DNA/RNA1 fragments from the reads based on mapped oligos
//

params.options = [:]

include { BIN_ENCODE as BIN_OLIGOS } from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fasta']] ) // Bin input oligos from fasta file

workflow FRAGMENTS_RETRIEVE {
    take:
        TableChunks // Channel with loaded TableChunks
        AlignedOligos // channel: [ val(meta), [ aligned ] ]
        HitsSmallOligos // channel: [ val(meta), [ hits ] ]
        HitsComplementary // channel: [ val(meta), [ hits ] ]

    main:

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
