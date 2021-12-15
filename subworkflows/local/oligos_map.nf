//
// Map oligos and check strict substrings in these oligos for RedC.
//

params.options = [:]

include { BIN_ENCODE as BIN_OLIGOS } from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fasta']] ) // Bin input oligos from fasta file
include { BIN_ENCODE as BIN_FASTQ }  from '../../modules/local/bin_encode/main' addParams( options: [args: [mode: 'fastq']] ) // Bin input fastq

include { OLIGOS_ALIGN } from '../../modules/local/align_encoded/main' addParams( options: [args: [:]] ) // Align oligos
include { OLIGOS_CHECK as OLIGOS_CHECK_GA } from '../../modules/local/oligos_check/main' addParams( options: [args: [oligo: 'GA', position: 35, orientation: 'F']] ) // Check oligos presence

workflow OLIGOS_MAP {
    take:
        FastqChunks // Channel with loaded FastqChunks
        TableChunks // Channel with loaded TableChunks

    main:

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

        MappedOligos =  IndexFastqs.bin.combine(IndexOligos.bin).multiMap { it ->
                                                                                reads: update_meta( [it[0], it[1]], [oligo: it[2].id, side: it[2].side] )
                                                                                oligos: [it[2], it[3]]
                                                                           } | OLIGOS_ALIGN

        aligned = MappedOligos.aligned

        /* Check presence of GA dinucleotide in the bridge: */
         hits = OLIGOS_CHECK_GA(TableChunks, MappedOligos.aligned.filter { it[0].oligo=='bridge_forward' } ).hits

    emit:
    aligned
    hits // channel: [ val(meta), [ hits ] ]
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