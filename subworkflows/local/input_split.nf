//
// Split fastq into chunks and return channel with additional meta
//

params.options = [:]

def chunksize = params.options.getOrDefault('chunksize', 100000000000)*4

include { FASTQ_SPLIT } from '../../modules/local/fastq_split/main' addParams( options: [args: [chunksize : chunksize]], outdir: "${params.outdir}" )

workflow INPUT_SPLIT {
    take:
        fastq

    main:

        chunks = FASTQ_SPLIT (
            fastq
        ).fastq

        fastq_chunks = chunks.transpose()
            .map{ update_meta(it) }

    emit:
    fastq_chunks // channel: [ val(meta_updated), [ reads ] ]
}

/* Useful Groovy methods */

// Parse forward (left) and reverse (right) pieces of the file,
// assert that they correspond each other and return index
String parseChunkPair(left_file, right_file) {
    def chunk_left_idx  =  left_file.toString().tokenize('.')[-3]
    def chunk_right_idx = right_file.toString().tokenize('.')[-3]
    def chunk_left_side  =  left_file.toString().tokenize('.')[-2]
    def chunk_right_side = right_file.toString().tokenize('.')[-2]
    assert chunk_left_idx == chunk_right_idx
    assert chunk_left_side != chunk_right_side
    return chunk_left_idx
}
String parseChunkSingle(file) {
    def chunk_idx  =  file.toString().tokenize('.')[-3]
    return chunk_idx
}

// Update metadata:
def update_meta( it ) {
    def meta = [:]
    def keys = it[0].keySet() as String[]
    for( def key in keys ) {
        meta[key] = it[0][key]
    }
    def file1 = it[1]
    def file2 = ""
    def chunk_index = ""

    if (meta.single_end) {
        chunk_index = parseChunkSingle(file1)
    } else {
        file2 = it[2]
        chunk_index = parseChunkPair(file1, file2)
    }

    meta.original_id = meta.id
    meta.id          = meta.id + "_" + chunk_index
    meta.chunk       = chunk_index

    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(file1) ] ]
    } else {
        array = [ meta, [ file(file1), file(file2) ] ]
    }
    return array
}