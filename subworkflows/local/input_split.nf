//
// Split fastq into chunks and return channel with additional meta
//

params.options = [:]

def chunksize = params.options.get('chunksize', 100000000000)*4

include { FASTQ_SPLIT } from '../../modules/local/fastq_split/main' addParams( options: [args: [chunksize : chunksize]], outdir: 'fastq_chunks' )

/* Useful Groovy method */
// Parse forward (left) and reverse (right) pieces of the file,
// assert that they correspond each other and return index
String parseChunkPair(left_file, right_file) {
    chunk_left_idx  =  left_file.toString().tokenize('.')[-3]
    chunk_right_idx = right_file.toString().tokenize('.')[-3]
    chunk_left_side  =  left_file.toString().tokenize('.')[-2]
    chunk_right_side = right_file.toString().tokenize('.')[-2]
    assert chunk_left_idx == chunk_right_idx
    assert chunk_left_side != chunk_right_side
    return chunk_left_idx
}
String parseChunkSingle(file) {
    chunk_idx  =  file.toString().tokenize('.')[-3]
    return chunk_idx
}

workflow INPUT_SPLIT {
    take:
        fastq

    main:

        FASTQ_SPLIT (
            fastq
        )

        FASTQ_SPLIT.out.fastqs.transpose()
            .map{ update_meta(it) }
            .set { output }

    emit:
    output // channel: [ val(meta_updated), [ reads ] ]
}

def update_meta( it ) {
    def meta = [:]
    keys = it[0].keySet() as String[]
    for( key in keys ) {
        meta[key] = it[0][key]
    }
    def file1 = it[1]
    def file2 = ""

    if (meta.single_end) {
        chunk_index = parseChunkSingle(file1)
    } else {
        file2 = it[2]
        chunk_index = parseChunkPair(file1, file2)
    }

    meta.original_id = meta.id
    meta.id          = meta.id + "_" + chunk_index
    meta.chunk       = chunk_index

   if (meta.single_end) {
        array = [ meta, [ file(file1) ] ]
    } else {
        array = [ meta, [ file(file1), file(file2) ] ]
    }
    return array
}