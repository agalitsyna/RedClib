//
// Trim reads with TRIMMOMATIC and create trimtable
//

params.options = [:]

def trimmomatic_params = params.get('trimmomatic_params', '')

include { TRIMMOMATIC as FASTQ_TRIM } from '../../modules/local/trimmomatic/main' addParams( options: [args: trimmomatic_params, suffix:'.trim', args2: [gzip: false]])
include { TABLE_TRIM } from '../../modules/local/table_trimmomatic' addParams( options: [:] )

workflow TRIMTABLE_CREATE {
    take:
        Fastq // Channel with full-size FASTQ files
        Table // Channel with input table
        
    main:
        // Deduplication of input sequences:
        Trimmed = FASTQ_TRIM(Fastq).fastq
        // TODO: deprecate join_2_channels
        trimtable = TABLE_TRIM( join_2_channels(Table, Trimmed, 'id') ).table

    emit:
        trimtable // channel: [ val(meta), [ nudup_readids ] ]
}

// TODO: deprecate
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
