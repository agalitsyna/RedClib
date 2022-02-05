//
// Trim reads with TRIMMOMATIC and create trimtable
//

params.options = [:]

def trimmomatic_params = params.getOrDefault('protocol', [:]).getOrDefault('trimmomatic_params', '')

include { TRIMMOMATIC as FASTQ_TRIM } from '../../modules/local/trimmomatic/main' addParams( options: [args: trimmomatic_params, suffix:'.trim', args2: [gzip: false]])
include { TABLE_TRIM } from '../../modules/local/table_trimmomatic' addParams( options: [:] )

workflow TRIMTABLE_CREATE {
    take:
        Fastq // Channel with full-size FASTQ files
        Table // Channel with input table
        
    main:
        // Deduplication of input sequences:
        Trimmed = FASTQ_TRIM(Fastq).fastq
        TrimmingInput = Table.combine(Trimmed)
                             .filter{ it[0].id==it[2].id }
                             .multiMap{meta, table, meta_trim, trimout ->
                                input: [meta, table]
                                input_trim: [meta_trim, trimout]
                             }
        trimtable = TABLE_TRIM( TrimmingInput ).table

    emit:
        trimtable // channel: [ val(meta), [ nudup_readids ] ]
}

