//
// Deduplication with fastuniq
//

params.options = [:]

def dedup_crop = params.getOrDefault('protocol', [:]).getOrDefault('dedup_crop', 50) // Crop 50 bp by default

include { TRIMMOMATIC as FASTQ_CROP } from '../../modules/local/trimmomatic/main' addParams( options: [args: 'CROP:'+dedup_crop, suffix:'.crop', args2: [gzip: false]])
include { FASTUNIQ } from '../../modules/local/fastuniq/main' addParams( options: [:] )

workflow DEDUP_FASTUNIQ {
    take:
        Fastq // Channel with full-size FASTQ files

    main:
        // Deduplication of input sequences:
        FASTQ_CROP( Fastq )
        FASTUNIQ( FASTQ_CROP.out.fastq )
        Nodup = FASTUNIQ.out.nodup_readids

    emit:
        Nodup // channel: [ val(meta), [ nudup_readids ] ]
}