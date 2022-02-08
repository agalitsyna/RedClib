//
// Check input samplesheet and get read channels
//

params.options = [:]

include { FASTQ_DOWNLOAD as FASTQDOWNLOAD } from '../../modules/local/fastq_download/main' addParams( options: [:] )

/* Useful Groovy method */
// Return True if the file is sra
Boolean isSRA(line) {
    return (line.startsWith("sra:"))
}

workflow INPUT_CHECK_DOWNLOAD {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:

        reads = Channel
            .fromPath( samplesheet )
            .splitCsv ( header:true, sep:',' )
            .map { format_fastq_channels(it) }
            .branch {
                meta, fastq_path ->
                    for_download  : isSRA(fastq_path.flatten()[0])
                        return [ meta, fastq_path.flatten() ]
                    ready : !isSRA(fastq_path.flatten()[0])
                        return [ meta, fastq_path.flatten() ]
            }

        fastq_downloaded = FASTQDOWNLOAD(
                reads.for_download
            ).fastq

        fastq = fastq_downloaded
            .mix(reads.ready)
            .map { create_fastq_channels(it) }

    emit:
        fastq // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def format_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.group        = row.group
    meta.single_end   = row.single_end.toBoolean()
    meta.rlen         = row.rlen
    //meta.strandedness = row.strandedness

    def array = []
    if (meta.single_end) {
        array = [ meta, [ row.fastq_1 ] ]
    } else {
        array = [ meta, [ row.fastq_1, row.fastq_2 ] ]
    }
    return array
}

def create_fastq_channels( it ) {
    def meta = it[0]
    def files = it[1]

    def array = []
    if (meta.single_end) {
        if (!file(files).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${files[0]}"
        }
        array = [ meta, [ file(files) ] ]
    } else {
        if (!file(files[0]).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${files[0]}"
        }
        if (!file(files[1]).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${files[1]}"
        }
        array = [ meta, [ file(files[0]), file(files[1]) ] ]
    }
    return array
}