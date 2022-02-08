// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process RKLIB_SEQ2HASH {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
//    cache "${params.cache}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bin"), emit: bin
    path  "*.version.txt"         , emit: version

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def mode   = options.args.getOrDefault('mode', 'fasta')

    // If the read length is not the same in the file, the hashing will work, but the rklib will fail to produce correct answers,
    // because it assumes equal read lengths. Thus we may force the reads to equal length (filling with N):
    def rlen = options.args.getOrDefault('force_rlen', 0)
    def extraN = "N"*(rlen+1)

    def Cmd = ""
    def software = ""

    if (mode=='fasta') {
        software = "fasta2bin"
        if (rlen) { // Extend reads to equal length and then convert:
            if (meta.single_end) {
                Cmd += "awk '{if (NR%2==0) print substr(\$0\"${extraN}\", 1, ${rlen}); else print \$0}' ${reads[0]} > tmp_r1.fa\n"
                Cmd += "fasta2hash tmp_r1.fa ${prefix}.bin\n"
            } else {
                Cmd += "awk '{if (NR%2==0) print substr(\$0\"${extraN}\", 1, ${rlen}); else print \$0}' ${reads[0]} > tmp_r1.fa\n"
                Cmd += "awk '{if (NR%2==0) print substr(\$0\"${extraN}\", 1, ${rlen}); else print \$0}' ${reads[1]} > tmp_r2.fa\n"
                Cmd += "fasta2hash tmp_r1.fa ${prefix}.R1.bin\n"
                Cmd += "fasta2hash tmp_r2.fa ${prefix}.R2.bin\n"
            }
        } else { // No need to extend reads to equal length:
            if (meta.single_end) {
                Cmd += "fasta2hash ${reads[0]} ${prefix}.bin\n"
            } else {
                Cmd += "fasta2hash ${reads[0]} ${prefix}.R1.bin\n"
                Cmd += "fasta2hash ${reads[1]} ${prefix}.R2.bin\n"
            }
        }
    } else {
        software = "fastq2bin"
        if (rlen) { // Extend reads to equal length and then convert:
            if (meta.single_end) {
                Cmd += "awk '{if (NR%2==0) print substr(\$0\"${extraN}\", 1, ${rlen}); else print \$0}' ${reads[0]} > tmp_r1.fq\n"
                Cmd += "fastq2hash tmp_r1.fq ${prefix}.bin\n"
            } else {
                Cmd += "awk '{if (NR%2==0) print substr(\$0\"${extraN}\", 1, ${rlen}); else print \$0}' ${reads[0]} > tmp_r1.fq\n"
                Cmd += "awk '{if (NR%2==0) print substr(\$0\"${extraN}\", 1, ${rlen}); else print \$0}' ${reads[1]} > tmp_r2.fq\n"
                Cmd += "fastq2hash tmp_r1.fq ${prefix}.R1.bin\n"
                Cmd += "fastq2hash tmp_r2.fq ${prefix}.R2.bin\n"
            }
        } else { // No need to extend reads to equal length:
            if (meta.single_end) {
                Cmd += "fastq2hash ${reads[0]} ${prefix}.bin\n"
            } else {
                Cmd += "fastq2hash ${reads[0]} ${prefix}.R1.bin\n"
                Cmd += "fastq2hash ${reads[1]} ${prefix}.R2.bin\n"
            }
        }
    }

    """

    ${Cmd}

    echo $VERSION > ${software}.version.txt

    """
}
