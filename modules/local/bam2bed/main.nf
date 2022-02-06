// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAM2BED {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    cache "${cacheMode}"
    conda (params.enable_conda ? "bioconda::bedtools=2.30 bioconda::samtools=1.14" : null)

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def extraTags = options.args.getOrDefault('extraTags', [])
    def filter = options.args.getOrDefault('samFilter', '')
    def tmpFiles = ["tmp.txt"]

    Cmd = "samtools view -h ${bam}"
    if (filter) {Cmd += " | " + filter + " | "}

    if (extraTags.size()>0) {
        Cmd += " | tee >(bedtools bamtobed -i - -cigar > tmp.txt) "
        extraTags.eachWithIndex{ tag, i ->
            Cmd += ">(bedtools bamtobed -i - -tag ${tag} | awk '{print \$5}' > tmp${i}.txt) "
            tmpFiles += "tmp${i}.txt"
        }
    }
    Cmd += " | cat "
    """

    ${Cmd}

    paste ${tmpFiles.join(' ')} > ${prefix}.bed

    bedtools --version &> ${software}.version.txt
    samtools --version | head -n 2 > ${software}.version.txt
    """
}

//    def cmd = ""
//    def filter1 = params.options.getOrDefault("filter", "")
//    def filter2 = params.options.getOrDefault("filter2", "")
//    def filter3 = params.options.getOrDefault("filter3", "")
//
//    if (filter1) {
//        cmd += "${filter1} ${bam} "
//    }
//    if (filter2) {
//        cmd += " | ${filter2} "
//    }
//    if (filter3) {
//        cmd += " | ${filter3} "
//    }
//
//    cmd += " | bedtools bamtobed -cigar -i stdin "

//     // Filter unmapped, retain reads with up to 2 substitutions:
//     filter: 'samtools view -h -F 4 -d XM:0 -d XM:1 -d XM:2',
//     // Filter uniquely mapped:
//     filter2: 'samtools view -h -d NH:1 -'
