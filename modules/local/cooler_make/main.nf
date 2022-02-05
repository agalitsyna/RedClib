// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = '0.0'

process COOLER_MAKE {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

//    cache "${params.cache}"
    conda (params.enable_conda ? "bioconda::cooler=0.8.11" : null)

    input:
    tuple val(meta), path(table)
    path(chromsizes)

    output:
    tuple val(meta), path("*.cool"), emit: cool
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    // Define the parameters. Pick from meta first, then check process options:
    def assembly = meta.getOrDefault("assembly", options.getOrDefault("assembly", "assembly"))
    def resolution = meta.getOrDefault("resolution", options.getOrDefault("resolution", 1000000))
    // By default, we assume regular pairs file:
    def c1 = meta.getOrDefault("c1", options.getOrDefault("c1", 1))
    def c2 = meta.getOrDefault("c2", options.getOrDefault("c2", 2))
    def p1 = meta.getOrDefault("p1", options.getOrDefault("p1", 3))
    def p2 = meta.getOrDefault("p2", options.getOrDefault("p2", 4))

    """
    cooler makebins ${chromsizes} ${resolution} > ${assembly}.${resolution}.bins.txt
    cooler cload pairs ${options.args} \
      -c1 ${c1} -c2 ${c2} \
      -p1 ${p1} -p2 ${p2} \
      ${assembly}.${resolution}.bins.txt <(tail -n +2 ${table}) "${prefix}.${assembly}.${resolution}.cool"
    """
}
