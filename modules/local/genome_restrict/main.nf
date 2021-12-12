// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GENOME_RESTRICT {
    tag "$meta.assembly $meta.renzyme"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:'') }

//    conda (params.enable_conda ? "anaconda::python=3.7 conda-forge::numpy=1.21.0 conda-forge::pandas=1.2.5 conda-forge::biopython=1.79 bioconda::coreutils=8.25" : null)
    conda (params.enable_conda ? "conda-forge::numpy=1.21.0 conda-forge::pandas=1.2.5 conda-forge::biopython=1.79 bioconda::coreutils=8.25" : null)

    input:
    val(meta)
    file(genome_fasta)

    output:
    tuple val(meta), path("${meta.assembly}.${meta.renzyme}.bed"), emit: genome_restricted
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def assembly = meta['assembly']
    def renzyme = meta['renzyme']
    """
    detect_restriction_sites.py ${genome_fasta} ${renzyme} ${assembly}.${renzyme}.nonsorted.bed
    sort -k1,1 -k2,2n --parallel=${task.cpus} ${assembly}.${renzyme}.nonsorted.bed > ${assembly}.${renzyme}.bed
    rm ${assembly}.${renzyme}.nonsorted.bed

    python -c "import Bio; print('BioPython', Bio.__version__)" > ${software}.version.txt
    """
    }
