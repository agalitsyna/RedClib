// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

/* Useful Groovy methods */
// Return True if the file is URL
Boolean isURL(line) {
    return (line.startsWith("http://") || line.startsWith("https://") || line.startsWith("ftp://"))
}
// Return True if the file is gzipped
Boolean isGZ(line) {
    return (line.endsWith(".gz"))
}

params.options = [:]
options        = initOptions(params.options)

process GENOME_PREPARE_RNA_ANNOTATIONS {
    tag "$assembly"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:assembly, publish_by_meta:'') }

//    cache "${params.cache}"
    conda (params.enable_conda ? "anaconda::python=3.7 bioconda::tabix=1.11" : null)

    input:
    val assembly

    output:
    path("*.spliced_genes.txt"), emit: genome_splicesites
    path("*.gtf")              , emit: genome_rna_annotation
    path  "*.version.txt"      , emit: version

    script:
    def software = getSoftwareName(task.process)
    def genes_gtf = options.args.genes_gtf
    def rna_annot_name = assembly
    def rna_annotation_suffix = options.args.getOrDefault('rna_annotation_suffix', '')
    def suffix = isGZ(genes_gtf) ? ".gz" : ""
    def getRNAAnnot = ""
    if (isURL(genes_gtf)) {
        getRNAAnnot = """
        wget ${genes_gtf} -O ${rna_annot_name}${rna_annotation_suffix}.gtf${suffix}
        """
    }
    else {
        getRNAAnnot = """
        cp ${genes_gtf} > ${rna_annot_name}${rna_annotation_suffix}.gtf${suffix}
        """
    }
    def unpackRNAAnnot = ""
    if (isGZ(genes_gtf)){
        unpackRNAAnnot = """
            bgzip -d -@ ${task.cpus} ${rna_annot_name}${rna_annotation_suffix}.gtf.gz
        """
    }

    """
    ${getRNAAnnot}
    ${unpackRNAAnnot}
    hisat2_extract_splice_sites.py ${rna_annot_name}${rna_annotation_suffix}.gtf > ${rna_annot_name}${rna_annotation_suffix}.spliced_genes.txt

    echo 'Python=3.7' > ${software}.version.txt
    """
}
