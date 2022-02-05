// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

/* Useful Groovy methods */
// Return True if the file is URL
Boolean isURL(line) {
    return (line.startsWith("http://") || line.startsWith("https://") || line.startsWith("ftp://"))
}
// Return True if the file is gzipped
Boolean isGZ(line) {
    return (line.endsWith(".gz"))
}

process GENOME_PREPARE {
    tag "$assembly"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:assembly, publish_by_meta:'') }

    // Without lenient cache, long list of index files can occasionally produce wrong hash
    cache "lenient"
    conda (params.enable_conda ? "bioconda::hisat2=2.2.1 bioconda::pyfaidx=0.5.9.5 bioconda::tabix=1.11" : null)

    input:
    tuple val(assembly), path(samplesheet)

    output:
    path "${assembly}.fa"            , emit: genome_fasta
    path "${assembly}.chromsizes.txt", emit: genome_chromsizes
    path "${assembly}.fa.*ht2"          , emit: genome_index //".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"
    path  "*.version.txt"  , emit: version

    script:

    def fasta_file = ""
    def suffix = ""
    def getGenomeCmd = ""
    def unpackGenomeCmd = ""
    def chromsizes_file = options.args.genome.getOrDefault("chromsizes", "")
    def getChromsizesCmd = ""
    def getIndexCmd = ""

    def software = getSoftwareName(task.process)
    if (options.args.auto_download_genome) {
        """
        wget http://hgdownload.cse.ucsc.edu/goldenPath/${assembly}/bigZips/${assembly}.fa.gz -O ${assembly}.fa.gz
        bgzip -d -@ ${task.cpus} ${assembly}.fa.gz
        faidx ${assembly}.fa -i chromsizes > ${assembly}.chromsizes.txt
        hisat2-build -p ${task.cpus} ${assembly}.fa ${assembly}.fa

        echo 'hisat 2.2.0' > ${software}.version.txt
        """
    }
    else {
        // Download/copy genomic fasta
        fasta_file = options.args.genome.getOrDefault("fasta", "")
        suffix = isGZ(fasta_file) ? ".gz" : ""

        // Get FASTA from URL or copy
        if (isURL(fasta_file)){
            getGenomeCmd = """
                wget ${fasta_file} -O ${assembly}.fa${suffix}
            """
        }
        else {
            getGenomeCmd = """
                cp ${fasta_file} ${assembly}.fa${suffix}
            """
        }

        // Unpack genomic fasta
        if (isGZ(fasta_file)){
            unpackGenomeCmd = """
                bgzip -d -@ ${task.cpus} ${assembly}.fa.gz
            """
        }

        // Download/copy/build chromosome sizes
        if (isURL(chromsizes_file)){
            getChromsizesCmd = """
                wget ${chromsizes_file} -O ${assembly}.chromsizes.txt
            """
        }
        else if (chromsizes_file.size()>0) {
            getChromsizesCmd = """
                cp ${chromsizes_file} ${assembly}.chromsizes.txt
            """
        }
        else {
            getChromsizesCmd = """
                faidx ${assembly}.fa -i chromsizes > ${assembly}.chromsizes.txt
            """
        }

        // Run index if it does not exist or fasta file was downloaded from URL
        suffix_list = [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
        if ((options.args.genome.getOrDefault('index_prefix', '').size()==0) | (isURL(fasta_file))) {
            getIndexCmd = "hisat2-build -p ${task.cpus} ${assembly}.fa ${assembly}.fa"
        } else {
            for (suf in suffix_list) {
                getIndexCmd = getIndexCmd + "cp ${options.args.genome.index_prefix}${suf} ${assembly}.fa${suf}; "
            }
        }

        """
        ${getGenomeCmd}
        ${unpackGenomeCmd}
        ${getChromsizesCmd}
        ${getIndexCmd}

        echo 'hisat 2.2.1' > ${software}.version.txt
        """
    }
}
