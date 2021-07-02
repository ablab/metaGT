// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PROKKA {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    conda (params.enable_conda ? "bioconda::blast=2.9 bioconda::prokka" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.gff'), emit: genome_gff
    tuple val(meta), path('*.ffn'), emit: ffn

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    """
    [ ! -f  ${prefix}.fasta ] && ln -s $fasta ${prefix}.fasta
    prokka  ${prefix}.fasta --outdir ./ --force --prefix ${prefix} --metagenome 
    
    """
}


