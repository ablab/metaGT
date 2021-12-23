// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRANSDECODER {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    //conda (params.enable_conda ? "bioconda::transdecoder" : null)

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.fasta'), emit: cds

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    
    """
    TransDecoder.LongOrfs -t $fasta --output_dir ./ 
    mv ./longest_orfs.cds ${prefix}_cds_from_all_transcripts.fasta
    
    """
}


