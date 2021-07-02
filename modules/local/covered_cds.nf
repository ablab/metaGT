// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COVERED_CDS {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    input:
    tuple val(meta), path(sam)
    tuple val(meta), path (gff)
    tuple val(meta), path (genes)


    output:
    tuple val(meta), path('*.csv'), emit: table
    tuple val(meta), path('*.fasta'), emit: fasta

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """

    extract_covered_cds.py --gff $gff --sam $sam --genes $genes --output ${prefix}_covered_cds
    """
}

