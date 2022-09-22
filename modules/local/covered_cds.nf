// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COVERED_CDS {
    tag "$meta.id"
    time { 40.hour }

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    input:
    tuple val(meta), path(bam)
    tuple val(meta), path (gff)
    tuple val(meta), path (genome)


    output:
    tuple val(meta), path('*.csv'), emit: table
    tuple val(meta), path('*.fasta'), emit: fasta

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def treshhold = params.threshold ? "-t ${params.threshold}" : ""

    """
    samtools index $bam

    extract_covered_cds.py --gff $gff --bam $bam --genome $genome --output ${prefix}_covered_cds 
    """
}

