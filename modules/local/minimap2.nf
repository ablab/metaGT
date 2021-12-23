// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MINIMAP2 {
    tag "$meta.id"

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

   // conda (params.enable_conda ? "bioconda::minimap2" : null)

    input:
    tuple val(meta), path(genome)
    tuple val(meta_t), path (transcriptome)


    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta_t), path('*.all_transcripts.fasta'), emit: transcripts

    script:
    def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """

    minimap2 -aY --MD $genome $transcriptome > ${prefix}.align.sam

    samtools sort ${prefix}.align.sam -o ${prefix}.align.sorted.bam

    change_name.py  $transcriptome ${meta_t.id}.all_transcripts.fasta

    """
  //  TransDecoder.LongOrfs -t ${meta_t.id}.unaligned_transcripts.fasta --output_dir ./ 
// mv ./longest_orfs.cds ${meta_t.id}_cds_from_unaligned_transcripts.fasta
}


