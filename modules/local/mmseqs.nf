// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MMSEQS_CLUSTER {
    tag "$meta.id"

    input:

    tuple val(meta), path(cov_transcripts)
    tuple val(meta_t), path(cds_from_unaligned)

    output:
    path '*.rep_seq.fasta', emit: rep_seq

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options,  publish_dir:'rep_seq', publish_id:meta.id) }
    
 //  conda (params.enable_conda ? "bioconda::mmseqs2" : null)
    
    


    script:
    def prefix  = params.prefix ? "${params.prefix}" : "${meta.id}"
    """
    cat $cov_transcripts $cds_from_unaligned > all.fasta
    mmseqs easy-linclust $cov_transcripts res tmp --min-seq-id ${params.cluster_idy}
    mv res_rep_seq.fasta ${prefix}.rep_seq.fasta
    """
}
