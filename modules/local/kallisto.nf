// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

//conda (params.enable_conda ? "bioconda::kallisto" : null)

process KALLISTO_INDEX {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process)) }

    
input:
    path fasta

output:
    path 'index', emit: index

script:
    //def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """

    kallisto index -i index $fasta 
    """
}

process KALLISTO_QUANT {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process),) }

input:
    path index
    tuple val(meta), path(reads)

output:
    path 'abudance.tsv', emit: report

script:
    //def prefix  = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = params.yaml ? "`parse_yaml.py $reads`" : "${reads[0]} ${reads[1]}"
    """

    kallisto quant -i $index $input_reads -t $task.cpus -o ./
    cp ./abundance.tsv abudance.tsv   
    """
}
