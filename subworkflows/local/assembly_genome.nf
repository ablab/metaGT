include { SPADES } from '../../modules/nf-core/software/spades/main.nf' addParams(spades_hmm: false, options: ['args': '--meta --only-assembler'])
include { FASTQC } from '../../modules/nf-core/software/fastqc/main.nf' addParams(options: [:])


workflow ASSEMBLY_GENOME {
    take:
    reads         // channel: [ val(meta), [ reads ] ]

    main:
    // Cannot do --meta with single-ends
    reads
        .filter { meta, fastq -> !meta.single_end }
        .set { ch_reads }

    FASTQC(ch_reads)
    SPADES(ch_reads, [], false)

    SPADES
        .out
        .scaffolds
        .filter { meta, scaffold -> scaffold.size() > 0 }
        .set { ch_scaffolds }

    emit:
    scaffolds          = SPADES.out.scaffolds               // channel: [ val(meta), [ scaffolds ] ]
    
}
