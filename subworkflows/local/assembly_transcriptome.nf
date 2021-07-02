include { SPADES } from '../../modules/nf-core/software/spades/main.nf' addParams(spades_hmm: false, options: ['args': '--rna'])
include { FASTQC } from '../../modules/nf-core/software/fastqc/main.nf' addParams(options: [:])


workflow ASSEMBLY_TRANSCRIPTOME {
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
        .transcripts
        .filter { meta, transcripts -> transcripts.size() > 0 }
        .set { ch_scaffolds }
    emit:
    transcripts          = SPADES.out.transcripts               // channel: [ val(meta), [ transcripts ] ]

    
}
