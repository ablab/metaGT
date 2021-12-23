#!/usr/bin/env nextflow

/*
========================================================================================
    metaGT
========================================================================================
    Github : https://github.com/ablab/metaGT
----------------------------------------------------------------------------------------
 */

nextflow.enable.dsl=2

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run metaGT --dna_reads 'input_dna_reads_{1,2}.fastq.gz' --rna_reads 'input_rna_reads_{1,2}.fastq.gz' "
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
// FIXME: add more params
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    abricate --version > v_abricate.txt
    echo \$(mmseqs 2>&1) > v_varscan.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

include { PROKKA } from './modules/local/prokka' addParams(options: [:])
include { COVERED_CDS } from './modules/local/covered_cds' addParams(options: [:])
include { MMSEQS_CLUSTER } from './modules/local/mmseqs' addParams(options: [:])
include { TRANSDECODER } from './modules/local/transdecoder' addParams(options: [:])


workflow {
    

    if (params.rna_reads) {
        include { ASSEMBLY_TRANSCRIPTOME } from './subworkflows/local/assembly_transcriptome.nf' addParams(options : [:])
        input_rna_reads =
            Channel.fromFilePairs(params.rna_reads,size: -1)
            .map{ item -> [ [id : item[0], yaml : false], item[1] ] }
    }
    if (params.transcriptome) {
        ch_transcriptome = Channel.fromPath(params.transcriptome, checkIfExists: true)
            .map{ item -> [ [id : file(item).getBaseName(), single_end : false], item ] }
    } else if (params.rna_reads) {
        
        ASSEMBLY_TRANSCRIPTOME(input_rna_reads)
        ch_transcriptome = ASSEMBLY_TRANSCRIPTOME.out.transcripts
    }
    else {
        exit 1, "ERROR: Please check input RNA reads and transcriptome do not exist"
    }


    if (params.dna_reads) {
        include { ASSEMBLY_GENOME } from './subworkflows/local/assembly_genome.nf' addParams(options : [:])
        input_reads =
            Channel.fromFilePairs(params.dna_reads,size: -1)
            .map{ item -> [ [id : item[0], yaml : false], item[1] ] }
        ASSEMBLY_GENOME(input_reads)
        ch_genome = ASSEMBLY_GENOME.out.scaffolds
    } 
    
    else if (params.genome) {
        ch_genome = Channel.fromPath(params.genome, checkIfExists: true)
            .map{ item -> [ [id : file(item).getBaseName(), single_end : false], item ] }
    }
    else {
        exit 1, "ERROR: Please check input DNA reads and genome do not exist"
    }

    if (params.sam) {
        ch_align = Channel.fromPath(params.sam, checkIfExists: true)
            .map{ item -> [ [id : file(item).getBaseName(), single_end : false], item ] }
    } else { 
        include { MINIMAP2 } from './modules/local/minimap2' addParams(options: [:])
        MINIMAP2(ch_genome, ch_transcriptome)
        ch_align = MINIMAP2.out.bam
    }
    
    if (params.gff) {
        ch_genome_gff = Channel.fromPath(params.gff, checkIfExists: true)
           .map{ item -> [ [id : file(item).getBaseName(), single_end : false], item ] }
    } else {
        PROKKA(ch_genome) 
        ch_genome_gff = PROKKA.out.genome_gff
           }

    
    
    COVERED_CDS(ch_align, ch_genome_gff, ch_genome)
    TRANSDECODER(ch_transcriptome)
    MMSEQS_CLUSTER(COVERED_CDS.out.fasta, TRANSDECODER.out.cds) 
    
    if (params.rna_reads) {
        include { KALLISTO_INDEX; KALLISTO_QUANT } from './modules/local/kallisto' addParams(options: [:])
        KALLISTO_INDEX(MMSEQS_CLUSTER.out.rep_seq)
        KALLISTO_QUANT(KALLISTO_INDEX.out.index, input_rna_reads)

    }
}

workflow.onComplete {
    // Set up the e-mail variables
    def subject = "[metaGT] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[metaGT] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir" ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[metaGT] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            mail_cmd.execute() << email_html
            log.info "[metaGT] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[metaGT]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[metaGT]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}
