//
// Subworkflow with functionality specific to the IntGenomicsLab/lrsomatic pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    // Parse the input samplesheet CSV and build a per-sample BAM channel
    // Each samplesheet row describes one tumor (+ optional normal) sample
    // Columns: sample_id, bam_tumor, bam_normal, method, sex, fiber,
    //          clair3_model, clairSTO_model, clairS_model, tumor_replicate, normal_replicate
    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        // Step 1: build a combined meta map from the samplesheet columns
        // paired_data = true if a normal BAM is present; false for tumor-only
        .map { meta, bam_tumor, bam_normal, method, sex, fiber, clair3_model, clairSTO_model, clairS_model, tumor_replicate, normal_replicate ->
            def real_clair3_model = (clair3_model == null ) ? null : clair3_model
            def real_clairS_model = (clairS_model == null ) ? null : clairS_model
            def real_clairSTO_model = (clairSTO_model == null ) ? null : clairSTO_model
            def paired_data = bam_normal ? true : false
            def meta_info = meta + [ paired_data: paired_data,
                                     platform: method,         // 'ont' or 'pb'
                                     sex: sex,                 // 'XX', 'XY', or null (for ASCAT)
                                     fiber: fiber,             // 'y' or 'n' (fiber-seq data flag)
                                     clair3_model: real_clair3_model,
                                     clairS_model: real_clairS_model,
                                     clairSTO_model: real_clairSTO_model,
                                     tumor_replicate: tumor_replicate,
                                     normal_replicate: normal_replicate]
            return [ meta_info, [ bam_tumor ], [ bam_normal ?: [] ] ]
        }
        // Flatten BAM lists (handles multi-run entries where bam_tumor/bam_normal are lists)
        .map { meta, bam_tumor, bam_normal ->
           [ meta, bam_tumor.flatten(), bam_normal.flatten() ]
        }
        // Step 2: split each row into separate tumor and normal items
        // flatMap emits 1 item (tumor-only) or 2 items (tumor + normal) per samplesheet row
        // Each item gets type='tumor' or type='normal' and the appropriate replicate ID
        .flatMap { meta, tumor_bam, normal_bam ->
            def meta_tumor = meta.clone()
            meta_tumor.type = 'tumor'
            meta_tumor.replicate = meta_tumor.tumor_replicate
            meta_tumor = meta_tumor.subMap('id',
                                           'paired_data',
                                           'type',
                                           'platform',
                                           'sex',
                                           'fiber',
                                           'clair3_model',
                                           'clairS_model',
                                           'clairSTO_model',
                                           'replicate')
            def result = [[meta_tumor, tumor_bam]]
            // result so far: [[meta_tumor, [tumor_bam_path...]]]

            if (normal_bam) {
                def meta_normal = meta.clone()
                meta_normal.type = 'normal'
                meta_normal.replicate = meta_normal.normal_replicate
                meta_normal = meta_normal.subMap('id',
                                                 'paired_data',
                                                 'type',
                                                 'platform',
                                                 'sex',
                                                 'fiber',
                                                 'clair3_model',
                                                 'clairS_model',
                                                 'clairSTO_model',
                                                 'replicate')
                result << [meta_normal, normal_bam]
                // result now: [[meta_tumor, [tumor_bams]], [meta_normal, [normal_bams]]]
            }

            return result
        }
        .set { ch_samplesheet }

    // Count replicates per sample+type and embed the count in meta as n_replicates.
    // This allows downstream groupTuple() to use groupKey() for eager per-sample release
    // instead of waiting for ALL samples to finish (global synchronization barrier).
    // This groupTuple is safe: the source is a fully-materialized list from
    // samplesheetToList(), so the channel closes immediately without blocking any process.
    ch_samplesheet
        .map { meta, bams -> [[meta.id, meta.type], meta, bams] }
        .groupTuple(by: 0)
        .flatMap { key, metas, bams_list ->
            def n = metas.size()
            [metas, bams_list].transpose().collect { m, b ->
                [m + [n_replicates: n], b]
            }
        }
        .set { ch_samplesheet }

    // ch_samplesheet: [meta, [bam...]]
    //   meta fields: id, paired_data, type ('tumor'|'normal'), platform ('ont'|'pb'),
    //                sex, fiber ('y'|'n'), clair3_model, clairS_model, clairSTO_model,
    //                replicate, n_replicates
    //   paired_data: true for both items in a T/N pair (same value for tumor AND normal rows)
    //   n_replicates: total number of replicates for this sample+type combination
    //   bam: list of paths (multiple runs for same sample remain as a list until SAMTOOLS_CAT)
    //
    // NOTE: tumor-only rows emit ONE item (type='tumor', paired_data=false)
    //       paired rows emit TWO items — tumor (paired_data=true) + normal (paired_data=true)
    //       Both share the same 'id' to allow downstream joins

    emit:
    samplesheet = ch_samplesheet  // [meta, [bam...]]  -- see channel structure above
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, bams) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], bams ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Resolve a VEP plugin resource: an explicit --vep_* wins, otherwise the
// per-assembly default from conf/igenomes.config. Read where it is used rather
// than assigned back onto params, since these params are declared in
// nextflow.config and a runtime assignment to a declared param is dropped --
// the same reason vep_plugin_args is not declared there at all.
//
def vepPluginResource(name) {
    return params[name] ?: getGenomeAttribute(name)
}

//
// True when no plugin annotation should happen at all.
//
def vepPluginsSkipped() {
    return params.skip_vep || params.skip_vep_plugins
}

//
// file(), checking existence for local paths only. A remote release is checked
// when Nextflow stages it, and some endpoints reject the HEAD request that
// checkIfExists relies on -- EVE's bulk download answers 405 to HEAD and 200
// to GET.
//
def vepPluginFile(value) {
    return value.toString().contains('://') ? file(value) : file(value, checkIfExists: true)
}

//
// VEP plugin data files, keyed by the param that supplies them. The value is
// the param holding the tabix index, or null when the resource needs none.
//
def vepPluginIndexParams() {
    return [
        'vep_alphamissense'   : 'vep_alphamissense_tbi',
        'vep_alphamissense_aa': 'vep_alphamissense_aa_tbi',
        'vep_polyphen_sift_db': null,
        'vep_clinvar'         : 'vep_clinvar_tbi',
        'vep_cadd_snv'        : 'vep_cadd_snv_tbi',
        'vep_cadd_indel'      : 'vep_cadd_indel_tbi',
        'vep_revel'           : 'vep_revel_tbi',
        'vep_eve'             : 'vep_eve_tbi'
    ]
}

//
// Whether a resource still has to be reshaped before VEP can read it. Three of
// them cannot be used as published: AlphaMissense ships without a tabix index,
// and REVEL and EVE ship as zip archives. The pipeline reshapes those itself in
// a prep task, so the shape of what was supplied is what decides -- a data file
// with no index, or a .zip, is the raw release; anything else is already
// prepared and is used as it stands.
//
def vepPluginNeedsPrep(data_param) {
    def value = vepPluginResource(data_param)
    if (!value) {
        return false
    }
    if (['vep_revel', 'vep_eve'].contains(data_param)) {
        return value.toString().toLowerCase().endsWith('.zip')
    }
    if (['vep_alphamissense', 'vep_alphamissense_aa'].contains(data_param)) {
        return !vepPluginResource(vepPluginIndexParams()[data_param])
    }
    return false
}

//
// The filename a prep task writes. Also what the VEP argument has to reference,
// since the module stages every plugin file into the task workdir root.
// Deterministic, so the argument string can still be assembled up front, before
// any prep task has run.
//
def vepPluginPreparedName(data_param) {
    // AlphaMissense is indexed where it lies, so it keeps the name it was
    // published under. The rest are written under a fixed name by their task.
    if (data_param == 'vep_alphamissense') {
        return vepPluginFile(vepPluginResource(data_param)).name
    }
    return [
        'vep_alphamissense_aa': 'alphamissense_protein.tsv.gz',
        'vep_revel'           : 'revel_grch38.tsv.gz',
        'vep_eve'             : 'eve_merged.vcf.gz'
    ][data_param]
}

//
// Exit pipeline if VEP plugin params are inconsistent with each other or with
// the target assembly. Checked up front so a run does not fail hours later
// inside VEP, or after an 80 GB download.
//
def validateVepPluginParams() {
    if (vepPluginsSkipped()) {
        return
    }

    def errors = []

    // An already-prepared resource needs its index supplied explicitly: the
    // index may not sit next to the data file when the data file is a URL. The
    // resources the pipeline can prepare itself are exempt, since for them a
    // missing index is the signal to prepare one.
    def index_advice = [
        'vep_clinvar'   : 'ClinVar publishes a .tbi alongside every VCF.',
        'vep_cadd_snv'  : 'CADD publishes a .tbi alongside every score file.',
        'vep_cadd_indel': 'CADD publishes a .tbi alongside every score file.',
        'vep_revel'     : 'Either supply the index, or pass the published revel-v1.3_all_chromosomes.zip and the pipeline will prepare both.',
        'vep_eve'       : 'Either supply the index, or pass the published EVE_all_data.zip and the pipeline will prepare both.'
    ]

    vepPluginIndexParams().each { data_param, index_param ->
        if (!index_param || !vepPluginResource(data_param) || vepPluginNeedsPrep(data_param)) {
            return
        }
        if (!vepPluginResource(index_param)) {
            errors << "  --${data_param} is set but --${index_param} is not. ${index_advice[data_param] ?: 'Both are required.'}"
        }
    }

    // Resources published only in GRCh37/GRCh38 coordinates cannot be used
    // against the T2T-CHM13 cache.
    def grch38_only = [
        'vep_alphamissense': 'Use --vep_alphamissense_aa instead, which is keyed in protein space.',
        'vep_cadd_snv'     : 'CADD scores non-coding positions and has no protein-space form, so it is unavailable on CHM13.',
        'vep_cadd_indel'   : 'CADD scores non-coding positions and has no protein-space form, so it is unavailable on CHM13.',
        'vep_revel'        : 'REVEL is published for GRCh37/GRCh38 only.',
        'vep_eve'          : 'EVE is published for GRCh38 only.'
    ]

    if (params.vep_genome == 'T2T-CHM13v2.0') {
        grch38_only.each { data_param, advice ->
            if (vepPluginResource(data_param)) {
                errors << "  --${data_param} is a GRCh38-only resource and cannot be used with --vep_genome T2T-CHM13v2.0. ${advice}"
            }
        }
    }
    else {
        if (vepPluginResource('vep_alphamissense_aa')) {
            errors << "  --vep_alphamissense_aa is the CHM13 route to AlphaMissense. On ${params.vep_genome} use --vep_alphamissense instead."
        }
    }

    if (errors) {
        error(
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Invalid VEP plugin configuration:\n" +
            errors.join("\n") + "\n" +
            "  See the VEP plugins section of docs/usage.md.\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        )
    }
}

//
// Say what a run is about to fetch, and under what terms. The plugin data is on
// by default, so both facts have to be visible without reading the docs first.
//
// Only resources still resolving to a URL are reported: a local path means the
// user has already downloaded it and nothing will be fetched.
//
def warnVepPluginDownloads() {
    if (vepPluginsSkipped()) {
        return
    }

    // Approximate published sizes, for the resources big enough that this is
    // worth knowing before the download starts rather than after.
    def sizes = [
        'vep_cadd_snv'        : '81 GB',
        'vep_polyphen_sift_db': '13 GB',
        'vep_eve'             : '9.6 GB',
        'vep_cadd_indel'      : '1.2 GB',
        'vep_alphamissense_aa': '1.1 GB',
        'vep_revel'           : '667 MB',
        'vep_alphamissense'   : '613 MB'
    ]

    def fetching = sizes.findAll { data_param, _size ->
        vepPluginResource(data_param)?.toString()?.contains('://')
    }

    if (fetching) {
        log.warn(
            "VEP plugin data will be downloaded: ${fetching.collect { data_param, size -> "${data_param} (${size})" }.join(', ')}. " +
            "Nothing prepared is published, so a fresh work directory downloads again -- pass a local path to the matching " +
            "--vep_* param to reuse a copy you already have. See the VEP plugins section of docs/usage.md."
        )
    }

    // Non-commercial resources are enabled by default, so say so rather than
    // leaving the user to discover it in a licence file.
    def non_commercial = [
        'vep_cadd_snv'  : 'CADD',
        'vep_cadd_indel': 'CADD',
        'vep_revel'     : 'REVEL',
        'vep_eve'       : 'EVE'
    ]

    def in_use = non_commercial
        .findAll { data_param, _tool -> vepPluginResource(data_param) }
        .values()
        .toList()
        .unique()

    if (in_use) {
        log.warn(
            "VEP annotation includes ${in_use.join(', ')}, which are free for non-commercial use only. " +
            "Observing those terms is your responsibility. --skip_vep_plugins turns the plugin annotation off."
        )
    }
}

//
// Register an already-prepared VEP plugin data file, and its index when there
// is one, on the `staged` accumulator. Returns the basename that VEP should
// reference, since the module stages these into the task workdir root.
//
def stageVepPluginFile(staged, data_param, index_param) {
    def data_file = vepPluginFile(vepPluginResource(data_param))
    staged << data_file
    if (index_param && vepPluginResource(index_param)) {
        staged << vepPluginFile(vepPluginResource(index_param))
    }
    return data_file.name
}

//
// Register one resource on the accumulators, returning the basename VEP should
// reference. Either the resource is usable as supplied -- staged under its own
// name -- or it is a raw release, in which case a prep task is recorded and the
// name that task will write is returned instead.
//
def registerVepPlugin(staged, prepare, data_param) {
    if (vepPluginNeedsPrep(data_param)) {
        prepare[data_param] = vepPluginFile(vepPluginResource(data_param))
        return vepPluginPreparedName(data_param)
    }
    return stageVepPluginFile(staged, data_param, vepPluginIndexParams()[data_param])
}

//
// Resolve the configured VEP plugins into three things: the VEP argument string,
// the already-prepared files to stage into the VEP task directory, and the raw
// releases still to be reshaped by a prep task.
//
// The module stages plugin files through its `extra_files` input, which lands
// them in the task workdir root, so every argument references a bare basename
// rather than the original path or URL.
//
def resolveVepPlugins() {
    if (vepPluginsSkipped()) {
        return [ args: '', ready_files: [], prepare: [:] ]
    }

    def staged = []
    def prepare = [:]
    def args = []

    if (vepPluginResource('vep_alphamissense')) {
        args << "--plugin AlphaMissense,file=${registerVepPlugin(staged, prepare, 'vep_alphamissense')}"
    }

    if (vepPluginResource('vep_alphamissense_aa')) {
        // Our own plugin, so the .pm has to travel alongside its data. Adding
        // --dir_plugins only prepends to @INC, leaving the plugins bundled in
        // the container reachable.
        staged << file("${projectDir}/assets/vep_plugins/AlphaMissenseProtein.pm", checkIfExists: true)
        args << "--dir_plugins ."
        args << "--plugin AlphaMissenseProtein,file=${registerVepPlugin(staged, prepare, 'vep_alphamissense_aa')}"
    }

    if (vepPluginResource('vep_polyphen_sift_db')) {
        args << "--plugin PolyPhen_SIFT,db=${registerVepPlugin(staged, prepare, 'vep_polyphen_sift_db')}"
    }

    if (vepPluginResource('vep_clinvar')) {
        // --custom takes a %-separated field list, unlike the comma-separated
        // form used everywhere else.
        def fields = (params.vep_clinvar_fields ?: '').tokenize(',').collect { it.trim() }.findAll().join('%')
        def clinvar = "--custom file=${registerVepPlugin(staged, prepare, 'vep_clinvar')},short_name=ClinVar,format=vcf,type=exact,coords=0"
        args << (fields ? "${clinvar},fields=${fields}" : clinvar)
    }

    if (vepPluginResource('vep_cadd_snv') || vepPluginResource('vep_cadd_indel')) {
        def cadd = []
        if (vepPluginResource('vep_cadd_snv')) {
            cadd << "snv=${registerVepPlugin(staged, prepare, 'vep_cadd_snv')}"
        }
        if (vepPluginResource('vep_cadd_indel')) {
            cadd << "indels=${registerVepPlugin(staged, prepare, 'vep_cadd_indel')}"
        }
        args << "--plugin CADD,${cadd.join(',')}"
    }

    if (vepPluginResource('vep_revel')) {
        args << "--plugin REVEL,file=${registerVepPlugin(staged, prepare, 'vep_revel')}"
    }

    if (vepPluginResource('vep_eve')) {
        args << "--plugin EVE,file=${registerVepPlugin(staged, prepare, 'vep_eve')}"
    }

    // --vep_custom keeps its existing mechanism: it is passed to the module as
    // its own input, and the user supplies the matching --custom skeleton in
    // --vep_args which the module rewrites to the staged path. Nothing to do
    // here beyond leaving that first --custom entry alone.

    return [ args: args.join(' '), ready_files: staged, prepare: prepare ]
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016),",
            "Samtools (Li et al. 2009),",
            "Mosdepth (Pedersen and Quinlan 2018)."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
            "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. doi: 10.1093/bioinformatics/btp352</li>",
            "<li>Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867-868. doi: 10.1093/bioinformatics/btx699</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
