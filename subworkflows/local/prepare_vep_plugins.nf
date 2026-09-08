//
// Turn the resolved VEP plugin configuration into the files VEP needs staged.
//
// Most resources are published in a form VEP can read and are staged straight
// from their release. Three cannot be: AlphaMissense ships without a tabix
// index, and REVEL and EVE ship as zip archives. Those are reshaped here, so a
// default run does not require the user to prepare anything by hand first.
//
// Nothing is published. The prepared files live in the work directory, are
// reused by every VEP task in the run, and survive a -resume against the same
// work directory. Passing a prepared file to the matching --vep_* param skips
// both the download and the prep task.
//

include { UNZIP as UNZIP_REVEL            } from '../../modules/nf-core/unzip/main.nf'
include { UNZIP as UNZIP_EVE              } from '../../modules/nf-core/unzip/main.nf'
include { VEPPLUGIN_ALPHAMISSENSE_INDEX   } from '../../modules/local/vepplugin/alphamissense_index/main.nf'
include { VEPPLUGIN_ALPHAMISSENSE_PROTEIN } from '../../modules/local/vepplugin/alphamissense_protein/main.nf'
include { VEPPLUGIN_REVEL                 } from '../../modules/local/vepplugin/revel/main.nf'
include { VEPPLUGIN_EVE                   } from '../../modules/local/vepplugin/eve/main.nf'

include { resolveVepPlugins               } from './utils_nfcore_lrsomatic_pipeline'

workflow PREPARE_VEP_PLUGINS {

    main:

    // Resolved here rather than taken as an input: which prep tasks to invoke
    // is a plain Groovy decision, and a workflow input would arrive wrapped in
    // a channel. resolveVepPlugins() is pure, so calling it again is free.
    def plugins = resolveVepPlugins()
    def prepare = plugins.prepare

    ch_versions = channel.empty()
    // Each entry is a channel of files to stage; mixed into one list below.
    def staged = [ channel.fromList(plugins.ready_files) ]

    //
    // MODULE: VEPPLUGIN_ALPHAMISSENSE_INDEX (label: process_low)
    // Input:  the AlphaMissense score file as published, which has no index
    // Output: .tbi -- the index, alongside the score file it belongs to
    //
    if (prepare.containsKey('vep_alphamissense')) {
        VEPPLUGIN_ALPHAMISSENSE_INDEX(
            prepare['vep_alphamissense']
        )
        // The score file is staged as published, so it travels with its index.
        staged << channel.value(prepare['vep_alphamissense'])
        staged << VEPPLUGIN_ALPHAMISSENSE_INDEX.out.tbi
    }

    //
    // MODULE: VEPPLUGIN_ALPHAMISSENSE_PROTEIN (label: process_medium)
    // Input:  the AlphaMissense protein-space release + the UniProt ID mapping
    // Output: .files -- alphamissense_protein.tsv.gz and its index
    //
    // The ID mapping is passed in rather than downloaded by the script, so the
    // task needs no network access and Nextflow caches the file.
    //
    if (prepare.containsKey('vep_alphamissense_aa')) {
        VEPPLUGIN_ALPHAMISSENSE_PROTEIN(
            prepare['vep_alphamissense_aa'],
            file(params.vep_uniprot_idmapping, checkIfExists: true)
        )
        staged << VEPPLUGIN_ALPHAMISSENSE_PROTEIN.out.files
    }

    //
    // MODULES: UNZIP_REVEL -> VEPPLUGIN_REVEL (labels: process_single, process_medium)
    // Input:  revel-v1.3_all_chromosomes.zip
    // Output: .files -- revel_grch38.tsv.gz and its index
    //
    // Unzipped by the p7zip module rather than inside the prep container, so no
    // container needs both htslib and unzip.
    //
    if (prepare.containsKey('vep_revel')) {
        UNZIP_REVEL(
            channel.value([ [ id: 'revel' ], prepare['vep_revel'] ])
        )
        VEPPLUGIN_REVEL(
            UNZIP_REVEL.out.unzipped_archive.map { _meta, dir -> dir }
        )
        staged << VEPPLUGIN_REVEL.out.files
        ch_versions = ch_versions.mix(UNZIP_REVEL.out.versions)
    }

    //
    // MODULES: UNZIP_EVE -> VEPPLUGIN_EVE (labels: process_single, process_medium)
    // Input:  EVE_all_data.zip -- one VCF per protein
    // Output: .files -- eve_merged.vcf.gz and its index
    //
    if (prepare.containsKey('vep_eve')) {
        UNZIP_EVE(
            channel.value([ [ id: 'eve' ], prepare['vep_eve'] ])
        )
        VEPPLUGIN_EVE(
            UNZIP_EVE.out.unzipped_archive.map { _meta, dir -> dir }
        )
        staged << VEPPLUGIN_EVE.out.files
        ch_versions = ch_versions.mix(UNZIP_EVE.out.versions)
    }

    // One list of every plugin file, consumed by both the germline and the
    // somatic VEP task -- so it has to be a value channel, which can be read
    // more than once. collect() gives one, but emits nothing at all on an empty
    // upstream, which would leave the VEP tasks waiting on a list that never
    // arrives; hence the ifEmpty, which is what carries the no-plugins case.
    ch_extra_files = staged
        .inject(channel.empty()) { acc, ch -> acc.mix(ch) }
        .flatten()
        .collect()
        .ifEmpty([])

    emit:
    extra_files = ch_extra_files // channel: value list of plugin .pm and data files
    versions    = ch_versions    // channel: versions.yml files
}
