process VEPPLUGIN_ALPHAMISSENSE_INDEX {
    tag "${tsv}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    path tsv

    output:
    // Only the index: the score file itself is already staged into the VEP task
    // straight from the release, under this same name.
    path "*.tbi", emit: tbi
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prepare_vep_plugin_data.sh alphamissense ${tsv} ${args}
    """

    stub:
    """
    touch ${tsv}.tbi
    """
}
