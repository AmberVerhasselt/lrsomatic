process VEPPLUGIN_ALPHAMISSENSE_PROTEIN {
    tag "${aa_substitutions}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    path aa_substitutions
    path idmapping

    output:
    // The table and its index together: both are staged into the VEP task.
    path "alphamissense_protein.tsv.gz{,.tbi}", emit: files
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    build_alphamissense_protein_table.sh \\
        -a ${aa_substitutions} \\
        -m ${idmapping} \\
        -o alphamissense_protein.tsv.gz \\
        -t . \\
        ${args}
    """

    stub:
    """
    touch alphamissense_protein.tsv.gz
    touch alphamissense_protein.tsv.gz.tbi
    """
}
