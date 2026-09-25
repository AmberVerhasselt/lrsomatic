process VEPPLUGIN_CLINVAR {
    tag "${vcf_url.toString().tokenize('/').last()}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data'
        : 'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e'}"

    input:
    tuple val(vcf_url), val(tbi_url), val(md5)

    output:
    path "${vcf_name}{,.tbi}", emit: files
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -1 | cut -d ' ' -f 3"), topic: versions, emit: versions_wget

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    vcf_name = vcf_url.toString().tokenize('/').last()
    def check = md5 ? "echo '${md5}  ${vcf_name}' | md5sum -c -" : ''
    """
    # NCBI answers bursts with 503, so retry those rather than fail on the first one
    wget \\
        --no-verbose \\
        --tries=5 \\
        --waitretry=30 \\
        --retry-on-http-error=429,500,502,503,504 \\
        ${args} \\
        -O ${vcf_name} \\
        ${vcf_url}

    # Saved next to the VCF under the name VEP looks for, whatever the host calls it
    wget \\
        --no-verbose \\
        --tries=5 \\
        --waitretry=30 \\
        --retry-on-http-error=429,500,502,503,504 \\
        ${args} \\
        -O ${vcf_name}.tbi \\
        ${tbi_url}

    # A pinned checksum keeps the release fixed: a host that re-publishes under the same name fails here
    ${check}
    """

    stub:
    vcf_name = vcf_url.toString().tokenize('/').last()
    """
    echo "" | gzip > ${vcf_name}
    touch ${vcf_name}.tbi
    """
}
