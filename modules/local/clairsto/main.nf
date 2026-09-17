process CLAIRSTO {
    tag "$meta.id"
    label 'process_very_high'

    // Conda is not supported: the image is the fork that lets Verdict resolve its CNA resources
    // from --cna_resource_dir instead of hardcoded GRCh38 names, so it can run against another
    // assembly (e.g. T2T-CHM13). It also checks the loci against --ref_fn, makes the replication
    // timing file optional, and fixes the GC window selection. Return to docker.io/hkubal/clairs-to
    // once HKU-BAL/ClairS-TO carries these changes.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/clairs-to-sif:0.5.1-verdict-chm13'
        : 'ghcr.io/ljwharbers/clairs-to:0.5.1-verdict-chm13'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(model), path(pon_vcfs), val(pon_flags)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)
    // Verdict's ASCAT loci/allele/GC set for the assembly the BAM was aligned to, or [] to use the
    // GRCh38 set shipped inside the image. Built by CLAIRSTO_CNA_RESOURCES.
    tuple val(meta4), path(cna_resources)

    output:
    tuple val(meta), path("indel.vcf.gz"),      emit: indel_vcf
    tuple val(meta), path("indel.vcf.gz.tbi"),  emit: indel_tbi
    tuple val(meta), path("snv.vcf.gz"),        emit: snv_vcf
    tuple val(meta), path("snv.vcf.gz.tbi"),    emit: snv_tbi
    tuple val("${task.process}"), val('clairsto'), eval("run_clairs_to  --version |& sed '1!d ; s/run_clairs_to //'"), topic: versions, emit: versions_clairsto

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // The launcher activates its own conda environment, which it cannot locate for itself inside a
    // read-only image; apptainer needs this as much as singularity does.
    def conda_prefix = (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') ? '--conda_prefix /opt/micromamba/envs/clairs-to' : ''
    // Omitted for GRCh38, where the image's own resource set is the right one. A set that cannot
    // belong to --ref_fn makes ClairS-TO disable Verdict with a warning rather than mistag.
    def cna_resource_dir = cna_resources ? "--cna_resource_dir ${cna_resources}" : ''
    def pon_string   = pon_vcfs.join(',')
    def flags_string = pon_flags.join(',')

    """
    /opt/bin/run_clairs_to \
        --tumor_bam_fn $tumor_bam \\
        --ref_fn $reference \\
        --platform $model \\
        --threads $task.cpus \\
        --output_dir . \\
        --sample_name ${prefix} \\
        --panel_of_normals ${pon_string} \\
        --panel_of_normals_require_allele_matching ${flags_string} \\
        $conda_prefix \\
        $cna_resource_dir \\
        $args
    """

    stub:
    """
    mkdir -p output
    echo "" | gzip > snv.vcf.gz
    touch snv.vcf.gz.tbi
    echo "" | gzip > indel.vcf.gz
    touch indel.vcf.gz.tbi
    """
}
