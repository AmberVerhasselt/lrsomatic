process CLAIRSTO {
    tag "$meta.id"
    label 'process_very_high'

    // Fork of ClairS-TO 0.5.1 that resolves Verdict's CNA resources from --cna_resource_dir
    // instead of hardcoded GRCh38 names. No conda build; revert once upstream carries it.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/clairs-to-sif:0.5.1-verdict-chm13-c0687e8'
        : 'ghcr.io/ljwharbers/clairs-to:0.5.1-verdict-chm13-c0687e8'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), val(model), path(pon_vcfs), val(pon_flags)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(index)
    // Verdict's ASCAT set for this assembly, from CLAIRSTO_CNA_RESOURCES; [] uses the image's own
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
    // The launcher cannot locate its own conda environment inside a read-only image
    def conda_prefix = (workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') ? '--conda_prefix /opt/micromamba/envs/clairs-to' : ''
    // Omitted for GRCh38; a set that cannot belong to --ref_fn disables Verdict with a warning
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

    # 0.5.1 names the VCFs after --sample_name and ignores --snv_output_prefix while it holds its
    # default. Rename back, globbed because ClairS-TO sanitises the sample name inside the path.
    mv -- snv_*.vcf.gz snv.vcf.gz
    mv -- snv_*.vcf.gz.tbi snv.vcf.gz.tbi
    mv -- indel_*.vcf.gz indel.vcf.gz
    mv -- indel_*.vcf.gz.tbi indel.vcf.gz.tbi
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
