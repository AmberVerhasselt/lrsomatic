process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    // Conda is not supported: the image carries the lrsomatic_report tool itself, not just its
    // dependencies, so an environment.yml would install the R/Quarto stack without render_report.R
    // (the guard in `script:` stops conda/mamba runs). Switch to a bioconda `lrsomatic-report`
    // package once the recipe at github.com/ljwharbers/lrsomatic_report/tree/main/recipe is merged.
    // Updating the tool is a tag bump here: tag upstream, let its container workflow build, edit these two lines.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/lrsomatic-report-sif:1.5.0'
        : 'ghcr.io/ljwharbers/lrsomatic-report:1.5.0'}"

    input:
    // Every path input is optional (`[]` when skipped); tumor/normal QC stage into separate dirs because a matched pair shares meta.id
    tuple val(meta), path(vep_somatic), path(sv_vep), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*'), path(wakhan_files, stageAs: 'wakhan/*')
    // Builtin gene panel TSVs, owned by this pipeline rather than by the tool, so the panel set
    // can change without a tool release; reaches the tool as --gene-lists-dir
    path(gene_lists, stageAs: 'gene_lists')
    // User-supplied gene panel TSVs (`[]` for builtins); the matching `--gene-panel gene_panels/<base>` args are built in conf/modules.config
    path(gene_panels, stageAs: 'gene_panels/*')

    output:
    tuple val(meta), path("*_report.html"), emit: report
    tuple val("${task.process}"), val('lrsomatic_report'), eval('render_report.R --version'), topic: versions, emit: versions_lrsomaticreport

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "LRSOMATICREPORT does not support Conda: the report tool ships only inside its container. Use Docker / Singularity / Apptainer, or --skip_report."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'male'

    // Discovery is recursive and matches on base name, so suffix-distinct files can be linked flat
    def flat_inputs = [vep_somatic, sv_vep, severus_vcf, ascat_files].flatten().findAll { f -> f }
    def link_flat = flat_inputs ? """
    for f in ${flat_inputs.collect { f -> "\"${f}\"" }.join(' ')}; do ln -s "\$PWD/\$f" "sample_dir/\$f"; done
    """ : ''

    // The VAF/depth/phasing source is looked up at a literal path
    def link_somatic = somatic_vcf ? """
    mkdir -p sample_dir/variants/phased
    ln -s "\$PWD/${somatic_vcf}" sample_dir/variants/phased/somatic_smallvariants.vcf.gz
    """ : ''

    """
    # Quarto/Deno write under \$HOME and \$TMPDIR, which clusters may mount read-only
    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp TMP=\$PWD/tmp TEMP=\$PWD/tmp
    mkdir -p "\$TMPDIR"

    mkdir -p sample_dir
    ${link_flat}
    ${link_somatic}

    # Link file by file: R's list.files(recursive = TRUE) does not descend into symlinked dirs
    if [ -d qc_tumor ]; then
        mkdir -p sample_dir/qc/tumor
        for f in qc_tumor/*; do ln -s "\$PWD/\$f" "sample_dir/qc/tumor/\$(basename "\$f")"; done
    fi
    if [ -d qc_normal ]; then
        mkdir -p sample_dir/qc/normal
        for f in qc_normal/*; do ln -s "\$PWD/\$f" "sample_dir/qc/normal/\$(basename "\$f")"; done
    fi

    # Wakhan is addressed by fixed path: sample_dir/wakhan must hold solutions_ranks.tsv, the heatmap and solution_<rank>/
    if [ -d wakhan ]; then
        mkdir -p sample_dir/wakhan
        for f in wakhan/*; do ln -s "\$PWD/\$f" "sample_dir/wakhan/\$(basename "\$f")"; done
    fi

    render_report.R \\
        --sample-dir sample_dir \\
        --sample-id "${prefix}" \\
        --sex "${sex}" \\
        --reference auto \\
        --gene-lists-dir gene_lists \\
        --output "${prefix}_report.html" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
