#!/usr/bin/env bash
#
# Reshape the REVEL and EVE releases into the form VEP expects. The pipeline runs
# this itself; run it by hand to pre-build files for --vep_revel / --vep_eve.
#
# Usage:
#   bin/prepare_vep_plugin_data.sh revel <zip-or-unpacked-dir> <outdir>
#   bin/prepare_vep_plugin_data.sh eve <eve-vcf-dir> <outdir>
#
# Requires: bgzip and tabix (htslib), awk, sort, and unzip for the REVEL zip.

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

need() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null || die "$tool not found on PATH"
    done
}

usage() {
    sed -n '2,10p' "$0" | sed 's/^#\{1,2\} \{0,1\}//'
    exit "${1:-1}"
}

# REVEL ships comma-separated, sorted on its GRCh37 column (2); GRCh38 (column 3) needs a re-sort and its own index.
prep_revel() {
    local src=${1:-}
    local outdir=${2:-.}
    [[ -n $src ]] || die "usage: $0 revel <revel-v1.3_all_chromosomes.zip|unpacked-dir> [outdir]"
    [[ -r $src ]] || die "cannot read $src"
    need bgzip tabix awk sort

    mkdir -p "$outdir"
    local work
    work=$(mktemp -d "${TMPDIR:-/tmp}/revel.XXXXXX")
    # EXIT rather than RETURN: a RETURN trap does not fire when set -e aborts mid-function.
    # The path is expanded now, not when the trap fires -- by then this function has returned
    # and its local $work is gone, which under set -u is a fatal "unbound variable".
    trap "rm -rf '$work'" EXIT

    # Either the release zip or a directory it has already been unpacked into
    local searchdir=$src
    if [[ -f $src ]]; then
        need unzip
        echo ">> unpacking $src" >&2
        unzip -o -q -d "$work" "$src"
        searchdir=$work
    elif [[ ! -d $src ]]; then
        die "$src is neither a zip file nor a directory"
    fi

    # -L because the pipeline passes a Nextflow-staged directory, which is a symlink;
    # plain find reports it and never descends.
    #
    # Collected with mapfile rather than piped to `head -n1`: under `set -o pipefail` an
    # early-exiting reader leaves the producer with SIGPIPE, and the pipeline then returns
    # 141, which `set -e` turns into an abort. See prep_eve, where that is not merely a risk.
    local raw
    mapfile -t revel_candidates < <(find -L "$searchdir" \( -name 'revel_with_transcript_ids' -o -name 'revel_all_chromosomes.csv' \) | sort)
    raw=${revel_candidates[0]:-}
    [[ -n $raw ]] || die "could not find the REVEL table in $src"
    echo "   using $(basename "$raw")" >&2

    local out="$outdir/revel_grch38.tsv.gz"
    echo ">> converting to tab-separated and re-sorting on the GRCh38 column" >&2
    {
        # header row, tab-separated and '#'-prefixed so tabix treats it as a comment
        head -n1 "$raw" | tr ',' '\t' | sed '1s/^/#/'
        # data rows: drop those with no GRCh38 position, then sort on it
        tail -n +2 "$raw" | tr ',' '\t' | awk -F'\t' '$3 != "." && $3 != ""' \
          | sort -T "$work" -k1,1 -k3,3n
    } | bgzip -c > "$out"

    tabix -f -s 1 -b 3 -e 3 -c '#' "$out"
    echo ">> wrote $out and $out.tbi" >&2
    echo >&2
    echo "Pass these to the pipeline with:" >&2
    echo "  --vep_revel $out --vep_revel_tbi $out.tbi" >&2
}

# EVE ships one VCF per protein; the plugin expects a single sorted, indexed VCF.
prep_eve() {
    local vcfdir=${1:-}
    local outdir=${2:-.}
    [[ -n $vcfdir ]] || die "usage: $0 eve <dir-of-per-protein-vcfs> [outdir]"
    [[ -d $vcfdir ]] || die "$vcfdir is not a directory"
    need bgzip tabix awk sort

    mkdir -p "$outdir"
    local work
    work=$(mktemp -d "${TMPDIR:-/tmp}/eve.XXXXXX")
    # EXIT rather than RETURN: a RETURN trap does not fire when set -e aborts mid-function.
    # The path is expanded now, not when the trap fires -- by then this function has returned
    # and its local $work is gone, which under set -u is a fatal "unbound variable".
    trap "rm -rf '$work'" EXIT

    # The bulk zip nests the per-protein files under vcf_files_missense_mutations.
    local src=$vcfdir
    if [[ -d "$vcfdir/vcf_files_missense_mutations" ]]; then
        src="$vcfdir/vcf_files_missense_mutations"
    fi

    # -L throughout because the pipeline passes a Nextflow-staged directory, which is a
    # symlink; plain find reports it and never descends.
    #
    # Collected once with mapfile rather than piped to `head -n1`. With thousands of files
    # `head` exits long before `sort` has finished writing, `sort` takes SIGPIPE, and under
    # `set -o pipefail` the pipeline returns 141 -- which `set -e` turns into an abort. The
    # release has ~3000 VCFs, so that fired every time.
    local vcfs count first
    mapfile -t vcfs < <(find -L "$src" -name '*.vcf' | sort)
    count=${#vcfs[@]}
    [[ $count -gt 0 ]] || die "no .vcf files found under $src"
    echo ">> merging $count per-protein VCFs from $src" >&2

    first=${vcfs[0]}

    local out="$outdir/eve_merged.vcf.gz"
    {
        grep '^#' "$first"
        find -L "$src" -name '*.vcf' -exec grep -hv '^#' {} + \
          | sort -T "$work" -k1,1V -k2,2n
    } | bgzip -c > "$out"

    tabix -f -p vcf "$out"
    echo ">> wrote $out and $out.tbi" >&2
    echo >&2
    echo "Pass these to the pipeline with:" >&2
    echo "  --vep_eve $out --vep_eve_tbi $out.tbi" >&2
}

case "${1:-}" in
    revel)         shift; prep_revel "$@" ;;
    eve)           shift; prep_eve "$@" ;;
    -h|--help|help) usage 0 ;;
    "")            die "no resource given. Try: $0 --help" ;;
    *)             die "unknown resource '$1'. Expected revel or eve." ;;
esac
