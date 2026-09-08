#!/usr/bin/env bash
#
# Preparation of VEP plugin data files.
#
# AlphaMissense, REVEL and EVE cannot be handed to VEP exactly as published:
# AlphaMissense needs a tabix index, REVEL ships comma-separated and sorted on
# its GRCh37 column, and EVE ships as thousands of per-protein VCFs.
#
# The pipeline calls this script itself, from modules/local/vepplugin/*, so a
# default run needs no manual preparation. Run it by hand to produce files you
# can then pass to --vep_alphamissense / --vep_revel / --vep_eve, which skips
# both the download and the prep task on every subsequent run.
#
# Nothing is redistributed: the inputs are fetched from their original source
# and the outputs stay on your filesystem. See docs/usage.md for the download
# URLs and the licence terms of each resource -- CADD, REVEL and EVE are free
# for non-commercial use only, and observing that is the user's responsibility.
#
# Usage:
#   bin/prepare_vep_plugin_data.sh alphamissense AlphaMissense_hg38.tsv.gz
#   bin/prepare_vep_plugin_data.sh revel <zip-or-unpacked-dir> <outdir>
#   bin/prepare_vep_plugin_data.sh eve <eve-vcf-dir> <outdir>
#
# ClinVar, CADD and the Ensembl pangenome PolyPhen/SIFT database need no
# preparation: ClinVar and CADD ship their own .tbi, and the pangenome file is
# an SQLite database.
#
# Requires: bgzip and tabix (htslib), awk, sort, and unzip when REVEL is given
# the release zip rather than an already-unpacked directory.

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

need() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null || die "$tool not found on PATH"
    done
}

usage() {
    sed -n '2,29p' "$0" | sed 's/^#\{1,2\} \{0,1\}//'
    exit "${1:-1}"
}

# AlphaMissense needs only an index. Count the leading non-data lines rather than hardcoding -S 1:
# tabix applies -S before comment detection, so a wrong count silently parses the header row as data.
prep_alphamissense() {
    local tsv=${1:-}
    [[ -n $tsv ]] || die "usage: $0 alphamissense <AlphaMissense_hg38.tsv.gz>"
    [[ -r $tsv ]] || die "cannot read $tsv"
    need bgzip tabix

    local skip
    # `|| true` because awk exits at the first data line, so gzip takes a SIGPIPE that
    # pipefail would otherwise turn into a silent abort of the whole script.
    skip=$( { gzip -dc "$tsv" || true; } | awk -F'\t' '
        $2 ~ /^[0-9]+$/ { print NR - 1; found = 1; exit }
        NR > 100        { exit }
        END             { if (!found) print "NONE" }
    ')
    [[ $skip != "NONE" && -n $skip ]] \
        || die "found no data line in the first 100 lines of $tsv; is this an AlphaMissense coordinate file?"
    echo ">> indexing $tsv (skipping $skip header line(s))" >&2

    tabix -s 1 -b 2 -e 2 -f -S "$skip" "$tsv"
    echo ">> wrote $tsv.tbi" >&2
    echo >&2
    echo "Pass these to the pipeline with:" >&2
    echo "  --vep_alphamissense $tsv --vep_alphamissense_tbi $tsv.tbi" >&2
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
    trap 'rm -rf "$work"' RETURN

    # Either the release zip or a directory it has already been unpacked into,
    # so the pipeline can leave the unpacking to a container that has unzip.
    local searchdir=$src
    if [[ -f $src ]]; then
        need unzip
        echo ">> unpacking $src" >&2
        unzip -o -q -d "$work" "$src"
        searchdir=$work
    elif [[ ! -d $src ]]; then
        die "$src is neither a zip file nor a directory"
    fi

    local raw
    raw=$(find "$searchdir" -name 'revel_with_transcript_ids' -o -name 'revel_all_chromosomes.csv' | head -n1)
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
    trap 'rm -rf "$work"' RETURN

    # The bulk zip nests the per-protein files under vcf_files_missense_mutations.
    local src=$vcfdir
    if [[ -d "$vcfdir/vcf_files_missense_mutations" ]]; then
        src="$vcfdir/vcf_files_missense_mutations"
    fi

    local count
    count=$(find "$src" -name '*.vcf' | wc -l)
    [[ $count -gt 0 ]] || die "no .vcf files found under $src"
    echo ">> merging $count per-protein VCFs from $src" >&2

    local first
    first=$(find "$src" -name '*.vcf' | sort | head -n1)

    local out="$outdir/eve_merged.vcf.gz"
    {
        grep '^#' "$first"
        find "$src" -name '*.vcf' -exec grep -hv '^#' {} + \
          | sort -T "$work" -k1,1V -k2,2n
    } | bgzip -c > "$out"

    tabix -f -p vcf "$out"
    echo ">> wrote $out and $out.tbi" >&2
    echo >&2
    echo "Pass these to the pipeline with:" >&2
    echo "  --vep_eve $out --vep_eve_tbi $out.tbi" >&2
}

case "${1:-}" in
    alphamissense) shift; prep_alphamissense "$@" ;;
    revel)         shift; prep_revel "$@" ;;
    eve)           shift; prep_eve "$@" ;;
    -h|--help|help) usage 0 ;;
    "")            die "no resource given. Try: $0 --help" ;;
    *)             die "unknown resource '$1'. Expected alphamissense, revel or eve." ;;
esac
