#!/usr/bin/env bash
#
# Build the protein-space AlphaMissense lookup table used by the
# AlphaMissenseProtein VEP plugin (assets/vep_plugins/AlphaMissenseProtein.pm),
# which is how AlphaMissense scores reach a T2T-CHM13 run.
#
# AlphaMissense_aa_substitutions.tsv.gz is already keyed in protein space, by
# UniProt accession and amino-acid substitution. This script re-keys it to gene
# symbol -- the join key the plugin uses, because CHM13 rapid-release stable IDs
# are unrelated to GRCh38 Ensembl IDs -- and writes a tabix-indexed table.
#
# The pipeline calls this script itself, from
# modules/local/vepplugin/alphamissense_protein, so a CHM13 run needs no manual
# preparation. Run it by hand to produce a table you can then pass to
# --vep_alphamissense_aa, which skips both the download and the prep task.
#
# Nothing here is redistributed: the inputs are fetched from their original
# source and the output stays on your filesystem.
#
# DATA LICENCE
#   AlphaMissense Database, Copyright (2023) DeepMind Technologies Limited,
#   licensed CC BY 4.0 <https://creativecommons.org/licenses/by/4.0/legalcode>.
#   CC BY permits this transformation and its redistribution, provided credit is
#   given and changes are indicated. CHANGES MADE BY THIS SCRIPT: values are
#   re-keyed from UniProt accession to gene symbol and reshaped into a
#   tabix-indexed table. No score is altered.
#
#   The UniProt ID mapping is from UniProtKB, licensed CC BY 4.0.
#
# Usage:
#   bin/build_alphamissense_protein_table.sh \
#       -a AlphaMissense_aa_substitutions.tsv.gz \
#       -o alphamissense_protein.tsv.gz \
#       [-m HUMAN_9606_idmapping.dat.gz] \
#       [-t tmpdir]
#
# Inputs, both public and unauthenticated:
#   AlphaMissense  https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz
#   ID mapping     https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
#                  (downloaded automatically if -m is not given)
#
# Requires: bgzip and tabix (htslib), awk, sort, and curl when -m is omitted.

set -euo pipefail

ALPHAMISSENSE=""
IDMAPPING=""
OUTPUT=""
TMPDIR_ARG=""

usage() {
    sed -n '2,42p' "$0" | sed 's/^#\{1,2\} \{0,1\}//'
    exit "${1:-1}"
}

while getopts ":a:m:o:t:h" opt; do
    case "$opt" in
        a) ALPHAMISSENSE=$OPTARG ;;
        m) IDMAPPING=$OPTARG ;;
        o) OUTPUT=$OPTARG ;;
        t) TMPDIR_ARG=$OPTARG ;;
        h) usage 0 ;;
        *) echo "ERROR: unknown option -$OPTARG" >&2; usage ;;
    esac
done

[[ -n $ALPHAMISSENSE ]] || { echo "ERROR: -a (AlphaMissense_aa_substitutions.tsv.gz) is required" >&2; usage; }
[[ -n $OUTPUT ]] || { echo "ERROR: -o (output path) is required" >&2; usage; }
[[ -r $ALPHAMISSENSE ]] || { echo "ERROR: cannot read $ALPHAMISSENSE" >&2; exit 1; }

for tool in bgzip tabix awk sort; do
    command -v "$tool" >/dev/null || { echo "ERROR: $tool not found on PATH" >&2; exit 1; }
done

# Sorting this table needs real scratch space, so honour $TMPDIR rather than
# defaulting to /tmp.
WORK=$(mktemp -d "${TMPDIR_ARG:-${TMPDIR:-/tmp}}/amprot.XXXXXX")
trap 'rm -rf "$WORK"' EXIT

if [[ -z $IDMAPPING ]]; then
    command -v curl >/dev/null || { echo "ERROR: curl not found; supply -m instead" >&2; exit 1; }
    IDMAPPING="$WORK/HUMAN_9606_idmapping.dat.gz"
    echo ">> downloading UniProt ID mapping" >&2
    curl -fsSL -o "$IDMAPPING" \
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
fi
[[ -r $IDMAPPING ]] || { echo "ERROR: cannot read $IDMAPPING" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 1. UniProt accession -> gene symbol.
#
# HUMAN_9606_idmapping.dat.gz is three columns: accession, id type, value. The
# Gene_Name rows carry the symbol. An accession can list several; keep the first,
# which is the primary name.
# ---------------------------------------------------------------------------
echo ">> extracting UniProt accession to gene symbol map" >&2
gzip -dc "$IDMAPPING" \
  | awk -F'\t' '$2 == "Gene_Name" && !($1 in seen) { seen[$1] = 1; print $1 "\t" $3 }' \
  > "$WORK/acc2sym.tsv"

MAPPED=$(wc -l < "$WORK/acc2sym.tsv")
echo "   $MAPPED accessions mapped to a gene symbol" >&2
[[ $MAPPED -gt 1000 ]] || { echo "ERROR: ID mapping looks empty or malformed" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 2. Re-key AlphaMissense onto gene symbols.
#
# Input columns (after the leading '#' comment lines and a header row):
#   uniprot_id  protein_variant  am_pathogenicity  am_class
# protein_variant is a one-letter substitution such as V2L, matching the form of
# VEP's pep_allele_string.
#
# Rows whose accession has no symbol are dropped: without a join key the plugin
# could never retrieve them.
# ---------------------------------------------------------------------------
echo ">> re-keying AlphaMissense onto gene symbols" >&2
gzip -dc "$ALPHAMISSENSE" \
  | awk -F'\t' -v mapfile="$WORK/acc2sym.tsv" '
      BEGIN {
        OFS = "\t"
        while ((getline line < mapfile) > 0) {
          split(line, f, "\t")
          sym[f[1]] = f[2]
        }
        close(mapfile)
      }
      # skip the release comment lines and the column-name row
      /^#/      { next }
      $1 == "uniprot_id" { next }
      NF < 4    { next }
      {
        acc = $1
        # AlphaMissense lists some accessions with an isoform suffix (P12345-2);
        # the mapping is keyed on the bare accession.
        base = acc
        sub(/-[0-9]+$/, "", base)
        gene = (acc in sym) ? sym[acc] : (base in sym ? sym[base] : "")
        if (gene == "") { unmapped++; next }

        # V2L -> aaref V, aapos 2, aaalt L
        if ($2 !~ /^[A-Z][0-9]+[A-Z]$/) { malformed++; next }
        variant = $2
        aaref = substr(variant, 1, 1)
        aaalt = substr(variant, length(variant), 1)
        aapos = substr(variant, 2, length(variant) - 2)
        if (aapos + 0 <= 0) { malformed++; next }

        print gene, aapos, aaref, aaalt, $3, $4, acc
        kept++
      }
      END {
        printf("   kept %d rows; dropped %d unmapped accessions, %d malformed variants\n",
               kept, unmapped, malformed) > "/dev/stderr"
      }
    ' \
  > "$WORK/rekeyed.tsv"

[[ -s $WORK/rekeyed.tsv ]] || { echo "ERROR: no rows survived re-keying" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 3. Sort, deduplicate, index.
#
# Sorted by gene then numeric amino-acid position, which is what tabix needs.
# Where several UniProt accessions share a gene symbol the same substitution can
# appear more than once; keep the first and note how many were dropped, since
# the plugin takes the first exact match anyway.
# ---------------------------------------------------------------------------
echo ">> sorting and indexing" >&2
mkdir -p "$(dirname "$OUTPUT")"

{
    printf '#gene\taapos\taaref\taaalt\tam_pathogenicity\tam_class\tuniprot_acc\n'
    sort -T "$WORK" -k1,1 -k2,2n -k3,3 -k4,4 "$WORK/rekeyed.tsv" \
      | awk -F'\t' '
          { key = $1 "\t" $2 "\t" $3 "\t" $4 }
          key != prev { print; prev = key; next }
          { dup++ }
          END { if (dup) printf("   dropped %d duplicate gene/substitution rows\n", dup) > "/dev/stderr" }
        '
} | bgzip -c > "$OUTPUT"

# -s 1 -b 2 -e 2: the "sequence" column is the gene symbol, a plain string, and
# the "position" columns are the amino-acid position. -c '#' keeps the header
# readable via `tabix -H`, which is how the plugin discovers the columns.
tabix -f -s 1 -b 2 -e 2 -c '#' "$OUTPUT"

echo ">> wrote $OUTPUT and $OUTPUT.tbi" >&2
echo ">> genes indexed: $(tabix -l "$OUTPUT" | wc -l)" >&2
echo >&2
echo "Pass these to the pipeline with:" >&2
echo "  --vep_alphamissense_aa $OUTPUT --vep_alphamissense_aa_tbi $OUTPUT.tbi" >&2
