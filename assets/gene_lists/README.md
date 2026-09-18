# Gene Panel Lists

These panels belong to the pipeline, not to the report tool. `LRSOMATICREPORT` stages this
directory into the task and passes it as `--gene-lists-dir`, which replaces the set bundled
inside `ghcr.io/ljwharbers/lrsomatic-report`. Adding, editing or removing a panel here is a
pipeline change and needs no release of
[lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) — drop a TSV in, and it
becomes a builtin that `--report_gene_panel` accepts and the report offers as a tickbox.

The files originated in that tool (MIT licensed) and were moved here when the pipeline
stopped vendoring its source.

Each file is a TSV with a required `gene` column (HGNC symbol). Coordinate columns
`chrom` (or `chr`), `start` and `end` are optional but **all-or-nothing** — a file with
some but not all three is rejected rather than quietly falling back to symbol matching.

| Panel columns             | Small-variant filter | SV filter                                                                                |
| ------------------------- | -------------------- | ---------------------------------------------------------------------------------------- |
| `gene` only               | symbol match         | direct-hit symbol match on either breakend's VEP gene — no windows                       |
| `gene, chrom, start, end` | symbol match         | coordinate match: within 1 Mb of a breakend (BND) or 100 kb of the SV span (other types) |

Coordinate matching is what makes breakend filtering reliable: whether a BND carries a
VEP gene symbol at all depends on the sample's VEP invocation (1.6%–90% of breakends
across the samples measured), so a symbol-only panel can hide the very translocations it
exists to find. Matching on coordinates needs no annotation on the row.

Optional metadata columns (`panel`, `notes`) are ignored by the loader and kept for the
reader.

## Scoping a gene to one table: `applies_to`

By default every gene filters **both** tables. An optional `applies_to` column scopes a
row to one of them:

| `applies_to`           | Small-variant table | SV table |
| ---------------------- | ------------------- | -------- |
| _blank_, `both`, `all` | ✓                   | ✓        |
| `snv`, `small`         | ✓                   | —        |
| `sv`, `structural`     | —                   | ✓        |

Values are case-insensitive; **anything else is a hard error**, so a typo cannot quietly
change what is filtered. A file with no `applies_to` column behaves exactly as before.

This exists because the two tables want different genes. A rearrangement partner — `IGH`,
`IGK`, `IGL`, `TRA/D`, `TRB`, `TRG` — is the point of an SV panel and noise in a
small-variant one: it is not a coding-mutation target, and VEP will annotate ordinary
variants across a 1.3 Mb locus. A coding-mutation gene such as `MYD88` or `NOTCH1` is the
reverse. Putting both in one flat panel filters each table on genes that do not belong
in it.

Coordinates are therefore required **only for rows that can match an SV** (blank or `sv`).
An `snv` row may leave `chrom`/`start`/`end` empty, because symbol matching is all it ever
uses. The column-level all-or-nothing rule above is unchanged — `chrom` without `end` is
still rejected — and a blank coordinate on a blank-or-`sv` row is still a hard error.

A panel whose every row is `snv`-scoped has nothing to match SVs on, so it reports no
coordinate mode however many coordinate columns it carries.

The column reaches the report through the tool, so a panel that uses it needs
`lrsomatic_report` ≥ 1.6.0 — the tag pinned in
[`modules/local/lrsomaticreport/main.nf`](../../modules/local/lrsomaticreport/main.nf).

## Reference declaration

Panel coordinates are only valid for the reference they were built on — matching an hg38
panel against a T2T sample produces wrong hits with no error anywhere. A
coordinate-carrying panel must therefore declare its reference, either as a leading
comment line:

```
# reference: hg38
gene	chrom	start	end
MYC	chr8	127735434	127742951
```

or as a `reference` column. A panel whose declared reference differs from the one the
report is rendered against is a **hard error**. A panel that declares none loads, but the
SV section footnote says "reference unverified". Symbol-only panels are
reference-agnostic and need no declaration.

Builtin panels ship one file per reference (`lymphoid.hg38.tsv`, `lymphoid.t2t.tsv`) and
are presented as a single selectable `lymphoid` entry, resolved against the detected
reference.

## Supplying a custom panel

A panel does not have to live here. `--report_gene_panel` takes file paths too, and the
pipeline stages them into the task alongside these:

```bash
nextflow run . --report_gene_panel /path/to/my_genes.tsv
nextflow run . --report_gene_panel lymphoid,/path/to/my_genes.tsv
```

The difference between the two routes: a TSV in this directory is a _builtin_, so its
`.hg38`/`.t2t` pair collapses to one entry resolved against the rendered reference, and it
is embedded in every report whether or not it was asked for. A panel given by path is
registered under its filename and embedded only because it was named.

A one-column file of symbols (with or without a `gene` header) is accepted, and gives
symbol-only matching. The report's "Custom…" textarea takes bare symbols, so it is
symbol-only too.

Several panels can be applied at once — a builtin and your own list together, say. A variant
or SV is kept if it hits any of them; `conf/modules.config` turns each comma-separated entry
into one `--gene-panel` flag.

Each panel is registered under its filename stem; two files sharing a basename both stay
selectable, the second as `<name>-custom`. Every registered panel is a checkbox in the
report, so the reader can retick them without re-rendering.

## Bundled panels

| File                      | Contents                                                                                                                                  |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `lymphoid.{hg38,t2t}.tsv` | 234 non-Hodgkin lymphoma genes, **scoped per table** — see below                                                                          |
| `sarcoma.hg38.tsv`        | 140 soft-tissue and bone sarcoma genes — tumour suppressors, amplification targets and recurrent fusion partners — GENCODE v46 gene spans |
| `sarcoma.t2t.tsv`         | the same 140 genes, spans from the CHM13v2.0 RefSeq Liftoff v5.1 annotation                                                               |

`sarcoma` carries no `applies_to` column, so every one of its genes filters both tables.

### `lymphoid`

Merged from two curated NHL lists that are deliberately not interchangeable:

- **SV list** — `NHL_Genes_{T2T,hg38}.bed`, 128 genes with per-reference spans, including
  the loci that only make sense as rearrangement partners (`IGH`, `IGK`, `IGL`, `TRA/D`,
  `TRB`, `TRG`, `DUSP22`).
- **Small-variant list** — `twist_genes.tsv`, 149 symbols with no coordinates, the panel's
  coding-mutation content; its `TSG`/`OG` class and remarks are kept in `notes`.

42 genes are on both lists and carry a blank `applies_to`; the 86 bed-only genes are `sv`
and the 106 twist-only genes are `snv`, giving 234 rows. The beds are 4-column,
headerless, CRLF and carry no trailing newline (`wc -l` reports 127 for 128 genes), and
their coordinates are already 1-based inclusive gene spans in the same convention as the
other builtins — no BED half-open conversion is applied. `TRA/D` is kept verbatim: it is
not an HGNC symbol, but it is an `sv` row and so is matched positionally.

One consequence worth stating: the 106 twist-only genes are `snv`-scoped, so a deletion
that removes `MYD88` or `NOTCH1` entirely does not appear in the panel-filtered SV table.
That is deliberate, and it is a one-cell edit per gene to change.

The two source lists live in the report tool's repository, so this pair is re-synced from
`assets/gene_lists/` in [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report)
rather than rebuilt here. Editing the panel in place is still a pipeline change; a rebuild
from the beds belongs upstream.

#### Corrections applied on top of the upstream lists

Four rows were repaired here rather than taken verbatim, so this pair is **not** a byte-for-byte
copy of the upstream files. Re-syncing means re-applying these, or fixing them upstream first:

- `EWSR1`, `KDM6B` and `SIK3` carried their **CHM13 coordinates in the hg38 file** — the only
  three of 128 coordinate rows whose hg38 and t2t cells were identical, off by +463 kb, −96 kb
  and +15 kb against GENCODE v46. All three are `sv`-scoped, so they are matched on position
  alone: an hg38 `EWSR1` rearrangement sat well outside the 100 kb SV-span rule and was filtered
  out of a `lymphoid`-ticked SV table. Replaced with the GENCODE v46 gene spans.
- `PRKBC` was a transposition of `PRKCB`. With no coordinates on an `snv` row, symbol matching is
  all it has, so the row matched nothing and `PRKCB` — in the panel before this rebuild — was
  silently absent.
- `RCK` is an obsolete alias for `DDX6`, and its coordinates are DDX6's. Renamed, so the
  `panel_hit` cell names a symbol that VEP output can actually carry.
- `FAM46C` was dropped: it is the former name of `TENT5C`, which is already in the list, so the
  row was a dead duplicate. This is what takes the count from 235 rows to 234.

Regenerating the `sarcoma` spans is mechanical — gene spans keyed on `gene_name`, taken from
`gene` features (GENCODE) or the min/max of `transcript` features (Liftoff, which has no
`gene` feature), restricted to `chr1`–`chr22`, `chrX`, `chrY`:

- hg38: `references/GRCh38.alt-masked-V2/annotation/gencode.v46.basic.annotation.gtf.gz`
- t2t: `references/chm13_v2.0_maskedY.rCRS/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gtf`

Neither annotation is keyed on current HGNC symbols throughout, so aliases in the source
gene list were mapped by hand. For `sarcoma`, from an input list of 151 lines:
`VEGFR2`→`KDR`, `VEGFR3`→`FLT4`, `MKL2`→`MRTFB`, `MGEA5`→`OGA`, plus `HER2`→`ERBB2`,
`SYT`→`SS18`, `H3F3A`→`H3-3A` and `H3F3B`→`H3-3B`, whose targets were already listed —
those merged, as did seven verbatim duplicates, leaving 140 genes.

Two `sarcoma` symbols need the T2T annotation handled specially, and both are recorded in
that file's comment header:

- `POU2AF3` — the Liftoff annotation predates the rename and carries it as `COLCA2`.
- `DUX4L10` — the Liftoff annotation has no D4Z4 paralogs at all (only `DUX4` itself), so
  this span comes from `GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz`, whose chromosomes
  are NCBI accessions (`NC_060925.1` = `chr1` … `NC_060948.1` = `chrY`). Its 10q26
  position agrees with the hg38 locus, so this is a genuine match rather than a guess.

On a row that can match an SV (blank or `sv` in `applies_to`), a missing coordinate is a
**hard error** in `load_gene_panel()`, not a per-row downgrade to symbol matching — so a
gene that resolves in one annotation and not the other has to be either mapped or dropped
from that reference's file. Only an `snv`-scoped row may be left blank, and that is a
declaration that the gene is never matched on position, not a fallback for a lookup that
failed.
