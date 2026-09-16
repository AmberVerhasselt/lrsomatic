# Phasing section (WhatsHap stats from qc/whatshap_stats/); these are germline statistics

register_section(list(
  id = "whatshap",
  title = "Phasing",

  locate = function(sample_dir, sample_id) {
    d = sample_dir

    find1 = function(pattern) {
      hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) hits[1] else NULL
    }

    # qc/whatshap_stats/ is not tumour/normal-scoped, so a plain recursive match is correct
    list(stats_tsv = find1("_whatshap_stats\\.tsv$"))
  },

  parse = function(inputs, section_data) {
    f = inputs$stats_tsv
    if (is.null(f) || !file.exists(f)) return(NULL)

    dt = tryCatch(
      fread(f, sep = "\t", header = TRUE),
      error = function(e) {
        message("Failed to parse WhatsHap stats: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(dt) || nrow(dt) == 0) return(NULL)

    # The header line is "#sample\tchromosome\t..." — fread keeps the leading "#".
    setnames(dt, sub("^#", "", names(dt)))
    if (!"chromosome" %in% names(dt)) {
      message("WhatsHap stats has no 'chromosome' column; skipping section")
      return(NULL)
    }

    # bp_per_block_sum reads as integer64, which DT renders badly; widen to double
    for (col in names(dt)) {
      if (inherits(dt[[col]], "integer64")) dt[, (col) := as.numeric(get(col))]
    }

    list(
      per_chrom = dt[chromosome != "ALL"],
      all       = if (any(dt$chromosome == "ALL")) as.list(dt[chromosome == "ALL"][1]) else NULL,
      vcf       = if ("file_name" %in% names(dt)) dt$file_name[1] else NA_character_
    )
  }
))

# Genome-wide totals for the header card: whatshap's own ALL row when it wrote one, otherwise
# the per-chromosome rows summed. NULL when there is nothing to total, which is what lets the
# card be omitted rather than render a permanent "N/A". Every field is read through a names()
# guard: `w$phased_fraction` is NULL when the column is absent, and `is.na(NULL)` is logical(0),
# which errors rather than falling back.
whatshap_totals = function(whatshap) {
  if (is.null(whatshap)) return(NULL)

  num1 = function(x, key) {
    if (is.null(x) || !key %in% names(x)) return(NA_real_)
    v = suppressWarnings(as.numeric(x[[key]]))
    if (length(v) == 0) NA_real_ else v[1]
  }

  all_row = whatshap$all
  if (!is.null(all_row)) {
    phased   = num1(all_row, "phased")
    variants = num1(all_row, "variants")
    fraction = num1(all_row, "phased_fraction")
  } else {
    # No ALL row: sum the per-chromosome rows rather than go quiet under a healthy table
    pc = whatshap$per_chrom
    if (is.null(pc) || nrow(pc) == 0) return(NULL)
    # na.rm = TRUE makes an all-missing column sum to 0, which then survives the
    # is.na(phased) && is.na(fraction) guard below and renders a confident "0" in the header
    # card. whatshap leaves the field empty for a chromosome with no het variants, so a file
    # where that is true of every row is the case this distinguishes.
    sum1 = function(key) {
      if (!key %in% names(pc)) return(NA_real_)
      v = suppressWarnings(as.numeric(pc[[key]]))
      if (all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)
    }
    phased   = sum1("phased")
    variants = sum1("variants")
    fraction = NA_real_
  }

  # phased_fraction is a 0-1 fraction; recompute when the column is absent or the rows were summed
  if (is.na(fraction) && !is.na(phased) && !is.na(variants) && variants > 0)
    fraction = phased / variants
  if (is.na(phased) && is.na(fraction)) return(NULL)

  list(phased = phased, variants = variants, fraction = fraction)
}
