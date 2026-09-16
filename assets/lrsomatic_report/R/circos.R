suppressPackageStartupMessages({
  library(circlize)
  library(data.table)
})

# Colour palettes — keep in sync with --circos-* in assets/styles/report.scss

# SBS-6 SNV palette (SigProfiler/COSMIC standard, softened slightly toward report ink/paper)
SNV_COLOURS = c(
  "C>A" = "#2EBAED",
  "C>G" = "#1b1e22",
  "C>T" = "#b3402f",
  "T>A" = "#c7c2b8",
  "T>C" = "#ADCC54",
  "T>G" = "#F0D0CE"
)

# SV_COLOURS and SV_YPOS are defined in R/parse_severus.R, next to the code that stamps them
# onto each row — see the note there on why there is only one copy.

# CNV colours — tied to the report spine
CNV_COLOURS = c(
  major = "#b3402f",  # brick = "more"
  minor = "#0d5c75",  # teal  = "less"
  total = "#1b1e22"   # ink
)

# BND/translocation link colour
BND_COLOUR = "#8a5fa3"

# Classify SNV into 6 SBS categories (C/T-ref normalised)
.classify_mut = function(ref, alt) {
  comp = c(A = "T", T = "A", C = "G", G = "C")
  ref = toupper(ref); alt = toupper(alt)
  use_comp = !(ref %in% c("C", "T"))
  norm_ref = ifelse(use_comp, comp[ref], ref)
  norm_alt = ifelse(use_comp, comp[alt], alt)
  paste0(norm_ref, ">", norm_alt)
}

# Draw a circos plot to output_path (SVG or PNG by extension)
draw_circos = function(snv_data = NULL,
                       sv_nontrans = NULL,
                       sv_trans = NULL,
                       cnv_data = NULL,
                       cytobands,
                       chrom_lengths,
                       chromosomes,
                       output_path) {

  # Filter cytobands and lengths to displayed chromosomes
  cyto_filt = cytobands[cytobands$chrom %in% chromosomes, ]
  lens_filt  = chrom_lengths[names(chrom_lengths) %in% chromosomes]
  lens_filt  = lens_filt[chromosomes[chromosomes %in% names(lens_filt)]]

  # Prepare SNV data
  if (!is.null(snv_data) && nrow(snv_data) > 0) {
    snv = as.data.table(snv_data)[nchar(ref) == 1 & nchar(alt) == 1]
    snv = snv[chrom %in% chromosomes]
    snv[, mut_cat := .classify_mut(ref, alt)]
    snv[, circos_col := SNV_COLOURS[mut_cat]]
    snv[is.na(circos_col), circos_col := "#AAAAAA"]
  } else {
    snv = data.table(chrom = character(), pos = integer(),
                     mut_cat = character(), circos_col = character())
  }

  # Prepare SV (non-BND) data
  if (!is.null(sv_nontrans) && nrow(sv_nontrans) > 0) {
    sv_nt = as.data.table(sv_nontrans)[chrom %in% chromosomes]
  } else {
    sv_nt = data.table(chrom = character(), pos = integer(), end = integer(),
                       svtype = character(), circos_pos = numeric(), circos_col = character())
  }

  # Prepare translocation (BND) data
  if (!is.null(sv_trans) && nrow(sv_trans) > 0) {
    sv_tr = as.data.table(sv_trans)[chrom %in% chromosomes & chrom2 %in% chromosomes]
  } else {
    sv_tr = data.table(chrom = character(), pos = integer(),
                       chrom2 = character(), pos2 = integer())
  }

  # Prepare CNV data
  if (!is.null(cnv_data) && nrow(cnv_data) > 0) {
    cnv = as.data.table(cnv_data)[chr %in% chromosomes]
    cnv = cnv[order(chr, startpos)]
  } else {
    cnv = data.table(chr = character(), startpos = integer(), endpos = integer(),
                     major_cn = numeric(), minor_cn = numeric(), total_cn = numeric())
  }

  # Open device
  ext = tolower(tools::file_ext(output_path))
  if (ext == "svg") {
    svglite::svglite(output_path, width = 8, height = 8)
  } else {
    png(output_path, width = 2400, height = 2400, res = 300)
  }

  plot.new()
  circos.clear()

  n_chr = length(chromosomes)
  gap_degrees = c(rep(1.5, n_chr - 1), 7)

  # Single quiet track surface; colours match the report.scss border/surface tokens
  track_bg     = "#fbfaf6"
  track_border = "#e4e0d6"

  circos.par(
    "start.degree" = 90,
    "gap.degree"   = gap_degrees,
    "track.margin" = c(0.006, 0.006),
    "cell.padding" = c(0, 0, 0, 0)
  )

  # Build cytobands list as expected by circos.initializeWithIdeogram
  cyto_list = list(
    df          = cyto_filt,
    chromosome  = chromosomes[chromosomes %in% unique(cyto_filt$chrom)],
    chr.len     = lens_filt
  )

  circos.initializeWithIdeogram(cyto_list$df,
                                chromosome.index = cyto_list$chromosome,
                                plotType         = c("ideogram", "labels"),
                                labels.cex       = 0.8)

  # Pre-compute jitter once so it varies per chromosome but stays reproducible
  set.seed(42)

  # ---- Track 1: SNV dots (coloured by mutation category) ------------------
  circos.trackPlotRegion(
    factors      = chromosomes,
    ylim         = c(0, 1),
    bg.border    = track_border,
    bg.col       = track_bg,
    track.height = 0.16,
    panel.fun    = function(region, value, ...) {
      chr = get.cell.meta.data("sector.index")
      sub_snv = snv[chrom == chr]
      if (nrow(sub_snv) == 0) return(invisible(NULL))
      y_jitter = runif(nrow(sub_snv), 0.05, 0.95)
      # Translucent so a dense cloud reads as a tint, not confetti.
      circos.points(
        x   = sub_snv$pos,
        y   = y_jitter,
        col = adjustcolor(sub_snv$circos_col, alpha.f = 0.65),
        pch = 19,
        cex = 0.18
      )
    }
  )

  # ---- Track 2: Non-BND SVs (DEL/DUP/INV/INS as horizontal segments) -----
  circos.trackPlotRegion(
    factors      = chromosomes,
    ylim         = c(0, 1),
    bg.border    = track_border,
    bg.col       = track_bg,
    track.height = 0.10,
    panel.fun    = function(region, value, ...) {
      chr = get.cell.meta.data("sector.index")
      sub_sv = sv_nt[chrom == chr & !is.na(circos_pos)]
      if (nrow(sub_sv) == 0) return(invisible(NULL))
      for (i in seq_len(nrow(sub_sv))) {
        x1 = sub_sv$pos[i]
        x2 = if (!is.na(sub_sv$end[i]) && sub_sv$end[i] > x1) sub_sv$end[i] else x1 + 1L
        circos.segments(
          x0  = x1, x1  = x2,
          y0  = sub_sv$circos_pos[i], y1 = sub_sv$circos_pos[i],
          col = sub_sv$circos_col[i],
          lwd = 2.5
        )
      }
    }
  )

  # Ask circlize which track that was rather than counting.
  # circos.initializeWithIdeogram() creates a track per plotType group *before* any of
  # ours -- one for axis/labels and one for the ideogram -- so the SV track is index 4,
  # not 3. A wrong track.index does not error, it silently draws on the neighbouring
  # track, which is how these labels ended up against the SNV dots.
  sv_track_index = get.current.track.index()

  # Y-axis labels for SV track. `at`/`labels` come from SV_YPOS (R/parse_severus.R) — the
  # same object severus_circos_tracks() positioned the segments from, not a second copy of
  # its values, so the two cannot drift apart.
  sv_at = sort(SV_YPOS)
  tryCatch(
    circos.yaxis(
      side              = "left",
      at                = unname(sv_at),
      labels            = names(sv_at),
      track.index       = sv_track_index,
      sector.index      = chromosomes[1],
      labels.niceFacing = TRUE,
      labels.cex        = 0.45
    ),
    error = function(e) NULL
  )

  # ---- Track 3: ASCAT copy-number -----------------------------------------
  circos.trackPlotRegion(
    factors      = chromosomes,
    ylim         = c(0, 4),
    bg.border    = track_border,
    bg.col       = track_bg,
    track.height = 0.18,
    panel.fun    = function(region, value, ...) {
      chr = get.cell.meta.data("sector.index")
      sub_cnv = cnv[chr == get.cell.meta.data("sector.index")]
      if (nrow(sub_cnv) == 0) return(invisible(NULL))

      xmax = lens_filt[chr]
      if (!is.na(xmax)) {
        for (y_ref in c(1, 2, 3)) {
          circos.lines(c(0, xmax), c(y_ref, y_ref),
                       col = track_border, lwd = 0.4, lty = "dotted")
        }
      }

      for (i in seq_len(nrow(sub_cnv))) {
        xl = sub_cnv$startpos[i]; xr = sub_cnv$endpos[i]
        maj = sub_cnv$major_cn[i]
        circos.rect(xl, maj + 0.02, xr, maj + 0.12,
                    col = CNV_COLOURS["major"], border = CNV_COLOURS["major"], lwd = 0.05)
        min_cn = sub_cnv$minor_cn[i]
        circos.rect(xl, min_cn - 0.12, xr, min_cn - 0.02,
                    col = CNV_COLOURS["minor"], border = CNV_COLOURS["minor"], lwd = 0.05)
        tot = sub_cnv$total_cn[i]
        circos.rect(xl, tot - 0.03, xr, tot + 0.03,
                    col = CNV_COLOURS["total"], border = CNV_COLOURS["total"], lwd = 0.05)
      }
    }
  )

  # Y-axis labels for the copy-number track. Drawn once, outside panel.fun: in there it
  # ran per sector with a fixed sector.index, so it was re-drawn once for every
  # chromosome carrying ASCAT segments (~23 overlapping copies, which is why it rendered
  # heavier than the SV axis). The empty case keeps its old behaviour deliberately: the
  # track is created unconditionally, but parse_ascat_segments() returns NULL when there is
  # no ASCAT output at all, and a 0..4+ scale against a blank ring reads as "copy number is
  # zero everywhere" rather than "no copy-number data".
  cnv_track_index = get.current.track.index()
  if (nrow(cnv) > 0) tryCatch(
    circos.yaxis(
      side              = "left",
      at                = c(0, 1, 2, 3, 4),
      labels            = c("0", "1", "2", "3", "4+"),
      track.index       = cnv_track_index,
      sector.index      = chromosomes[1],
      labels.niceFacing = TRUE,
      labels.cex        = 0.40
    ),
    error = function(e) NULL
  )

  # ---- Translocation links (BND): one arc per mate-collapsed rearrangement ----
  if (nrow(sv_tr) > 0) {
    for (i in seq_len(nrow(sv_tr))) {
      tryCatch(
        circos.link(
          sector.index1 = sv_tr$chrom[i],  point1 = sv_tr$pos[i],
          sector.index2 = sv_tr$chrom2[i], point2 = sv_tr$pos2[i],
          col = adjustcolor(BND_COLOUR, alpha.f = 0.45),
          lwd = 0.9
        ),
        error = function(e) NULL
      )
    }
  }

  circos.clear()
  dev.off()
  invisible(output_path)
}
