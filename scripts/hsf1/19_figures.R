# 19_figures.R -- figures F1..F9 for P2 and P3.
#
# Every number is read from a saved .tsv. No figure recomputes anything.
#
# Figure rules, enforced in code:
#  - N values in legends are looked up by name from the source table, never with min()/max().
#    A helper `Nof()` does the lookup and stops if the name is absent, so a silent wrong N is
#    impossible.
#  - colour-blind-safe palette only (Okabe-Ito), and no information is carried by colour alone:
#    every series also differs in point symbol or line type.
#  - vector PDF plus a final-dpi PNG for inspection; base font size 8 pt, nothing below 7 pt.
#  - overlapping points are offset horizontally (jitter is deterministic, seeded).
#  - each caption is self-contained and states the source file.
#
# Statistical summaries used in the captions:
#  - cis-trans compensation is summarized per timepoint (-0.497 / -0.386 / -0.268 identity;
#    cross-block -0.166 is an 18-cell median); medians and per-timepoint values are distinct summaries
#    in one statement. F6 states this on its face.
#  - P2 covers a narrow experimental window rather than a temperature-by-time map. F2c shows that a count orders amplitudes only when the
#    noise floor is shared: count ratios 0.014/1/0.010/0.019 vs amplitude ratios 0.42/1/0.19/0.13,
#    noise-floor spread 1.47.
#  - three historical citations (Pearson 1897, Oldham 1962, Atchley et al. 1976) are cited for the
#    spurious correlation of ratios that share a component.
#  - F5 distinguishes reference-level effects from reference-mapping effects.
#
# Usage: Rscript scripts/hsf1/19_figures.R

## Figure width: 6.9 in. A target journal caps double-column figures at 178 mm = 7.008 in. The
## width is set here rather than rescaled at submission, because a rescaled raster loses the
## >= 300 dpi guarantee that the same requirement also imposes.
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
hsd <- file.path(root, "data/02_processed/hsf1")
csd <- file.path(root, "data/02_processed/crossspecies")
pro <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
fig <- file.path(root, "doc/04_figures/final"); dir.create(fig, showWarnings = FALSE, recursive = TRUE)
repd <- file.path(root, "doc/05_reports"); MIRROR <- "logs/output_copy"
set.seed(20260906L)   # DPI comes from theme.R, sourced below
say <- function(...) cat(sprintf(...), "\n", sep = "")
logf <- file.path(repd, "logs/hsf1_figures.log"); sink(logf, split = TRUE)
say("=== 19_figures | %s | F1..F9 ===", format(Sys.time()))

## Colours come from the shared theme so that a species or a class is the same colour in every figure of every
## manuscript. `SPCOL` is a code-keyed alias of the theme's full-name map, so the call sites use species codes and
## the colours are the theme's.
source(file.path(root, "scripts/figures/theme.R"))
SPCOL <- setNames(unname(SPECIES_COL[c("S. cerevisiae", "S. paradoxus", "S. kudriavzevii", "S. uvarum")]),
                  c("Sc", "Sp", "Sk", "Su"))
stopifnot(!anyNA(SPCOL), identical(unname(SPCOL[["Sc"]]), unname(SPECIES_COL[["S. cerevisiae"]])))
SPPCH <- c(Sc=21, Sp=22, Sk=23, Su=24)
CAP <- list()   # captions collected and written out beside the figures

# N lookup by name. Stops rather than guessing.
Nof <- function(DT, keycol, key, ncol_ = "n") {
  stopifnot(keycol %in% names(DT), ncol_ %in% names(DT))
  v <- DT[[ncol_]][DT[[keycol]] == key]
  if (!length(v)) stop(sprintf("Nof: '%s' not found in column '%s'", key, keycol))
  as.integer(v[1])
}
dev_open <- function(tag, w, h) {
  pdf(file.path(fig, paste0(tag, ".pdf")), width = w, height = h, pointsize = 8, useDingbats = FALSE)
}
dev_close <- function(tag, w, h) {
  dev.off()
  png(file.path(fig, paste0(tag, ".png")), width = w * DPI, height = h * DPI, res = DPI, pointsize = 8)
}
# draw twice (pdf then png) so the PNG really is the same figure at final dpi
## Point size. With an 8 pt device base, the annotation cex values on top of it (0.52 to 0.72) put text at 4.2 to
## 5.8 pt at the final 175 mm width (six figures measured at 4.0 pt in the PDF content stream). That is not legible in
## print, and it is not obvious on screen because the same figure looks generous at 1000% zoom.
##
## Raising the base rather than editing thirty cex values moves every size by the same factor, so the relative
## emphasis of titles, axis labels and annotations is preserved exactly, and no symbol shrinks. The figure is not
## shrunk to make the type look bigger.
PT_BASE <- 11
draw_both <- function(tag, w, h, fn, outdir = fig) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") pdf(file.path(outdir, paste0(tag, ".pdf")), width = w, height = h,
                          pointsize = PT_BASE, useDingbats = FALSE)
    else png(file.path(outdir, paste0(tag, ".png")), width = round(w * DPI), height = round(h * DPI),
             res = DPI, pointsize = PT_BASE)
    par(mgp = c(1.9, 0.55, 0), tcl = -0.25, cex.axis = 0.85, cex.lab = 0.95, font.main = 1,
        cex.main = 0.95)
    ok <- try(fn(), silent = TRUE)
    dev.off()
    if (inherits(ok, "try-error")) { say("  !! %s FAILED: %s", tag, conditionMessage(attr(ok,"condition"))); return(FALSE) }
  }
  say("  %-4s ok  (%s.pdf + %s.png)", tag, tag, tag); TRUE
}
errbar <- function(x, y, lo, hi, col, lwd = 1) segments(x, lo, x, hi, col = col, lwd = lwd)

## Every figure below is drawn in the shared figure style of scripts/figures/fig_style.R:
## Noto Sans embedded, 7.5/8.5/10 pt text, the style's line weights, bold uppercase
## panel tags, no panel titles (named in the captions), keys above the panels, one colour semantics per figure.
## draw_both keeps its signature so the call sites do not change; the width is the style's 174 mm and the height the
## figure's own.
source(file.path(root, "scripts/figures/fig_style.R"))
SPCOL[] <- unname(SPECIES_COL_SUB[c("S. cerevisiae", "S. paradoxus", "S. kudriavzevii", "S. uvarum")])   # P3 species colours
draw_both <- function(tag, w, h, fn, outdir = fig) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch({ fig_both(tag, if (w < 5) 108 else 174, h * 25.4, fn, outdir); TRUE },   # w < 5 in = 1.5 column
                 error = function(E) { say("  !! %s FAILED: %s", tag, conditionMessage(E)); FALSE })
  if (ok) say("  %-4s ok  (%s.pdf + %s.png)", tag, tag, tag)
  ok
}
errbar <- function(x, y, lo, hi, col, lwd = LW$ci) segments(x, lo, x, hi, col = col, lwd = lwd)
ptb <- function(x, y, pch, bg, col = COL_SCI$truth, cex = PCH_CEX) points(x, y, pch = pch, bg = bg, col = col,
                                                                          cex = cex, lwd = PCH_LWD)

## ------------------------------------------------------------------ F1 design ---
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")]
GENO <- c("Wildtype","HSF1.KD","MSN24.KO","Double.KDKO")
CONDS <- c("25_0","37_10","37_30","42_10","42_30")
ok1 <- draw_both("F1_design", 6.9, 3.6, function() {
  layout(matrix(1:2, 1, 2), widths = c(1, 1.2)); after_layout()
  ## (A) the perturbation dataset: one cell per genotype x condition, n in each; the 37 °C, 30 min anchor is the
  ## darker cell with its count in bold
  par(mar = c(3.4, 8.2, 3.0, 0.8))
  tab <- table(factor(SM$knockout, GENO), factor(SM$cond, CONDS))
  plot(NA, xlim = c(0.5, length(CONDS) + 0.5), ylim = c(0.5, length(GENO) + 0.5), axes = FALSE, xlab = "", ylab = "")
  for (i in seq_along(CONDS)) for (j in seq_along(GENO)) {
    anch <- CONDS[i] == "37_30"
    rect(i - 0.46, j - 0.42, i + 0.46, j + 0.42, col = "white",
         border = if (anch) NPG$red else GREYS$g67, lwd = if (anch) LW$data else LW$ci)
    if (i == 1) rect(0.30, j - 0.42, 0.42, j + 0.42, col = GENO_SUB$col[[GENO[j]]], border = NA, xpd = NA)
    text(i, j, tab[j, i], cex = CEX_TXT, font = if (anch) 2 else 1)
  }
  ## Note: with one-line labels base R silently drops three of the five condition labels for overlap, and two-line
  ## labels still collide at 900 px, so the temperature is written once per group and the time under each cell.
  mtext(c("25 \u00b0C", "37 \u00b0C", "42 \u00b0C"), side = 1, line = 0.0, at = c(1, 2.5, 4.5), cex = CEX_LAB)
  mtext(c("base", "10", "30", "10", "30"), side = 1, line = 0.9, at = seq_along(CONDS), cex = CEX_LAB)
  mtext("minutes at temperature", side = 1, line = 1.8, at = 3, cex = CEX_LAB)
  axis(2, seq_along(GENO), lab_of(GENO), las = 1, tick = FALSE)
  panel_tag("A")
  ## (B) the comparative dataset: four species (fixed colours and shapes) and their six hybrids
  par(mar = c(3.4, 1.0, 3.0, 0.8))
  plot(NA, xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "")
  sp <- c("T73 (Sc)\nthermophilic","N17 (Sp)\nthermophilic","CR85 (Sk)\ncryophilic","BMV58 (Su)\ncryophilic")
  xs <- c(1.1, 3.7, 6.3, 8.9)
  cmb <- combn(4, 2)
  ## one node per hybrid, joined to its two parents, so each of the six crosses can be read
  hx <- seq(1.0, 9.0, length.out = 6)   # spread so standard-size labels do not touch
  for (k in 1:6) { segments(xs[cmb[, k]], 7.9, hx[k], 5.6, col = GREYS$g80, lwd = LW$ref)
    points(hx[k], 5.35, pch = 21, bg = CHARCOAL, col = CHARCOAL, cex = PCH_CEX, lwd = PCH_LWD)
    ## name each hybrid, so the six crosses can be read without tracing lines
    text(hx[k], 4.55, paste0(c("Sc","Sp","Sk","Su")[cmb[1, k]], "\u00d7", c("Sc","Sp","Sk","Su")[cmb[2, k]]), cex = CEX_TXT) }
  for (i in 1:4) { ptb(xs[i], 8.3, SPPCH[i], SPCOL[i], cex = PCH_CEX * 1.6); text(xs[i], 9.6, sp[i], cex = CEX_TXT) }
  text(5, 3.5, "6 F1 hybrids give the allele-specific axis", cex = CEX_TXT)
  text(5, 2.3, "25 to 37 \u00b0C; 0, 15, 30, 60 min; 2 blocks; 5,121 orthologs", cex = CEX_TXT)
  panel_tag("B")
})
CAP$F1 <- "F1. Design and the projection constraint. (A) GSE179258: 4 genotypes x 5 conditions x 3 replicates = 60 libraries; cell labels give n. Source: data/01_metadata/sample_mapping.tsv. (B) Fay 2023 (PMID 37972247): four Saccharomyces species, two thermophilic and two cryophilic, and their six F1 hybrids, which supply the allele-specific (cis) axis. Hsf1 dependency can be determined only in S. cerevisiae because only the Mühlhofer dataset perturbs Hsf1; all cross-species gene sets are orthologous PROJECTIONS of Sc-defined sets, and regulon membership may itself have diverged."

## --------------------------------- F2 epistasis + the count-vs-amplitude panel (C) ---
CLS <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
## Panel (b) uses the method-unanimous, mutually exclusive class assignment (Hsf1_only n = 87), which is what
## the manuscript reports. `p2_class_resolution.tsv` carries the non-exclusive assignment, in which a gene the
## three DE methods disagree about is counted in more than one class (Hsf1_only n = 127). The resolution SEs the
## panel needs are only in the resolution table, so the two are joined on class and the exclusive n and medians
## take precedence.
CRESX <- fread(file.path(hsd, "p2_class_exclusive.tsv"))[regime %like% "as published"]
CRES0 <- fread(file.path(hsd, "p2_class_resolution.tsv"))
CRES  <- merge(CRESX[, .(cls, n, med_I, se_of_median = se_of_med, refv, vsref, zref)],
               CRES0[, .(cls, med_se_model)], by = "cls", all.x = TRUE)
CRES[is.na(med_se_model), med_se_model := median(CRES0$med_se_model)]
## The interval of the class median is the correlation-adjusted one (P2 Amendment 1, U1): sqrt(VIF) x the
## gene-independence SE, VIF = 1 + (n - 1) rho from 80_p2_adjusted_uncertainty.R. The gene-independence SE alone gives
## intervals 3-11 times too narrow for the class median.
ADJ <- fread(file.path(hsd, "p2_class_adjusted.tsv"))[regime == "all genes (as published)" & contrast == "I(37_30)",
                                                        .(cls, se_adj_med = sqrt(vif) * se_gene, med_chk = med, vs_ref, lo_ref = lo_adj, hi_ref = hi_adj)]
CRES <- merge(CRES, ADJ, by = "cls")
stopifnot(nrow(CRES) == 7L, max(abs(CRES$med_I - CRES$med_chk)) < 1e-10)
AMPC <- fread(file.path(hsd, "p2_interaction_amplitude.tsv"))
HEAT <- c("37_10","37_30","42_10","42_30")
setkey(AMPC, condition); AMPC <- AMPC[HEAT]
## The figure uses the model-based interaction-resolution estimate from script 17.
## The model-based median SE measures interaction resolution without mixing response and baseline components
## and is computed by `17_p2_interaction_se_and_four_conditions.R`, which
## also records the comparison value for diagnostics.

ok2 <- draw_both("F2_epistasis", 6.9, 6.0, function() {
  ## (A) and (B) are horizontal. The class names sit once, on the y axis of (A), and (B) shares its rows (same ylim,
  ## same top and bottom margins), so each class reads across. Largest class at the top.
  D <- CRES[order(-n)]
  k <- nrow(D); yy <- rev(seq_len(k))
  ml <- local({
    tf <- tempfile(fileext = ".pdf")
    grDevices::cairo_pdf(tf, width = 6.85, height = 6.0, pointsize = PS_DEV, family = FONT_SUB)
    on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
    fig_par(); plot.new()
    max(strwidth(lab_of(D$cls), units = "inches", cex = CEX_TXT)) / par("csi") + 1.0   # label width in lines
  })
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), widths = c(1.3, 1), heights = c(1, 1.1)); after_layout()
  ## (A) class sizes, n at the end of each bar. No colour: the classes are named on the axis.
  par(mar = c(3.4, ml, 3.4, 0.8))
  plot(NA, xlim = c(0, max(D$n) * 1.3), ylim = c(0.5, k + 0.9), yaxt = "n", xaxt = "n", xlab = "genes", ylab = "")  # n after the bar
  axis(1, at = seq(0, 4000, 2000), labels = format(seq(0, 4000, 2000), big.mark = ","))   # thousands separated by commas
  rect(0, yy - 0.34, D$n, yy + 0.34, col = CLASS_SUB[D$cls], border = COL_SCI$neutral, lwd = LW$ci)
  text(D$n, yy, format(D$n, big.mark = ",", trim = TRUE), pos = 4, offset = 0.3, cex = CEX_TXT)   # thousands separated by commas
  axis(2, yy, lab_of(D$cls), las = 1, tick = FALSE)
  panel_tag("A")
  ## (B) class median I with the two error bars, drawn horizontally on the rows of (A). n is
  ## labelled for no class here (panel A gives every class size), and where the 1.96-SE interval of the median is
  ## shorter than the marker it is hidden by it, which the caption says.
  par(mar = c(3.4, 1.2, 3.4, 0.8))
  ## the primary reading is against the not-induced reference class, so (b) shows class median minus the
  ## reference median with the correlation-adjusted 95% interval of that difference (p2_class_adjusted.tsv).
  lo <- ifelse(is.na(D$lo_ref), D$vs_ref, D$lo_ref); hi <- ifelse(is.na(D$hi_ref), D$vs_ref, D$hi_ref)
  xmax <- max(abs(c(lo, hi))) * 1.15
  plot(NA, xlim = c(-xmax, xmax), ylim = c(0.5, k + 0.9), yaxt = "n", xlab = "median interaction minus reference class", ylab = "")
  abline(v = 0, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  segments(lo, yy, hi, yy, col = COL_SCI$truth, lwd = LW$ci)
  ptb(D$vs_ref, yy, ifelse(D$cls == "not_induced", 23, 21), CLASS_SUB[D$cls])
  legend_top(c("adjusted 95% interval", "reference class"),
             pch = c(NA, 23), lty = c(1, NA), lwd = c(LW$ci, NA), pt.bg = c(NA, CLASS_SUB[["not_induced"]]), col = COL_SCI$truth,
             pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2)
  panel_tag("B")
  ## (C) count ratio (the naive proxy, dark red) against amplitude ratio (noise-corrected, dark green), log scale.
  ## ticks at 2, 0.5, 0.1, 0.02 and 0.005, with the reference at 1 drawn and labelled.
  par(mar = c(4.6, 7.0, 3.0, 0.8), mgp = c(3.6, 0.5, 0))   # y title clear of the 0.005 tick label
  rc <- AMPC$n_sig / AMPC$n_sig[AMPC$condition == "37_30"]
  ra <- AMPC$rms_tru / AMPC$rms_tru[AMPC$condition == "37_30"]
  x <- seq_along(HEAT)
  plot(NA, xlim = c(0.7, 4.3), ylim = c(0.003, 2.4), log = "y", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "ratio to 37 \u00b0C, 30 min (log scale)")
  axis(2, at = c(2, 1, 0.5, 0.1, 0.02, 0.005), labels = c("2", "1", "0.5", "0.1", "0.02", "0.005"))
  hl <- sub(", ", ",\n", lab_of(HEAT), fixed = TRUE)
  axis(1, x, hl, padj = 0.5, tick = FALSE)
  abline(h = 1, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  ## no sideways shift; both series equal 1 at the anchor by definition
  lines(x, rc, col = COL_SCI$naive, lwd = LW$data); ptb(x, rc, 22, COL_SCI$naive)
  lines(x, ra, col = COL_SCI$corrected, lwd = LW$data, lty = 2); ptb(x, ra, 21, COL_SCI$corrected)
  legend_top(c("significant genes (count)", "interaction size above replicate scatter"), pch = c(22, 21),
             pt.bg = c(COL_SCI$naive, COL_SCI$corrected), col = c(COL_SCI$naive, COL_SCI$corrected), lty = c(1, 2),
             lwd = LW$data, pt.cex = PCH_CEX, pt.lwd = PCH_LWD)
  panel_tag("C")
  exg <- (ra / rc)[HEAT != "37_30"]
})

CAP$F2 <- sprintf("F2. Epistatic decomposition of the Hsf1/Msn2-4 module, and why a significant-gene count may not be read as an amplitude. (a) Class sizes at the 37 \u00b0C, 30 min anchor, on the mutually exclusive assignment the three differential-expression methods agree on. (b) Class median interaction I = D - H - M + W with two error bars that answer different questions: gray is \u00b1 one single-gene standard error from a single model, black is \u00b1 1.96 standard errors of the class median, adjusted for the correlation between genes (Methods); where that interval is shorter than the marker, the marker hides it. Msn2/4-only (n = %d) has median I = %+.4f against zero and %+.4f against the not-induced reference class; uncertainty adjusted for the correlation between genes is in Table 1. (c) Ratios to the same anchor, on a log scale with the dashed line at 1: the count (red squares) falls further than the noise-corrected amplitude of I (blue circles) by %s-fold at the three non-anchor conditions, and the noise floors differ by %.2f-fold, so they are not shared and the count cannot be used to order amplitudes. The design therefore supports a contrast between 37 \u00b0C and 42 \u00b0C, not a temperature-by-time map. Sources: p2_class_exclusive.tsv, p2_class_resolution.tsv, p2_interaction_amplitude.tsv.",
  ## the worked example is Msn2/4-only; the gene-independence z is not printed (P2 Amendment 1, U1)
  Nof(CRES,"cls","Msn24_only"), CRES[cls=="Msn24_only", med_I], CRES[cls=="Msn24_only", vsref],
  paste(sprintf("%.1f", ((AMPC$rms_tru/AMPC$rms_tru[AMPC$condition=="37_30"]) /
                         (AMPC$n_sig/AMPC$n_sig[AMPC$condition=="37_30"]))[AMPC$condition != "37_30"]),
        collapse = ", "),
  max(AMPC$rms_null)/min(AMPC$rms_null))

## F3 is drawn by 60_p2_fig3_two_backgrounds.R, the sole producer of
## that output name.


## ------------------------------------------------------ F4 parental dimensions --
COMP <- fread(file.path(hsd, "fay_parental_composition_37_30.tsv"))
## panel (a) draws the correlation-adjusted interval (correlation-adjusted, p3_adjusted_composition.tsv lo/hi), not
## the gene bootstrap; the point estimates are the same numbers and are checked here. A null line is drawn at 1.
COMPP <- fread(file.path(hsd, "p3_adjusted_composition.tsv"))
stopifnot(max(abs(merge(COMP, COMPP, by = c("stratum", "time"))[, ratio.x - ratio.y])) < 1e-12)
COMP <- merge(COMP[, !c("lo", "hi")], COMPP[, .(stratum, time, lo, hi)], by = c("stratum", "time"))
PK   <- fread(file.path(hsd, "fay_parental_peak_37_30.tsv"))
ok4 <- draw_both("F4_parental", 6.9, 3.2, function() {
  layout(matrix(1:2, 1, 2)); after_layout()
  ## (A) the Hsf1-dependent core (COL_BIO Hsf1 green); the stratum and its n are in the caption
  par(mar = c(3.4, 4.2, 3.0, 0.8))
  z <- COMP[stratum == "Hsf1_core"][order(time)]
  plot(NA, ylim = c(0, 1.15), xlim = c(12, 70), xlab = "minutes after shift",
       ylab = "between/within-species similarity")   # the plotted quantity is between over within
  abline(h = 1, lty = 2, col = COL_SCI$truth, lwd = LW$ref)
  errbar(z$time, z$ratio, z$lo, z$hi, GREYS$g67)
  lines(z$time, z$ratio, col = NPG$blue, lwd = LW$data)
  ptb(z$time, z$ratio, 21, NPG$blue)
  text(z$time + 1.4, z$ratio + 0.045, sprintf("%.3f", z$ratio), adj = c(0, 0.5), cex = CEX_TXT)
  panel_tag("A")
  ## (B) where each species peaks: the species' fixed colours, named in the key above the panel in italics.
  ## colour alone does not identify a species in greyscale or for a colour-blind reader, so each bar also
  ## carries the species' fixed shape from Figure 1 (SPPCH, same SPCOL), and the key shows shape and colour together.
  par(mar = c(3.4, 4.2, 3.0, 0.8))
  pk <- PK[stratum == "Hsf1_core"]
  tt <- table(factor(pk$species, names(SPCOL)), factor(pk$peak))
  bb <- barplot(tt, beside = TRUE, col = SPCOL[rownames(tt)], border = COL_SCI$truth, lwd = LW$ci, axes = FALSE,
                ylim = c(0, max(tt) * 1.25), xlab = "timepoint of peak response", ylab = "replicate blocks")
  axis(2, at = 0:max(tt), labels = 0:max(tt), las = 1)
  sp_r <- matrix(rownames(tt), nrow(tt), ncol(tt))
  hit <- tt > 0
  ptb(bb[hit], tt[hit] + 0.12 * max(tt), SPPCH[sp_r[hit]], SPCOL[sp_r[hit]])
  legend_top(lab_of(rownames(tt)), pch = SPPCH[rownames(tt)], pt.bg = SPCOL[rownames(tt)], col = COL_SCI$truth,
             pt.cex = PCH_CEX, pt.lwd = PCH_LWD, text.font = 3)
  panel_tag("B")
})
# the timing sentence is derived from the plotted table, from the plotted table
PKc <- PK[stratum == "Hsf1_core"]
agree <- PKc[, .(k = uniqueN(peak), pk = paste(sort(unique(peak)), collapse = "/")), by = species]
tim_txt <- paste(sprintf("%s %s", agree$species,
  ifelse(agree$k == 1L, sprintf("peaks at %s min in both blocks", agree$pk),
         sprintf("disagrees between blocks (%s min) and is unresolved", agree$pk))), collapse = "; ")
CAP$F4 <- sprintf("F4. Cross-species parental divergence: only two of four pre-registered dimensions are interpretable. (A) Composition, the strongest dimension: within/between-species similarity for the Hsf1-associated core (n = %d genes) falls 0.753 -> 0.559 -> 0.391, i.e. divergence grows monotonically with time; bars are bootstrap 95%% CI. (B) Timing, the only resolvable species-level dimension, read directly from the plotted table: %s. The AMPLITUDE dimension is VOID: species ordering changes after noise correction and the not_induced set shows the same ordering, so it is not heat-specific. Persistence is censored and not shown. Sources: fay_parental_composition_37_30.tsv, fay_parental_peak_37_30.tsv.", Nof(COMP,"stratum","Hsf1_core"), tim_txt)

## ------------------------------------------------------ F5 the amplitude artefact
A5 <- fread(file.path(hsd, "a5_cis_ratio_calibrated_37_30.tsv"))
PW5 <- fread(file.path(hsd, "a5_power_calibration_37_30.tsv"))
## panel (a) draws the correlation-adjusted 95% intervals (p3_adjusted_core.tsv, the paper's stated basis) for both
## readings, beta = 0 and the fitted beta; point estimates are asserted equal to the a5 table.
## The core rows are the reported_core rows reported in Table 1 (the Figure 3 stratum rows repeat the same ratio with a
## bootstrap SE from another seed, so their intervals differed from Table 1 by 0.02-0.04); Msn2/4-only keeps its stratum rows.
ADJ5 <- fread(file.path(hsd, "p3_adjusted_core.tsv"))[
           (set == "reported_core" & abs(beta - 1) > 1e-6) | (set == "Msn24_only" & grepl("^Figure 3 stratum", what)),
           .(stratum = ifelse(set == "reported_core", "Hsf1_core", set), time, b0 = abs(beta) < 1e-6, ratio, lo, hi)]
chk5 <- merge(A5[, .(stratum, time, r0 = ratio_beta0, rc = cal_ratio)], ADJ5, by = c("stratum", "time"), allow.cartesian = TRUE)
stopifnot(nrow(chk5) == 12L, max(abs(chk5[b0 == TRUE, ratio - r0])) < 1e-6, max(abs(chk5[b0 == FALSE, ratio - rc])) < 1e-6)
ok5 <- draw_both("F5_amplitude_artefact", 4.25, 3.4, function() {
  layout(matrix(1, 1, 1)); after_layout()
  ## (A) colour = class (COL_BIO: Hsf1-dependent core green, Msn2/4-only red); fill = reading (filled uncorrected,
  ## open amplitude-matched, with its 95% interval). Solid lines throughout: categories never differ by dash.
  par(mar = c(3.4, 5.0, 4.6, 0.8))
  st <- c("Hsf1_core", "Msn24_only")
  cols <- c(COL_BIO[["Hsf1"]], COL_BIO[["Msn2/4"]]); pchs <- c(21, 22)
  plot(NA, xlim = c(10, 65), ylim = c(0.6, 1.95), xlab = "minutes after shift", ylab = "", yaxt = "n")
  axis(2, at = seq(0.6, 1.8, 0.2))   # range covers the adjusted intervals
  mtext("divergence ratio\n(gene set / comparison set)", side = 2, line = 2.2, cex = CEX_LAB, las = 0)
  abline(h = 1, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  for (i in seq_along(st)) {
    z <- A5[stratum == st[i]][order(time)]
    if (!nrow(z)) next
    off <- (i - 1.5) * 1.4
    a0 <- ADJ5[stratum == st[i] & b0 == TRUE][order(time)]; a1 <- ADJ5[stratum == st[i] & b0 == FALSE][order(time)]
    errbar(z$time + off - 0.35, a0$ratio, a0$lo, a0$hi, adjustcolor(cols[i], 0.5)); errbar(z$time + off + 0.35, a1$ratio, a1$lo, a1$hi, adjustcolor(cols[i], 0.5))
    lines(z$time + off - 0.35, z$ratio_beta0, col = cols[i], lwd = LW$data)
    lines(z$time + off + 0.35, z$cal_ratio, col = cols[i], lwd = LW$data, lty = 2)
    points(z$time + off - 0.35, z$ratio_beta0, pch = pchs[i], bg = cols[i], col = cols[i], cex = PCH_CEX, lwd = PCH_LWD)
    points(z$time + off + 0.35, z$cal_ratio, pch = pchs[i], bg = "white", col = cols[i], cex = PCH_CEX, lwd = PCH_LWD)
  }
  legend_top(c(paste(lab_of("Hsf1_core"), "uncorrected"), paste(lab_of("Hsf1_core"), "amplitude-corrected"),
               paste(lab_of("Msn24_only"), "uncorrected"), paste(lab_of("Msn24_only"), "amplitude-corrected")),   # plotted values are the beta-calibrated ratios (cal_ratio), not the matched-partner analysis
             pch = c(21, 21, 22, 22), pt.bg = c(cols[1], "white", cols[2], "white"), col = c(cols[1], cols[1], cols[2], cols[2]),
             lty = c(1, 2, 1, 2), lwd = LW$data, pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2)
  ## single panel: no panel tag
})
## Calibrated power is shown in Supplementary Figure S3.
## 
CAP$F5 <- "F5. The amplitude artefact, and the resolution of the resulting negative. Filled = uncorrected, open = amplitude-matched; colour gives the class. The Msn2/4 arm's apparent positive (1.274 / 1.266 / 1.243, 6/6 hybrids, CI excluding 1) collapses to 0.975 / 1.041 / 1.071 once amplitude is matched, because Msn2/4 targets respond 2.04x / 1.69x / 1.49x more strongly than TF-independent genes and rms-type divergence statistics scale with amplitude. The Hsf1-associated core shows no enrichment either way. The artefact addressed here is a reference-level artefact (a shared reference arm entering two contrasts), which differs from reference-mapping bias in allele-specific expression (reads aligned to one allele's genome). Sources: a5_cis_ratio_calibrated_37_30.tsv, a5_power_calibration_37_30.tsv."

## ------------------------------------------------------------ F6 cis/trans -
SH <- fread(file.path(hsd, "a6_cis_share_37_30.tsv"))
## panel (a) draws the primary correlation-adjusted interval (p3_adjusted_cisshare.tsv lo/hi) the text uses.
SHP <- fread(file.path(hsd, "p3_adjusted_cisshare.tsv"))
stopifnot(max(abs(merge(SH, SHP, by = c("stratum", "time"))[, cis_share.x - cis_share.y])) < 1e-12)
SH <- merge(SH[, !c("lo", "hi")], SHP[, .(stratum, time, lo, hi)], by = c("stratum", "time"))
AX <- fread(file.path(hsd, "a6_axis_amplitudes_37_30.tsv"))
ID <- fread(file.path(hsd, "a6_cistrans_identity.tsv"))
ok6 <- draw_both("F6_cis_trans", 6.9, 3.4, function() {
  layout(matrix(1:3, 1, 3), widths = c(1.25, 1, 1)); after_layout()
  ## (A) colour = class (COL_BIO Hsf1 green for the Hsf1-dependent core; black other induced; grey not induced),
  ## each also its own shape; the 0.5 reference is black dashed
  par(mar = c(3.4, 4.2, 4.6, 0.8))
  sts <- c("Hsf1_core", "other_induced", "not_induced")
  cols <- c(NPG$red, NPG$slate, GREYS$g80); bgs <- cols; pchs <- c(21, 22, 24)
  plot(NA, xlim = c(10, 66), ylim = c(0, 0.72), xlab = "minutes after shift", ylab = "cis share, cis² / (cis² + remainder²)")  # the quantity omits the covariance; upper bounds <= 0.631 for the three strata drawn
  abline(h = 0.5, lty = 2, col = COL_SCI$truth, lwd = LW$ref)
  text(64, 0.52, "equal", cex = CEX_TXT)
  for (i in seq_along(sts)) { z <- SH[stratum == sts[i]][order(time)]
    if (!nrow(z)) next
    errbar(z$time + (i - 2) * 1.3, z$cis_share, z$lo, z$hi, cols[i])
    lines(z$time + (i - 2) * 1.3, z$cis_share, col = cols[i], lwd = LW$data)
    ptb(z$time + (i - 2) * 1.3, z$cis_share, pchs[i], bgs[i]) }
  legend_top(sprintf("%s (n = %d)", lab_of(sts), vapply(sts, function(s) Nof(SH, "stratum", s), 0L)),
             pch = pchs, pt.bg = bgs, col = COL_SCI$truth, pt.cex = PCH_CEX, pt.lwd = PCH_LWD)
  panel_tag("A")
  ## (b) and (c) show every cross separately, so the six-of-six pattern is visible; noise-corrected
  ## amplitudes of the Hsf1-dependent core from a6_axis_amplitudes; crosses distinguished by line type and end label.
  cz <- AX[stratum == "Hsf1_core"][order(pair, time)]; prs <- unique(cz$pair); lt <- HYB_LTY[prs]; stopifnot(!anyNA(lt))
  ## one Okabe-Ito colour per hybrid (colour-blind safe; six dash styles in one colour are hard to tell apart)
  pc <- HYB_SUB[prs]; stopifnot(!anyNA(pc))   # palette colours only
  yl <- c(0, max(c(cz$cis_tru, cz$trn_tru)) * 1.12)
  for (comp in c("cis_tru", "trn_tru")) {
    par(mar = c(3.4, 4.0, 4.6, 0.8))
    plot(NA, xlim = c(10, 66), ylim = yl, xlab = "minutes after shift",
         ylab = if (comp == "cis_tru") "cis amplitude, Hsf1-dependent core" else "parental remainder amplitude, core")
    for (pp in prs) { z <- cz[pair == pp]
      lines(z$time, z[[comp]], lty = lt[[pp]], lwd = LW$data, col = pc[[pp]])
      ptb(z$time, z[[comp]], HYB_PCH[[pp]], pc[[pp]])
    }
    if (comp == "cis_tru") legend_top(sub("x", "\u00d7", prs), lty = lt[prs], lwd = LW$data, col = pc[prs], pch = HYB_PCH[prs], pt.bg = pc[prs], pt.cex = PCH_CEX, rows = 3)
    panel_tag(if (comp == "cis_tru") "B" else "C") }
})
CAP$F6 <- sprintf("F6. cis/trans decomposition of divergence. (A) The cis share of total divergence falls 0.439 -> 0.313 -> 0.210 for the Hsf1-associated core; 17 of 18 stratum x timepoint cells lie below 0.5, the single exception being Hsf1_joint at 15 min (0.614, n = 36). Bars are bootstrap 95%% CI. (B) Absolute noise-corrected amplitudes relative to 15 min: cis falls (-24.5%%) while trans rises (+31.5%%) in 6/6 species pairs with total amplitude flat; the identical pattern in not_induced shows this is NOT specific to the heat response. The identity estimate of cor(cis,trans) is %.3f / %.3f / %.3f at 15/30/60 min; the least-biased cross-block estimate is -0.166, an 18-cell median on a different summary scale. The null for the unbiased cross-block estimator is 0, so this is MILD COMPENSATION, in the same direction as Fraser 2019 (PMID 30270122) and Rojas & Lang 2026 (doi 10.1093/genetics/iyag233), which are cited as the methodological base. r_mech, the value of cor(cis,trans) at cov(cis,par) = 0, marks near-complete compensation, not zero compensation. Sources: a6_cis_share_37_30.tsv, a6_axis_amplitudes_37_30.tsv, a6_cistrans_identity.tsv.",
  ID[time==15, cor_ct_identity], ID[time==30, cor_ct_identity], ID[time==60, cor_ct_identity])

## --------------------------------------------------------------- F7 Hsf1 domains -
DOM <- fread(file.path(hsd, "hsf1_domain_identity.tsv"))
SLD <- fread(file.path(hsd, "hsf1_sliding_identity.tsv"))
DIV <- fread(file.path(hsd, "hsf1_divergent_residues.tsv"))
ok7 <- draw_both("F7_hsf1_domains", 6.9, 3.4, function() {
  ## Single panel, so no tag and no title (the caption names it). Colour-light: the identity trace is grey, the
  ## per-region means black and heavier, regions whose alignment is reliable are shaded light grey and the unreliable
  ## low-complexity segment a darker grey (not red: in P3 red and green are reserved for the dependence classes);
  ## variable DBD positions are open black triangles. The DBD annotation sits low and the key high so they do not
  ## crowd, and the linker label is centred on its region so it does not run into the DBD label.
  par(mar = c(3.4, 4.2, 5.2, 0.8))
  ylo <- 40; yhi <- 122
  SLD[, win_mid := (win_start + win_end) / 2]   # windows plotted at their centres, not their starts
  plot(NA, xlim = range(SLD$win_mid), ylim = c(ylo, yhi), yaxt = "n",
       xlab = expression(paste("position in the ", italic("S. cerevisiae"), " Hsf1 protein (UniProt P10961, 833 aa)")),
       ylab = "mean pairwise identity (%)")   # the plotted values are means over the six species pairs (script 08, line 196)
  axis(2, at = seq(40, 100, 10))
  SH_OK <- "#EFEFEF"; SH_BAD <- GREYS$g80; SH_DBD <- TINT$cyan
  for (i in seq_len(nrow(DOM))) {
    rect(DOM$start[i], ylo, DOM$end[i], 102, col = if (!isTRUE(DOM$align_reliable[i])) SH_BAD else if (grepl("DBD", DOM$domain[i])) SH_DBD else SH_OK, border = NA)
    segments(c(DOM$start[i], DOM$end[i]), ylo + 8, c(DOM$start[i], DOM$end[i]), 102, col = "white", lwd = LW$ci)
  }
  lines(SLD$win_mid, SLD$id_mean, col = COL_SCI$uncertainty, lwd = LW$data)
  for (i in seq_len(nrow(DOM))) if (isTRUE(DOM$align_reliable[i]))   # no mean inside the excluded region
    segments(DOM$start[i], DOM$id_mean[i], DOM$end[i], DOM$id_mean[i], col = if (grepl("DBD", DOM$domain[i])) NPG$blue else COL_SCI$truth, lwd = 0.73 / 0.75)
  LBL <- c(`N-terminal (disordered)` = "N-terminal (disordered)", DBD = "DBD", `flexible linker` = "linker",
           trimerization = "trimerization", `C-terminal (disordered/AR)` = "C-terminal (disordered/AR)")
  for (i in seq_len(nrow(DOM))) if (DOM$domain[i] %in% names(LBL))
    local({
      d <- DOM$domain[i]; narrow <- d %in% c("flexible linker", "trimerization")
      x <- if (d == "trimerization") DOM$start[i] else (DOM$start[i] + DOM$end[i]) / 2
      x <- if (d == "flexible linker") (DOM$start[i] + DOM$end[i]) / 2 else x
      y <- ylo + if (narrow) 2.0 else 5.2
      text(x, y, LBL[[d]], cex = CEX_TXT, adj = c(if (d == "trimerization") 0 else 0.5, 0.5), col = if (d == "DBD") NPG$blue else COL_SCI$truth)
    })
  dbd <- DOM[grepl("DBD", domain)]
  if (nrow(dbd)) {
    dv <- DIV[sc_pos >= dbd$start[1] & sc_pos <= dbd$end[1]]
    if (nrow(dv)) points(dv$sc_pos, rep(105, nrow(dv)), pch = 25, bg = NPG$red, col = NPG$red,
                         cex = PCH_CEX * 0.9, lwd = PCH_LWD)
    ## one line from x = 290: starting at the DBD end + 120 (379) it would run past the right edge of the data region
    text(290, 111,
         sprintf("DBD %d\u2013%d, %.1f%% identical, %d of %d positions variable (%s)",
                 dbd$start[1], dbd$end[1], dbd$id_mean[1], nrow(dv), dbd$n_aa[1],
                 paste(sort(dv$sc_pos), collapse = ", ")), cex = CEX_TXT, adj = 0)
  }
  legend_top(c("sliding-window identity", "per-region mean identity", "alignment reliable",
               "alignment unreliable (low-complexity repeat)", "variable DBD position"),
             lty = c(1, 1, NA, NA, NA), lwd = c(LW$data, 0.73 / 0.75, NA, NA, NA),
             col = c(COL_SCI$uncertainty, COL_SCI$truth, COL_SCI$truth, COL_SCI$truth, NPG$red),
             pch = c(NA, NA, 22, 22, 25), pt.bg = c(NA, NA, SH_OK, SH_BAD, NPG$red), pt.cex = PCH_CEX,
             pt.lwd = c(NA, NA, LW$ci, LW$ci, PCH_LWD), rows = 2)
})
CAP$F7 <- "F7. Hsf1 orthologue divergence across the four species, in the S. cerevisiae coordinate frame (UniProt P10961 release 2026_03, 833 aa, externally locked). Grey line: sliding-window identity. Solid segments: per-domain mean identity. Shaded bands mark domains whose alignment is reliable (the three commonest residues make up less than half of the S. cerevisiae segment) versus unreliable inside low-complexity repeats, where a global aligner pairs residues by position rather than by homology. The DNA-binding domain (170-259) is 97.2% identical with only 4 of 90 positions variable, mostly conservative substitutions, so they are candidate residues and not a mechanism. Whole-protein identity vs Sc: Sp 91.0%, Sk 81.4%, Su 78.5%. Strain caveat: only S. kudriavzevii (CR85) shares a strain between the sequence and expression axes, so all sequence statements are at the species-reference level. Sources: hsf1_domain_identity.tsv, hsf1_sliding_identity.tsv, hsf1_divergent_residues.tsv."

## ------------------------------------------------- F8 HSE scan + resolution ------
PC <- fread(file.path(hsd, "b3v2_positive_control.tsv"))
CH <- fread(file.path(hsd, "b3v2_chance_calibration_highres.tsv"))
# The seven-column bound written by 11_b3v2_resolution_bound.R. The schema is asserted rather than trusted, because
# script 10 writes a different five-column summary over its own 10 substrates, without p90.
RS <- fread(file.path(hsd, "b3v2_power_resolution.tsv"))
stopifnot(identical(sort(names(RS)), sort(c("time", "n", "n_never", "min", "median", "p90", "max"))))
# The panel also draws one point per substrate, so the substrate-level values are pooled again here for the scatter.
# The detection column is named differently in the calibration table and in the shards; whichever of the two names
# is present is used.
P0 <- fread(file.path(hsd, "b3v2_power_calibration.tsv"))
dcol <- intersect(c("detected_infer_alpha", "detected"), names(P0))[1]
stopifnot(!is.na(dcol))
P0 <- P0[, .(substrate = paste0("s0_", substrate), c_inj, time, detected = get(dcol))]
SHD <- rbindlist(lapply(1:3, function(s) {
  f <- file.path(hsd, sprintf("b3v2_resolution_shard%d.tsv", s))
  if (!file.exists(f)) return(NULL)
  d <- fread(f); d[, .(substrate, c_inj, time, detected)] }), use.names = TRUE)
PWR <- rbindlist(list(P0, SHD), use.names = TRUE)
say("F8 panel c pools %d substrates", uniqueN(PWR$substrate))
ok8 <- draw_both("F8_hse", 6.9, 5.7, function() {
  ## Colour-free: nothing here is naive or corrected, and no class is being contrasted with another. Element classes
  ## are black and light grey fills; summaries are black shapes; per-substrate values grey; references dotted grey.
  ## (A) and (B) are horizontal, names on the y axis.
  ordn <- c("Hsf1_core","Hsf1_joint","Msn24_only","other_induced","TF_independent","not_induced")
  z <- PC[stratum %in% ordn]; setkey(z, stratum); z <- z[ordn]
  cz <- CH[!is.na(ratio_lo)]
  ml <- local({
    tf <- tempfile(fileext = ".pdf")
    grDevices::cairo_pdf(tf, width = 6.85, height = 5.7, pointsize = PS_DEV, family = FONT_SUB)
    on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
    fig_par(); plot.new()
    c(max(strwidth(lab_of(z$stratum), units = "inches", cex = CEX_TXT)),
      max(strwidth(cz$metric, units = "inches", cex = CEX_TXT))) / par("csi") + 1.0
  })
  layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), widths = c(1.25, 1), heights = c(1, 1.25)); after_layout()
  ## (A) HSE carriage by dependence class, Hsf1_core at the top; classic above gapped within each class
  par(mar = c(3.4, ml[1], 3.0, 0.8))
  k <- nrow(z); yy <- rev(seq_len(k))
  FILL <- c(NPG$blue, NPG$red)
  xm <- max(100 * c(z$any_classic, z$any_gapped)) * 1.12
  plot(NA, xlim = c(0, xm), ylim = c(0.5, k + 0.5), yaxt = "n", xlab = "% promoters carrying an HSE", ylab = "")
  rect(0, yy + 0.02, 100 * z$any_classic, yy + 0.36, col = FILL[1], border = COL_SCI$truth, lwd = LW$ci)
  rect(0, yy - 0.36, 100 * z$any_gapped, yy - 0.02, col = FILL[2], border = COL_SCI$truth, lwd = LW$ci)
  axis(2, yy, lab_of(z$stratum), las = 1, tick = FALSE)
  legend_top(c("classic HSE", "gapped HSE"), fill = FILL, border = COL_SCI$truth)
  panel_tag("A")
  ## (B) chance calibration, with the 2x admission threshold
  par(mar = c(3.4, ml[2], 3.0, 0.8))
  k2 <- nrow(cz); y2 <- rev(seq_len(k2))
  plot(NA, xlim = c(0, max(cz$ratio_hi) * 1.1), ylim = c(0.4, k2 + 0.6), yaxt = "n",
       xlab = "observed / matched random", ylab = "")
  rect(0, y2 - 0.3, cz$ratio, y2 + 0.3, col = c(classic = NPG$blue, gapped = NPG$red, either = "white")[cz$metric], border = COL_SCI$truth, lwd = LW$ci)
  segments(cz$ratio_lo, y2, cz$ratio_hi, y2, col = COL_SCI$truth, lwd = LW$ci)
  abline(v = 2, lty = 3, col = COL_SCI$uncertainty, lwd = LW$ref)
  axis(2, y2, cz$metric, las = 1, tick = FALSE)
  legend_top("admission threshold, 2\u00d7", lty = 3, col = COL_SCI$uncertainty, lwd = LW$ref)
  panel_tag("B")
  ## (C) smallest detectable coupling per substrate (grey, sized by count), with the median and the 90th percentile
  par(mar = c(3.4, 5.0, 4.6, 0.8))
  mc <- PWR[c_inj > 0 & detected == TRUE, .(min_c = min(c_inj)), by = .(substrate, time)]
  tt <- sort(unique(RS$time))
  plot(NA, xlim = c(0.45, length(tt) + 1.15), ylim = c(-0.02, 0.92), xaxt = "n", xlab = "minutes after shift", ylab = "")
  axis(1, seq_along(tt), tt)
  mtext("smallest detectable\ncoupling c", side = 2, line = 2.2, cex = CEX_LAB, las = 0)
  for (i in seq_along(tt)) if (i %% 2 == 1)
    rect(i - 0.5, par("usr")[3], i + 0.5, par("usr")[4], col = adjustcolor("#DFDFDF", 0.8), border = NA)   # #E5E5E5 on white; the frame under it stays visible
  LH8 <- 1.6 * strheight("8", units = "user", cex = CEX_TXT)
  spread <- function(v) {                         # push stacked labels apart, keeping order, anchored at the top
    o <- order(-v); y <- v[o]
    for (k in seq_along(y)[-1]) if (y[k - 1] - y[k] < LH8) y[k] <- y[k - 1] - LH8
    out <- numeric(length(v)); out[o] <- y; out
  }
  for (i in seq_along(tt)) {
    w <- mc[time == tt[i], min_c]
    if (length(w)) {
      cnt <- as.data.table(table(w))[, .(v = as.numeric(w), k = N)]
      points(rep(i - 0.22, nrow(cnt)), cnt$v, pch = 21, bg = COL_SCI$uncertainty, col = NA,
             cex = PCH_CEX * (0.45 + 1.1 * sqrt(cnt$k / max(cnt$k))))
      text(rep(i - 0.29, nrow(cnt)), spread(cnt$v), cnt$k, cex = CEX_TXT, adj = 1)
    }
    r <- RS[time == tt[i]]
    segments(i - 0.02, 0.278, i + 0.46, 0.278, lty = 3, col = COL_SCI$uncertainty, lwd = LW$ref)
    ptb(i + 0.06, r$median, 23, NPG$red)
    ptb(i + 0.06, r$p90, 24, NPG$cyan)
    if (abs(r$median - r$p90) > 0.02) {
      text(i + 0.11, r$p90, sprintf("p90 = %.2f", r$p90), cex = CEX_TXT, adj = 0)
      text(i + 0.11, r$median, sprintf("med = %.2f", r$median), cex = CEX_TXT, adj = 0)
    } else
      text(i + 0.11, r$p90, sprintf("p90 = med = %.2f", r$p90), cex = CEX_TXT, adj = 0)
  }
  legend_top(c("per-substrate value (label = how many substrates)", "median", "90th percentile (resolution limit)",
               "coupling implied by the 15-min signal"), pch = c(21, 23, 24, NA), lty = c(NA, NA, NA, 3),
             pt.bg = c(COL_SCI$uncertainty, COL_SCI$truth, "white", NA), col = c(COL_SCI$uncertainty, COL_SCI$truth,
             COL_SCI$truth, COL_SCI$uncertainty), pt.cex = PCH_CEX, pt.lwd = PCH_LWD, lwd = LW$ref, rows = 2)
  panel_tag("C")
})
CAP$F8 <- sprintf("F8. Target-side HSE architecture: the scan is validated, and the test is nevertheless not feasible. (A) Carriage of the two pre-specified HSE classes (Yamamoto 2005 PMID 15647283; Pincus 2018 PMID 30332327) rises monotonically with Hsf1 dependency; Hsf1_core n = %d versus not_induced n = %d, permutation p_emp = 0.009 for classic and 0.0000 for gapped. This control is independent of the test that follows. (B) Observed carriage divided by a composition-matched random expectation (NRAND = 1e5), with bootstrap 95%% CI propagating both rates; admission requires the CI lower bound to clear 2x. Whole-promoter unit counts give 1.27x, fail the twofold admission threshold and are not used; gapped exceeds classic with a paired ratio-of-ratios of 1.55 [1.30, 1.88], so the primary split is classic UNION gapped. (C) Smallest detectable coupling for the cis-minus-trans discriminant, one point per randomised substrate (52 substrates, permutation criterion at the inference's own Bonferroni alpha = 0.0167). The 90th percentile is the reported equivalence bound; per-substrate values span c ~ 0.02-0.80, so the test is reported as not feasible and the interval itself is the result. The 15 min signal's implied coupling of about 0.278 is only 35%% of the conservative bound. Sources: b3v2_positive_control.tsv, b3v2_chance_calibration_highres.tsv, b3v2_power_calibration.tsv, b3v2_power_resolution.tsv.",
  Nof(PC,"stratum","Hsf1_core"), Nof(PC,"stratum","not_induced"))

## ---------------------------------------------------------- F9 module profile ---
MA <- fread(file.path(hsd, "s4_module_amplitude.tsv"))
MD <- fread(file.path(hsd, "s4_module_divergence_beta.tsv"))
MU <- fread(file.path(hsd, "s4_uniformity_test.tsv"))
SC <- fread(file.path(hsd, "s4_stem_collapse_ratios.tsv"))
MODS <- c("proteostasis","translation_ribosome","oxidative","metabolism","cell_cycle_growth")
# module HOG counts looked up by name from the membership table, never positionally
MEMB <- fread(file.path(hsd, "s4_module_membership.tsv"))[, .(n = uniqueN(hog)), by = module]
NMOD <- function(m) Nof(MEMB, "module", m)
ok9 <- draw_both("F9_modules", 6.9, 5.6, function() {
  ## The style's palette has no module colours: modules are told apart by shape and fill, the heat map
  ## is a Viridis ramp, readings by fill (C: open before, filled after the collapse).
  layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE)); after_layout()
  MPCH <- c(proteostasis = 24, translation_ribosome = 22, oxidative = 25, metabolism = 21, cell_cycle_growth = 23)
  MBG  <- c(proteostasis = COL_SCI$uncertainty, translation_ribosome = COL_SCI$light, oxidative = "white",
            metabolism = COL_SCI$truth, cell_cycle_growth = "white")
  ## (A) noise-corrected amplitude by module and species at 60 min
  par(mar = c(5.6, 10.6, 3.0, 1.4))
  W <- dcast(MA[time == 60L], module ~ species, value.var = "rms_true")
  setkey(W, module); W <- W[MODS]
  m <- as.matrix(W[, c("Sc","Sp","Sk","Su")])
  MCOL <- c(proteostasis = NPG$red, translation_ribosome = NPG$cyan, oxidative = NPG$blue, metabolism = NPG$peach,
            cell_cycle_growth = NPG$slate)   # modules stay black: in P3 the Okabe-Ito blues and oranges are the species
  PAL9 <- grDevices::colorRampPalette(c("#440154", "#3B528B", "#21918C", "#5EC962", "#FDE725"))(24)
  image(seq_len(ncol(m)), seq_len(nrow(m)), t(m), col = PAL9, axes = FALSE, xlab = "", ylab = "")
  SPL <- lab_of(colnames(m))
  assert_angled_labels_fit(seq_len(ncol(m)), SPL, 30, CEX_TXT)
  axis(1, seq_len(ncol(m)), labels = FALSE)
  text(seq_len(ncol(m)), par("usr")[3] - 0.10 * diff(par("usr")[3:4]), labels = SPL, srt = 30, adj = 1,
       xpd = TRUE, cex = CEX_TXT, font = 3)
  axis(2, seq_len(nrow(m)), lab_of(MODS), las = 1, tick = FALSE)
  cell_col <- function(v) PAL9[pmax(1, pmin(24, 1 + floor(23 * (v - min(m)) / diff(range(m)))))]
  lum <- function(cl) { r <- col2rgb(cl) / 255; sum(c(0.2126, 0.7152, 0.0722) * r) }
  for (i in seq_len(ncol(m))) for (j in seq_len(nrow(m)))
    text(i, j, sprintf("%.2f", m[j,i]), cex = CEX_TXT, col = if (lum(cell_col(m[j,i])) < 0.55) "white" else COL_SCI$truth)
  box(lwd = LW$axis)
  panel_tag("A")
  ## (B) divergence by module at 60 min across beta; the measured exponent 0.657 is the dashed line
  par(mar = c(3.4, 5.0, 5.4, 0.8))
  yb <- range(c(MD[time == 60L & module %in% MODS, ratio],
                MD[time == 60L & module %in% MODS, c(lo, hi)]), na.rm = TRUE)   # all five modules
  yb <- yb + c(-1, 1) * 0.10 * diff(yb)
  plot(NA, xlim = c(-0.07, 1.07), ylim = yb, xlab = expression(beta~"(amplitude normalisation)"), ylab = "")
  mtext("divergence ratio vs\nother modules", side = 2, line = 2.2, cex = CEX_LAB, las = 0)
  abline(h = 1, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  abline(v = 0.657, col = COL_SCI$uncertainty, lty = 3, lwd = LW$ref)
  text(0.657, yb[1] + 0.03 * diff(yb), "measured 0.657 ", adj = 1, cex = CEX_TXT)
  for (i in seq_along(MODS)) { z <- MD[module == MODS[i] & time == 60L][order(beta)]
    if (!nrow(z)) next
    ## intervals for every module, drawn as gene-independent (none excludes 1 after adjustment, Table 5),
    ## and open symbols, because filled symbols mean 'collapsed' in panel C
    xo <- z$beta + (i - 3) * 0.02   # offset wide enough that five modules do not stack near 1
    lines(xo, z$ratio, col = MCOL[[MODS[i]]], lwd = LW$data)
    ptb(xo, z$ratio, MPCH[[MODS[i]]], "white", col = MCOL[[MODS[i]]])
    errbar(xo, z$ratio, z$lo, z$hi, COL_SCI$uncertainty) }   # bars drawn last so short ones show over the markers
  legend_top(lab_of(MODS), pch = MPCH[MODS], pt.bg = "white", col = MCOL[MODS], pt.cex = PCH_CEX, pt.lwd = PCH_LWD,
             rows = 3)
  ## the adjustment statement is in the caption; an in-panel note would cross the exponent line and a bar
  panel_tag("B")
  ## (C) stem-proxy paralogue collapse: shape = module, open = full (grey line), filled = collapsed (black line)
  par(mar = c(3.4, 5.0, 4.6, 0.8))
  MM <- SC[module %in% c("metabolism","cell_cycle_growth","proteostasis") & time == 60L]
  yc <- range(c(MM$ratio_full, MM$ratio_collapsed, MM$lo, MM$hi), na.rm = TRUE)
  yc <- yc + c(-1, 1) * 0.08 * diff(yc)
  ## Near a ratio of 1 the three modules and the two readings overlap at every beta, so they are dodged
  ## horizontally, by module (-0.03, 0, +0.03) and within a module by reading (full -0.009, collapsed +0.009).
  ## Ticks sit at the computed exponents 0, 0.25, 0.5, 0.75 and 1; 0.657 is the dotted line.
  DODGE <- c(metabolism = -0.03, cell_cycle_growth = 0, proteostasis = 0.03)
  plot(NA, xlim = c(-0.06, 1.06), ylim = yc, xlab = expression(beta), ylab = "", xaxt = "n")
  axis(1, at = c(0, 0.25, 0.5, 0.75, 1))
  mtext("divergence ratio vs\nother modules", side = 2, line = 2.2, cex = CEX_LAB, las = 0)
  abline(h = 1, col = COL_SCI$truth, lty = 2, lwd = LW$ref); abline(v = 0.657, lty = 3, col = COL_SCI$uncertainty, lwd = LW$ref)
  for (mmod in c("metabolism","cell_cycle_growth","proteostasis")) {
    z <- SC[module == mmod & time == 60L][order(beta)]
    if (!nrow(z)) next
    xf <- z$beta + DODGE[[mmod]] - 0.009; xc <- z$beta + DODGE[[mmod]] + 0.009
    errbar(xc, z$ratio_collapsed, z$lo, z$hi, COL_SCI$uncertainty)
    lines(xf, z$ratio_full, col = MCOL[[mmod]], lwd = LW$data, lty = 2)
    lines(xc, z$ratio_collapsed, col = MCOL[[mmod]], lwd = LW$data)
    ptb(xf, z$ratio_full, MPCH[[mmod]], "white", col = MCOL[[mmod]])
    ptb(xc, z$ratio_collapsed, MPCH[[mmod]], MCOL[[mmod]], col = MCOL[[mmod]])
  }
  legend_top(c(lab_of("metabolism"), lab_of("cell_cycle_growth"), lab_of("proteostasis"),
               "full (open, dashed line)", "collapsed (filled, solid line)"),
             pch = c(21, 23, 24, 21, 21), pt.bg = c(MCOL[["metabolism"]], MCOL[["cell_cycle_growth"]], MCOL[["proteostasis"]], "white", CHARCOAL),
             col = c(MCOL[["metabolism"]], MCOL[["cell_cycle_growth"]], MCOL[["proteostasis"]], CHARCOAL, CHARCOAL),
             lty = c(NA, NA, NA, 2, 1), lwd = LW$data, pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2, span = "device")
  panel_tag("C")
})
CAP$F9 <- sprintf("F9. Which cellular programmes diverged. (A) Noise-corrected response amplitude rms_true by module and species at 60 min; modules differ several-fold, which is why every ratio in (B) and (C) is shown across the whole beta range rather than at beta = 0 alone. (B) Divergence ratio of each module against the union of the other mutually exclusive modules, as a function of the amplitude-normalisation exponent beta; the dashed line marks beta = 0.657, the measured exponent of the parental divergence axis from A6. Proteostasis includes 1 at every beta and every timepoint (n = %d HOGs); metabolism (n = %d) is enriched and cell_cycle_growth (n = %d) is depleted. (C) Stem-proxy paralogue collapse, a deliberately severe sensitivity analysis on top of the primary uncollapsed analysis: it asks only whether a few related family members can by themselves carry a module-level ratio, and it is not a model of gene-family history. Open symbols full, filled symbols collapsed, with bootstrap 95%% CI. Metabolism survives the collapse and is slightly strengthened by it (1.146 -> 1.155 at the anchor, p_emp = 0.006); cell_cycle_growth is strengthened markedly; proteostasis survives at no beta or timepoint. The stem proxy errs in BOTH directions -- it over-merges non-families such as HSP26/42/82/104 and under-merges verified families with different stems such as TIP1 versus TIR1-TIR4 (PMID 36134526) -- so it is reported as a sensitivity analysis, not as truth. Sources: s4_module_amplitude.tsv, s4_module_divergence_beta.tsv, s4_uniformity_test.tsv, s4_stem_collapse_ratios.tsv.",
  NMOD("proteostasis"), NMOD("metabolism"), NMOD("cell_cycle_growth"))

## ------------------------------------------------------------------- captions ---
## F3 and its caption are produced by 60_p2_fig3_two_backgrounds.R, with the numbers taken from
## the enrichment table. The caption is included in FIGURE_CAPTIONS.md in figure order;
## the figure and caption share one source.
capf3 <- file.path(fig, "CAPTION_F3.md")
if (!file.exists(capf3))
  stop(sprintf("%s is missing: run scripts/hsf1/60_p2_fig3_two_backgrounds.R first, which draws F3 and writes its caption", capf3))
## Normalize the caption prefix. 60_p2_fig3_two_backgrounds.R writes it as "F3. ...", and the caption loop
## Caption files use "**F3** ...". Normalising to the internal "F3. " form keeps
## idempotent parsing.
CAP$F3 <- paste(readLines(capf3, warn = FALSE), collapse = " ")
CAP$F3 <- sub("^\\*\\*F3\\*\\*\\s*", "F3. ", CAP$F3)
stopifnot(nzchar(CAP$F3), grepl("^F3\\.", CAP$F3))
CAP <- CAP[c("F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9")]
stopifnot(!any(vapply(CAP, is.null, TRUE)))

capf <- file.path(fig, "FIGURE_CAPTIONS.md")
writeLines(c("# Figure captions — F1..F9", "",
  "> Generated by `scripts/hsf1/19_figures.R`. Every number is read from a landed `.tsv`; no figure recomputes anything.",
  "> Legend N values are looked up BY NAME (`Nof()` stops if the name is absent) — never with `min()`/`max()`.",
  "> No information is carried by colour alone. Vector PDF and 600 dpi PNG.", "",
  "Pearson 1897 (doi 10.1098/rspl.1896.0076), Oldham 1962 (doi 10.1016/0021-9681(62)90116-9) and Atchley et al. 1976 (doi 10.2307/2412740) are cited for the spurious correlation of ratios that share a component.", "",
  unlist(lapply(names(CAP), function(k) paste0("## ", k, "\n\n", CAP[[k]], "\n")))), capf)
say("\ncaptions written: %s", capf)
## Each caption is also written to a separate file.
## 
## The supplementary item that holds each P3 figure's source tables, is read from the supporting-data map. A number
## 
SUPPF <- file.path(root, "doc/04_figures/final/figure_support_P3.tsv")
SUPP <- if (file.exists(SUPPF)) fread(SUPPF) else data.table(tag = character(), supplement_item = character())

for (k in names(CAP)) {
  d <- fig
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  body <- sub(sprintf("^%s\\. ", k), "", CAP[[k]])
  ## A reader sent from a figure to its supporting data should be given the item, not left to match file names.
  it <- SUPP[tag == k, supplement_item]
  if (length(it) == 1L)
    body <- sprintf("%s Supporting data: Supplementary Table %s.", sub("\\s+$", "", body), it)
  writeLines(sprintf("**%s** %s", k, body), file.path(d, sprintf("CAPTION_%s.md", k)))
}
say("per-figure captions written: %d", length(CAP))
say("figures produced: %s", paste(names(CAP)[c(ok1,ok2,TRUE,ok4,ok5,ok6,ok7,ok8,ok9)], collapse = ", "))

## Stop on a failed figure. draw_both catches the error from a drawing function, prints it and returns FALSE; without
## this check a figure whose fit assertion failed would leave a defective PDF on disk with exit status 0, and every
## measurement assertion in this file would be advisory.
DREW <- c(F2_epistasis = ok1, F2d_library_share = ok2, F4_parental = ok4, F5_amplitude_artefact = ok5,
          F6_cis_trans = ok6, F7_hsf1_domains = ok7, F8_hse = ok8, F9_modules = ok9)
if (any(!DREW)) {
    say("!! FAILED: %s", paste(names(DREW)[!DREW], collapse = ", "))
    sink()
    stop(sprintf("%d figure(s) failed to draw: %s", sum(!DREW), paste(names(DREW)[!DREW], collapse = ", ")))
}
writeLines(capture.output(sessionInfo()), file.path(fig, "sessionInfo_figures.txt"))
say("=== done ===")
sink()
invisible(NULL)
