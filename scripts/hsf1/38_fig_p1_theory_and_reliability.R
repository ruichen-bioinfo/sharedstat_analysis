#!/usr/bin/env Rscript
## 38_fig_p1_theory_and_reliability.R -- the P1 figures for the analytical result and its corollaries
## (Figure 1) and for the measured reliability (Figure 3, Supplementary Figure S4).
##
## No new analysis. Nothing here is computed that is not already computed and saved: Figure P1.1 draws the
## closed form the manuscript states and the saved replicate-axis table, and Figure P1.3 draws the saved
## reliability table. No number produced here is new.
##
## The replicate panel (Figure 1B) is drawn from p1_replicate_axis.tsv, not from the design grid's n_rep
## axis. That axis holds h constant by construction (the grid sets the per-observation error to
## sqrt(var * n_rep) and then divides by sqrt(n_rep)), so a panel drawn along it can only be flat.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1"); figd <- file.path(root, "doc/04_figures/final")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
## Colours come from the shared theme, so a colour means the same thing in every figure of every
## manuscript. A local copy of the Okabe-Ito vector would drift from the theme unnoticed.
source(file.path(root, "scripts/figures/theme.R"))

## Device base point size 11. With the annotation cex values used here, an 8 pt base puts text at
## 4 to 6 pt at the final figure width, which is not legible in print. Raising the base moves every
## size by one factor and preserves the relative emphasis; the physical width does not change, so
## this is not a figure shrunk to make its type look bigger.
PT_BASE <- 11
OI <- as.list(OI)   # this file accesses colours as OI$name, so the theme vector is taken as a list
draw_both <- function(stem, w, h, f) {          # 6.9 in wide, as every other figure in this project
  for (dev in c("pdf", "png")) {
    p <- file.path(figd, paste0(stem, ".", dev))
    ## cairo_pdf, not pdf(): the default PDF device cannot encode the Greek gamma these panels need and
    ## drops it with a conversion warning, so the axis labels would have arrived at the journal incomplete.
    ## Both devices get an explicit pointsize; otherwise R's default of 12 applies and the cex values here sit on a
    ## size this script never chose.
    if (dev == "pdf") grDevices::cairo_pdf(p, width = w, height = h, pointsize = PT_BASE)
    else png(p, width = w, height = h, units = "in", res = DPI, type = "cairo", pointsize = PT_BASE)
    ## No global `par(cex = 0.72)`: a global multiplier scales every call-site cex, so an annotation at cex = 0.64
    ## comes out at 0.72 x 0.64 x base, 5.1 pt at an 11 pt base, and raising the base cannot fix it because the
    ## multiplier scales with it. Explicit axis, label and title expansions keep the visual proportions of 0.72
    ## without shrinking annotations.
    par(cex.axis = 0.72, cex.lab = 0.80, cex.main = 0.85, mgp = c(2.1, 0.7, 0), tcl = -0.25); f(); dev.off()
  }
  say("%s ok (pdf + png)", stem)
}

G  <- fread(file.path(hsd, "p1_design_grid.tsv"))
RA <- fread(file.path(hsd, "p1_replicate_axis.tsv"))
## The figure must use the same h as the text. `p1_measured_reliability.tsv` holds the uncentred form; the
## manuscript reports the centred one, because the gain estimators are built from centred moments. B5's table holds
## both, the centred value computed by `reliability_from_blocks(..., center = TRUE)` -- the package function itself,
## not a copy of its arithmetic -- together with its bootstrap interval, so one table serves both panels.
B5 <- fread(file.path(hsd, "p1x_b5_reliability_centred_ci.tsv"))[dataset == "Fay 2023"]
B5[, c("stratum", "hybrid", "time") := tstrsplit(arm, "|", fixed = TRUE)]
B5[, time := as.integer(time)]
stopifnot(nrow(B5) == 108L, all(c("h_centred", "h_uncentred", "lo", "hi") %in% names(B5)))
## Fidelity: the uncentred column must reproduce the reference table for the same 108 cells.
MRU <- fread(file.path(hsd, "p1_measured_reliability.tsv"))
chk <- merge(MRU[, .(stratum, hybrid, time, h_old = h)], B5[, .(stratum, hybrid, time, h_uncentred)],
             by = c("stratum", "hybrid", "time"))
stopifnot(nrow(chk) == 108L, max(abs(chk$h_old - chk$h_uncentred)) < 1e-12)
## Figure 3 shows the reference-arm reliability, the h of the gamma*h identity (scripts 91/93, same 108 cells).
## B5 above is the DoD-contrast reliability; it is kept only for the fidelity check and is not plotted.
MR <- fread(file.path(hsd, "p1_reference_arm_h_adjusted.tsv"))
stopifnot(nrow(MR) == 108L, all(c("h", "hi_adj2") %in% names(MR)))

## ============================ FIGURE P1.1 -- the analytical result ============================
GAM <- c(1.0, 1.1, 1.3, 1.6)
hh  <- seq(0.30, 1.00, by = 0.005)
## gamma is drawn in Noto Sans Italic everywhere in P1's figures. Upright Noto Sans draws it almost as a "y"
## (the two glyphs overlap with IoU 0.82 once scaled to one box); in italics the overlap is 0.20 against an italic y and
## 0.30 against an upright one, and a variable in italics is the manuscript's own convention. Plotmath's `gamma` symbol
## is set in the device's symbol font, which fontconfig resolves to upright Noto Sans here, so the character is given
## explicitly inside italic().
G_IT <- "\u03b3"
## Drawn in the shared figure style (scripts/figures/fig_style.R). No panel titles; the corollary each panel shows is
## named in the caption. Colour: every naive reading is dark red; black is the model and the reference.
## Categories differ by shape, fill or a direct label, never by dash pattern.
source(file.path(root, "scripts/figures/fig_style.R"))
gl <- function(v) bquote(italic(.(G_IT)) == .(sprintf("%.1f", v)))
fig_both("FP1_1_analytical_result", 174, 76, function() {
  layout(matrix(1:3, 1, 3)); after_layout()
  ## (A) the closed form: gamma_hat = gamma * h; the heavier line is the no-divergence case gamma = 1
  ## The end labels start 1.5 mm outside the plot region (which runs to h = 1.028), in a right margin widened to
  ## hold them; at h = 1.01 they would sit on the region's right edge.
  par(mar = c(3.4, 3.6, 3.4, 3.4))
  plot(NA, xlim = range(hh), ylim = c(0.3, 1.7), xlab = "reference-arm reliability, h", ylab = "expected naive gain")
  x_lab <- grconvertX(grconvertX(par("usr")[2], "user", "inches") + 1.5 / 25.4, "inches", "user")
  abline(h = 1, col = GREYS$g67, lty = 3, lwd = LW$ref)
  GCOL_A <- function(g) if (g == 1) CHARCOAL else GAMMA_COL[[sprintf("%.1f", g)]]
  ## the no-divergence line is black, so it cannot be confused with the red naive-estimate lines
  for (i in seq_along(GAM)) lines(hh, GAM[i] * hh, col = GCOL_A(GAM[i]), lty = if (GAM[i] == 1) 2 else 1, lwd = LW$data)
  points(c(0.5, 0.7, 0.9), c(0.5, 0.7, 0.9), pch = 22, bg = CHARCOAL, cex = PCH_CEX, lwd = PCH_LWD)
  op <- par(xpd = NA)
  for (i in seq_along(GAM)) text(x_lab, GAM[i], as.expression(gl(GAM[i])), adj = c(0, 0.5), cex = CEX_TXT,
                                 col = GCOL_A(GAM[i]))
  par(op)
  ## each true gain is labelled at its line end in its own colour; the key names only the dashed no-divergence line
  legend_top("true gain = 1", col = CHARCOAL, lty = 2, lwd = LW$data, rows = 1)
  panel_tag("A")
  ## (B) Corollary 2: the bias is O(1) in features and O(1/n) in replicates.
  ## One logarithmic x axis, the fold increase in design size over each series' smallest design (500 genes;
  ## 2 replicates per arm): each x is n / min(n), so both series share one axis. Genes are open squares,
  ## replicates filled circles, both in the naive colour; the key is the same two-row key as panels A and C,
  ## so the top margin matches theirs.
  par(mar = c(3.4, 3.6, 3.4, 0.8))
  FA <- RA[axis == "features"][order(n_gene)]; RB <- RA[axis == "replicates"][order(n_rep)]
  xf <- FA$n_gene / min(FA$n_gene); xr <- RB$n_rep / min(RB$n_rep)
  stopifnot(max(xf) == 256, min(xf) == 1, min(xr) == 1, max(xr) <= 256)
  ## drawn on log2(fold) with labelled ticks rather than log = "x": legend_top() converts usr[1] as a data value, which
  ## on a log axis starting below 1 is the log of a negative number, and the key silently comes out at NA
  xf <- log2(xf); xr <- log2(xr)
  plot(NA, xlim = c(0, 8), ylim = c(0, max(abs(c(RB$naive_emp - RB$gamma, FA$naive_emp - FA$gamma))) * 1.12),
       xaxt = "n", xlab = "fold increase in design size", ylab = bquote(paste("|bias| at true ", italic(.(G_IT)), " = 1.3")))
  axis(1, at = 2 * (0:4), labels = 4^(0:4))
  lines(xf, abs(FA$naive_emp - FA$gamma), col = NPG$cyan, lwd = LW$data)
  points(xf, abs(FA$naive_emp - FA$gamma), pch = 22, bg = NPG$cyan, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  lines(xr, abs(RB$naive_emp - RB$gamma), col = NPG$red, lwd = LW$data)
  points(xr, abs(RB$naive_emp - RB$gamma), pch = 21, bg = NPG$red, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  legend_top(c("genes", "replicates per arm"), pch = c(22, 21), pt.bg = c(NPG$cyan, NPG$red),
             col = c(NPG$cyan, NPG$red), pt.cex = PCH_CEX, lwd = LW$data, rows = 2)
  panel_tag("B")
  ## (C) Corollary 3: the composed bias crosses 1, classical dilution alone does not
  par(mar = c(3.4, 3.6, 3.4, 0.8))
  g <- 1.3
  plot(NA, xlim = range(hh), ylim = c(0.35, 1.45), xlab = "reference-arm reliability, h",
       ylab = bquote(paste("expected estimate (true ", italic(.(G_IT)), " = 1.3)")))
  polygon(c(hh[hh < 1/g], rev(hh[hh < 1/g])), c(g * hh[hh < 1/g], rep(1, sum(hh < 1/g))),
          col = TINT$red, border = NA)
  abline(h = 1, col = CHARCOAL, lty = 2, lwd = LW$ref)
  abline(v = 1/g, col = COL_SCI$uncertainty, lty = 3, lwd = LW$ref)
  lines(hh, g * hh, col = COL_SCI$naive, lwd = LW$data)
  lines(hh, g * hh + (1 - hh), col = COL_SCI$corrected, lwd = LW$data)
  text(1/g - 0.035, 0.42, bquote(paste("1/", italic(.(G_IT)), " = ", .(sprintf("%.3f", 1/g)))),
       cex = CEX_TXT, srt = 90, adj = 0)
  legend_top(c("naive gain", "dilution only"), col = c(COL_SCI$naive, COL_SCI$corrected),
             lwd = LW$data, rows = 2)
  panel_tag("C")
}, figd)
say("FP1_1_analytical_result ok (pdf + png)")

## caption facts derived from the plotted tables, never typed
FA <- RA[axis == "features"]; RB <- RA[axis == "replicates"][order(n_rep)]
cap1 <- sprintf(paste0(
  "Figure 1. The analytical result and its two consequences, drawn from the closed form and from the ",
  "landed replicate-axis table. (A) The closed form: the naive estimator returns the true gain multiplied ",
  "by the reference arm's reliability. The heavy line is the no-divergence case, true gain 1; squares mark ",
  "h = 0.5, 0.7 and 0.9, where it returns exactly those values. (B) Corollary 2: the bias is O(1) in features ",
  "and O(1/n_rep) in replicates. Over a %.0f-fold change in feature count it moves by %.4f (open ",
  "squares); over %d to %d replicates it falls from %.4f to %.4f, a %.0f-fold reduction (filled ",
  "circles). Both series are drawn against one logarithmic axis, the fold increase in design size over the ",
  "smallest design of each (500 genes; 2 replicates per arm). Both series are at a true gain of 1.3. (C) Corollary 3 at a true gain ",
  "of 1.3. Classical regression dilution alone gives gamma*h + (1 - h) (blue), which lies between 1 and ",
  "gamma and cannot cross 1. The two effects composed give gamma*h (red), which crosses 1 at h = %.3f; in the ",
  "shaded region a true increase is reported as a decrease."),
  max(FA$n_gene)/min(FA$n_gene), diff(range(abs(FA$naive_emp - FA$gamma))),
  min(RB$n_rep), max(RB$n_rep), abs(RB$naive_emp - RB$gamma)[1],
  abs(RB$naive_emp - RB$gamma)[nrow(RB)],
  abs(RB$naive_emp - RB$gamma)[1] / abs(RB$naive_emp - RB$gamma)[nrow(RB)], 1/1.3)

## ==================== FIGURE P1.3 -- the reliability measured in the data ====================
GS <- c(1.1, 1.2, 1.3, 1.5, 2.0)
## h per gene class is Supplementary Figure S4; the interval split by true gain is drawn separately below
fig_both("FS4_p1_reliability_by_class", 108, 100, function() {
  ST <- sort(unique(MR$stratum))
  SRT <- 40
  SLAB <- lab_of(ST)
  mb <- local({
    ## measured in the same font as the figure; pdf(NULL) would measure Helvetica, about 6% narrower
    tf <- tempfile(fileext = ".pdf")
    grDevices::cairo_pdf(tf, width = 6.85, height = 3.94, pointsize = PS_DEV, family = FONT_SUB)
    on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
    fig_par(); plot.new()
    mar_lines_for_angled(SLAB, SRT, CEX_TXT, pad_lines = 2.0)   # 2.0 lines keeps the labels clear of the edge
  })
  ## (A) measured h per stratum; the three 1/gamma thresholds are grey dotted and labelled where they end
  ## the first angled label runs left of its tick by its width x cos(40 degrees); the left margin holds it
  par(mar = c(mb / 0.75, 5.4, 2.6, 4.6))
  bx <- lapply(ST, function(s) MR[stratum == s, h])
  BCOL <- unname(CLASS_SUB[ST]); stopifnot(!anyNA(BCOL))   # class colour on outline and median only
  BCOL[ST == "not_induced"] <- "#7A7A7A"   # the light-grey class colour would vanish as an outline in print
  boxplot(bx, names = rep("", length(ST)), col = "white", border = BCOL, outpch = 21, outcol = GREYS$g53,
          outbg = GREYS$g53, outcex = PCH_CEX * 0.8, lwd = LW$ci, medlwd = LW$data,
          ylim = c(0, 1.02), yaxt = "n",
          ylab = "reference-arm reliability, h", frame.plot = FALSE)
  axis(2, at = seq(0, 1, 0.2))
  axis(1, at = seq_along(ST), labels = FALSE)
  assert_angled_labels_fit(seq_along(ST), SLAB, SRT, CEX_TXT)
  text(seq_along(ST), par("usr")[3] - 0.07, labels = SLAB, srt = SRT, adj = 1, xpd = TRUE, cex = CEX_TXT)
  op <- par(xpd = NA)
  ## the three thresholds are 0.06-0.08 apart, less than a line of text: labels are pushed apart by a measured
  ## line height, and a short leader joins each to its line
  thr <- 1 / GS[1:3]; lh <- 1.7 * strheight(as.expression(gl(1.1)), units = "user", cex = CEX_TXT)
  ly <- thr; for (k in 2:length(ly)) if (ly[k - 1] - ly[k] < lh) ly[k] <- ly[k - 1] - lh
  xr <- par("usr")[2]
  for (i in seq_along(thr)) {
    segments(par("usr")[1], thr[i], xr, thr[i], col = COL_SCI$uncertainty, lty = 3, lwd = LW$ref)
    segments(xr, thr[i], xr + 0.12, ly[i], col = COL_SCI$uncertainty, lwd = LW$ref)
    text(xr + 0.16, ly[i], as.expression(gl(GS[i])), adj = c(0, 0.5), cex = CEX_TXT)
  }
  par(op)
}, figd)
say("FS4_p1_reliability_by_class ok (pdf + png)")
fig_both("FP1_3_measured_reliability", 108, 100, function() {
  ## (B) what the text reports. For each true gain, the 108 cells are split by where the whole rho^2-widened 95%
  ## interval lies relative to 1/gamma: wholly below (confidently in the reversal regime), across it (unresolved), or wholly
  ## above (confidently outside). The point-estimate share below the threshold is drawn as a marker on each bar.
  par(mar = c(4.0, 4.6, 3.8, 0.8))
  lo2 <- pmax(0, MR$h - sqrt(MR$vif2) * (MR$h - MR$lo)); hi2 <- MR$hi_adj2
  SPL <- sapply(GS, function(g) { t <- 1/g; c(below = 100 * mean(hi2 < t), across = 100 * mean(hi2 >= t & lo2 <= t),
                                              above = 100 * mean(lo2 > t)) })
  stopifnot(all(abs(colSums(SPL) - 100) < 1e-9))
  frac_pt <- vapply(GS, function(g) 100 * mean(MR$h < 1/g), 0)
  FILL <- c(NPG$red, GREYS$g90, NPG$cyan)   # mid grey so "spans" and "above" are distinct
  bp <- barplot(SPL, names.arg = sprintf("%.1f", GS), ylim = c(-5, 114), yaxt = "n",   # room for diamonds near 0 and counts
                col = FILL, border = COL_SCI$truth, lwd = LW$ci,
                xlab = bquote(paste("true gain, ", italic(.(G_IT)))), ylab = "stratum-timepoints (% of 108)")
  axis(2, at = seq(0, 100, 20))
  ## the cell counts behind each bar (below / spans / above), since the "below" share can be too thin to see
  CNT <- round(SPL * nrow(MR) / 100)
  text(bp, 107, sprintf("%d / %d / %d", CNT[1, ], CNT[2, ], CNT[3, ]), cex = CEX_TXT)
  points(bp, frac_pt, pch = 23, bg = NPG$blue, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  legend_top(c("interval of h wholly below 1/γ", "interval of h spans 1/γ", "interval of h wholly above 1/γ", "point estimate of h below 1/γ"),
             fill = c(FILL, NA), border = c(rep(COL_SCI$truth, 3), NA), pch = c(NA, NA, NA, 23),
             pt.bg = c(NA, NA, NA, NPG$blue), col = c(NA, NA, NA, COL_SCI$truth), pt.cex = PCH_CEX, rows = 2)
}, figd)
say("FP1_3_measured_reliability ok (pdf + png)")
HA_cap <- copy(MR)   # reference-arm table
frac_iv_cap <- vapply(GS, function(g) 100 * mean(HA_cap$hi_adj2 < 1/g), 0)
cap3 <- sprintf(paste0(
  "Figure 3. The reliability that controls the bias, measured in the data being analysed rather than ",
  "assumed. (A) h = 1 - var(null)/var(signal) for the reference arm, the first-named allele's heat response in the two replicate blocks, per stratum over %d cells; ",
  "nothing is fitted and all %d cells are estimable. Medians run from %.3f to %.3f and the full range is ",
  "%.3f to %.3f. Horizontal lines are the inversion thresholds 1/gamma. (B) For each true gain, the 108 cells split by where ",
  "the whole 95%% interval of h, widened for the correlation between genes, lies relative to 1/gamma: wholly below (red), ",
  "across it (grey) or wholly above (blue); diamonds give the share whose point estimate lies below (%s). ",
  "At gamma = 1.1, 2 cells are confidently below, 10 confidently above and 96 unresolved. gamma is not known in these data; ",
  "the panel states what the measured reliability implies for a hypothetical gain."),
  nrow(MR), nrow(MR), min(MR[, median(h), by = stratum]$V1), max(MR[, median(h), by = stratum]$V1),
  min(MR$h), max(MR$h),
  paste(sprintf("%.0f%% at gamma = %.1f", vapply(GS, function(g) 100 * mean(MR$h < 1/g), 0), GS),
        collapse = ", "))

writeLines(c(cap1, "", cap3), file.path(root, "doc/04_figures/final/FP1_captions.txt"))
say("\n--- Figure 1 caption ---"); say("%s", cap1)
say("\n--- Figure 3 caption ---"); say("%s", cap3)
