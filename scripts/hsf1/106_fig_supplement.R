#!/usr/bin/env Rscript
## 106_fig_supplement.R -- supplementary figures for scripts 103-105b, in the shared figure
## style (scripts/figures/fig_style.R).
##  FS_P2_torc1_split: P2 Supplementary Figure S1. (a) rapamycin readout at 25 C (script 103), (b) held-out class interaction (script 104b).
##  FS_P3_calibration: P3 Supplementary Figure S3. (a) exclusion rate of three intervals for random sets (script 105),
##                     (b) power after injecting a known divergence excess (script 105b).
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
source(file.path(root, "scripts/figures/fig_style.R")); out <- file.path(root, "doc/04_figures/final")
ptb <- function(x, y, pch, bg, col = COL_SCI$truth) points(x, y, pch = pch, bg = bg, col = col, cex = PCH_CEX, lwd = PCH_LWD)
errbar <- function(x, lo, hi, col = COL_SCI$truth) segments(x, lo, x, hi, col = col, lwd = LW$ci)
## ---------------- P2 ----------------
T1 <- fread(file.path(hsd, "p2_torc1_readout.tsv")); SP <- fread(file.path(hsd, "p2_class_split_b.tsv"))
fig_both("FS_P2_torc1_split", WIDTH_MM[["full"]], 70, function() {
  layout(matrix(1:2, 1, 2), widths = c(1.15, 1)); after_layout()
  par(mar = c(3.4, 4.4, 3.6, 0.8))
  sets <- c("ribosomal protein (RPL/RPS)", "ribosome biogenesis (GO:0042254)", "Msn2/4-only class")
  labs <- c("ribosomal\nprotein", "ribosome\nbiogenesis", "Msn2/4-only\nclass")
  cts <- c("H-W 25", "D-M 25", "WT 37/10 vs 25"); pch <- c(21, 22, 24)
  bg <- c(NPG$red, NPG$cyan, NPG$blue); off <- c(-0.22, 0, 0.22)
  plot(NA, xlim = c(0.5, 3.5), ylim = c(-2.9, 3.2), xaxt = "n", xlab = "", ylab = "median log2 fold change")
  axis(1, at = 1:3, labels = FALSE); mtext(labs, side = 1, at = 1:3, line = 1.6, cex = CEX_TXT)
  abline(h = 0, lty = 2, lwd = LW$ref, col = CHARCOAL)
  for (i in seq_along(sets)) for (j in seq_along(cts)) { z <- T1[set == sets[i] & contrast == cts[j]]
    errbar(i + off[j], z$lo, z$hi); ptb(i + off[j], z$median_log2fc, pch[j], bg[j]) }
  legend_top(c("Hsf1-depleted vs wild type, 25 °C", "double vs msn2Δ msn4Δ, 25 °C", "wild type, 37 °C/10 min vs 25 °C"),
             pch = pch, pt.bg = bg, col = COL_SCI$truth, pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 3)
  panel_tag("A")
  par(mar = c(3.4, 4.4, 3.6, 0.8))
  cl <- c("Msn24_only", "redundant", "TF_independent", "Hsf1_only"); cll <- c("Msn2/4-\nonly", "redun-\ndant", "TF-\nindep.", "Hsf1-\nonly")
  plot(NA, xlim = c(0.5, 4.5), ylim = c(-3.2, 0.6), xaxt = "n", xlab = "", ylab = "class minus reference, log2")
  axis(1, at = 1:4, labels = FALSE); mtext(cll, side = 1, at = 1:4, line = 1.6, cex = CEX_TXT)
  abline(h = 0, lty = 2, lwd = LW$ref, col = COL_SCI$truth)
  for (i in seq_along(cl)) { z <- SP[cls == cl[i]][order(rep)]
    points(i - 0.12 + c(-0.05, 0, 0.05), z$held_out, pch = 22, bg = COL_SCI$corrected, col = COL_SCI$truth, cex = PCH_CEX * 0.8, lwd = PCH_LWD)
    points(i + 0.12 + c(-0.05, 0, 0.05), z$with_k, pch = 21, bg = "white", col = COL_SCI$naive, cex = PCH_CEX * 0.8, lwd = PCH_LWD)
    segments(i - 0.24, mean(z$held_out), i, mean(z$held_out), lwd = LW$data, col = COL_SCI$corrected)
    segments(i, mean(z$with_k), i + 0.24, mean(z$with_k), lwd = LW$data, col = COL_SCI$naive) }
  legend_top(c("held-out replicate", "same replicate in definition"), pch = c(22, 21),
             pt.bg = c(COL_SCI$corrected, "white"), col = c(COL_SCI$truth, COL_SCI$naive), pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2)
  panel_tag("B")
}, out)
## ---------------- P3 ----------------
if (file.exists(file.path(hsd, "p3_power_injection.tsv"))) {
NC <- fread(file.path(hsd, "p3_null_calibration.tsv")); PW <- fread(file.path(hsd, "p3_power_injection.tsv"))
## the two interval methods: correlation-adjusted deep blue, gene bootstrap red
F11_A <- NPG$blue; F11_B <- NPG$red
fig_both("F11_calibration", WIDTH_MM[["full"]], 66, function() {
  layout(matrix(1:2, 1, 2)); after_layout()
  par(mar = c(3.4, 4.4, 3.6, 0.8))
  x <- seq_len(nrow(NC)); NC <- NC[order(size, time)]
  plot(NA, xlim = c(0.5, nrow(NC) + 0.5), ylim = c(-0.8, 12), xaxt = "n", xlab = "", ylab = "random sets excluding 1 (%)")
  axis(1, at = x, labels = FALSE); mtext(sprintf("%d min", NC$time), side = 1, at = x, line = 0.7, cex = CEX_TXT)
  mtext(c("83 groups", "134 groups"), side = 1, at = c(2, 5), line = 1.8, cex = CEX_TXT)
  abline(h = 5, lty = 2, lwd = LW$ref, col = CHARCOAL)
  ptb(x - 0.15, 100 * NC$excl_adjusted, 21, F11_A); ptb(x + 0.15, 100 * NC$excl_boot, 22, "white", F11_B)
  legend_top(c("correlation-adjusted (primary)", "gene bootstrap", "nominal 5%"), pch = c(21, 22, NA),
             lty = c(NA, NA, 2), pt.bg = c(F11_A, "white", NA), col = c(COL_SCI$truth, F11_B, COL_SCI$truth), pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 3)
  panel_tag("A")
  par(mar = c(3.4, 4.4, 3.6, 1.6))
  plot(NA, xlim = c(1.05, 1.6), ylim = c(-4, 106), xlab = "median ratio after injection", ylab = "power (% of sets excluding 1)")
  abline(h = 80, lty = 2, lwd = LW$ref, col = COL_SCI$truth)
  for (s in c(83, 134)) { z <- PW[size == s, .(r = mean(median_ratio), p = 100 * mean(power_adjusted), pb = 100 * mean(power_boot)), by = .(rr = r)][order(r)]
    pc <- if (s == 83) 21 else 22
    lines(z$r, z$p, lwd = LW$data, col = F11_A); ptb(z$r, z$p, pc, F11_A)
    lines(z$r, z$pb, lwd = LW$data, col = F11_B, lty = 1); ptb(z$r, z$pb, pc, "white", F11_B) }
  legend_top(c("adjusted, 83", "adjusted, 134", "bootstrap, 83", "bootstrap, 134"), pch = c(21, 22, 21, 22),
             pt.bg = c(F11_A, F11_A, "white", "white"), col = c(COL_SCI$truth, COL_SCI$truth, F11_B, F11_B),
             pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2)
  panel_tag("B")
}, out) }
cat("done\n")
