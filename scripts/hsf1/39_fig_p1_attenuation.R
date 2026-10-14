#!/usr/bin/env Rscript
## 39_fig_p1_attenuation.R -- P1 Figure 4, a single panel: the attenuation that makes a cross-block
## estimate a lower bound.
##
## Scope. The figure shows only what is P1's own: the closed-form attenuation, computed from P1's saved
## table, which is what the manuscript's section 5 cites the figure for. Counts against noise-corrected
## amplitudes and the compositional offset before and after re-normalisation are, in content and in source
## table, the companion manuscript's F2c and F2d (p2_interaction_amplitude.tsv, p2_library_share.tsv), and
## the cis/trans estimator ladder over the 18 hybrid x timepoint cells belongs to the third manuscript. The
## related-manuscripts disclosure states "No figure panel is reused", so none of those is drawn here.
##
## Source: data/02_processed/hsf1/p1_a2_scale_invariance.tsv

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
fig  <- file.path(root, "doc/04_figures/final")
DPI  <- 300
say  <- function(...) cat(sprintf(...), "\n")

f <- file.path(root, "data/02_processed/hsf1/p1_a2_scale_invariance.tsv")
if (!file.exists(f)) stop("missing landed table: p1_a2_scale_invariance.tsv")
ATT <- fread(f)

## ---- fidelity gate: the two values the manuscript prints ----
gate <- function(what, got, want, tol = 5e-4) {
  ok <- is.finite(got) && abs(got - want) <= tol
  say("  %-46s %9.4f vs %8.4f  %s", what, got, want, if (ok) "OK" else "MISMATCH")
  if (!ok) stop(sprintf("fidelity gate failed: %s", what))
}
say("=== fidelity gate ===")
gate("attenuated correlation at var ratio 1",  ATT[var_ratio == 1]$cor_ct,  -0.3259)
gate("attenuated correlation at var ratio 64", ATT[var_ratio == 64]$cor_ct, -0.4266)
gate("true correlation used in the sweep",     ATT$cor_true[1],             -0.5000)
## the closed form and the simulation must agree, which is the point of the panel
d <- max(abs(ATT$cor_ct - ATT$cor_pred))
say("  %-46s %9.2e  %s", "max |simulated - closed form|", d, if (d < 2e-3) "OK" else "TOO LARGE")
if (d >= 2e-3) stop("simulated and predicted attenuation disagree")

## Drawn in the shared figure style (scripts/figures/fig_style.R) at the 84 mm single-column width.
## Colour: black is the model and the truth, grey the reference. The simulated estimate is a small filled dot inside the
## closed form's open circle, because the two agree to under 2e-3.
source(file.path(root, "scripts/figures/fig_style.R"))
fig_both("FP1_4_attenuation", 84, 80, function() {
  par(mar = c(3.4, 5.0, 4.0, 0.8))
  x <- seq_along(ATT$var_ratio)
  plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = c(-0.54, -0.30), xaxt = "n", yaxt = "n",
       xlab = expression(paste("var(", italic(trans), ") / var(", italic(cis), ")")), ylab = "")
  axis(2, at = seq(-0.50, -0.30, 0.05))
  ## mtext() counts lines at the full 10 pt line height (12 pt), not at mex; 2.3 lines sits outside the ticks
  mtext("cis/trans correlation", side = 2, line = 2.7, las = 0, cex = CEX_LAB)
  abline(h = ATT$cor_true[1], lty = 2, col = CHARCOAL, lwd = LW$ref)
  text(length(x) + 0.15, ATT$cor_true[1] + 0.012, "true value, \u22120.5", cex = CEX_TXT, pos = 2)
  lines(x, ATT$cor_ct, col = COL_SCI$naive, lwd = LW$data)
  points(x, ATT$cor_pred, pch = 21, bg = "white", col = COL_SCI$naive, cex = PCH_CEX, lwd = PCH_LWD)
  points(x, ATT$cor_ct, pch = 16, cex = 0.8, col = COL_SCI$naive)
  axis(1, at = x, labels = ATT$var_ratio)
  y15 <- mean(ATT$cor_ct[1:2])
  arrows(1.5, ATT$cor_true[1], 1.5, y15, length = 0.05, code = 3, col = COL_SCI$uncertainty, lwd = LW$arrow)
  text(1.54, ATT$cor_true[1] + 0.3 * (y15 - ATT$cor_true[1]), "attenuation", cex = CEX_TXT,
       col = COL_SCI$neutral, pos = 4)
  legend_top(c("replicate-split estimate, simulated", "closed form"), pch = c(16, 21), pt.bg = c(COL_SCI$naive, "white"),
             col = c(COL_SCI$naive, COL_SCI$naive), pt.cex = c(0.8, PCH_CEX), pt.lwd = PCH_LWD, rows = 2)
}, fig)

W <- 84 / 25.4
say("\n=== written ===")
for (nm in c("FP1_4_attenuation.pdf", "FP1_4_attenuation.png")) {
  p <- file.path(fig, nm); say("  %-30s %8d bytes", nm, file.info(p)$size)
}
say("  width %.2f in = %.1f mm (single-column limit 84 mm)", W, W * 25.4)
