# 99_fig_p2_schematic.R -- P2 Figure 1, schematic of the four-genotype interaction contrast and of the
# compositional offset. Numbers are the library shares at 37 C, 30 min already reported in Results (p2_library_share.tsv).
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); fig <- file.path(root, "doc/04_figures/final")
source(file.path(root, "scripts/figures/theme.R")); source(file.path(root, "scripts/figures/fig_style.R"))
SH <- c(W = 0.1551, H = 0.2785, M = 0.1041, D = 0.0703)
stopifnot(abs(SH[["D"]] - SH[["H"]] - SH[["M"]] + SH[["W"]] - (-0.1572)) < 1e-4)
box <- function(x0, y0, x1, y1, lab, fill = "white", txt = COL_SCI$truth, edge = COL_SCI$truth) {
  rect(x0, y0, x1, y1, col = fill, border = edge, lwd = LW$axis)
  text((x0 + x1) / 2, (y0 + y1) / 2, lab, cex = CEX_LAB, col = txt)
}
## Drawn at the full 174 mm width with two equal panels. Every text item is at least 7.5 pt (CEX_TXT); the cell
## letters and the contrast are at the 8.5 pt axis-title size so the grid reads at the size of panel b.
draw <- function() {
  layout(matrix(1:2, 1, 2), widths = c(1.05, 1)); after_layout()
  ## (a) the 2 x 2 design and the contrast
  par(mar = c(1, 0.5, 2.6, 0.5)); plot.new(); plot.window(c(0, 10), c(0, 12))
  cx <- c(3.7, 7.7); hw <- 1.9                     # cell centres and half-width
  ry <- c(9.0, 5.2); hh <- 1.8                      # row centres and half-height
  text(cx[1], 11.55, "Msn2/4 +", cex = CEX_TXT); text(cx[2], 11.55, "Msn2/4 \u2212", cex = CEX_TXT)
  text(cx[1] - hw - 0.3, ry[1], "Hsf1 +", adj = c(1, 0.5), cex = CEX_TXT)
  text(cx[1] - hw - 0.3, ry[2], "Hsf1 \u2212", adj = c(1, 0.5), cex = CEX_TXT)
  ## cells filled in the genotype colour; black text on pale fills, white on dark ones
  lum <- function(cl) { r <- col2rgb(cl) / 255; sum(c(0.2126, 0.7152, 0.0722) * r) }
  cell <- function(i, j, lab, g) { f <- GENO_SUB$col[[g]]
    box(cx[j] - hw, ry[i] - hh, cx[j] + hw, ry[i] + hh, lab, f, if (lum(f) > 0.4) "#000000" else "white", COL_SCI$truth) }
  cell(1, 1, "W\n+1", "Wildtype"); cell(1, 2, "M\n\u22121", "MSN24.KO")
  cell(2, 1, "H\n\u22121", "HSF1.KD"); cell(2, 2, "D\n+1", "Double.KDKO")
  text(mean(cx), 2.05, "I = D \u2212 H \u2212 M + W", cex = CEX_LAB, font = 2)
  text(mean(cx), 0.7, "I = 0 if the effects add", cex = CEX_TXT)
  panel_tag("A")
  ## The library-share arithmetic is shown in Figure 3b and is not repeated here.
  ## (b) library-share dilution as four stacked bars, not to scale; dark segments ordered as the
  ## measured shares at 37 C, 30 min (H largest, then W, M, D); no numbers on the axis
  par(mar = c(3.0, 3.6, 2.6, 7.6))
  sh <- SH[c("W", "H", "M", "D")] / max(SH) * 0.6
  M2 <- rbind(sh, 1 - sh)
  bp <- barplot(M2, col = c(COL_SCI$truth, GREYS$g90), border = COL_SCI$truth, lwd = LW$ci, axes = FALSE,
                names.arg = rep("", 4), ylab = "share of library (fixed total)", space = 0.35, width = 1)
  rect(bp - 0.5, 0, bp + 0.5, sh, col = GENO_SUB$col[c("Wildtype","HSF1.KD","MSN24.KO","Double.KDKO")], border = COL_SCI$truth, lwd = LW$ci)
  axis(2, at = c(0, 1), labels = FALSE)
  mtext(c("W", "H", "M", "D"), side = 1, at = bp, line = 0.7, cex = CEX_LAB)
  text(bp[4] + 0.7, sh[4] / 2, "responding\ngenes", adj = 0, cex = CEX_TXT, xpd = NA)
  text(bp[4] + 0.7, sh[4] + (1 - sh[4]) / 2, "other genes", adj = 0, cex = CEX_TXT, xpd = NA)
  panel_tag("B")
}
fig_both("F0s_schematic", 174, 74, draw, fig)
