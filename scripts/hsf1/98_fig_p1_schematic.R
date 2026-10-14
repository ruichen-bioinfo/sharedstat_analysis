# 98_fig_p1_schematic.R -- P1 Figure 1, a schematic of the problem, the threshold, and the sharedstat workflow.
# Uses only numbers stated in the manuscript (gain 1.1; threshold 1/1.1 = 0.91; median reliability 0.89).
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); fig <- file.path(root, "doc/04_figures/final")
source(file.path(root, "scripts/figures/theme.R")); source(file.path(root, "scripts/figures/fig_style.R"))
box <- function(x0, y0, x1, y1, lab, fill = "white", font = 1, edge = GREYS$g67, txt = COL_SCI$truth) {
  rect(x0, y0, x1, y1, col = fill, border = edge, lwd = LW$axis)
  text((x0 + x1) / 2, (y0 + y1) / 2, lab, cex = CEX_TXT, font = font, col = txt)
}
arr <- function(x0, y0, x1, y1, col = CHARCOAL) arrows(x0, y0, x1, y1, length = 0.05, lwd = LW$arrow, col = col)
draw <- function() {
  layout(matrix(1:2, 1, 2), widths = c(1, 1.1)); after_layout()
  ## The box diagram of the model is in Box 1 of the manuscript, not in this figure.
  ## (A) the threshold, with the paper's numbers
  par(mar = c(3.6, 3.8, 3.2, 0.8))
  h <- seq(0.8, 1, length.out = 50)
  plot(h, 1.1 * h, type = "l", lwd = LW$data, col = COL_SCI$naive, xlim = c(0.8, 1), ylim = c(0.86, 1.12),
       xlab = "reference-arm reliability, h", ylab = expression(paste("expected naive gain (true ", gamma, " = 1.1)")))
  abline(h = 1, lty = 2, lwd = LW$ref, col = CHARCOAL)
  abline(v = 1 / 1.1, lty = 3, lwd = LW$ref, col = COL_SCI$neutral)
  points(0.89, 1.1 * 0.89, pch = 21, bg = COL_SCI$naive, col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  text(0.905, 0.875, "threshold\n1/1.1 = 0.91", adj = 1, cex = CEX_TXT, col = COL_SCI$truth)
  text(0.925, 0.955, "median\nh = 0.89", adj = 0, cex = CEX_TXT, col = COL_SCI$truth)
  segments(0.893, 0.976, 0.923, 0.96, lwd = LW$ref, col = COL_SCI$neutral)
  text(0.805, 1.015, "no divergence", adj = 0, cex = CEX_TXT, col = COL_SCI$truth)
  panel_tag("A")
  ## (B) sharedstat workflow
  par(mar = c(1, 0.5, 3.2, 0.5)); plot.new(); plot.window(c(0, 10), c(0, 10))
  box(0.3, 7.6, 9.7, 9.6, "measure h\nreliability_from_blocks()", fill = "white")
  box(0.3, 4.0, 9.7, 6.6, "diagnose the bias\nnaive_gain_bias()\nadmission_checklist()", fill = "white")
  ## The workflow forks: without a credible instrument it stops at sensitivity; with one it gives an instrumented
  ## gain that is conditional on relevance and exclusion
  box(0.1, 0.2, 4.8, 2.9, "no instrument:\nreport\nsensitivity", fill = "#FFFFFF")
  box(5.2, 0.2, 9.9, 2.9, "instrument:\niv_gain(),\nconditional", fill = "#FFFFFF", edge = COL_SCI$corrected, txt = COL_SCI$corrected)
  arr(5, 7.6, 5, 6.6); arr(3.6, 4.0, 2.45, 2.9); arr(6.4, 4.0, 7.55, 2.9, col = COL_SCI$corrected)
  panel_tag("B")
}
fig_both("FP1_0_schematic", 174, 70, draw, fig)
