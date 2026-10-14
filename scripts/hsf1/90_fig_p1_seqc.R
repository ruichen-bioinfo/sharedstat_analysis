#!/usr/bin/env Rscript
## 90_fig_p1_seqc.R -- P1 Supplementary Figure S3: the SEQC known-truth test. Draws only what
## 89_p1_seqc_known_truth.R wrote; nothing is estimated here.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
source(file.path(root, "scripts/figures/theme.R"))
PT_BASE <- 11
source(file.path(root, "scripts/figures/fig_style.R"))
fig <- file.path(root, "doc/04_figures/final"); hsd <- file.path(root, "data/02_processed/hsf1")
LV <- rbindlist(lapply(c("bgi", "may"), function(s) fread(file.path(hsd, sprintf("p1_seqc_%s_levels.tsv", s)))))
VD <- rbindlist(lapply(c("bgi", "may"), function(s) fread(file.path(hsd, sprintf("p1_seqc_%s_verdicts.tsv", s)))))
draw <- function() {
  layout(matrix(1:2, 1, 2)); after_layout()
  yl <- c(0, max(c(LV$gamma_naive, LV$gamma_iv, VD$gamma_truth)) * 1.08)
  for (k in 1:2) {
    s <- c("BGI", "MAY")[k]; z <- LV[site == s][order(h)]; g <- VD[site == s, gamma_truth]
    par(mar = c(3.4, 3.8, 4.2, 0.8))
    plot(NA, xlim = c(0, 1), ylim = yl, xlab = expression(paste("measured reliability ", italic(h))),
         ylab = "estimated gain")
    abline(h = 1, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
    abline(h = g, col = CHARCOAL, lty = 3, lwd = LW$ref)
    abline(v = 1 / g, col = COL_SCI$uncertainty, lty = 3, lwd = LW$ref)
    lines(c(0, 1), c(0, g), col = COL_SCI$naive, lty = 2, lwd = LW$ref)
    text(1 / g, yl[2] * 0.97, "1/\u03b3", pos = 4, offset = 0.15, cex = CEX_TXT, col = COL_SCI$neutral)
    jn <- z$h >= 0.10; ji <- jn & z$first_stage >= 0.3
    points(z$h, z$gamma_naive, pch = 21, bg = ifelse(jn, COL_SCI$naive, "white"), col = COL_SCI$naive,
           cex = PCH_CEX, lwd = PCH_LWD)
    points(z$h, z$gamma_iv, pch = 22, bg = ifelse(ji, COL_SCI$corrected, "white"), col = COL_SCI$corrected,
           cex = PCH_CEX, lwd = PCH_LWD)
    legend_top(c("naive", "instrumented", "\u03b3h", "\u03b3 from ERCC"), pch = c(21, 22, NA, NA),
               lty = c(NA, NA, 2, 3), col = c(COL_SCI$naive, COL_SCI$corrected, COL_SCI$naive, CHARCOAL),
               pt.bg = c(COL_SCI$naive, COL_SCI$corrected, NA, NA), pt.cex = PCH_CEX, lwd = LW$ref, rows = 2)
    mtext(s, side = 3, line = 0.2, adj = 1, cex = CEX_TXT)
    panel_tag(c("A", "B")[k])
  }
}
fig_both("FS3_seqc_known_truth", 174, 80, draw, fig)
cat("FS3_seqc_known_truth written\n")
