# 97_fig_p2_time_course.R -- P2 Figure 6, from 96_p2_time_contrasts.R, with the class-definition point marked:
# Msn2/4-only induction over time, and peak classes.
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
fig <- file.path(root, "doc/04_figures/final")
source(file.path(root, "scripts/figures/theme.R")); source(file.path(root, "scripts/figures/fig_style.R"))
S <- fread(file.path(hsd, "p2_time_class_medians.tsv"))[cls == "Msn24_only"]
P <- fread(file.path(hsd, "p2_time_peak_classes.tsv"))[cls == "Msn24_only" & test == "BH+LFC" & t == 1 & genotype %in% c("Wildtype", "HSF1.KD")]
GENO <- c("Wildtype", "HSF1.KD", "MSN24.KO", "Double.KDKO")
GL <- c(Wildtype = "wild type", HSF1.KD = "Hsf1-depleted", MSN24.KO = "msn2Δ msn4Δ", Double.KDKO = "double perturbation")
traj <- function(T) rbindlist(lapply(GENO, function(g) data.table(g = g, x = c(0, 10, 30),
  med = c(0, S[contrast == paste0("E", T, "|", g), med], S[contrast == paste0("F", T, "|", g), med]),
  lo  = c(0, S[contrast == paste0("E", T, "|", g), lo],  S[contrast == paste0("F", T, "|", g), lo]),
  hi  = c(0, S[contrast == paste0("E", T, "|", g), hi],  S[contrast == paste0("F", T, "|", g), hi]))))
A37 <- traj("37"); A42 <- traj("42"); YL <- range(c(A37$lo, A37$hi, A42$lo, A42$hi, 0)) * c(1.1, 1.25)
off <- c(Wildtype = -1.2, HSF1.KD = -0.4, MSN24.KO = 0.4, Double.KDKO = 1.2)
panel_traj <- function(D, T, tag, leg) {
  par(mar = c(3.6, 3.8, 5.0, 0.8))
  plot(NA, xlim = c(-3, 33), ylim = YL, xaxt = "n", xlab = "minutes after shift", ylab = "median log2 fold change vs 25 °C")
  axis(1, c(0, 10, 30)); abline(h = 0, col = COL_SCI$truth, lwd = LW$ref, lty = 3)
  for (gg in GENO) { g <- gg; d <- D[D$g == gg]; xx <- d$x + off[[g]]
    segments(xx, d$lo, xx, d$hi, col = COL_SCI$uncertainty, lwd = LW$ci)
    lines(xx, d$med, col = GENO_SUB$col[[g]], lwd = LW$data, lty = GENO_SUB$lty[[g]])
    points(xx, d$med, pch = GENO_SUB$pch[[g]], bg = GENO_SUB$bg[[g]], col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD) }
  text(-2, YL[2], paste0(T, " °C"), adj = c(0, 1), cex = CEX_LAB, font = 2, col = COL_SCI$truth)
  # the classes were selected at 37 °C, 30 min; mark that point so the 10-to-30 min change is read with it in mind
  if (T == "37") { text(33, YL[2], "classes\ndefined here", adj = c(1, 1), cex = CEX_TXT, col = COL_SCI$uncertainty, xpd = NA)
    arrows(30, YL[2] - 0.62, 30, max(D[D$x == 30]$hi) + 0.12, length = 0.05, col = COL_SCI$uncertainty, lwd = LW$ref) }
  if (leg) legend_top(GL[GENO], pch = GENO_SUB$pch[GENO], pt.bg = GENO_SUB$bg[GENO], col = COL_SCI$truth,
                      pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2, span = "device")
  panel_tag(tag)
}
draw <- function() {
  layout(matrix(1:3, 1, 3), widths = c(1, 1, 1.05)); after_layout()
  panel_traj(A37, "37", "A", TRUE); panel_traj(A42, "42", "B", FALSE)
  par(mar = c(4.6, 3.8, 5.0, 0.8))
  P[, key := paste(genotype, temp)]; ord <- c("Wildtype 37", "HSF1.KD 37", "Wildtype 42", "HSF1.KD 42"); setkey(P, key); Q <- P[ord]
  M <- rbind(Q$transient, Q$sustained, Q$rising) / rep(Q$n_induced, each = 3) * 100
  cols <- c(PAL$sign_neg, PAL$sign_zero, PAL$sign_pos)
  bp <- barplot(M, col = cols, border = COL_SCI$truth, lwd = LW$ci, ylim = c(0, 100), axes = FALSE, names.arg = rep("", 4),
                ylab = "% of induced Msn2/4-only genes", space = c(0.3, 0.3, 0.9, 0.3))
  axis(2)
  for (i in seq_along(bp)) mtext(rep(c("wild\ntype", "anchor-\naway"), 2)[i], side = 1, at = bp[i], line = 1.3, cex = CEX_TXT, col = COL_SCI$truth)
  mtext(c("37 °C", "42 °C"), side = 1, at = c(mean(bp[1:2]), mean(bp[3:4])), line = 3.1, cex = CEX_LAB, font = 2, col = COL_SCI$truth)
  legend_top(c("lower at 30 min", "no change", "higher at 30 min"), fill = cols, border = COL_SCI$truth, rows = 3)
  panel_tag("C")
}
fig_both("F5t_time_course", 174, 88, draw, fig)
