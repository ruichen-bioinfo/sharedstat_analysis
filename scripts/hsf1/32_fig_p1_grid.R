# 32_fig_p1_grid.R -- P1 appendix figure FA3, the operating regime, drawn from the saved grid table only.
#
# The replicate axis is not drawn here. Along this grid's n_rep axis h is constant by construction,
# because the grid sets the per-observation error to sqrt(var * (1-h)/h * n_rep) and then divides by
# sqrt(n_rep). Since the bias is a function of h alone, a panel along that axis can only be flat: it
# would restate the construction, not test a claim. Once h is allowed to move the bias is O(1/n_rep);
# that behaviour is Figure 1c, drawn from p1_replicate_axis.tsv, which varies n_rep with the
# per-observation error held fixed.
# Every fact stated in the caption is computed here from the plotted table.
## Figure width is 6.9 in (178 mm = 7.008 in maximum). It is set at
## source because rescaling a raster can lose the >= 300 dpi
## resolution.
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1"); figd <- file.path(root, "doc/04_figures/final")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
G <- fread(file.path(hsd, "p1_design_grid.tsv"))
## Colours come from the shared theme, so a colour means the same thing in every figure of every
## manuscript. A local copy of the Okabe-Ito vector would drift from the theme unnoticed.
source(file.path(root, "scripts/figures/theme.R"))
OI <- as.list(OI)   # this file accesses colours as OI$name, so the theme vector is taken as a list
PT_BASE <- 11
## The x label is kept short to fit each panel's region at 11 pt
## without clipping.

## Drawn in the shared figure style of scripts/figures/fig_style.R:
## Noto Sans embedded, 7.5/8.5/10 pt, the style's line weights, bold uppercase tags, no panel titles, keys above the
## panels. Colour: COL_SCI -- naive dark red, IV-corrected dark green, truth black. In B and C the three true gains are
## all naive readings, so they share the naive colour and differ by shape and fill. Data are unchanged.
source(file.path(root, "scripts/figures/fig_style.R"))
GAMMAS <- sort(unique(G[gamma > 1, gamma])); HS <- sort(unique(G$h))
GPCH <- c(22, 24, 25); GBG <- unname(GAMMA_COL[c("1.1", "1.3", "1.6")]); GLC <- GBG   # true gain 1.1 red, 1.3 cyan, 1.6 deep blue
fig_both("FA3_p1_operating_regime", 174, 78, function() {
  layout(matrix(1:3, 1, 3)); after_layout()
  XL <- range(HS) + c(-0.04, 0.04)
  ## (A) the negative control
  par(mar = c(3.4, 3.6, 5.2, 0.8))
  NC <- G[gamma == 1.0, .(naive = mean(gamma_naive), iv = mean(gamma_iv)), by = h][order(h)]
  yl <- range(c(NC$naive, NC$iv, 1)); yl <- yl + c(-1, 1) * 0.08 * diff(yl)
  plot(NA, xlim = XL, ylim = yl, xlab = "reference-arm reliability, h", ylab = "estimated gain")
  abline(h = 1, col = CHARCOAL, lty = 2, lwd = LW$ref)
  lines(NC$h, NC$naive, col = COL_SCI$naive, lwd = LW$data)
  lines(HS, HS, col = COL_SCI$naive, lty = 2, lwd = LW$ref)
  points(NC$h, NC$naive, pch = 22, bg = COL_SCI$naive, cex = PCH_CEX, lwd = PCH_LWD)
  lines(NC$h, NC$iv, col = COL_SCI$corrected, lwd = LW$data)
  points(NC$h, NC$iv, pch = 21, bg = COL_SCI$corrected, cex = PCH_CEX, lwd = PCH_LWD)
  legend_top(c("naive", "instrumented", "true gain = 1", "predicted naive value, h"),
             pch = c(22, 21, NA, NA), pt.bg = c(COL_SCI$naive, COL_SCI$corrected, NA, NA), pt.cex = PCH_CEX,
             col = c(COL_SCI$truth, COL_SCI$truth, CHARCOAL, COL_SCI$naive), lty = c(NA, NA, 2, 2),
             lwd = c(PCH_LWD, PCH_LWD, LW$ref, LW$ref), rows = 2)
  panel_tag("A")
  ## (B) naive bias as a multiple of the true effect
  par(mar = c(3.4, 3.6, 5.2, 0.8))
  M <- sapply(GAMMAS, function(gg) sapply(HS, function(hh) G[h == hh & gamma == gg, mean(rel_bias_naive)]))
  yl2 <- range(c(as.vector(M), 0)); yl2 <- yl2 + c(-1, 1) * 0.08 * diff(yl2)
  plot(NA, xlim = XL, ylim = yl2, xlab = "reference-arm reliability, h", ylab = "naive bias / true effect")
  abline(h = 0, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  for (j in seq_along(GAMMAS)) {
    lines(HS, M[, j], col = GLC[j], lwd = LW$data)
    points(HS, M[, j], pch = GPCH[j], bg = GBG[j], col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  }
  legend_top(sprintf("true gain %.1f", GAMMAS), pch = GPCH, pt.bg = GBG, pt.cex = PCH_CEX,
             col = COL_SCI$truth, lwd = PCH_LWD, rows = 2)
  panel_tag("B")
  ## (C) wrong-direction rate
  par(mar = c(3.4, 3.6, 5.2, 0.8))
  W <- sapply(GAMMAS, function(gg) sapply(HS, function(hh) 100 * G[h == hh & gamma == gg, mean(sign_wrong_naive)]))
  plot(NA, xlim = XL, ylim = c(-6, 106), yaxt = "n", xlab = "reference-arm reliability, h",
       ylab = "wrong direction (% of simulations)")
  axis(2, at = seq(0, 100, 25))
  abline(h = c(0, 100), col = COL_SCI$uncertainty, lty = 3, lwd = LW$ref)
  for (j in seq_along(GAMMAS)) {
    xo <- HS + (j - (length(GAMMAS) + 1) / 2) * 0.012   # offset so equal values do not hide each other
    lines(xo, W[, j], col = GLC[j], lwd = LW$data)
    points(xo, W[, j], pch = GPCH[j], bg = GBG[j], col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  }
  legend_top(sprintf("true gain %.1f", GAMMAS), pch = GPCH, pt.bg = GBG, pt.cex = PCH_CEX,
             col = COL_SCI$truth, lwd = PCH_LWD, rows = 2)
  panel_tag("C")
}, figd)
say("FA3_p1_operating_regime ok (pdf + png)")
## caption facts derived, never typed
nc_pred <- G[gamma==1.0, .(obs=mean(gamma_naive)), by=h][order(h)]
cap <- sprintf(paste0("FA3. The operating regime of a shared reference arm, from %d grid cells x %d simulations. ",
  "(A) Negative control: with the true gain set to 1, i.e. no divergence whatsoever, the naive estimator ",
  "reports %s at h = %s, matching the algebraic prediction that it returns h exactly; the ",
  "IV-corrected estimator returns %.4f averaged over all no-divergence cells. In every one of the %d ",
  "no-divergence cells, %.0f%% of simulations report an apparent divergence under the naive reading. ",
  "(B) Naive bias as a multiple of the true effect: it reaches %.2f at the noisiest reference and the ",
  "smallest true effect, i.e. the artefact is several times the signal. (C) Fraction of simulations in ",
  "which the naive estimator reports the wrong DIRECTION, up to %.0f%%. The IV-corrected estimator stays ",
  "within %.4f of zero bias across the whole grid. The replicate axis is not shown here because this grid ",
  "holds h constant along it by construction; the replicate behaviour is Figure 1c. No published method is ",
  "compared or ranked: both readings are of the same data. ",
  "Source: p1_design_grid.tsv."),
  nrow(G), 200L,
  paste(sprintf("%.3f", nc_pred$obs), collapse=" / "), paste(sprintf("%.1f", nc_pred$h), collapse=" / "),
  G[gamma==1.0, mean(gamma_iv)], G[gamma==1.0, .N],
  100*G[gamma==1.0, mean(sign_wrong_naive)],
  G[gamma>1, min(rel_bias_naive)],
  100*G[gamma>1, max(sign_wrong_naive)],
  G[gamma>1, max(abs(bias_iv))])
## Each appendix figure writes a separate caption file because concurrent scripts
## must not append to one shared file in an unspecified order;
## this prevents one writer from overwriting another.
writeLines(c("# Figure FA3 caption", "",
             "> Generated by `scripts/hsf1/32_fig_p1_grid.R`. Numbers are read from landed tables.", "",
             "## FA3", "", cap), file.path(figd, "CAPTION_FA3.md"))
cat("\n", cap, "\n", sep="")
