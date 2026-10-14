## Figure width 6.9 in: a target journal caps double-column figures at 178 mm = 7.008 in. Set at
## source rather than rescaled at submission, because rescaling a raster forfeits the >= 300 dpi
## the same requirement imposes.
# 24_fig_a2_ladder.R -- P1 appendix figure FA2, the estimator ladder.
# Shows, on simulated data with a known true correlation, what each of the three estimators returns.
# It shows why the r_mech inference does not hold: the naive estimator is biased toward
# compensation, the cross-block estimator is unbiased, and cor(cis,par) is positive at zero
# compensation, so the condition cov(cis,par)=0 is near-complete compensation, not zero compensation.
# Panel b shows that the true correlation is scale-invariant, which rules out a second prediction.
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
fig <- file.path(root, "doc/04_figures/final"); hsd <- file.path(root, "data/02_processed/hsf1")
MIRROR <- "logs/output_copy"; DPI <- 300; set.seed(20260906L)
## Colours come from the shared theme, so a colour means the same thing in every figure of every
## manuscript. A local copy of the Okabe-Ito vector would drift from the theme unnoticed.
source(file.path(root, "scripts/figures/theme.R"))

## Device base point size 11. With the annotation cex values used here, an 8 pt base puts text at
## 4 to 6 pt at the final figure width, which is not legible in print. Raising the base moves every
## size by one factor and preserves the relative emphasis; the physical width does not change, so
## this is not a figure shrunk to make its type look bigger.
PT_BASE <- 11
N <- 200000L; SEC <- 0.6; SEP <- 0.6
sim <- function(rho, sc = 1, st = 1) {
  Z <- matrix(rnorm(2 * N), N, 2)
  c_ <- sc * Z[, 1]; t_ <- st * (rho * Z[, 1] + sqrt(max(0, 1 - rho^2)) * Z[, 2]); p_ <- c_ + t_
  ch1 <- c_ + rnorm(N, 0, SEC); ch2 <- c_ + rnorm(N, 0, SEC)
  ph1 <- p_ + rnorm(N, 0, SEP); ph2 <- p_ + rnorm(N, 0, SEP)
  c(naive = cor(ch1, ph1 - ch1), cross = cor(ch1, ph2 - ch2), cp = cor(ch1, ph2))
}
RHO <- seq(0, -0.99, by = -0.11)
S <- rbindlist(lapply(RHO, function(r) as.list(c(rho = r, sim(r)))))
fwrite(S, file.path(hsd, "p1_a2_estimator_ladder.tsv"), sep = "\t")
# scale-invariance check
# Two curves are needed here and must not be conflated: the true correlation is scale-invariant
# (flat by construction), but the observed cross-block estimate is attenuated by measurement error,
# and the attenuation depends on the variance ratio:
#   cor(chat1, that2) = rho * sd(c)*sd(t) / [ sqrt(var(c)+se_c^2) * sqrt(var(t)+se_p^2+se_c^2) ]
# As var(trans) grows at fixed noise, the t-side attenuation shrinks and |cor| rises toward rho.
V <- rbindlist(lapply(c(1, 4, 16, 64), function(v) {
  st <- sqrt(v)
  att <- (1 * st) / (sqrt(1 + SEC^2) * sqrt(v + SEP^2 + SEC^2))    # analytic attenuation factor
  data.table(var_ratio = v, cor_true = -0.5, cor_ct = sim(-0.5, sc = 1, st = st)[["cross"]],
             cor_pred = -0.5 * att)
}))
fwrite(V, file.path(hsd, "p1_a2_scale_invariance.tsv"), sep = "\t")
cat("scale invariance: var(trans)/var(cis) =", paste(V$var_ratio, collapse = " / "),
    "-> cross-block cor =", paste(sprintf("%.3f", V$cor_ct), collapse = " / "), "\n")
draw <- function() {
  layout(matrix(1:2, 1, 2), widths = c(1.25, 1))
  ## layout() silently sets cex to 0.66, the same reduction par(mfrow=) documents. Reset it so the cex values
    ## below sit on the device base this script chose.
  par(cex = 1)
  ## An outer margin keeps content clear of the bottom trim at the larger type size.
  par(oma = c(0.7, 0.2, 0.2, 0.3))
  par(mgp = c(1.9, 0.55, 0), tcl = -0.25, cex.axis = 0.85, cex.lab = 0.95, font.main = 1,
      cex.main = 0.95, mar = c(4.4, 3.9, 2.4, 0.6))
  plot(NA, xlim = c(0, -1), ylim = c(-1, 0.75), xlab = "true cor(cis, trans)",
       ylab = "estimated value", main = "a  three estimators against a known truth")
  abline(h = 0, col = OI[["grey"]], lty = 3); abline(0, 1, col = OI[["grey"]], lty = 2)
  lines(S$rho, S$naive, type = "b", pch = 22, bg = OI[["vermillion"]], col = OI[["vermillion"]], cex = 1.1)
  lines(S$rho, S$cross, type = "b", pch = 21, bg = OI[["blue"]], col = OI[["blue"]], cex = 1.1)
  lines(S$rho, S$cp,    type = "b", pch = 24, bg = OI[["green"]], col = OI[["green"]], cex = 1.1)
  r0 <- S[which.min(abs(rho))]
  points(r0$rho, r0$naive, pch = 22, bg = "white", col = OI[["vermillion"]], cex = 2.0, lwd = 2)
  text(-0.02, r0$naive - 0.10, sprintf("naive at true 0:\n%.3f (biased)", r0$naive), cex = 0.64, adj = 0)
  text(-0.02, r0$cp + 0.09, sprintf("cor(cis,par) at true 0:\n%+.3f, and positive", r0$cp), cex = 0.64, adj = 0)
  z <- S[which.min(abs(cp))]
  arrows(z$rho, 0.26, z$rho, 0.05, length = 0.06, col = OI[["purple"]], lwd = 2)
  text(z$rho + 0.30, 0.40, sprintf("cor(cis,par) = 0 at true rho = %.2f:\nnear-complete compensation,\nnot an absence of it", z$rho),
       cex = 0.64, col = OI[["purple"]], adj = 0.5)
  legend("bottomright", bty = "n", cex = 0.68,
         legend = c("naive (same block both sides)", "cross-block (independent blocks)",
                    "cor(cis, par)", "identity line"),
         pch = c(22, 21, 24, NA), pt.bg = c(OI[["vermillion"]], OI[["blue"]], OI[["green"]], NA),
         col = c(OI[["vermillion"]], OI[["blue"]], OI[["green"]], OI[["grey"]]), lty = c(1, 1, 1, 2))
  par(mar = c(4.4, 3.9, 2.4, 0.6))
  plot(V$var_ratio, V$cor_true, log = "x", type = "b", pch = 24, bg = OI[["vermillion"]],
       col = OI[["vermillion"]], ylim = c(-0.62, -0.25), xlab = "var(trans) / var(cis)",
       ylab = "cor(cis, trans)", main = "b  truth invariant, estimate attenuated",
       cex = 1.2, lty = 2)
  lines(V$var_ratio, V$cor_ct, type = "b", pch = 21, bg = OI[["blue"]], col = OI[["blue"]], cex = 1.3)
  lines(V$var_ratio, V$cor_pred, type = "l", col = OI[["green"]], lwd = 2, lty = 3)
  legend("bottomright", bty = "n", cex = 0.64,
         legend = c("true rho (invariant by construction)", "cross-block estimate (simulated)",
                    "analytic attenuation prediction"),
         pch = c(24, 21, NA), pt.bg = c(OI[["vermillion"]], OI[["blue"]], NA),
         col = c(OI[["vermillion"]], OI[["blue"]], OI[["green"]]), lty = c(2, 1, 3), lwd = c(1,1,2))
  mtext("the ratio does not amplify the true |cor|;", 1, 3.0, cex = 0.64)
  mtext("it changes how much error attenuates the estimate", 1, 3.7, cex = 0.64)
}
for (d in c("pdf","png")) {
  if (d == "pdf") pdf(file.path(fig, "FA2_estimator_ladder.pdf"), width = 6.9, height = 3.4,
                      pointsize = PT_BASE, useDingbats = FALSE)
  else png(file.path(fig, "FA2_estimator_ladder.png"), width = round(6.9*DPI), height = round(3.4*DPI),
           res = DPI, pointsize = PT_BASE)
  draw(); dev.off()
}
cap <- sprintf("FA2. What each estimator returns when the truth is known (simulation, n = 2 x 10^5 per point, measurement SD 0.6 on both axes). (a) At true cor(cis,trans) = 0 the naive estimator, which uses the same replicate block on both sides, returns %.3f -- biased toward apparent compensation, because trans is computed as par - cis and inherits -var(e_cis). The cross-block estimator, using independent blocks for the two sides, returns %.3f and is unbiased, so its null hypothesis is 0. cor(cis,par) is %+.3f at true zero compensation, i.e. positive; it reaches 0 only near a true correlation of %.2f, which is near-complete compensation. The value of cor(cis,trans) at cov(cis,par) = 0 is near-complete compensation, not 'zero compensation', which is the opposite end of the range. (b) Two distinct facts. The true correlation is scale-invariant: holding it at -0.5 while raising var(trans)/var(cis) from 1 to 64 leaves it at -0.5 by construction, so the variance ratio cannot amplify the true |cor|, . The observed cross-block estimate is attenuated by measurement error and the attenuation depends on the variance ratio: it runs %.3f to %.3f across the same range, tracking the analytic factor rho * sd(c)sd(t) / [sqrt(var(c)+se_c^2) sqrt(var(t)+se_p^2+se_c^2)]. Consequence for any reported compensation strength: a cross-block correlation is a lower bound on |rho|, and designs with more replicates attenuate less, so compensation magnitudes are not directly comparable across studies with different replication. Sources: p1_a2_estimator_ladder.tsv, p1_a2_scale_invariance.tsv.",
  S[which.min(abs(rho)), naive], S[which.min(abs(rho)), cross], S[which.min(abs(rho)), cp],
  S[which.min(abs(cp)), rho],
  ## Ordered by variance ratio, not by min and max. These are negative numbers, so min() returns the
  ## value at the largest ratio and max() the value at the smallest, which would print the range in the
  ## reverse of the direction the sentence describes and the mechanism gives.
  V[var_ratio == min(var_ratio), cor_ct], V[var_ratio == max(var_ratio), cor_ct])
## The shared caption file FIGURE_CAPTIONS.md is written only by 19_figures.R. Several scripts
## appending to one file in an unspecified order let the last writer win and erase the others, so
## each appendix figure writes its own caption file.
writeLines(c("# Figure FA2 caption", "",
             "> Generated by `scripts/hsf1/24_fig_a2_ladder.R`. Numbers are read from landed tables, and the\n             two ends of a signed range are taken by sorting on the independent variable, never with min()/max().", "",
             "## FA2", "", cap), file.path(fig, "CAPTION_FA2.md"))
for (f in list.files(fig, pattern = "^(FA2|CAPTION_FA2)", full.names = TRUE))
  file.copy(f, file.path(MIRROR, "doc/04_figures/final", basename(f)), overwrite = TRUE)
for (f in list.files(hsd, pattern = "^p1_a2", full.names = TRUE))
  file.copy(f, file.path(MIRROR, "data/02_processed/hsf1", basename(f)), overwrite = TRUE)
file.copy(file.path(root,"scripts/hsf1/24_fig_a2_ladder.R"), file.path(MIRROR,"scripts/hsf1/"), overwrite=TRUE)
cat("FA2 written
")
