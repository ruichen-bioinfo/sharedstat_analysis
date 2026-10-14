#!/usr/bin/env Rscript
## 43_fig_p1_robustness.R -- Supplementary figure: the two estimators when an assumption is violated.
##
## Scope. The closed form uses only second moments and three orthogonalities, and the Lean proofs settle the
## algebra given that structure. The algebra cannot say how far each estimator drifts when one assumption is
## false in a way this design plausibly makes it false. One violation at a time, at a single operating point
## taken from the real data (gamma = 1.3, h = 0.75, near the median measured reliability), reporting only the
## bias of the naive and the instrumented estimator. No new estimand, no claim about the real data, and no
## existing estimate is altered.
##
## The four violations:
##   cov(u,v) != 0    two arms sharing a library share an error component. Section 1 predicts the slope
##                    shifts by cov(u,v)/var(Rhat), so a positive covariance pushes the naive estimate
##                    toward the truth. Both signs are shown, because the real design sits on the forgiving
##                    side and the figure should show that directly.
##   heavy tails      t on 3 df at matched variance. The algebra is distribution-free; this asks whether the
##                    finite-sample behaviour is.
##   unequal blocks   error variances in the two replicate blocks in a 1:3 ratio, violating the equal-
##                    variance assumption behind h = 1 - var(null)/var(signal).
##   block offset     an additive shift between blocks, which the null contrast absorbs as error.

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
## Sourced for DPI and the shared palette, so this figure is exported at 600 dpi like the others. Note: at
## 300 dpi a 175 mm figure looks like an 88 mm one when the width is inferred from the pixel count.
source(file.path(root, "scripts/figures/theme.R"))

## Device base point size. With the annotation cex values used here, an 8 pt base puts text at 4 to
## 6 pt at the final figure width, which is not legible in print. The base scales every size by one
## factor and preserves the relative emphasis; the physical width does not change.
PT_BASE <- 11
fig  <- file.path(root, "doc/04_figures/final")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

GAMMA <- 1.3; HTRUE <- 0.75; NGENE <- 3000L; NSIM <- 400L; SEED <- 20260912L
vR <- 1; vu <- vR * (1 - HTRUE) / HTRUE

one_sim <- function(kind, strength, rho_su = 0) {
  R <- rnorm(NGENE, sd = sqrt(vR))
  if (kind == "heavy") {
    z <- rt(NGENE, df = 3); u <- z / sd(z) * sqrt(vu)
    z2 <- rt(NGENE, df = 3); v <- z2 / sd(z2) * sqrt(vu)
  } else if (kind == "cov") {
    e1 <- rnorm(NGENE); e2 <- strength * e1 + sqrt(max(0, 1 - strength^2)) * rnorm(NGENE)
    u <- e1 * sqrt(vu); v <- e2 * sqrt(vu)
  } else if (kind == "blocks") {
    s1 <- sqrt(2 * vu / (1 + strength)); s2 <- sqrt(2 * vu * strength / (1 + strength))
    u <- (rnorm(NGENE, sd = s1) + rnorm(NGENE, sd = s2)) / sqrt(2)
    v <- rnorm(NGENE, sd = sqrt(vu))
  } else if (kind == "offset") {
    u <- rnorm(NGENE, sd = sqrt(vu)) + strength
    v <- rnorm(NGENE, sd = sqrt(vu))
  } else {
    u <- rnorm(NGENE, sd = sqrt(vu)); v <- rnorm(NGENE, sd = sqrt(vu))
  }
  Rh <- R + u; A <- GAMMA * R + v; e <- A - Rh
  ## The instrument. rho_su is the only violation that can touch the instrumented estimator, because that
  ## estimator assumes S is orthogonal to both errors and uses nothing else. It is the formal version of the
  ## allele-specific mapping bias named in the limitations: a per-gene systematic error that is shared
  ## between hybrids with a common parent, and so is carried by any stratum-level instrument built from them.
  S0 <- 0.8 * R + rnorm(NGENE, sd = 0.6)
  S  <- if (rho_su == 0) S0 else {
    ## inject cor(S, u) = rho_su while holding var(S) fixed, so the panel varies one quantity
    Sr <- resid(lm(S0 ~ u))
    rho_su * sd(S0) * (u - mean(u)) / sd(u) + sqrt(max(0, 1 - rho_su^2)) * sd(S0) * Sr / sd(Sr)
  }
  c(naive = 1 + cov(e, Rh) / var(Rh),
    iv    = 1 + cov(e, S)  / cov(Rh, S),
    pred  = GAMMA * (var(R) / var(Rh)))
}

CONDS <- list(list("none",   0,    "assumptions hold"),
              list("cov",   +0.2,  "cov(u,v) positive"),
              list("cov",   -0.2,  "cov(u,v) negative"),
              list("heavy",  0,    "errors t(3)"),
              list("blocks", 3,    "block variances 1:3"),
              list("offset", 0.5,  "block offset 0.5"))

set.seed(SEED)
RES <- rbindlist(lapply(CONDS, function(cd) {
  M <- t(vapply(seq_len(NSIM), function(i) one_sim(cd[[1]], cd[[2]]), numeric(3)))
  data.table(lab = cd[[3]],
             naive_bias = mean(M[, "naive"]) - GAMMA,
             iv_bias    = mean(M[, "iv"])    - GAMMA,
             naive_sd   = sd(M[, "naive"]),
             iv_sd      = sd(M[, "iv"]))
}))

say("=== bias, gamma = %.1f, h = %.2f, %d genes, %d simulations ===", GAMMA, HTRUE, NGENE, NSIM)
print(RES[, .(lab, naive_bias = round(naive_bias, 4), iv_bias = round(iv_bias, 4),
              naive_sd = round(naive_sd, 4), iv_sd = round(iv_sd, 4))], row.names = FALSE)

base <- RES[lab == "assumptions hold"]
pos  <- RES[lab == "cov(u,v) positive"]
neg  <- RES[lab == "cov(u,v) negative"]
say("\n=== validation checks ===")
say("  closed form predicts a naive bias of gamma*h - gamma = %.4f ; observed %.4f",
    GAMMA * HTRUE - GAMMA, base$naive_bias)
stopifnot(abs(base$naive_bias - (GAMMA * HTRUE - GAMMA)) < 0.02)
say("  the instrumented estimator stays within 0.03 of unbiased in every condition : %s",
    all(abs(RES$iv_bias) < 0.03))
stopifnot(all(abs(RES$iv_bias) < 0.03))
say("  a positive cov(u,v) moves the naive estimate toward the truth (%.4f > %.4f) : %s",
    pos$naive_bias, base$naive_bias, pos$naive_bias > base$naive_bias)
stopifnot(pos$naive_bias > base$naive_bias, neg$naive_bias < base$naive_bias)

## ---- panel (b): the one violation the instrumented estimator cannot absorb -------------------------------
## Everything in panel (a) leaves the instrumented estimate alone, which is expected rather than reassuring:
## those are assumptions of the naive model. The question a reader should be able to answer from the figure is
## how large a correlation between the instrument and the reference-arm error has to be before the correction
## stops being worth making. This sweeps it and reports the fraction of the correction that survives.
RHO <- c(0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4)
set.seed(SEED + 1L)
SW <- rbindlist(lapply(RHO, function(r) {
  M <- t(vapply(seq_len(NSIM), function(i) one_sim("none", 0, rho_su = r), numeric(3)))
  data.table(rho = r, iv_bias = mean(M[, "iv"]) - GAMMA, iv_sd = sd(M[, "iv"]),
             naive_bias = mean(M[, "naive"]) - GAMMA)
}))
SW[, surviving := 1 - abs(iv_bias) / abs(naive_bias)]

say("\n=== cor(S, u) sweep: the only violation that reaches the instrumented estimator ===")
print(SW[, .(rho, iv_bias = round(iv_bias, 4), iv_sd = round(iv_sd, 4),
             fraction_of_correction_surviving = round(surviving, 3))], row.names = FALSE)

say("\n=== validation checks on the sweep ===")
stopifnot(abs(SW[rho == 0, iv_bias]) < 0.03)
say("  at rho = 0 the instrumented estimator is unbiased (%.4f)", SW[rho == 0, iv_bias])
stopifnot(all(diff(abs(SW$iv_bias)) > -0.005))
say("  |bias| is monotone in rho, so the panel reads as a single direction : TRUE")
half <- SW[surviving < 0.5, min(rho)]
say("  more than half the correction is lost once cor(S, u) exceeds %.3f", half)
stopifnot(is.finite(half))

fwrite(RES, file.path(root, "data/02_processed/hsf1/p1_robustness_grid.tsv"), sep = "\t")
fwrite(SW,  file.path(root, "data/02_processed/hsf1/p1_robustness_rho_sweep.tsv"), sep = "\t")

## Note: gamma in plotmath is set in the symbol font, which is upright Noto Sans here and reads as "y"; it is drawn
## as an italic Noto Sans gamma instead (see 38_fig_p1_theory_and_reliability.R for the measurement).
## Drawn in the shared figure style (scripts/figures/fig_style.R). No panel titles (panels are named in the caption);
## colour is COL_SCI: naive dark red, instrumented (the corrected estimator) dark green, the null black.
source(file.path(root, "scripts/figures/fig_style.R"))
draw <- function() {
  layout(matrix(1:2, 1, 2), widths = c(1.32, 1)); after_layout()
  ## (A) one violation at a time
  par(mar = c(3.4, 11.0, 3.4, 0.8))
  y  <- rev(seq_len(nrow(RES)))
  xl <- range(c(RES$naive_bias, RES$iv_bias)) + c(-0.04, 0.06)
  plot(NA, xlim = xl, ylim = c(0.5, nrow(RES) + 1.2), yaxt = "n", ylab = "", xlab = "bias of the estimate")
  segments(0, 0.5 - 1, 0, nrow(RES) + 0.6, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  segments(GAMMA * HTRUE - GAMMA, 0.5 - 1, GAMMA * HTRUE - GAMMA, nrow(RES) + 0.6, col = COL_SCI$uncertainty, lty = 3,
           lwd = LW$ref)
  ## "γh − γ" as one italic term, a true minus (U+2212), and no gap between γ and h (otherwise read as "γ h − γ")
  text(GAMMA * HTRUE - GAMMA, nrow(RES) + 0.85, expression(paste(italic("\u03b3h"), " \u2212 ", italic("\u03b3"))),
       cex = CEX_TXT, col = COL_SCI$neutral, pos = 4, offset = 0.1)
  segments(RES$naive_bias, y, RES$iv_bias, y, col = COL_SCI$uncertainty, lwd = LW$ci)
  points(RES$naive_bias, y, pch = 21, bg = COL_SCI$naive, cex = PCH_CEX, lwd = PCH_LWD)
  points(RES$iv_bias, y, pch = 22, bg = COL_SCI$corrected, cex = PCH_CEX, lwd = PCH_LWD)
  axis(2, at = y, labels = RES$lab, tick = FALSE)
  legend_top(c("naive", "instrumented"), pch = c(21, 22), pt.bg = c(COL_SCI$naive, COL_SCI$corrected),
             pt.cex = PCH_CEX, pt.lwd = PCH_LWD)
  panel_tag("A")
  ## (B) the instrumented estimator against its correlation with the reference-arm error
  par(mar = c(3.4, 3.8, 3.4, 0.8))
  plot(NA, xlim = range(SW$rho), ylim = range(c(SW$iv_bias, base$naive_bias, 0)),
       xlab = expression(paste("cor(", italic(S), ", ", italic(u), ")")), ylab = "bias of the estimate")
  abline(h = 0, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
  abline(h = base$naive_bias, col = COL_SCI$naive, lwd = LW$data)
  lines(SW$rho, SW$iv_bias, col = COL_SCI$corrected, lwd = LW$data)
  points(SW$rho, SW$iv_bias, pch = 22, bg = COL_SCI$corrected, cex = PCH_CEX, lwd = PCH_LWD)
  legend_top(c("instrumented", "naive"), pch = c(22, NA), lty = c(1, 1), col = c(COL_SCI$corrected, COL_SCI$naive),
             pt.bg = COL_SCI$corrected, pt.cex = PCH_CEX, lwd = LW$data)
  panel_tag("B")
}
W <- 174 / 25.4
fig_both("FS2_robustness", 174, 74, draw, fig)
say("\nwritten: FS2_robustness.{pdf,png}  (%.1f mm wide)", W * 25.4)
