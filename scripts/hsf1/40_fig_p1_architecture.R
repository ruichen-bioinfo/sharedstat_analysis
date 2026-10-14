#!/usr/bin/env Rscript
## 40_fig_p1_architecture.R -- the shared-component architecture, as a concept diagram.
##
## Purpose. Every result in P1 follows from one arrangement: a single measured quantity is used twice, once
## inside the contrast and once as the variable the contrast is examined against. The paper states this
## in prose and then works algebraically. A reader meeting it for the first time has to hold the
## arrangement in their head while following the algebra. This draws the arrangement.
##
## Supplementary Figure S1. It is generated as supplementary material.
##  It is
## 
##
## No new numbers. The diagram contains no data. The only quantities on it are the two symbolic results
## already in the manuscript, gamma*h and h < 1/gamma, and they are drawn as symbols, not as values.

root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
fig  <- file.path(root, "doc/04_figures/final")
DPI  <- 300
say  <- function(...) cat(sprintf(...), "\n")

## Drawn in the shared figure style (scripts/figures/fig_style.R). Headings and box labels 8.5 pt, all
## other text 7.5 pt; boxes 0.55 pt, arrows 0.45 pt, guides 0.35 pt; greys are the style's (COL_SCI).
source(file.path(root, "scripts/figures/fig_style.R"))
GREY <- COL_SCI$neutral; LGREY <- COL_SCI$uncertainty; ACC <- COL_SCI$truth

box <- function(x, y, w, h, lab, sub = NULL, lwd = LW$data, border = ACC, bg = "white", cex = CEX_LAB, tcol = ACC, font = 1) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2, col = bg, border = border, lwd = lwd)
  if (is.null(sub)) {
    text(x, y, lab, cex = cex, col = tcol, font = font)
  } else {
    text(x, y + h * 0.16, lab, cex = cex, col = tcol, font = font)
    ## The sub-label is set at CEX_TXT, not as a multiplier on the caller's cex (cex * 0.74 at cex = 0.68 gives
    ## 5.5 pt), so no caller can push it below the 7 pt print minimum. At 7 pt a sub-label can be wider than its
      ## box and spill over the arrow captions beside it, so its width is checked against the box and a sub-label
      ## that does not fit stops the script instead of being drawn across its neighbour.
      scex <- CEX_TXT
      if (strwidth(sub, units = "user", cex = scex) > w * 0.94)
          stop(sprintf(paste0("box sub-label is wider than its box and would overlap what is beside it: '%s' needs ",
                              "%.1f user units at cex %.2f, the box is %.1f wide"),
                       sub, strwidth(sub, units = "user", cex = scex), scex, w))
      text(x, y - h * 0.22, sub, cex = scex, col = GREY)
  }
}

arr <- function(x0, y0, x1, y1, lab = NULL, col = ACC, lty = 1, lwd = LW$arrow, pos = 3, cex = CEX_TXT,
                off = 0.3, curve = 0) {
  if (curve == 0) {
    arrows(x0, y0, x1, y1, length = 0.07, col = col, lty = lty, lwd = lwd)
    if (!is.null(lab)) text((x0 + x1)/2, (y0 + y1)/2, lab, cex = cex, col = col, pos = pos, offset = off)
  } else {
    ## a quadratic bend, so the two paths out of the reference do not overlap
    tt <- seq(0, 1, length.out = 60)
    mx <- (x0 + x1)/2 + curve; my <- (y0 + y1)/2
    bx <- (1-tt)^2 * x0 + 2*(1-tt)*tt*mx + tt^2 * x1
    by <- (1-tt)^2 * y0 + 2*(1-tt)*tt*my + tt^2 * y1
    lines(bx, by, col = col, lty = lty, lwd = lwd)
    arrows(bx[length(bx)-1], by[length(by)-1], x1, y1, length = 0.07, col = col, lty = lty, lwd = lwd)
    if (!is.null(lab)) text(mx, my, lab, cex = cex, col = col, pos = pos, offset = off)
  }
}

## gamma in plotmath is set in the symbol font, which is upright Noto Sans here and reads as "y"; it is drawn
## as an italic Noto Sans gamma instead (see 38_fig_p1_theory_and_reliability.R for the measurement).
draw <- function() {
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot(NA, xlim = c(0, 100), ylim = c(0, 58), axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")

  ## ---------- left: what is measured ----------
  text(15, 55, "what is measured", font = 2, cex = CEX_LAB)
  box(15, 45.5, 25, 9,  expression(hat(R) == R + u),
      "reference arm", bg = TINT$red, border = COL_SCI$naive)
  box(15, 31.5, 25, 9,  expression(A == italic("\u03b3") * R + v),
      "response arm", bg = "white")

  ## ---------- middle: the arrangement ----------
  text(50, 55, "the arrangement", font = 2, cex = CEX_LAB)
  box(50, 38.5, 24, 8.5, expression(e == A - hat(R)), "the contrast", bg = "white")
  box(50, 20, 24, 8.5, expression(paste("regress  ", italic(e), "  on  ", hat(R))), bg = "white")

  ## the same measured quantity travels along both paths: that is the whole architecture
  arr(27.5, 46.5, 38, 41, NULL, col = ACC)
  arr(27.5, 31.5, 38, 36.5, NULL, col = ACC)
  ## The reference is used a second time, as the regressor. Routed to the right of the response box so it
  ## crosses nothing: a straight drop from the reference would pass through the response box and read as a
  ## link between the two arms, which is not what the architecture says.
  arr(28.5, 43, 37.8, 21.5, NULL, col = COL_SCI$naive, curve = 5.5)
  ## The label sits in the empty band below the curve, ending just left of the regress box and below the end of the
  ## curved arrow it describes: on the curve it would collide with the response arm's arrow and the box edge, and at
  ## 7.5 pt two lines do not fit between the regress box and the identification rule below it.
  text(37.0, 16.4, "the same measurement,\nused a second time", cex = CEX_TXT, col = COL_SCI$naive, adj = c(1, 0.5))
  arr(50, 34.2, 50, 24.3, NULL, col = LGREY, lwd = LW$ref)

  ## the error enters twice with opposite signs, which is the mechanism
  ## Two lines, spaced by the measured line height rather than a fixed offset: at this type size a line is taller
  ## than 3.2 user units, so a fixed 3.2 spacing makes them run into each other and into the box sub-label to their
  ## left. The height is asserted positive.
  ARRH <- strheight("Mg", units = "user", cex = CEX_TXT)
  stopifnot(ARRH > 0)
  ay <- 50.0
  text(50, ay, expression(paste("u enters e with a minus sign and ", hat(R), " with a plus sign")),
       cex = CEX_TXT, col = GREY)
  text(50, ay - 1.75 * ARRH, "so the predictor is correlated with the error in the response", cex = CEX_TXT, col = GREY)

  ## ---------- right: the consequence ----------
  text(85, 55, "the consequence", font = 2, cex = CEX_LAB)
  box(85, 38.5, 26, 9, expression(paste(hat(italic("\u03b3")), " (naive) = 1 + slope ", "" == italic("\u03b3") * h)),   # shows the 1 + slope step
      expression(paste("h = var(R) / var(", hat(R), ")")), bg = TINT$red, border = COL_SCI$naive)
  arr(62.2, 38.5, 71.8, 38.5, NULL, col = ACC)

  box(85, 20, 26, 12, "", bg = "white", border = LGREY)
  text(85, 24.2, "two consequences", cex = CEX_TXT, col = ACC)
  text(85, 21.2, expression(paste(italic("\u03b3") == 1, "  gives  ", hat(italic("\u03b3")) == h)), cex = CEX_TXT)
  text(85, 17.6, expression(paste("direction inverts once  ", h < 1/italic("\u03b3"))), cex = CEX_TXT)
  arr(85, 33.8, 85, 26.3, NULL, col = LGREY, lwd = LW$ref)

  ## ---------- the identification route, kept visually separate ----------
  segments(4, 10.5, 96, 10.5, col = LGREY, lty = 3, lwd = LW$ref)
  text(6, 7.6, "identification", font = 2, cex = CEX_LAB, adj = 0, col = COL_SCI$corrected)
  text(6, 4.4, adj = 0, cex = CEX_TXT, col = COL_SCI$corrected,
       labels = expression(paste("an instrument ", italic(S), " correlated with ", italic(R),
                                 " but not with ", italic(u), " or ", italic(v), " recovers ", italic("\u03b3"),
                                 " without knowing ", italic(h))))
  ## plotmath sets a [subscript] at about 0.7 of its body, so `hat(gamma)[IV]` would put the IV at 0.72 x 11 x 0.7 =
  ## 5.5 pt, and raising the body to clear 7 pt makes the line wide enough to collide with the sentence to its left,
  ## which shares this baseline. The estimator is named in a full-size parenthetical instead, so nothing is below 7 pt.
  text(96, 3.5, adj = 1, cex = CEX_TXT, col = COL_SCI$corrected,
       labels = expression(paste(hat(italic("\u03b3")), " (IV) = 1 + cov(e, S) / cov(", hat(R), ", S) = ", italic("\u03b3"))))
}

W <- 174 / 25.4; H <- 3.4
fig_both("FS1_architecture", 174, 86, draw, fig)

say("=== written ===")
for (f in c("FS1_architecture.pdf", "FS1_architecture.png")) {
  p <- file.path(fig, f); say("  %-28s %8d bytes", f, file.info(p)$size)
}
say("  width %.2f in = %.1f mm; contains no data values", W, W * 25.4)
