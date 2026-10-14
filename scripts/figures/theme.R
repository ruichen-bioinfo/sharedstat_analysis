## theme.R -- the house style for every figure in the three manuscripts.
##
## Base-R defaults applied one script at a time give working plots rather than journal figures: sentences printed
## inside panels, panel titles that state conclusions, legends carrying analysis variable names (`Msn24_only`,
## `cis_not_induced`, `Sc`), tick labels around 5 pt at print size, and several ways of drawing a confidence
## interval. Each of these is fixed here once, so that the fix applies to all three manuscripts.
##
## Rules, taken from what GBE, G3 and NAR Genomics and Bioinformatics all ask for:
##   - draw at the final physical width: 88 mm single column, 180 mm double column. Figures are drawn at their final physical width.
##   - one sans-serif family; no text below 7 pt at print size; axis titles 8-9 pt; panel letters 10 pt bold.
##   - panel letters A, B, C at the top left. Panel titles use neutral phrases rather than conclusions, or
##     are omitted. Conclusions are stated in the caption.
##   - The plotting area contains no sentences, footnotes, "NOT", "=>", or sample-size asides.
##     Information needed to interpret a panel is provided in its caption.
##   - legend labels use reader-facing terminology. Species in italic binomials. Dependence classes spelled
##     as prose. No underscores.
##   - one colour-blind-safe palette (Okabe-Ito) with fixed assignments: each species always the same colour in
##     every figure of every manuscript; each dependence class likewise; cis and trans a fixed pair.
##   - continuous quantities use one sequential ramp, and that ramp is never mixed with the categorical colours.
##   - confidence intervals are drawn one way: a thin bar in the same colour as the point, or a translucent band
##     of the same hue. Intervals do not introduce a third colour.
##   - reference lines are thin grey dashes, one style everywhere.
##
## Shared settings are defined here because local settings
## would make the styles diverge.

## ------------------------------------------------------------------ device geometry ----
MM  <- 1 / 25.4                       # millimetres to inches
W1  <- 88  * MM                       # single column
W2  <- 180 * MM                       # double column (GBE allows 185; 180 is safe for all three journals)
DPI <- 600                            # raster export for journals that ask for TIFF/PNG
## Base pointsize; everything else is a multiple of it. With the annotation cex values used across the figure
## scripts (0.52 to 0.72), an 8 pt base puts text at 4 to 6 pt at the final 175 mm width (measured from the PDF
## content stream), so the base is 11. Changing the base moves every size by one factor, so the relative emphasis of
## titles, axis labels and annotations is unchanged, and the physical figure width does not move.
PT  <- 11

## Font sizes as multipliers of PT, chosen so that the smallest text is 7 pt at print size.
## Text sizes only. `point` below is a plotting-symbol size, not a font size, and is kept in a separate list so
## that a check on the smallest text does not have to reason about which entries are text.
## The comments give the resulting point size at the PT above and must be updated when PT changes.
CEX <- list(axis = 0.90,              # 9.9 pt tick labels
            lab  = 1.05,              # 11.6 pt axis titles
            panel = 1.25,             # 13.8 pt panel letters, drawn bold
            legend = 0.90)            # 9.9 pt
SYM <- list(point = 0.85)             # symbol expansion, not text

## Assert the rule this file exists to enforce, at source time: no text below 7 pt at print size.
stopifnot(all(PT * unlist(CEX) >= 7))

## ------------------------------------------------------------------ palette ----
OI <- c(black = "#000000", orange = "#FF7F0E", skyblue = "#17BECF", green = "#2CA02C",
        yellow = "#FF7F0E", blue = "#1F77B4", vermillion = "#D62728", purple = "#9467BD",
        grey = "#999999")

## Fixed assignments. A species is one colour in every figure in every manuscript; the same for a class.
SPECIES_COL <- c("S. cerevisiae"   = OI[["vermillion"]],
                 "S. paradoxus"    = OI[["orange"]],
                 "S. kudriavzevii" = OI[["blue"]],
                 "S. uvarum"       = OI[["skyblue"]])
SPECIES_PCH <- c("S. cerevisiae" = 21, "S. paradoxus" = 22,
                 "S. kudriavzevii" = 23, "S. uvarum" = 24)

CLASS_COL <- c("Hsf1-dependent core"            = "#9467BD",
               "Hsf1-joint"                     = "#B595D1",
               "Hsf1-only"                      = "#9467BD",
               "Msn2/4-only"                    = "#FF7F0E",
               "both required"                  = "#543005",
               "transcription-factor-independent" = "#636363",
               "other induced"                  = "#A6A6A6",
               "not induced"                    = "#FFFFFF")

MODULE_COL <- c("carbohydrate metabolism" = OI[["vermillion"]],
                "proteostasis"            = OI[["blue"]],
                "translation and ribosome" = OI[["green"]],
                "cell cycle and growth"   = OI[["orange"]],
                "oxidative stress"        = OI[["purple"]])

AXIS_COL <- c(cis = OI[["vermillion"]], trans = OI[["blue"]], total = OI[["grey"]])

## One sequential ramp for continuous quantities, never mixed with the categorical colours above.
ramp_seq <- function(n) grDevices::colorRampPalette(c("#F7F7F7", "#92C5DE", "#2166AC"))(n)

## The label a raw name should be printed as. A script that plots a name not in here gets an error, which is the
## point: a new class must be given a printable label before it can appear in a figure.
LABEL <- c(Sc = "S. cerevisiae", Sp = "S. paradoxus", Sk = "S. kudriavzevii", Su = "S. uvarum",
           Hsf1_core = "Hsf1-dependent core", Hsf1_joint = "Hsf1-joint", Hsf1_only = "Hsf1-only",
           Msn24_only = "Msn2/4-only", both_required = "both required",
           ## Printed as the manuscripts name them: P2 says "redundantly controlled", and "unassigned" is the
           ## group where the three DE methods do not agree on a class.
           redundant = "redundantly controlled", unassigned = "unassigned",
           TF_independent = "transcription-factor-independent",
           other_induced = "other induced", not_induced = "not induced",
           metabolism = "carbohydrate metabolism", proteostasis = "proteostasis",
           translation_ribosome = "translation and ribosome",
           cell_cycle_growth = "cell cycle and growth", oxidative = "oxidative stress",
           cis = "cis", trans = "trans", total = "total",
           ## Genotypes and conditions of the perturbation dataset. The raw names are the column keys of the
           ## count matrix and must stay as they are in the data; only what is printed changes.
           Wildtype = "wild type", HSF1.KD = "Hsf1-depleted", MSN24.KO = "msn2/4 deleted",
           Double.KDKO = "both perturbed",
           "25_0" = "25 \u00b0C baseline", "37_10" = "37 \u00b0C, 10 min", "37_30" = "37 \u00b0C, 30 min",
           "42_10" = "42 \u00b0C, 10 min", "42_30" = "42 \u00b0C, 30 min")

lab_of <- function(x) {
    miss <- setdiff(unique(x), names(LABEL))
    if (length(miss))
        stop(sprintf("no printable label for: %s -- add it to LABEL in theme.R rather than printing the raw name",
                     paste(miss, collapse = ", ")))
    unname(LABEL[x])
}

## Species names are italic in print. plotmath needs an expression, so build one per label.
italic_species <- function(x) {
    sapply(x, function(s) if (grepl("^S\\. ", s)) as.expression(bquote(italic(.(s)))) else as.expression(s))
}

## ------------------------------------------------------------------ device ----
## Draw one figure to PDF (vector, for submission) and PNG (raster, for reading copies), at the final width.
fig_open <- function(dir, tag, width_mm, height_mm, dev) {
    w <- width_mm * MM; h <- height_mm * MM
    if (dev == "pdf")
        grDevices::pdf(file.path(dir, paste0(tag, ".pdf")), width = w, height = h,
                       pointsize = PT, useDingbats = FALSE)
    else
        grDevices::png(file.path(dir, paste0(tag, ".png")), width = round(w * DPI), height = round(h * DPI),
                       res = DPI, pointsize = PT)
    graphics::par(mgp = c(1.8, 0.5, 0), tcl = -0.22, cex.axis = CEX$axis, cex.lab = CEX$lab,
                  font.main = 1, cex.main = CEX$lab, las = 1, xaxs = "r", yaxs = "r",
                  family = "sans", lend = "butt")
}

draw_figure <- function(dir, tag, width_mm, height_mm, fn, say = message) {
    for (dev in c("pdf", "png")) {
        fig_open(dir, tag, width_mm, height_mm, dev)
        ok <- try(fn(), silent = TRUE)
        grDevices::dev.off()
        if (inherits(ok, "try-error")) {
            say(sprintf("  !! %s FAILED: %s", tag, conditionMessage(attr(ok, "condition"))))
            return(FALSE)
        }
    }
    say(sprintf("  %-22s %.0f x %.0f mm", tag, width_mm, height_mm))
    TRUE
}

## ------------------------------------------------------------------ panel furniture ----
## A panel letter, top left, outside the plotting region. No title argument: a panel either carries a neutral
## phrase set with `sub_title()` or carries nothing.
## The letter goes at the top left of the panel, including its margins, not of the device. An offset computed from
## the plotting region in user coordinates becomes enormous with a wide left margin and is clipped, which puts both
## letters of a two-panel figure at the device's top left corner, one over the other. Converting from the figure
## region handles any margin width and any log axis.
panel_letter <- function(letter, dx_mm = 1.2, dy_mm = 1.2) {
    ## grconvert* live in `graphics`, not `grDevices`. Place the letter relative to the panel's figure region so
    ## that any margin width works and each letter of a multi-panel figure lands in its own panel.
    x0 <- graphics::grconvertX(0, from = "nfc", to = "user")
    y1 <- graphics::grconvertY(1, from = "nfc", to = "user")
    xin <- graphics::par("pin")[1]; yin <- graphics::par("pin")[2]
    xr <- diff(graphics::par("usr")[1:2]); yr <- diff(graphics::par("usr")[3:4])
    if (graphics::par("xlog")) xr <- diff(10^graphics::par("usr")[1:2])
    graphics::par(xpd = NA)
    graphics::text(x0, y1 - dy_mm / 25.4 / yin * yr, letter, cex = CEX$panel, font = 2, adj = c(0, 1))
    graphics::par(xpd = FALSE)
}

sub_title <- function(txt) graphics::mtext(txt, side = 3, line = 0.3, adj = 0, cex = CEX$axis)

## Does a panel title fit the panel it belongs to?
##
## mtext clips to the figure region, not to the page, so a title longer than its own panel silently loses its tail
## ("noise-corrected amplitude at 60 mi") while the PDF text layer still holds the whole string; the type-size,
## trim-edge and text-layer checks all pass it. Detecting this from the rendered page is unreliable: the word boxes
## pdftotext reports include side bearing, and an ordinary word's ink stops 2 to 5 pt short of its box, which is the
## same signal as a clipped glyph. Here the arithmetic is exact, because strwidth and the figure width are both known
## at draw time.
## Is a y position inside the figure region, so that what is drawn there will actually appear?
##
## A key placed below the axis with xpd = NA is clipped at the figure region, and clipping is silent: a key placed
## 5.4 lines below the axis can vanish entirely while the script reports success. Nothing in the type-size or edge
## measurements can see a legend that was never drawn.
assert_y_visible <- function(y, what = "annotation") {
    lo <- graphics::grconvertY(0, "nfc", "user")
    hi <- graphics::grconvertY(1, "nfc", "user")
    if (y < min(lo, hi) || y > max(lo, hi))
        stop(sprintf(paste0("%s is placed outside the figure region and would be clipped away without a warning: ",
                            "y = %.4g, the figure spans %.4g to %.4g"), what, y, min(lo, hi), max(lo, hi)))
    invisible(TRUE)
}

assert_title_fits <- function(txt, cex = CEX$lab) {
    need <- graphics::strwidth(txt, units = "inches", cex = cex)
    have <- graphics::par("fin")[1]
    if (need > have)
        stop(sprintf(paste0("panel title does not fit its panel and mtext would clip it silently: '%s' needs ",
                            "%.2f in at cex %.2f, the figure region is %.2f in. Shorten it or widen the panel."),
                     txt, need, cex, have))
    invisible(TRUE)
}

## One way to draw an interval: a bar in the point's own colour.
ci_bar <- function(x, lo, hi, col, lwd = 0.9, cap = 0) {
    graphics::segments(x, lo, x, hi, col = col, lwd = lwd)
    if (cap > 0) {
        graphics::segments(x - cap, lo, x + cap, lo, col = col, lwd = lwd)
        graphics::segments(x - cap, hi, x + cap, hi, col = col, lwd = lwd)
    }
}

## Or a translucent band of the same hue, for a continuous x.
ci_band <- function(x, lo, hi, col, alpha = 0.18) {
    graphics::polygon(c(x, rev(x)), c(lo, rev(hi)), border = NA,
                      col = grDevices::adjustcolor(col, alpha.f = alpha))
}

## One reference-line style.
ref_line <- function(h = NULL, v = NULL) {
    if (!is.null(h)) graphics::abline(h = h, lty = 2, lwd = 0.7, col = OI[["grey"]])
    if (!is.null(v)) graphics::abline(v = v, lty = 2, lwd = 0.7, col = OI[["grey"]])
}

## A legend, always outside the data or in a corner with no box, labels already passed through lab_of().
legend_plain <- function(pos, labels, col, pch = NULL, lty = NULL, italic = FALSE, ...) {
    graphics::legend(pos, legend = if (italic) italic_species(labels) else labels,
                     col = col, pch = pch, lty = lty, bty = "n", cex = CEX$legend,
                     pt.cex = SYM$point, seg.len = 1.2, ...)
}

## ------------------------------------------------------------------ guard ----
## Refuse to print a sentence inside a panel. Scripts call this on any string they are about to draw with text()
## or mtext() inside the plot region; it stops on anything that looks like prose or carries an ASCII arrow.
## The bottom margin a set of rotated tick labels actually needs, in lines, measured in the open device.
## Readable labels are long ("transcription-factor-independent" is 31 characters where the internal key is 14), so a
## margin guessed in lines can let the longest run out of the panel towards the cell below. `par("cin")[2]` is the
## character height in inches, which is what a margin line is measured in, so the conversion is a measurement rather
## than a constant.
mar_lines_for_rotated <- function(labels, cex, pad_lines = 1.4) {
    if (!length(labels)) return(pad_lines)
    w_in <- max(graphics::strwidth(labels, units = "inches", cex = cex))
    line_in <- graphics::par("cin")[2]
    ceiling(w_in / line_in * 10) / 10 + pad_lines
}

## Angled tick labels: the margin they need, and whether they will collide.
##
## Angled labels overlap when the distance between ticks, projected perpendicular to the text baseline, is smaller than
## a line of text. That projection is `spacing * sin(angle)`, which is why a shallow angle collides where a steep one
## does not. Wrapping each label onto two lines does not help: two stacked lines rotated 32 degrees run into each
## other.
mar_lines_for_angled <- function(labels, srt, cex, pad_lines = 1.2) {
    if (!length(labels)) return(pad_lines)
    w_in <- max(graphics::strwidth(labels, units = "inches", cex = cex))
    drop_in <- w_in * sin(srt * pi / 180)
    ceiling(drop_in / graphics::par("cin")[2] * 10) / 10 + pad_lines
}

assert_angled_labels_fit <- function(at, labels, srt, cex) {
    if (length(at) < 2L) return(invisible(TRUE))
    h <- graphics::strheight("Mg", units = "user", cex = cex)
    perp <- min(diff(sort(at))) * sin(srt * pi / 180)
    if (perp < h)
        stop(sprintf(paste0("angled labels at %g degrees would collide: ticks are %.3f user units apart, which gives ",
                            "%.3f perpendicular, and a line of text is %.3f. Steepen the angle or widen the panel; ",
                            "do NOT wrap them onto two lines, which makes the collision worse."),
                     srt, min(diff(sort(at))), perp, h))
    invisible(TRUE)
}

## Refuse in-panel text that will not fit inside the plot window. Annotating a line at its left end with pos = 3
## centres the string on that x, so it runs off the left edge and the word naming a number can be cut while the
## number survives. Clipping is silent, so it has to be asserted.
assert_text_inside <- function(x, labels, cex = CEX$axis, adj = 0) {
    u <- graphics::par("usr")
    xlo <- if (graphics::par("xlog")) 10^u[1] else u[1]
    xhi <- if (graphics::par("xlog")) 10^u[2] else u[2]
    w <- graphics::strwidth(labels, units = "user", cex = cex)
    left  <- x - adj * w
    right <- left + w
    bad <- which(left < xlo | right > xhi)
    if (length(bad))
        stop(sprintf(paste0("in-panel text would be clipped: %s. Window is [%.4g, %.4g] and the text spans ",
                            "[%.4g, %.4g]. Move it inside, or widen the window -- do not rely on xpd."),
                     paste(sprintf("'%s'", labels[bad]), collapse = ", "), xlo, xhi,
                     min(left[bad]), max(right[bad])))
    invisible(TRUE)
}

## Refuse axis labels that will overlap, measured in the device that is open. Text width is a device property: it
## depends on the font the device resolves and on the size it was opened at, so a width computed anywhere else is a
## guess. This covers the tick labels as well as the data points.
assert_axis_fits <- function(at, labels, cex = CEX$axis, gap_frac = 0.06, side = 1) {
    if (length(at) != length(labels))
        stop(sprintf("assert_axis_fits: %d positions for %d labels", length(at), length(labels)))
    if (length(at) < 2L) return(invisible(TRUE))
    w <- graphics::strwidth(labels, units = "user", cex = cex)
    ord <- order(at)
    at <- at[ord]; w <- w[ord]; labels <- labels[ord]
    span <- diff(range(at))
    need <- (w[-length(w)] + w[-1]) / 2 + span * gap_frac
    have <- diff(at)
    bad <- which(have < need)
    if (length(bad))
        stop(sprintf(paste0("axis labels would overlap on side %d: %s. ",
                            "Widest pair needs %.3f user units and has %.3f. ",
                            "Rotate them, shorten them, or give the panel more width -- do not shrink cex until they ",
                            "merely stop touching."),
                     side,
                     paste(sprintf("'%s' / '%s'", labels[bad], labels[bad + 1]), collapse = ", "),
                     max(need[bad]), min(have[bad])))
    invisible(TRUE)
}

no_prose <- function(s) {
    bad <- grepl("=>|->|\\bNOT\\b", s) | (vapply(strsplit(s, "\\s+"), length, 1L) > 6L)
    if (any(bad))
        stop(sprintf("in-panel prose annotation rejected: %s",
                     paste(utils::head(s[bad], 2), collapse = " | ")))
    invisible(s)
}

invisible(TRUE)
