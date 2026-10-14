## fig_style.R -- the submission figure style for base-R figures: sizes, line widths and colours.
##
## The standard is scripts/figures/base/theme_base.R. Its constants are sourced from that file, not copied,
## so the two cannot diverge; this file only converts them to base-R units.
##
## These conversions were measured, not assumed:
##   * the theme's LWMM values are drawn by ggplot2 at that many points (mm2lw() divides by 2.13, ggplot multiplies by
##     2.134), so 0.55 "mm" is 0.55 pt in the PDFs. In cairo, lwd = 1 is 0.75 pt, so lwd = pt / 0.75.
##   * a ggplot point of size 2.15, stroke 0.40 has an outer diameter of 5.76 pt; base-R pch 21 at pointsize 10 needs
##     cex = 1.125 for the same diameter (measured on the raster), and lwd = 0.76 for the same 0.57 pt outline.
##   * the device pointsize is 10, not 7.5: par("ps") is an integer, so a 7.5 pt device silently becomes 7 pt and every
##     size comes out wrong (7.93 and 9.33 pt measured). Sizes are therefore cex = pt / 10, and margins are counted in
##     7.5 pt lines (mex = 0.75).

local({
  here <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
  sys.source(file.path(here, "scripts/figures/base/theme_base.R"), envir = globalenv())
})
## Programme palette (P1, P2, P3): the six-colour Nature Classic set, one identity per colour within each paper, every
## colour paired with a shape, fill or line type; greys for references, uncertainty and background.
NPG <- list(red = "#E64B35", cyan = "#4DBBD5", teal = "#00A087", blue = "#3C5488", peach = "#F39B7F", slate = "#8491B4")
NPG$orange <- NPG$peach; NPG$sky <- NPG$cyan
CHARCOAL <- "#4D4D4D"
GREYS <- list(g20 = "#333333", g40 = "#666666", g53 = "#888888", g60 = "#999999", g67 = "#AAAAAA", g80 = "#CCCCCC",
              g86 = "#DCDCDC", g90 = "#E5E5E5", g94 = "#F0F0F0", g96 = "#F5F5F5")
TINT <- list(red = "#FAD2CC", cyan = "#D2EEF4", blue = "#CED4E1")   # pale versions for box fills, black text on them
## P1 estimators: naive red, naive/h cyan, instrumented deep blue, h + kappa peach; truth and nulls charcoal
PAL <- list(naive = NPG$red, naive_mid = NPG$red, naive_tint = TINT$red,
            corrected = NPG$blue, corrected_tint = TINT$blue, rel = NPG$cyan, hk = NPG$peach,
            sign_neg = NPG$blue, sign_zero = GREYS$g90, sign_pos = NPG$red, ref = GREYS$g40, fill = GREYS$g90)
COL_SCI$naive <- PAL$naive; COL_SCI$corrected <- PAL$corrected; COL_SCI$uncertainty <- PAL$ref
## P1 true-gain levels, used wherever gamma is a series
GAMMA_COL <- c("1.1" = NPG$red, "1.3" = NPG$cyan, "1.6" = NPG$blue)
## P2 genotypes and gene classes (also used for the class names in P1 S4 and P3)
COL_BIO[["Hsf1"]] <- NPG$red; COL_BIO[["Msn2/4"]] <- NPG$cyan
CLASS_SUB <- c(Hsf1_core = NPG$red, Hsf1_only = NPG$red, Hsf1_joint = NPG$blue, Msn24_only = NPG$cyan,
               both_required = NPG$blue, redundant = NPG$blue, TF_independent = NPG$slate,
               other_induced = NPG$peach, unassigned = GREYS$g60, not_induced = GREYS$g86)
## P3 hybrids: six colours, six shapes, six line types, one key for every panel that shows them
HYB_SUB <- c(ScxSk = NPG$red, ScxSp = NPG$cyan, ScxSu = NPG$teal, SpxSk = NPG$blue, SpxSu = NPG$peach, SkxSu = NPG$slate)
HYB_PCH <- c(ScxSk = 21, ScxSp = 22, ScxSu = 23, SpxSk = 24, SpxSu = 25, SkxSu = 21)
HYB_LTY <- c(ScxSk = 1, ScxSp = 2, ScxSu = 3, SpxSk = 4, SpxSu = 5, SkxSu = 6)
stopifnot(exists("COL_SCI"), exists("COL_BIO"), exists("PT"), exists("MM"), exists("LWMM"))

FONT_SUB <- FONT_FIG                                   # Noto Sans
PS_DEV   <- 10                                         # device pointsize (integer; see above)
CEX_TXT  <- PT$tick / PS_DEV                           # 7.5 pt ticks, legends, annotations
CEX_LAB  <- PT$axis_title / PS_DEV                     # 8.5 pt axis titles
CEX_TAG  <- PT$tag / PS_DEV                            # 10 pt bold panel tags
LW <- list(data  = LWMM$data  / 0.75,                  # 0.55 pt
           ci    = LWMM$ci    / 0.75,                  # 0.40 pt
           ref   = LWMM$ref   / 0.75,                  # 0.35 pt
           axis  = 0.35       / 0.75,                  # 0.35 pt, theme L55-56
           arrow = LWMM$arrow / 0.75)                  # 0.45 pt
PCH_CEX <- 1.125                                       # 5.76 pt outer diameter at pointsize 10, measured on the raster
PCH_LWD <- 0.57 / 0.75                                 # 0.57 pt outline, as measured in the reference PDFs
WIDTH_MM <- c(full = MM$full, onehalf = 108, single = MM$half)

## the species: one colour and one shape each
SPECIES_COL_SUB <- c("S. cerevisiae" = NPG$red, "S. paradoxus" = NPG$cyan,
                     "S. kudriavzevii" = NPG$blue, "S. uvarum" = NPG$teal)
if (exists("SPECIES_COL")) SPECIES_COL[names(SPECIES_COL_SUB)] <- SPECIES_COL_SUB
SPECIES_PCH_SUB <- c("S. cerevisiae" = 21, "S. paradoxus" = 22, "S. kudriavzevii" = 23, "S. uvarum" = 24)

## See PAL above. Genotypes keyed by the count-matrix column names, each also told apart by shape and line type.
GENO_SUB <- list(
  col = c(Wildtype = CHARCOAL, HSF1.KD = NPG$red, MSN24.KO = NPG$cyan, Double.KDKO = NPG$blue),
  bg  = c(Wildtype = CHARCOAL, HSF1.KD = NPG$red, MSN24.KO = NPG$cyan, Double.KDKO = NPG$blue),
  pch = c(Wildtype = 21, HSF1.KD = 22, MSN24.KO = 23, Double.KDKO = 24),
  lty = c(Wildtype = 1, HSF1.KD = 1, MSN24.KO = 1, Double.KDKO = 2))

## base-R graphical parameters equivalent to theme_classic() at the specified sizes
P2_FIG_TAGS <- c("F1_design", "F4_parental", "F5_amplitude_artefact", "F6_cis_trans", "F7_hsf1_domains", "F8_hse", "F9_modules",   # P3 lower-case for G3
                 "F0s_schematic", "F2_epistasis", "F2d_library_share", "F3_go_collapse", "FB4_reference_class_recovery", "F5t_time_course", "FS_P2_torc1_split", "FS_P3_calibration", "F11_calibration")
fig_par <- function() {
  line_pt <- PS_DEV * 1.2                              # one text line at cex 1, in pt (tcl is in these units)
  par(family = FONT_SUB, ps = PS_DEV, cex = 1, cex.axis = CEX_TXT, cex.lab = CEX_LAB, cex.main = CEX_TXT,
      font.main = 1, mex = 0.75, lheight = 1.25, lwd = LW$axis, las = 1, bty = "l", mgp = c(1.9, 0.5, 0),
      tcl = -(1.5 * 72 / 25.4) / line_pt,              # 1.5 mm ticks, theme L57
      xaxs = "r", yaxs = "r", col.axis = COL_SCI$truth, col.lab = COL_SCI$truth, fg = COL_SCI$truth)
}

## Open both devices at a journal width and draw the same function into each. Height is given in mm.
fig_both <- function(tag, width_mm, height_mm, fn, outdir) {
  ## G3 author guidelines require lower-case panel labels; P2 figures only (P1 NAR GAB and P3 keep upper case)
  op_pl <- options(pbl.panel_lower = tag %in% P2_FIG_TAGS); on.exit(options(op_pl), add = TRUE)
  stopifnot(any(abs(width_mm - WIDTH_MM) < 0.5))
  win <- width_mm / 25.4; hin <- height_mm / 25.4
  ok <- TRUE
  ## Plot regions are recorded as each panel is set up. cairo_pdf does not clip to the plot region and bty = "l" draws
  ## no frame, so without this record the geometry check finds no plot region at all and every per-region test passes
  ## by default. The record is written beside the PDF as <tag>.regions.tsv (points, top-left origin).
  REG <- new.env(); REG$rows <- list()
  rec <- function() {
    ## only the figure's own device: a measurement plot on a throwaway device fires the same hook
    if (grDevices::dev.cur() != REG$dev) return(invisible())
    x <- graphics::grconvertX(c(0, 1), "npc", "inches") * 72
    y <- graphics::grconvertY(c(0, 1), "npc", "inches") * 72
    ## bottom-origin points: the checker converts with the PDF's own MediaBox height, which cairo writes rounded
    ## (244 for a 244.8 pt device), so a top-origin record would drift by the rounding
    REG$rows[[length(REG$rows) + 1]] <- c(x[1], y[1], x[2], y[2])
  }
  for (dev in c("pdf", "png")) {
    if (dev == "pdf") grDevices::cairo_pdf(file.path(outdir, paste0(tag, ".pdf")), width = win, height = hin,
                                           pointsize = PS_DEV, family = FONT_SUB)
    else grDevices::png(file.path(outdir, paste0(tag, ".png")), width = round(win * 600), height = round(hin * 600),
                        res = 600, type = "cairo", pointsize = PS_DEV, family = FONT_SUB)
    REG$dev <- grDevices::dev.cur()
    fig_par()
    if (dev == "pdf") setHook("plot.new", rec)
    r <- try(fn(), silent = TRUE)
    if (dev == "pdf") {
      setHook("plot.new", NULL, "replace")
      if (length(REG$rows)) {
        m <- do.call(rbind, REG$rows)
        utils::write.table(data.frame(x0 = m[, 1], yb0 = m[, 2], x1 = m[, 3], yb1 = m[, 4]),
                           file.path(outdir, paste0(tag, ".regions.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
    grDevices::dev.off()
    if (inherits(r, "try-error")) { message(sprintf("  %s (%s) failed: %s", tag, dev, r)); ok <- FALSE }
  }
  ## A failed drawing stops the script. Reporting it and continuing would leave a partial or blank figure on disk
  ## that passes the style and geometry checks, which measure only what is drawn.
  if (!ok) stop(sprintf("fig_both: %s did not draw completely", tag))
  invisible(ok)
}

## after layout(), which resets cex to 0.66 (measured), restore the device base
after_layout <- function() par(cex = 1)

## Bold uppercase panel tag at the top-left of the panel's cell, 1.2 mm in, as tag_grob() in the theme (L75-77).
panel_tag <- function(letter) {
  x <- graphics::grconvertX(1.2 / 25.4, from = "inches", to = "user")
  x <- graphics::grconvertX(graphics::grconvertX(0, "nfc", "inches") + 1.2 / 25.4, "inches", "user")
  y <- graphics::grconvertY(graphics::grconvertY(1, "nfc", "inches") - 1.2 / 25.4, "inches", "user")
  op <- graphics::par(xpd = NA); on.exit(graphics::par(op))
  if (isTRUE(getOption("pbl.panel_lower"))) letter <- tolower(letter)   # G3 wants a, b, c
  graphics::text(x, y, letter, font = 2, cex = CEX_TAG, adj = c(0, 1), col = COL_SCI$truth)
}

## Legend above the plot region, horizontal, no title, no box (theme L58-60). Two rows when one does not fit the
## panel width. Everything else is passed to legend(). Needs a top margin of about 2.4 lines (one row) or 3.4 (two).
legend_top <- function(legend, ..., rows = NA, span = c("plot", "device")) {
  ## span = "plot": the key starts at the data region's left edge and must fit its width.
  ## span = "device": for a panel whose data region is narrow (a long label column on its left), the key starts
  ## 3 mm from the device's left edge and may use the device width less 6 mm.
  ## Widths are compared in inches, so log axes need no special case; legend() takes data coordinates.
  span <- match.arg(span)
  usr <- graphics::par("usr"); pin <- graphics::par("pin")
  a <- list(...)
  nm <- names(a)
  if (!is.null(nm)) a <- a[!nm %in% c("x", "y", "bty", "horiz", "ncol", "cex", "xpd", "inset")]
  dat <- function(u, log) if (log) 10^u else u
  x_in <- if (span == "plot") graphics::grconvertX(usr[1], "user", "inches") else 3 / 25.4
  avail_in <- if (span == "plot") pin[1] else graphics::par("din")[1] - 6 / 25.4
  x0 <- dat(graphics::grconvertX(x_in, "inches", "user"), graphics::par("xlog"))
  y0 <- dat(graphics::grconvertY(graphics::grconvertY(usr[4], "user", "inches") + 0.25 * graphics::par("csi"),
                                 "inches", "user"), graphics::par("ylog"))
  width_in <- function(nr) {
    pr <- do.call(graphics::legend, c(list(x = x0, y = y0, legend = legend, ncol = ceiling(length(legend) / nr),
                                           bty = "n", plot = FALSE, cex = CEX_TXT, x.intersp = 0.6, seg.len = 1.6), a))
    pr$rect$w / diff(usr[1:2]) * pin[1]      # on a log axis rect$w is in log10 units, as usr is (measured)
  }
  nrow <- if (is.na(rows)) 1L else rows
  while (nrow < length(legend) && width_in(nrow) > avail_in) nrow <- nrow + 1L
  if (width_in(nrow) > avail_in) stop(sprintf("legend_top: '%s' does not fit even as one column", legend[1]))
  invisible(do.call(graphics::legend, c(list(x = x0, y = y0, legend = legend, ncol = ceiling(length(legend) / nrow),
                                             bty = "n", xpd = NA, cex = CEX_TXT, yjust = 0, xjust = 0,
                                             x.intersp = 0.6, seg.len = 1.6), a)))
}

## Negative numbers are set with the minus sign U+2212, not the hyphen-minus base R writes. The complete Noto Sans
## (OFL, notofonts; the full Noto Sans, not the system variable font) carries U+2212; the system variable font does not.
## minus() converts labels; axis() below is used for every numeric axis.
MINUS <- "\u2212"
minus <- function(x) sub("^-", MINUS, x)
num_lab <- function(v, digits = NULL) {
  s <- if (is.null(digits)) format(v, trim = TRUE, drop0trailing = FALSE) else formatC(v, format = "f", digits = digits)
  minus(s)
}
## axis() with U+2212 on negative labels. When labels are not given, they are formatted as graphics::axis() would
## (axTicks) and then converted. Character labels are converted too.
axis <- function(side, at = NULL, labels = TRUE, ...) {
  if (isTRUE(labels) || is.null(labels)) {
    if (is.null(at)) at <- graphics::axTicks(side)
    labels <- minus(format(at, trim = TRUE))
  } else if (is.character(labels)) labels <- minus(labels)
  graphics::axis(side, at = at, labels = labels, ...)
}

## plot() draws its own axes with graphics' internal labels, which use the hyphen-minus. This wrapper suppresses them
## and draws the same axes through axis() above. Only the default axes change; everything else is passed through.
plot <- function(x, y = NULL, ..., axes = TRUE, xaxt = "s", yaxt = "s") {
  graphics::plot(x, y, ..., axes = FALSE)
  if (axes) {
    if (xaxt != "n") axis(1)
    if (yaxt != "n") axis(2)
    bt <- graphics::par("bty")
    if (bt != "n" && bt != "l") graphics::box(bty = bt)
    if (bt == "l") graphics::box(bty = "l")
  }
  invisible(NULL)
}

## In-panel annotation at 7.5 pt
annot <- function(x, y, labels, ...) graphics::text(x, y, labels, cex = CEX_TXT, ...)
