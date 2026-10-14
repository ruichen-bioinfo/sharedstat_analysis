## theme_base.R -- publication visual grammar (sizes, line widths, theme) shared by the figure scripts.
##
## It sets the palette, theme and line semantics. Default ggplot2 settings give an engineering look: bright blue
## and orange as the primary semantic pair, dashed lines carrying ordinary categories, and annotation text enlarged
## to fill whitespace. Using ggplot2 does not by itself produce a
## journal figure; the visual language has to be specified separately, which is what this file does.
##
## Three rules:
##   1. Data series are SOLID. Categories are distinguished by colour, shape and fill, never by dash pattern.
##   2. At most one reference line type per panel: black dashed for a null, grey dotted for a threshold.
##   3. Type is sized for the printed page. Nothing is enlarged to fill space, and nothing falls below 7 pt.

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(scales)
})

FONT_FIG <- "Noto Sans"

## ---- palette ------------------------------------------------------------------------------------------------
## Black carries truth, model and reference. Dark red carries the naive quantity and the problem. Dark green
## carries the corrected quantity and anything that survives a check. Grey carries uncertainty. Red and green are
## always accompanied by a shape or fill difference so the distinction survives deuteranopia and greyscale.
COL_SCI <- list(
  truth       = "#111111",
  naive       = "#D01C1F",
  corrected   = "#33A02C",
  uncertainty = "#8C8C8C",
  neutral     = "#4D4D4D",
  light       = "#F2F2F2"
)
## Differential-expression methods, each with its own shape.
COL_METHOD <- c(DESeq2 = "#111111", edgeR = "#D01C1F", limma = "#33A02C")
PCH_METHOD <- c(DESeq2 = 21, edgeR = 24, limma = 22)
## Biological pages.
COL_BIO <- c(Hsf1 = "#33A02C", `Msn2/4` = "#D01C1F", neutral = "#111111", double = "#111111")

## ---- type, in points at final printed size --------------------------------------------------------------------
PT <- list(tag = 10, axis_title = 8.5, strip = 8.5, tick = 7.5, legend = 7.5, annot = 7.5)
stopifnot(min(unlist(PT)) >= 7)

## ---- physical sizes and weights -------------------------------------------------------------------------------
MM <- list(full = 174, half = 84, gutter = 6)
## Line widths in mm; ggplot2 linewidth is in mm/2.13 units, so convert at the point of use.
LWMM <- list(data = 0.55, ci = 0.40, ref = 0.35, arrow = 0.45)
mm2lw <- function(mm) mm / 2.13
PSMM <- list(point = 2.15, stroke = 0.40)

theme_base <- function() {
  theme_classic(base_family = FONT_FIG, base_size = PT$tick) +
    theme(
      axis.title  = element_text(size = PT$axis_title, colour = COL_SCI$truth),
      axis.text   = element_text(size = PT$tick, colour = COL_SCI$truth),
      axis.line   = element_line(linewidth = mm2lw(0.35), colour = COL_SCI$truth),
      axis.ticks  = element_line(linewidth = mm2lw(0.35), colour = COL_SCI$truth),
      axis.ticks.length = unit(1.5, "mm"),
      legend.position   = "top",
      legend.direction  = "horizontal",
      legend.title      = element_blank(),
      legend.text       = element_text(size = PT$legend, colour = COL_SCI$truth),
      legend.key.width  = unit(4, "mm"),
      legend.key.height = unit(3.5, "mm"),
      legend.margin     = margin(0, 0, 1, 0),
      legend.box.spacing = unit(1, "mm"),
      strip.background  = element_blank(),
      strip.text        = element_text(size = PT$strip, face = "bold", colour = COL_SCI$truth),
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin   = margin(3, 4, 3, 3, unit = "mm")
    )
}

## Panel tag, outside the data region.
tag_grob <- function(letter) {
  textGrob(letter, x = unit(1.2, "mm"), y = unit(1, "npc") - unit(1.2, "mm"), hjust = 0, vjust = 1,
           gp = gpar(fontfamily = FONT_FIG, fontface = "bold", fontsize = PT$tag, col = COL_SCI$truth))
}

fmt_2dp <- function(x) formatC(x, format = "f", digits = 2)
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_pct <- function(x) paste0(formatC(x, format = "d"), "%")

## ---- export ---------------------------------------------------------------------------------------------------
## ---- side-by-side composition -----------------------------------------------------------------------------------
## Printing a ggplot into a hand-pushed viewport with print(newpage = FALSE) does NOT confine it: legend key glyphs
## escape and are drawn elsewhere on the page. Escaped glyphs can appear as stray red points and green squares inside the
## left-hand panel of a two-panel figure or as a stray mark inside another tile. Converting each plot to a
## gtable and laying the gtables out keeps every glyph inside its own cell. Panel tags are plain text grobs, which
## do not have this problem, so they are still drawn in viewports afterwards.
compose_row <- function(plots, widths, tags = LETTERS[seq_along(plots)]) {
  stopifnot(length(plots) == length(widths), length(tags) == length(plots), all(widths > 0))
  ## A negative legend margin makes ggplot2 draw the legend key glyphs outside the legend box, where they land inside
  ## the panel and read as data points. It is refused: use
  ## legend.box.spacing to tighten the gap instead.
  for (i in seq_along(plots)) {
    if (!inherits(plots[[i]], "ggplot")) next
    lm <- plots[[i]]$theme$legend.margin
    if (!is.null(lm) && any(as.numeric(lm) < 0))
      stop(sprintf("INVALID: plot %d has a negative legend.margin (%s); key glyphs would escape into the panel. ",
                   i, paste(as.numeric(lm), collapse = ", ")),
           "Use legend.box.spacing to close the gap.")
    ## A top or bottom legend justified left or right also puts key glyphs outside the legend box, where they land in
    ## the neighbouring cell and read as data. Centre justification
    ## is the only setting verified to keep the keys inside.
    lj <- plots[[i]]$theme$legend.justification
    if (!is.null(lj) && !identical(lj, "center") && !identical(lj, "centre"))
      stop(sprintf("INVALID: plot %d sets legend.justification = '%s'; ", i, paste(lj, collapse = ", ")),
           "key glyphs escape the legend box and are drawn over the adjacent panel. Leave it centred.")
  }
  ## Accepts a ggplot or an already-assembled grob, so a column that pairs a plot with a text strip can be passed in.
  as_cell <- function(x) if (inherits(x, "ggplot")) ggplot2::ggplotGrob(x) else x
  ## Each cell is drawn inside a CLIPPED viewport. ggplot2 can place a legend key glyph outside the legend box -- the
  ## trigger varies with justification and with which layers contribute keys
  ## In a laid-out figure, such a glyph can land in the neighbouring panel and read as a data point.
  ## Clipping at the cell boundary removes the whole class of problem instead of one instance of it.
  fr <- widths / sum(widths)
  x0 <- c(0, cumsum(fr)[-length(fr)])
  grid::grid.newpage()
  for (i in seq_along(plots)) {
    grid::pushViewport(grid::viewport(x = x0[i], width = fr[i], just = "left", clip = "on"))
    grid::grid.draw(as_cell(plots[[i]]))
    grid::popViewport()
  }
  for (i in seq_along(tags)) {
    grid::pushViewport(grid::viewport(x = x0[i], width = fr[i], just = "left"))
    grid::grid.draw(tag_grob(tags[i]))
    grid::popViewport()
  }
  invisible(NULL)
}

## ---- text measurement, on the device the figure is actually written with -----------------------------------------
## Measuring on the default device silently falls back to a PostScript metric for a font it does not have, so the
## number returned would describe the wrong font. Every caller gets the cairo metric.
text_mm <- function(txt, pt) {
  tmp <- tempfile(fileext = ".pdf")
  cairo_pdf(tmp, width = 8, height = 8); on.exit({ dev.off(); unlink(tmp) }, add = TRUE)
  g <- grid::textGrob(txt, gp = grid::gpar(fontfamily = "Noto Sans", fontsize = pt))
  grid::convertWidth(grid::grobWidth(g), "mm", valueOnly = TRUE)
}
## Refuses rather than warns: a label that does not fit is an overlap in the printed figure.
assert_text_fits <- function(txt, panel_mm, pt, what) {
  w <- text_mm(txt, pt)
  cat(sprintf("  text fit  %-24s %5.1f mm needed, %5.1f mm available\n", what, w, panel_mm))
  if (w > panel_mm) stop(sprintf("INVALID: %s needs %.1f mm but only %.1f mm is available", what, w, panel_mm))
  invisible(w)
}

export_fig <- function(tag, plot_or_fn, width_mm, height_mm, plot_data = NULL, caption = NULL) {
  outdir  <- file.path(Sys.getenv("PBL_ROOT", "<PBL_ROOT>"), "doc/04_figures/extra")
  datadir <- file.path(Sys.getenv("PBL_ROOT", "<PBL_ROOT>"), "data/02_processed/extra")
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
  win <- width_mm / 25.4; hin <- height_mm / 25.4
  draw <- function() if (is.function(plot_or_fn)) plot_or_fn() else print(plot_or_fn)

  grDevices::cairo_pdf(file.path(outdir, paste0(tag, ".pdf")), width = win, height = hin,
                       family = FONT_FIG, onefile = TRUE)
  draw(); grDevices::dev.off()
  grDevices::png(file.path(outdir, paste0(tag, ".png")), width = win, height = hin, units = "in",
                 res = 600, type = "cairo", family = FONT_FIG)
  draw(); grDevices::dev.off()
  ## A render at 75% of final size, where crowding shows first.
  grDevices::png(file.path(outdir, paste0(tag, "_75pct.png")), width = win * 0.75, height = hin * 0.75,
                 units = "in", res = 600, type = "cairo", family = FONT_FIG)
  draw(); grDevices::dev.off()

  if (!is.null(plot_data)) {
    stopifnot(is.data.frame(plot_data), nrow(plot_data) > 0)
    data.table::fwrite(plot_data, file.path(datadir, paste0(tag, "_plotdata.tsv")), sep = "\t")
  }
  if (!is.null(caption)) {
    stopifnot(nzchar(caption))
    ## A caption that still contains a doubled percent sign is a formatting bug, not a caption.
    if (grepl("%%", caption, fixed = TRUE))
      stop(sprintf("caption for %s contains '%%%%': fix the sprintf escaping", tag))
    if (grepl("\\band and\\b", caption)) stop(sprintf("caption for %s contains 'and and'", tag))
    writeLines(caption, file.path(outdir, paste0("CAPTION_", tag, ".md")))
  }
  message(sprintf("  %-22s %5.1f x %5.1f mm  pdf+png+75pct%s%s", tag, width_mm, height_mm,
                  if (!is.null(plot_data)) "+data" else "", if (!is.null(caption)) "+caption" else ""))
  invisible(TRUE)
}

message(sprintf("theme_base loaded: %s, classic, tags %.0f pt, axis %.1f pt, ticks %.1f pt",
                FONT_FIG, PT$tag, PT$axis_title, PT$tick))
