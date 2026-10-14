#!/usr/bin/env Rscript
## 60_p2_fig3_two_backgrounds.R -- P2 Figure 3, against both backgrounds, in the house style.
##
## Why both backgrounds. Against other induced genes, protein folding and chaperone binding survive the paralogue
## collapse and cell wall does not; against the whole genome, cell wall survives and protein folding does not.
## Both are true within their respective comparisons. The
## figure shows both columns side by side, with collapse survival marked per term, so each
## comparison can be interpreted separately:
##
##   cell wall enrichment holds against the whole genome and not against other induced genes, so it is a property
##   of being induced and not a property of the Hsf1-dependent core;
##   protein folding and chaperone binding stand out specifically against other induced genes.
##
## The figure uses reader-facing labels for comparison sets.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
source(file.path(root, "scripts/figures/theme.R"))
hsd <- file.path(root, "data/02_processed/hsf1")
fig <- file.path(root, "doc/04_figures/final")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
ALPHA <- 0.05

G <- fread(file.path(hsd, "hsf1_core_go_enrichment_37_30.tsv"))
SETS <- c(universe = "core_vs_universe", induced = "core_vs_other_induced")
COLL <- c(universe = "core_vs_universe_paralog_collapsed", induced = "core_vs_other_induced_paralog_collapsed")
stopifnot(all(c(SETS, COLL) %in% G$set))

## Terms to show: those significant in either background before collapse, so the figure does not select the
## comparison that suits the text. Ranked by the better of the two odds ratios.
sig <- G[set %in% SETS & padj < ALPHA, unique(go)]
say("=== terms significant against at least one background before collapse: %d ===", length(sig))
D <- G[go %in% sig & set %in% c(SETS, COLL)]
D[, bg := ifelse(set %in% c(SETS[["universe"]], COLL[["universe"]]), "universe", "induced")]
D[, collapsed := set %in% COLL]
W <- dcast(D, go + name ~ bg + collapsed, value.var = c("odds", "padj", "k", "m"))
setnames(W, gsub("_TRUE$", "_coll", gsub("_FALSE$", "_full", names(W))))
need <- c("odds_universe_full", "odds_induced_full", "padj_universe_full", "padj_induced_full",
          "padj_universe_coll", "padj_induced_coll")
stopifnot(all(need %in% names(W)))
W <- W[is.finite(odds_universe_full) & is.finite(odds_induced_full)]
W[, best := pmax(odds_universe_full, odds_induced_full)]
setorder(W, -best)
W <- W[1:min(12L, nrow(W))]
## Shorten by choosing a shorter name, not by cutting the string: "protein-folding chaperone bi..." is not a
## label. Anything not in this map keeps its full GO name and the width is set to fit it.
SHORT <- c("protein-folding chaperone binding" = "chaperone binding",
           "structural constituent of cell wall" = "cell wall structural constituent",
           "external encapsulating structure" = "external encapsulating structure",
           "fungal-type cell wall" = "cell wall, fungal-type")
W[, short := sub("^'de novo' ", "", name)]
W[, short := ifelse(short %in% names(SHORT), unname(SHORT[short]), short)]
stopifnot(max(nchar(W$short)) <= 34)
say("  showing the %d terms with the largest odds ratio in either background", nrow(W))
print(W[, .(term = short, or_genome = round(odds_universe_full, 1), p_genome = signif(padj_universe_full, 2),
            genome_survives = padj_universe_coll < ALPHA,
            or_induced = round(odds_induced_full, 1), p_induced = signif(padj_induced_full, 2),
            induced_survives = padj_induced_coll < ALPHA)], row.names = FALSE)

## ------------------------------------------------------------------ draw ----
## Drawn in the shared figure style (scripts/figures/fig_style.R). No panel titles: A is against the whole
## genome and B against the other induced genes, as the caption says. Colour follows the theme's semantics: dark green
## for what survives a check (significant and surviving the paralogue collapse), dark red for what does not, grey for not
## significant; the three also differ by shape (filled circle, open circle, cross).
source(file.path(root, "scripts/figures/fig_style.R"))
ST_COL <- c(survives = COL_SCI$corrected, fails = COL_SCI$naive, ns = COL_SCI$uncertainty)
ST_PCH <- c(survives = 16, fails = 21, ns = 4)
fig_both("F3_go_collapse", 174, 72, function() {
  layout(matrix(1:2, 1, 2), widths = c(1.55, 1)); after_layout()
  LAB_MAR <- mar_lines_for_rotated(W$short, cex = CEX_TXT, pad_lines = 1.0) / 0.75
  xr <- range(c(W$odds_universe_full, W$odds_induced_full), na.rm = TRUE)
  xr <- c(0.8, xr[2] * 1.18)
  for (p in list(list(k = "universe", tag = "A"), list(k = "induced", tag = "B"))) {
    par(mar = c(3.4, if (p$k == "universe") LAB_MAR else 1.2, 4.6, 0.8))   # room for the legend above the panel tags
    o <- W[[paste0("odds_", p$k, "_full")]]
    pf <- W[[paste0("padj_", p$k, "_full")]]
    pc <- W[[paste0("padj_", p$k, "_coll")]]
    y <- rev(seq_len(nrow(W)))
    ## each panel names its background in the axis title, so the figure reads without the caption
    plot(NA, xlim = xr, ylim = c(0.4, nrow(W) + 0.6), log = "x", axes = FALSE, ylab = "",
         xlab = if (p$k == "universe") "odds ratio vs whole genome" else "odds ratio vs other induced genes")
    abline(v = 1, col = COL_SCI$truth, lty = 2, lwd = LW$ref)
    ## one marker per term, chosen from three mutually exclusive states
    state <- ifelse(pf >= ALPHA, "ns", ifelse(pc < ALPHA, "survives", "fails"))
    segments(1, y, o, y, col = COL_SCI$uncertainty, lwd = LW$ci)
    points(o, y, pch = ST_PCH[state], col = ST_COL[state], bg = "white", cex = PCH_CEX, lwd = PCH_LWD)
    axis(1, at = c(1, 2, 5, 10, 20), labels = c("1", "2", "5", "10", "20"))
    if (p$k == "universe") axis(2, at = y, labels = W$short, tick = FALSE, las = 1)
    box(bty = "l", lwd = LW$axis)
    if (p$k == "universe")
      legend_top(c("significant, survives the paralog collapse", "significant, does not survive", "not significant"),
                 pch = ST_PCH, col = ST_COL, pt.cex = PCH_CEX, pt.lwd = PCH_LWD, span = "device", rows = 2)   # two rows: one row reaches the right edge
    panel_tag(p$tag)
  }
}, fig)

fwrite(W[, .(go, name, odds_genome = odds_universe_full, padj_genome = padj_universe_full,
             padj_genome_collapsed = padj_universe_coll,
             odds_induced = odds_induced_full, padj_induced = padj_induced_full,
             padj_induced_collapsed = padj_induced_coll)],
       file.path(hsd, "p2_fig3_two_backgrounds.tsv"), sep = "\t")

say("\n=== interpretation of the two backgrounds ===")
cw <- W[grepl("^cell wall$|^fungal-type cell wall$", name)]
pfd <- W[name == "protein folding"]
cb <- W[grepl("chaperone binding", name)]
for (r in list(list("cell wall", cw), list("protein folding", pfd), list("chaperone binding", cb))) {
  d <- r[[2]]
  if (!nrow(d)) { say("  %s not among the plotted terms", r[[1]]); next }
  d <- d[1]
  say("  %-20s genome %.1f-fold p %.3g, collapse %s | induced %.1f-fold p %.3g, collapse %s",
      r[[1]], d$odds_universe_full, d$padj_universe_full,
      ifelse(d$padj_universe_coll < ALPHA, "survives", "fails"),
      d$odds_induced_full, d$padj_induced_full,
      ifelse(d$padj_induced_coll < ALPHA, "survives", "fails"))
}

## ------------------------------------------------------- the caption for this figure ---
## The caption is written to a separate file that 19_figures.R reads into the
## caption collection.
capf3 <- file.path(fig, "CAPTION_F3.md")
## The row is taken by name and the name must match exactly one row. `[1]` on a multi-row match would silently
## caption the figure with another term's numbers.
CW <- W[name == "fungal-type cell wall"]
stopifnot(nrow(CW) == 1L)
writeLines(sprintf(paste0(
  "F3. Gene Ontology enrichment of the Hsf1-associated core at 37 \u00b0C, 30 min against two backgrounds, the whole ",
  "genome (a) and ",
  "the other induced genes (b), with the verified-family paralog collapse applied to each. Point shape and color give ",
  "the three outcomes: enriched and surviving the collapse (filled, blue), enriched but not surviving (open, red), ",
  "and not enriched (cross, gray). Which ",
  "enrichment survives depends on the background: cell-wall categories are enriched %.1f-fold against the genome ",
  "(BH-adjusted p %.2e) and survive the collapse (%.2e), but against the other induced genes they are enriched ",
  "%.1f-fold (%.2e) and do not survive (%.3f); protein folding behaves the other way round. Chaperone binding ",
  "survives against both. Source: %s."),
  CW$odds_universe_full,
  CW$padj_universe_full,
  CW$padj_universe_coll,
  CW$odds_induced_full,
  CW$padj_induced_full,
  CW$padj_induced_coll,
  "hsf1_core_go_enrichment_37_30.tsv"), capf3)
say("caption written: %s", capf3)
say("")
