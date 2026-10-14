## Figure width 6.9 in: a target journal caps double-column figures at 178 mm = 7.008 in. Set at
## source rather than rescaled at submission, because rescaling a raster forfeits the >= 300 dpi
## the same requirement imposes.
# 22_fig_f2d.R -- P2 figure F2d: library share, its non-additivity, and the offset.
#
# Purpose. P2 states that the non-zero genome-wide median of the epistasis contrast is
# compositional. The claim rests on three quantities, drawn here:
#   (a) the share of each library taken by the responding gene classes, per genotype and condition
#   (b) the same (+1,-1,-1,+1) contrast applied to that share, i.e. whether the dilution itself is
#       non-additive, because an additive dilution cancels exactly in the contrast
#   (c) the genome-wide median interaction before and after removing the responding classes from the
#       estimation of the normalisation factors
# If the claim is right, (b) and (c) must track each other across conditions. That is a falsifiable
# prediction and the panel is where a reader checks it.
#
# Class assignment is the exclusive one. The class file holds one assignment per gene per DE method
# and the three methods disagree for 312 genes; the union double-counts them. This script uses the
# method-unanimous, mutually exclusive assignment (script 21) so that the "responding classes" set
# is the same set the P2 class table is built on.
#
# Usage: Rscript scripts/hsf1/22_fig_f2d.R

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
pro <- file.path(root, "data/02_processed"); hsd <- file.path(pro, "hsf1")
meta <- file.path(root, "data/01_metadata"); repd <- file.path(root, "doc/05_reports")
fig <- file.path(root, "doc/04_figures/final"); MIRROR <- "logs/output_copy"
DPI <- 300
## Colours come from the shared theme, so a colour means the same thing in every figure of every
## manuscript. A local copy of the Okabe-Ito vector would drift from the theme unnoticed.
source(file.path(root, "scripts/figures/theme.R"))

## Device base point size 11. With the annotation cex values used here, an 8 pt base puts text at
## 4 to 6 pt at the final figure width, which is not legible in print. Raising the base moves every
## size by one factor and preserves the relative emphasis; the physical width does not change, so
## this is not a figure shrunk to make its type look bigger.
PT_BASE <- 11
GENO <- c("Wildtype","HSF1.KD","MSN24.KO","Double.KDKO")
GCOL <- c(Wildtype=OI[["blue"]], HSF1.KD=OI[["vermillion"]],
          MSN24.KO=OI[["green"]], Double.KDKO=OI[["orange"]])
GPCH <- c(Wildtype=21, HSF1.KD=22, MSN24.KO=23, Double.KDKO=24)
CONDS <- c("25_0","37_10","37_30","42_10","42_30")
HEAT  <- c("37_10","37_30","42_10","42_30")
say <- function(...) cat(sprintf(...), "\n", sep = "")
logf <- file.path(repd, "logs/hsf1_fig_f2d.log"); sink(logf, split = TRUE)
say("=== 22_fig_f2d | %s ===", format(Sys.time()))

## ------------------------------------------------------- the share table ------
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
stopifnot(identical(SM$sample, libs))
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]

CLS <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
CL <- unique(CLS[, .(gene, cls = class, method)])
U <- CL[, .(k = uniqueN(cls), nm = uniqueN(method), cls = cls[1]), by = gene]
U[, final := ifelse(k == 1L & nm == uniqueN(CL$method), cls, "unassigned")]
stopifnot(!any(duplicated(U$gene)))                 # one class per gene, checked
RESP <- U[final %in% c("Hsf1_only","Msn24_only","redundant","both_required"), gene]
say("responding classes (method-unanimous, exclusive): %d genes", length(RESP))

k <- rownames(M) %in% RESP
share <- colSums(M[k, , drop = FALSE]) / colSums(M)
SH <- data.table(genotype = SM$knockout, cond = SM$cond, share = share)[
      , .(share = mean(share), sd = sd(share)), by = .(cond, genotype)]
W <- dcast(SH, cond ~ genotype, value.var = "share"); setkey(W, cond); W <- W[CONDS]
W[, nonadd := Double.KDKO - HSF1.KD - MSN24.KO + Wildtype]
fwrite(merge(SH, W[, .(cond, nonadd)], by = "cond"),
       file.path(hsd, "p2_library_share.tsv"), sep = "\t")
say("\n%-8s %10s %10s %10s %12s %14s", "cond", "Wildtype", "HSF1.KD", "MSN24.KO", "Double.KDKO",
    "non-additivity")
for (i in seq_len(nrow(W))) say("%-8s %10.4f %10.4f %10.4f %12.4f %14.4f", W$cond[i],
  W$Wildtype[i], W$HSF1.KD[i], W$MSN24.KO[i], W$Double.KDKO[i], W$nonadd[i])

## --------------------------------------------- the offsets, both regimes ------
OFF <- fread(file.path(hsd, "p2_offset_origin_summary.tsv"))
OFF[, regime_short := ifelse(regime %like% "excluded", "corrected", "published")]
O <- dcast(OFF[condition %in% HEAT], condition ~ regime_short, value.var = "med_I")
setkey(O, condition); O <- O[HEAT]
O <- merge(O, W[, .(condition = cond, nonadd)], by = "condition"); setkey(O, condition); O <- O[HEAT]
say("\n%-8s %14s %14s %16s", "cond", "offset published", "offset corrected", "share non-additivity")
for (i in seq_len(nrow(O))) say("%-8s %+14.4f %+14.4f %+16.4f", O$condition[i],
  O$published[i], O$corrected[i], O$nonadd[i])
r_pub <- cor(O$nonadd, O$published); r_cor <- cor(O$nonadd, O$corrected)
say("\ncorrelation across the four heat conditions:")
say("  r(share non-additivity, offset PUBLISHED)  = %+.3f   <- the prediction", r_pub)
say("  r(share non-additivity, offset CORRECTED)  = %+.3f   <- should be weaker once the", r_cor)
say("     compositional term has been removed from the normalisation")
say("n = 4 conditions => DESCRIPTIVE ONLY, no significance is claimed for a 4-point correlation.")
fwrite(O, file.path(hsd, "p2_share_vs_offset.tsv"), sep = "\t")

## ----------------------------------------------------------------- the figure --
## Drawn in the shared figure style (scripts/figures/fig_style.R). No panel titles (panels are named in the caption).
## The published normalisation is filled black, responding classes excluded open black, the share contrast light
## grey bars.
source(file.path(root, "scripts/figures/fig_style.R"))
## Genotypes in black and greys (GENO_SUB in fig_style.R), each with its own shape; red and green are reserved in P2
## for naive against corrected. The end labels are black, so the grey series stay readable; each is joined to its
## line by a leader.
GCOL <- GENO_SUB$col; GBG <- GENO_SUB$bg; GPCH <- GENO_SUB$pch; GLTY <- GENO_SUB$lty
GLTY[] <- 1   # the style spec does not use dash pattern as a distinguisher; shape and fill carry it
draw <- function() {
  LW3 <- c(1.6, 0.9, 1.25)
  layout(matrix(1:3, 1, 3), widths = LW3); after_layout()   # B and C are separate panels
  cl <- lab_of(CONDS)
  par(mar = c(mar_lines_for_rotated(cl, CEX_TXT) / 0.75, 3.8, 5.0, 0.8))
  x <- seq_along(CONDS)
  ## The right end of the x range is solved from the measured width of the longest end label, so the labels
  ## finish 2 mm inside the plot region rather than at its edge.
  GLAB2 <- c(Wildtype = "wild type", HSF1.KD = "Hsf1-depleted", MSN24.KO = "msn2\u0394 msn4\u0394", Double.KDKO = "double perturbation")   # local labels
  lab_in <- max(strwidth(GLAB2[GENO], units = "inches", cex = CEX_TXT)) + 2 / 25.4
  ## before plot.new() par("pin") describes another panel; take panel A's width from the layout instead
  pinA <- par("din")[1] * LW3[1] / sum(LW3) - sum(par("mai")[c(2, 4)])
  a_end <- length(CONDS) + 0.32; frac <- lab_in / pinA
  stopifnot(frac < 0.6)
  xr <- 0.6 + (a_end - 0.6) / (1 - frac)
  plot(NA, xlim = c(0.6, xr), ylim = c(0, max(W[, ..GENO]) * 1.15), xaxs = "i",
       xaxt = "n", xlab = "", ylab = "library share taken by responding genes")
  axis(1, x, cl, las = 2)
  abline(v = which(CONDS == "37_30"), lty = 2, col = COL_SCI$truth, lwd = LW$ref)
  for (g in GENO) {
    lines(x, W[[g]], col = GCOL[[g]], lty = GLTY[[g]], lwd = LW$data)
    points(x, W[[g]], pch = GPCH[[g]], bg = GBG[[g]], col = COL_SCI$truth, cex = PCH_CEX, lwd = PCH_LWD)
  }
  glab <- GLAB2[GENO]
  gy <- vapply(GENO, function(g) W[[g]][length(x)], numeric(1))
  ## end labels pushed apart by a measured line height where two series end close together (wild type and
  ## Hsf1-depleted end 0.02 apart); a short leader joins each label to its line
  lh <- 1.5 * strheight("Mg", units = "user", cex = CEX_TXT)
  o <- order(-gy); ly <- gy[o]
  for (k in seq_along(ly)[-1]) if (ly[k - 1] - ly[k] < lh) ly[k] <- ly[k - 1] - lh
  gyl <- gy; gyl[o] <- ly
  segments(length(x) + 0.05, gy, length(x) + 0.28, gyl, col = COL_SCI$uncertainty, lwd = LW$ref)
  assert_text_inside(rep(length(x) + 0.32, length(GENO)), glab, cex = CEX_TXT, adj = 0)
  for (i in seq_along(GENO))
    text(length(x) + 0.32, gyl[i], glab[i], adj = c(0, 0.5), cex = CEX_TXT, col = COL_SCI$truth)
  panel_tag("A")
  ocl <- lab_of(O$condition)
  ## (B) the share contrast alone, on its own axis
  par(mar = c(mar_lines_for_rotated(ocl, CEX_TXT) / 0.75, 5.0, 5.0, 0.8), mgp = c(3.1, 0.5, 0))   # axis title clear of ticks
  yb <- range(c(O$nonadd, 0)) * 1.25
  bp <- barplot(O$nonadd, names.arg = rep("", length(ocl)), border = COL_SCI$neutral, lwd = LW$ci,
                col = ifelse(O$condition == "37_30", NPG$red, GREYS$g67), ylim = yb, axes = FALSE, ylab = "interaction contrast of library share")
  axis(2); axis(1, at = bp, labels = ocl, las = 2, tick = FALSE)
  abline(h = 0, col = COL_SCI$truth, lwd = LW$axis)
  panel_tag("B")
  ## (C) the genome-wide median interaction I under the two normalisations, on its own axis
  par(mar = c(mar_lines_for_rotated(ocl, CEX_TXT) / 0.75, 5.0, 5.0, 0.8), mgp = c(3.1, 0.5, 0))   # axis title clear of ticks
  xs <- seq_along(ocl); yc <- range(c(O$published, O$corrected, 0)); yc <- yc + c(-1, 1) * 0.15 * diff(yc)
  plot(NA, xlim = c(0.6, length(xs) + 0.4), ylim = yc, xaxt = "n", xlab = "", ylab = "genome-wide median interaction I")
  axis(1, at = xs, labels = ocl, las = 2, tick = FALSE)
  abline(h = 0, col = COL_SCI$truth, lwd = LW$axis)
  lines(xs, O$published, col = COL_SCI$naive, lwd = LW$data)
  points(xs, O$published, pch = 21, bg = COL_SCI$naive, cex = PCH_CEX, lwd = PCH_LWD)
  lines(xs, O$corrected, col = COL_SCI$corrected, lwd = LW$data)
  points(xs, O$corrected, pch = 24, bg = COL_SCI$corrected, cex = PCH_CEX, lwd = PCH_LWD)
  legend_top(c("published", "classes excluded"), pch = c(21, 24),
             pt.bg = c(COL_SCI$naive, COL_SCI$corrected), col = COL_SCI$truth, pt.cex = PCH_CEX, pt.lwd = PCH_LWD, rows = 2)
  panel_tag("C")
}
fig_both("F2d_library_share", 174, 112, draw, fig)
say("\nF2d written (pdf + png)")

CAP <- sprintf("F2d. The epistasis contrast carries a compositional offset. (a) Share of each library taken by the responding gene classes (method-unanimous, mutually exclusive assignment; %d genes), by genotype and condition. The share differs between genotypes, but an additive share difference cancels exactly in the (+1,-1,-1,+1) interaction contrast, so a share difference alone cannot produce an offset. (b) The same contrast applied to the share, that is, whether the dilution is itself non-additive. (c) The genome-wide median interaction I under the published normalization (red circles) and after excluding the responding classes from the estimation of the normalization factors (blue triangles). At 37 \u00b0C, 30 min the share non-additivity is %+.4f and the offset falls from %+.4f to %+.4f, that is, %.0f%% of it is removed; at 42 \u00b0C, 10 min and 42 \u00b0C, 30 min the share non-additivity is far smaller (%+.4f, %+.4f) and the offsets were already near zero. Across the four heat conditions r(share non-additivity, offset) = %+.3f under the published normalization and %+.3f after correction. With n = 4 conditions this correlation is descriptive only and no significance is claimed. Sources: p2_library_share.tsv, p2_share_vs_offset.tsv, p2_offset_origin_summary.tsv.",
  length(RESP), O$nonadd[O$condition=="37_30"], O$published[O$condition=="37_30"],
  O$corrected[O$condition=="37_30"],
  100*(1 - abs(O$corrected[O$condition=="37_30"])/abs(O$published[O$condition=="37_30"])),
  O$nonadd[O$condition=="42_10"], O$nonadd[O$condition=="42_30"], r_pub, r_cor)
## Each figure script writes its caption to a separate file. This avoids multiple scripts
## appending to one file in an unspecified order, where the last writer could erase earlier content, so
## caption outputs remain independent.
writeLines(c("# Figure F2d caption", "",
             "> Generated by `scripts/hsf1/22_fig_f2d.R`. Numbers are read from landed tables.", "",
             "## F2d", "", CAP), file.path(fig, "CAPTION_F2d.md"))
say("\ncaption written to CAPTION_F2d.md")
say("=== done ===")
sink()
invisible(NULL)
