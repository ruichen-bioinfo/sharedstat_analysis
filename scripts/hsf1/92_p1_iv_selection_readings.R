#!/usr/bin/env Rscript
## 92_p1_iv_selection_readings.R -- IV gain under three gene-selection readings; the primary reading was
## fixed before this analysis was run.
## The Table 1 IV keeps only genes with reference-arm padj<0.05
## (selection on the regressor R-hat). This runs three readings on the same data:
##   current: padj_t < 0.05           (as in 11_gain_iv.R)
##   (i)      no significance filter   (finite values only; the primary reading)
##   (ii)     instrument S significant (padj on the pooled instrument)
## Reuses 11_gain_iv.R's Rtab / instrument / IV construction verbatim. No existing table is overwritten.
suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
outd <- file.path(root, "data/02_processed/crossspecies")
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu")
METHODS <- c("DESeq2","edgeR","limma"); ALPHA <- 0.05
partner <- c(ScxSp="SkxSu",SkxSu="ScxSp",ScxSu="SpxSk",SpxSk="ScxSu",ScxSk="SpxSu",SpxSu="ScxSk")
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(20260905)

X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod.tsv.gz"))))
X <- X[condition == "heat" & group %in% HYB]
Rtab <- X[model == "A" & term == "time", .(gene, group, method, R = eff, padj_t = padj)]
instrument <- function(g, mt) {
  keep <- setdiff(HYB, c(g, partner[[g]]))
  # (ii) needs a significance handle on the instrument: use the min padj_t among contributing crosses
  s <- Rtab[group %in% keep & method == mt, .(S = mean(R, na.rm = TRUE),
            nS = sum(is.finite(R)), padj_S = suppressWarnings(min(padj_t, na.rm = TRUE))), by = gene]
  s[nS >= 3]
}
ivslope <- function(D) cov(D$ei, D$S) / cov(D$R, D$S)

readings <- c("current", "no_filter", "instrument")
out <- list()
for (rd in readings) for (g in HYB) for (mt in METHODS) {
  ii <- X[group == g & method == mt & model == "A" & term == "inter", .(gene, ei = eff)]
  D <- merge(merge(ii, Rtab[group == g & method == mt, .(gene, R, padj_t)], by = "gene"),
             instrument(g, mt), by = "gene")
  D <- D[is.finite(ei) & is.finite(R) & is.finite(S)]
  n_all <- nrow(D)
  if (rd == "current")    D <- D[!is.na(padj_t)  & padj_t  < ALPHA]
  if (rd == "instrument") D <- D[is.finite(padj_S) & padj_S < ALPHA]
  # no_filter: keep all finite
  if (nrow(D) < 200) next
  iv  <- ivslope(D)
  ols <- unname(coef(lm(ei ~ R, D))[2])
  corRS <- cor(D$R, D$S)
  out[[length(out)+1L]] <- data.table(reading = rd, cross = g, method = mt,
    n_kept = nrow(D), n_all = n_all, slope_ols = ols, slope_iv = iv,
    gain_1p = 1 + iv, artefact_pct = 100*abs(ols-iv)/abs(ols), cor_RS = corRS)
}
R <- rbindlist(out)
fwrite(R, file.path(outd, "p1_iv_selection_readings.tsv"), sep = "\t")
say("wrote p1_iv_selection_readings.tsv (%d rows)\n", nrow(R))

## median over methods per reading x cross, gain = 1 + slope_iv (the P1 estimand, as in 11_gain_iv.R)
say("=== corrected gain (1 + median slope_IV) by reading and hybrid ===")
say("%-8s %8s %8s %8s", "hybrid", "current", "no_filt", "instr")
med <- function(rd, g) { x <- R[reading==rd & cross==g]; if(nrow(x)) 1+median(x$slope_iv) else NA }
for (g in HYB) say("%-8s %8.3f %8.3f %8.3f", g, med("current",g), med("no_filter",g), med("instrument",g))

say("\n=== genes kept (median over methods), by reading x hybrid ===")
nk <- function(rd,g){x<-R[reading==rd&cross==g];if(nrow(x))median(x$n_kept) else NA}
say("%-8s %8s %8s %8s", "hybrid","current","no_filt","instr")
for (g in HYB) say("%-8s %8.0f %8.0f %8.0f", g, nk("current",g), nk("no_filter",g), nk("instrument",g))

say("\n=== cor(R, S) diagnostic (median over methods): does the filter move it? ===")
cr <- function(rd,g){x<-R[reading==rd&cross==g];if(nrow(x))median(x$cor_RS) else NA}
say("%-8s %8s %8s %8s", "hybrid","current","no_filt","instr")
for (g in HYB) say("%-8s %8.3f %8.3f %8.3f", g, cr("current",g), cr("no_filter",g), cr("instrument",g))

say("\n=== sign flip vs gain=1, per reading (how many of 6 hybrids have corrected gain>1) ===")
for (rd in readings) { v <- sapply(HYB, function(g) med(rd,g)); say("  %-11s gains: %s  | #>1: %d",
  rd, paste(sprintf("%.3f", v), collapse=" "), sum(v>1, na.rm=TRUE)) }
say("=== done (primary reading = no_filter) ===")
