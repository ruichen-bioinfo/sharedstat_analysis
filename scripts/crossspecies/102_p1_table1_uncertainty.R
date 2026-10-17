#!/usr/bin/env Rscript
## 102_p1_table1_uncertainty.R -- uncertainty and reading comparison for P1 Table 1. Per hybrid, the median over the
## three DE methods of the naive gain (1 + OLS slope) and of the IV gain with its interval (gain_iv.tsv, model A), the
## first-stage correlation of R with the instrument S, the reliability-corrected gain naive/h (h = h_nofilt from
## p1_anchor_diagnosis.tsv), the closure naive - IV x h, and the rank agreement of the three readings. These
## definitions were fixed before this analysis was run.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); outd <- file.path(root, "data/02_processed/crossspecies")
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu"); METHODS <- c("DESeq2","edgeR","limma")
partner <- c(ScxSp="SkxSu",SkxSu="ScxSp",ScxSu="SpxSk",SpxSk="ScxSu",ScxSk="SpxSu",SpxSu="ScxSk")
X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod_p1order.tsv.gz"))))[condition == "heat" & group %in% HYB]   # SpxSk refitted with the first-named allele as reference
Rtab <- X[model == "A" & term == "time", .(gene, group, method, R = eff)]
inst <- function(g, mt) { keep <- setdiff(HYB, c(g, partner[[g]]))
  Rtab[group %in% keep & method == mt, .(S = mean(R, na.rm = TRUE), nS = sum(is.finite(R))), by = gene][nS >= 3] }
FS <- rbindlist(lapply(HYB, function(g) rbindlist(lapply(METHODS, function(mt) {
  ii <- X[group == g & method == mt & model == "A" & term == "inter", .(gene, ei = eff)]
  D <- merge(merge(ii, Rtab[group == g & method == mt, .(gene, R)], by = "gene"), inst(g, mt), by = "gene")
  D <- D[is.finite(ei) & is.finite(R) & is.finite(S)]
  data.table(group = g, method = mt, n = nrow(D), first_stage_r = cor(D$R, D$S)) }))))
GI <- fread(file.path(outd, "gain_iv_p1order.tsv"))[model == "A"]   # first-named allele as reference
chk <- merge(GI, FS, by = c("group", "method")); stopifnot(all(chk$n.x == chk$n.y))
T <- GI[, .(naive = median(1 + slope_ols), iv = median(gain_iv), iv_lo = 1 + median(iv_lo), iv_hi = 1 + median(iv_hi)), by = group]
T <- merge(T, FS[, .(first_stage_r = median(first_stage_r)), by = group], by = "group")
H <- fread(file.path(root, "data/02_processed/hsf1/p1_anchor_diagnosis.tsv"))[, .(group = hybrid, h = h_nofilt)]
T <- merge(T, H, by = "group")[, `:=`(rel_corrected = naive / h, closure = naive - iv * h,
   iv_excludes_1 = iv_lo > 1 | iv_hi < 1)]
T[, `:=`(rank_naive = rank(-abs(naive - 1)), rank_iv = rank(-abs(iv - 1)), rank_rel = rank(-abs(rel_corrected - 1)))]
fwrite(T, file.path(outd, "p1_table1_uncertainty.tsv"), sep = "\t")
print(T[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
cat(sprintf("Spearman IV vs naive %.3f; rel vs naive %.3f; rel vs IV %.3f\n", cor(T$rank_iv, T$rank_naive), cor(T$rank_rel, T$rank_naive), cor(T$rank_rel, T$rank_iv)))
