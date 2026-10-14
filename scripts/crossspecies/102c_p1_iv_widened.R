#!/usr/bin/env Rscript
## 102c_p1_iv_widened.R -- widens the P1 Table 1 IV intervals (p1_table1_uncertainty.tsv) about the IV estimate by
## sqrt(VIF), VIF = 1 + (n - 1) rho, with rho the gene-weighted inter-gene correlation per hybrid from
## p1_intergene_correlation.tsv (primary rho, and the allele-level rho beside it), and counts the hybrids whose widened
## interval excludes 1. The widening was fixed before this analysis was run.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT"); cs <- file.path(root, "data/02_processed/crossspecies"); hs <- file.path(root, "data/02_processed/hsf1")
T <- fread(file.path(cs, "p1_table1_uncertainty.tsv"))
N <- fread(file.path(cs, "gain_iv.tsv"))[model == "A", .(n = median(n)), by = group]
R <- fread(file.path(hs, "p1_intergene_correlation.tsv"))[, .(rho = weighted.mean(rho, n_genes), rho_al = weighted.mean(rho_allele_level, n_genes)), by = .(group = hybrid)]
T <- merge(merge(T, N, by = "group"), R, by = "group")
w <- function(v) sqrt(1 + (T$n - 1) * v)
T[, `:=`(vif_sqrt = w(rho), lo_w = iv - w(rho) * (iv - iv_lo), hi_w = iv + w(rho) * (iv_hi - iv),
         lo_w_al = iv - w(rho_al) * (iv - iv_lo), hi_w_al = iv + w(rho_al) * (iv_hi - iv))]
T[, `:=`(excl_w = lo_w > 1 | hi_w < 1, excl_w_al = lo_w_al > 1 | hi_w_al < 1)]
fwrite(T, file.path(cs, "p1_table1_iv_widened.tsv"), sep = "\t")
print(T[, .(group, n, rho = signif(rho, 2), sqrtVIF = round(vif_sqrt, 2), iv = round(iv, 3), lo_w = round(lo_w, 3), hi_w = round(hi_w, 3), excl_w, lo_al = round(lo_w_al, 3), hi_al = round(hi_w_al, 3), excl_w_al)])
cat(sprintf("resolved below/above 1: primary widening %d of 6; allele-level widening %d of 6\n", sum(T$excl_w), sum(T$excl_w_al)))
