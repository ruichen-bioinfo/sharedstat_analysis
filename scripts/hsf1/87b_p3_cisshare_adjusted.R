#!/usr/bin/env Rscript
## 87b_p3_cisshare_adjusted.R -- Figure 4: the cis share of divergence variance per stratum and timepoint, script
## 05's statistic (median over the six hybrid cells of cis^2 / (cis^2 + trans^2), noise-corrected). A single gene set,
## so the adjustment reduces to one term: the gene-bootstrap SE on the logit scale times sqrt(VIF), rho-bar measured on per-library
## allele differences . Sensitivities: leave-one-hybrid-out jackknife with t(5); the mean-rho^2 form; the gene
## bootstrap as published. The hybrid-level counts (6 of 6 etc.) are design-level quantities and are not recomputed here.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
EX <- parse(file.path(root, "scripts/hsf1/05_cis_trans_decomposition.R"))
for (e in EX) { s <- paste(deparse(e), collapse = " "); if (grepl("^SH <- list\\(\\)", s)) break
  if (grepl("^sink\\(|^fwrite\\(|quit\\(|^say\\(", s)) next; eval(e, envir = globalenv()) }
stopifnot(is.function(cell_share), exists("SETS"), exists("tp_of"), exists("GENES"))
source(file.path(root, "scripts/hsf1/87_allele_rho.R"))
NB <- 1000L; T5 <- qt(0.975, 5); lg <- function(p) log(p / (1 - p)); il <- function(x) 1 / (1 + exp(-x))
hyb_of <- sub("_[0-9]+$", "", CELLS); HYBS <- unique(hyb_of); stopifnot(length(HYBS) == 6L)
st <- function(i, keep = rep(TRUE, length(CELLS))) { v <- cell_share(i); vapply(TP, function(t) median(v[tp_of == t & keep], na.rm = TRUE), 0) }
LS <- fread(file.path(root, "data/02_processed/hsf1/a6_cis_share_37_30.tsv"))
OUT <- rbindlist(lapply(SETS, function(s) { i <- idx_of(s); i <- i[!is.na(i)]; o <- st(i); set.seed(87)
  b <- t(replicate(NB, st(sample(i, length(i), TRUE)))); p <- allele_rho(GENES[i]); n <- length(i)
  jk <- sapply(HYBS, function(h) lg(st(i, hyb_of != h)))
  rbindlist(lapply(seq_along(TP), function(j) { sb <- sd(lg(b[, j]), na.rm = TRUE); x <- lg(o[j])
    se <- sb * sqrt(1 + (n - 1) * p[["rho"]]); se2 <- sb * sqrt(1 + (n - 1) * p[["rho2"]])
    jv <- jk[j, ]; sej <- if (all(is.finite(jv))) sqrt(5 / 6 * sum((jv - mean(jv))^2)) else NA_real_
    data.table(stratum = s, time = TP[j], n = n, cis_share = o[j], lo = il(x - 1.96 * se), hi = il(x + 1.96 * se), rho = p[["rho"]],
               lo_rho2 = il(x - 1.96 * se2), hi_rho2 = il(x + 1.96 * se2), lo_jk_t5 = il(x - T5 * sej), hi_jk_t5 = il(x + T5 * sej),
               lo_boot = quantile(b[, j], .025, na.rm = TRUE), hi_boot = quantile(b[, j], .975, na.rm = TRUE)) })) }))
ck <- merge(OUT[, .(stratum, time, cis_share)], LS[, .(stratum, time, landed = cis_share)], by = c("stratum", "time"))
stopifnot(nrow(ck) == nrow(LS), max(abs(ck$cis_share - ck$landed)) < 1e-9)
cat(sprintf("fidelity: %d landed cis shares reproduced\n", nrow(ck)))
fwrite(OUT, file.path(root, "data/02_processed/hsf1/p3_adjusted_cisshare.tsv"), sep = "\t")
## the primary P3 comparison: the Hsf1-core cis share falls from 15 to 60 min (0.44 to 0.21)
H <- OUT[stratum == "Hsf1_core"]
cat(sprintf("Hsf1 core: 15 min %.3f [%.3f, %.3f] vs 60 min %.3f [%.3f, %.3f]; adjusted intervals overlap: %s\n",
    H[time == 15, cis_share], H[time == 15, lo], H[time == 15, hi], H[time == 60, cis_share], H[time == 60, lo], H[time == 60, hi],
    H[time == 15, lo] < H[time == 60, hi]))
print(OUT[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
