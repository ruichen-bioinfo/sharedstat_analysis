#!/usr/bin/env Rscript
## 77_p2_sensitivity.R -- (i) class sizes and the redundant-class partition at alpha = 0.01, 0.05, 0.10, from the
## stored per-method differential-expression tables and 01_epistasis_decomposition.R's own classify() (read by
## parse, with ALPHA rebound), validated at 0.05 against the saved class file; (ii) class medians of I read against an
## alternative reference class, transcription-factor-independent (heat-induced, requiring neither regulator), beside
## the not-induced reference, on the level scale I(37_30) and on the response scale I(37_30) - I(25).
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
src <- parse(file.path(root, "scripts/hsf1/01_epistasis_decomposition.R"))
txt <- vapply(src, function(e) paste(deparse(e), collapse = " "), "")
want <- grepl("^(proc|outd|ref|METHODS|gtf|gl|sysn|stdn|nmap|lut|load_comp|classify) <-", txt); ANCHOR <- "37_30"
root_ <- root; for (e in src[want]) eval(e); root <- root_
ALPHA_SET <- c(0.01, 0.05, 0.10)
run <- function(a) { ALPHA <<- a; environment(classify)$ALPHA <- a
  CL <- rbindlist(lapply(METHODS, classify)); CL[, alpha := a]; CL }
ALL <- rbindlist(lapply(ALPHA_SET, run))
## validation at 0.05
LC <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
chk <- merge(ALL[alpha == 0.05, .N, by = .(method, class)], LC[, .N, by = .(method, class)], by = c("method", "class"))
stopifnot(nrow(chk) == 18L, all(chk$N.x == chk$N.y))
cat("fidelity: class counts at alpha 0.05 reproduce the saved class file in all 18 method x class cells\n")
excl <- function(CL) { U <- CL[, .(k = uniqueN(class), nm = uniqueN(method), cls = class[1]), by = gene]
  U[, cls := ifelse(k == 1L & nm == 3L, cls, "unassigned")]; U }
SZ <- rbindlist(lapply(ALPHA_SET, function(a) excl(ALL[alpha == a])[, .N, by = cls][, alpha := a]))
SZw <- dcast(SZ, cls ~ alpha, value.var = "N", fill = 0L)
PT <- ALL[class == "redundant", .N, by = .(alpha, method, sub)][, pct := 100 * N / sum(N), by = .(alpha, method)]
PTw <- dcast(PT, alpha + method ~ sub, value.var = "pct")
fwrite(SZ, file.path(hsd, "p2_alpha_class_sizes.tsv"), sep = "\t")
fwrite(PT, file.path(hsd, "p2_alpha_redundant_partition.tsv"), sep = "\t")
cat("\nexclusive class sizes by alpha:\n"); print(SZw)
cat("\nredundant-class partition (% of the class, per method) by alpha:\n"); print(PTw[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
## (ii) alternative reference class, published normalisation, exclusive classes at alpha 0.05
IX <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "p2_offset_origin_genewise.tsv.gz"))))[regime == "all genes (as published)"]
I25C <- fread(file.path(hsd, "p2_i25_classes.tsv"))
U <- excl(ALL[alpha == 0.05])
D <- merge(IX[, .(gene, I, se)], U[, .(gene, cls)], by = "gene")
lvl <- D[, .(n = .N, med = median(I), se_med = median(se) / sqrt(.N)), by = cls]
ref_alt <- lvl[cls == "TF_independent", med]; ref_ni <- lvl[cls == "not_induced", med]
lvl[, `:=`(scale = "level I(37_30)", vs_not_induced = med - ref_ni, vs_TF_independent = med - ref_alt)]
rsp <- I25C[regime == "all genes (as published)" & contrast == "I(37_30) - I(25)", .(cls, n, med = med_I, se_med = se_of_med)]
ra <- rsp[cls == "TF_independent", med]; rn <- rsp[cls == "not_induced", med]
rsp[, `:=`(scale = "response I(37_30) - I(25)", vs_not_induced = med - rn, vs_TF_independent = med - ra)]
RC <- rbind(lvl, rsp, fill = TRUE)[, `:=`(z_vs_not_induced = vs_not_induced / se_med, z_vs_TF_independent = vs_TF_independent / se_med)]
fwrite(RC, file.path(hsd, "p2_reference_class_sensitivity.tsv"), sep = "\t")
cat("\nclass medians against two reference classes (published normalisation):\n")
print(RC[order(scale, -n), .(scale, cls, n, vs_NI = round(vs_not_induced, 3), z_NI = round(z_vs_not_induced, 1),
                             vs_TFi = round(vs_TF_independent, 3), z_TFi = round(z_vs_TF_independent, 1))])
## (iii) how the classes were defined, and what that implies. Requirement for a regulator is a level difference at the
## anchor (HSF1.KD vs Wildtype at 37_30), not a difference in heat response; so a gene whose Hsf1-depleted level is
## already lower at 25 °C can be called Hsf1-only. For each exclusive class: the share of genes whose Hsf1-depleted
## level is already significantly lower than wild type at 25 °C (registered comparison HSF1.KD_vs_Wildtype_at_25_0),
## per method, and the median of that 25 °C log2 fold change. Reads the stored DE tables only.
B <- rbindlist(lapply(METHODS, function(m) { d <- load_comp(m, "HSF1.KD_vs_Wildtype_at_25_0"); d[, method := m]; d }))
B <- merge(B, U[, .(gene, cls)], by = "gene")
BS <- B[, .(n = .N, pct_lower_at_25 = 100 * mean(!is.na(padj) & padj < 0.05 & lfc < 0), med_lfc_25 = median(lfc, na.rm = TRUE)),
        by = .(cls, method)]
fwrite(BS, file.path(hsd, "p2_basal_hsf1_dependence.tsv"), sep = "\t")
cat("\nHsf1-depleted vs wild type at 25 °C, by exclusive class (per method):\n")
print(dcast(BS, cls ~ method, value.var = "pct_lower_at_25")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])
print(dcast(BS, cls ~ method, value.var = "med_lfc_25")[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
## (iv) I(25) split into its two halves, I = (D - M) - (H - W): the Hsf1-depletion effect in the msn2/4-deleted
## background minus the same effect in wild type, at 25 °C, per method, class median minus the not-induced median.
H2 <- rbindlist(lapply(METHODS, function(m) {
  hw <- load_comp(m, "HSF1.KD_vs_Wildtype_at_25_0")[, .(gene, hw = lfc)]
  md <- load_comp(m, "MSN24.KO_vs_Double.KDKO_at_25_0")[, .(gene, dm = -lfc)]
  merge(hw, md, by = "gene")[, method := m] }))
H2 <- merge(H2, U[, .(gene, cls)], by = "gene")
HS <- H2[, .(n = .N, hsf1_effect_in_wt = median(hw, na.rm = TRUE), hsf1_effect_in_msn24ko = median(dm, na.rm = TRUE)), by = .(cls, method)]
HS <- merge(HS, HS[cls == "not_induced", .(method, ref_wt = hsf1_effect_in_wt, ref_ko = hsf1_effect_in_msn24ko)], by = "method")
HS[, `:=`(wt_vs_ref = hsf1_effect_in_wt - ref_wt, ko_vs_ref = hsf1_effect_in_msn24ko - ref_ko)]
fwrite(HS, file.path(hsd, "p2_i25_halves.tsv"), sep = "\t")
cat("\nI(25) halves: Hsf1-depletion effect at 25 °C in wild type and in msn2/4-deleted background, class median minus not-induced median:\n")
print(HS[order(cls, method), .(cls, method, n, wt_vs_ref = round(wt_vs_ref, 3), ko_vs_ref = round(ko_vs_ref, 3))])
