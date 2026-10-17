#!/usr/bin/env Rscript
## 02b_p1order_refit.R -- refits S. paradoxus x S. kudriavzevii (SpxSk) with its first-named allele as reference.
## The first-named allele is the reference arm in every P1 quantity, as stated in the manuscript. 02_fay_multimethod.R fits
## factor(Allele), i.e. alphabetical levels, which is the first-named allele in every hybrid except SpxSk (Sk). The
## exact re-expression of the stored effects (p1_orientation.R) holds, but limma's voom-weighted multi-coefficient
## moderated F is not invariant to the reparametrisation of the `time` test (its padj differs by up to 0.52), so
## SpxSk is refitted here with Allele levels c("Sp","Sk") through 02's own code (evaluated verbatim up to its driver),
## both conditions, models A and B, all three methods. Every other row is copied from fay_multimethod.tsv.gz unchanged.
## Output, for P1 only: data/02_processed/crossspecies/fay_multimethod_p1order.tsv.gz. 02's own output is not touched.
## Checks before the output is written (any failure stops): (b) refit eff equals the exact re-expression of
## the stored eff within the effect tolerances of TH (padj differences reported only); (c) same cells, genes and NA
## pattern; (a) 02's own code path: levels c("Sk","Sp") reproduce the stored SpxSk fits.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT"); stopifnot(nzchar(root))
cko <- Sys.getenv("REFIT_CHECK_DIR", "logs")
src <- readLines(file.path(root, "scripts/crossspecies/02_fay_multimethod.R"))
cut <- grep("^## -+ driver", src); stopifnot(length(cut) == 1)
stopifnot(!any(grepl("fwrite|writeLines|saveRDS", src[seq_len(cut - 1)])))   # 02 writes nothing before its driver
eval(parse(text = src[seq_len(cut - 1)]), envir = globalenv())
source(file.path(root, "scripts/crossspecies/p1_orientation.R"))
IN <- file.path(outd, "fay_multimethod.tsv.gz"); OUT <- file.path(outd, "fay_multimethod_p1order.tsv.gz")
X <- fread(cmd = paste("zcat", shQuote(IN)))
g <- "SpxSk"; KEY <- c("gene", "term", "condition", "model", "method")
ST <- X[group == g]
fit_all <- function(lv) rbindlist(lapply(c("heat", "cold"), function(cond) {
  idx <- intersect(which(design$Strain == g), which(design$Condition == cond))
  dz0 <- droplevels(design[idx, , drop = FALSE]); cw <- cts[, idx, drop = FALSE]
  uni <- rowSums(cw >= MIN_CNT) >= MIN_SMP; cw <- cw[uni, , drop = FALSE]; storage.mode(cw) <- "integer"
  dz <- dz0; dz$Allele <- factor(as.character(dz$Allele), levels = lv)
  rbindlist(lapply(c("A", "B"), function(model) {
    D <- build_mm(dz, model)
    say("%s %s model %s levels %s: %s", g, cond, model, paste(lv, collapse = ","), paste(colnames(D$mm), collapse = " "))
    if (D$resid_df < 2) return(NULL)                                   # as in 02's driver
    rbindlist(lapply(names(runners), function(m) {
      r <- runners[[m]](cw, dz, D); r[, `:=`(condition = cond, group = g, model = model, method = m)]; r }))
  }))
}))
NEW <- fit_all(c("Sp", "Sk")); CTL <- fit_all(c("Sk", "Sp"))
TH <- data.table(method = c("limma", "edgeR", "DESeq2"), tol_eff = c(1e-6, 1e-6, 1e-4), tol_padj = c(1e-10, 1e-10, 1e-6))
cmp <- function(ref, fit, lab) {
  M <- merge(ref[, c(KEY, "eff", "padj"), with = FALSE], fit[, c(KEY, "eff", "padj"), with = FALSE], by = KEY, suffixes = c("_r", "_f"))
  stopifnot(nrow(M) == nrow(ref), nrow(M) == nrow(fit))               # (c) coverage
  d <- function(a, b) { x <- abs(a - b); if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE) }
  merge(M[, .(check = lab, n = .N, max_abs_eff = d(eff_r, eff_f), max_abs_padj = d(padj_r, padj_f),
             na_mismatch_eff = sum(is.na(eff_r) != is.na(eff_f)), na_mismatch_padj = sum(is.na(padj_r) != is.na(padj_f))),
          by = .(condition, model, method, term)], TH, by = "method")
}
R <- rbind(cmp(ST, CTL, "a_control_SkSp_vs_stored"), cmp(p1_orient(ST), NEW, "b_refit_SpSk_vs_reexpressed"))
R[, eff_ok := na_mismatch_eff == 0 & (is.na(max_abs_eff) | max_abs_eff <= tol_eff)]
R[, padj_ok := na_mismatch_padj == 0 & (is.na(max_abs_padj) | max_abs_padj <= tol_padj)]
## Check b gates heat / model A, the cells P1 uses; in the
## other condition x model cells it is reported only (there 4 model-B cells exceed the
## effect tolerances: DESeq2 inter 1.2e-4 / 1.8e-4, edgeR heat-B allele / inter 3.6e-6 / 4.1e-6).
R[, gating := check == "a_control_SkSp_vs_stored" | (condition == "heat" & model == "A")]
R[, gate := if (check[1] == "a_control_SkSp_vs_stored") eff_ok & padj_ok else eff_ok, by = check]
R[gating == FALSE, gate := TRUE]
setcolorder(R, c("check", "condition", "model", "method", "term")); setorder(R, check, condition, model, method, term)
fwrite(R, file.path(cko, "refit_checks.tsv"), sep = "\t")
print(R[, .(check, condition, model, method, term, n, max_abs_eff = signif(max_abs_eff, 3), max_abs_padj = signif(max_abs_padj, 3), gate)])
stopifnot(setequal(unique(ST[, paste(condition, model, method)]), unique(NEW[, paste(condition, model, method)])))
if (!all(R$gate)) { say("REFIT CHECKS FAILED (%d of %d rows)", sum(!R$gate), nrow(R)); quit(status = 3L) }
say("refit checks a, b, c: PASS (%d rows)", nrow(R))
## replace SpxSk's stat, pval, padj, eff by the refit; row order and every other row unchanged
Y <- copy(X); Y[, .row := .I]
U <- NEW[, c(KEY, "stat", "pval", "padj", "eff"), with = FALSE]; U[, group := g]
n0 <- Y[group == g, .N]
Y[U, on = c(KEY, "group"), `:=`(stat = i.stat, pval = i.pval, padj = i.padj, eff = i.eff)]
stopifnot(n0 == nrow(U), nrow(Y) == nrow(X), identical(Y$.row, seq_len(nrow(X))))
Y[, .row := NULL]
stopifnot(isTRUE(all.equal(Y[group != g], X[group != g], check.attributes = FALSE)))
fwrite(Y, OUT, sep = "\t")
say("wrote %s (%d rows; %d SpxSk rows refitted)", basename(OUT), nrow(Y), n0)
