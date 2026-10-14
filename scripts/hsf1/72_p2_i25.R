#!/usr/bin/env Rscript
## 72_p2_i25.R -- the (+1, -1, -1, +1) interaction at the 25 °C baseline, I(25), and the baseline-subtracted
## interaction I(c) - I(25).
##
## This analysis forms new contrasts without changing the three stored tables produced by 20_p2_offset_origin.R; it
## reads script 20's model construction by parse() (same count matrix, same filter, same TMM regimes, same voom
## group-means fit) and forms only new contrasts. It validates the model against the four saved heat-condition medians in
## p2_offset_origin_summary.tsv and stops if any differs. The amplitude uses
## 18_p2_amplitude_and_cistrans_identity.R's two-block estimator, also read directly, and is validated against the
## stored p2_interaction_amplitude.tsv before the baseline is added.
##
## Quantities, each under both normalisation regimes where a regime applies:
##   model-based  I(25) and I(c) - I(25) = [D - H - M + W]_c - [D - H - M + W]_25, one contrast of the same fit, for
##                every heat condition c; genome-wide median with SE of the median (the expression that gave +0.1319)
##                and the same medians per class, read against zero and against the not-induced class.
##   amplitude    rms_tru of I_r(25) and of I_r(c) - I_r(25), replicate by replicate, noise-corrected by the
##                signal/null decomposition over replicate pairs (script 18's estimator, published normalisation).
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
say <- function(...) cat(sprintf(...), "\n", sep = "")
## ---- read script 20 through RESP/I_of without running its logging or output expressions ----
src <- parse(file.path(root, "scripts/hsf1/20_p2_offset_origin.R"))
txt <- vapply(src, function(e) paste(deparse(e), collapse = " "), "")
stop_at <- which(grepl("^ALLG <- rownames\\(M\\)", txt))[1]; stopifnot(!is.na(stop_at))
for (k in seq_len(stop_at)) if (!grepl("^(sink|say|logf|fwrite)|^for \\(i in seq_len\\(nrow\\(W", txt[k])) eval(src[[k]])
stopifnot(exists("I_of"), exists("RESP"), exists("M"), exists("CL"))
## ---- validation against saved heat-condition medians ----
LS <- fread(file.path(hsd, "p2_offset_origin_summary.tsv"))
REG <- c(pub = "all genes (as published)", exc = "response classes excluded from normalisation")
NORM <- list(pub = ALLG, exc = setdiff(ALLG, RESP))
IXH <- rbindlist(lapply(names(REG), function(r) I_of(NORM[[r]], REG[[r]])))
chk <- merge(IXH[, .(med = median(I)), by = .(regime, condition)], LS, by = c("regime", "condition"))
stopifnot(nrow(chk) == 8L, max(abs(chk$med - chk$med_I)) < 1e-10)
say("validation: 8 saved heat-condition medians reproduced (max |diff| %.1e)", max(abs(chk$med - chk$med_I)))
## ---- the new contrasts: I(25) and I(c) - I(25), from the same fit as I_of ----
I_base <- function(norm_genes, tag) {
  dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  sub <- dge[rownames(dge) %in% norm_genes, , keep.lib.sizes = TRUE]
  dge$samples$norm.factors <- calcNormFactors(sub)$samples$norm.factors
  grp <- factor(SM$grp); design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  fit <- lmFit(voom(dge, design, plot = FALSE), design)
  mk <- function(c_) { ct <- setNames(rep(0, ncol(design)), colnames(design))
    ct[paste("Double.KDKO", c_, sep = "__")] <- 1; ct[paste("HSF1.KD", c_, sep = "__")] <- -1
    ct[paste("MSN24.KO", c_, sep = "__")] <- -1; ct[paste("Wildtype", c_, sep = "__")] <- 1; ct }
  HEAT <- c(ANCHOR, OTHER)
  CT <- cbind(`25_0` = mk("25_0"), sapply(HEAT, function(c_) mk(c_) - mk("25_0")))
  colnames(CT) <- c("I(25)", paste0("I(", HEAT, ") - I(25)"))
  stopifnot(all(abs(colSums(CT)) < 1e-12))
  cf <- contrasts.fit(fit, CT)
  rbindlist(lapply(seq_len(ncol(CT)), function(j) data.table(regime = tag, gene = rownames(cf$coefficients),
    contrast = colnames(CT)[j], I = cf$coefficients[, j], se = cf$stdev.unscaled[, j] * cf$sigma)))
}
IB <- rbindlist(lapply(names(REG), function(r) I_base(NORM[[r]], REG[[r]])))
SUMB <- IB[, .(n = .N, med_I = median(I), se_of_med = median(se) / sqrt(.N)), by = .(regime, contrast)]
SUMB[, z := med_I / se_of_med]
## Class medians on the method-unanimous, mutually exclusive assignment that P2 reports (21_p2_exclusive_classes.R:
## a gene takes a class only if all methods agree, otherwise "unassigned"), validated against the saved anchor
## medians in p2_class_exclusive.tsv.
CLS3 <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
C3 <- unique(CLS3[, .(gene = get(names(CLS3)[1]), cls = get(intersect(c("class","epistasis_class","cls"), names(CLS3))[1]),
                      method = get(intersect(c("method","meth"), names(CLS3))[1]))])
U <- C3[, .(k = uniqueN(cls), n_meth = uniqueN(method), cls = cls[1]), by = gene]
U[, cls := ifelse(k == 1L & n_meth == uniqueN(C3$method), cls, "unassigned")]; stopifnot(!any(duplicated(U$gene)))
LX <- fread(file.path(hsd, "p2_class_exclusive.tsv"))
ax <- merge(IXH[condition == ANCHOR], U[, .(gene, cls)], by = "gene")[, .(med = median(I)), by = .(regime, cls)]
ax <- merge(ax, LX[, .(regime, cls, med_I)], by = c("regime", "cls")); stopifnot(nrow(ax) == 14L, max(abs(ax$med - ax$med_I)) < 1e-10)
say("validation: 14 saved exclusive-class anchor medians reproduced")
IB <- rbind(IB, IXH[condition == ANCHOR, .(regime, gene, contrast = "I(37_30)", I, se)])
CLB <- merge(IB[contrast %in% c("I(25)", "I(37_30)", "I(37_30) - I(25)")], U[, .(gene, cls)], by = "gene")[
  , .(n = .N, med_I = median(I), se_of_med = median(se) / sqrt(.N)), by = .(regime, contrast, cls)]
CLB <- merge(CLB, CLB[cls == "not_induced", .(regime, contrast, ref = med_I)], by = c("regime", "contrast"))
CLB[, `:=`(vs_ref = med_I - ref, z_vs_0 = med_I / se_of_med, z_vs_ref = (med_I - ref) / se_of_med)]
## ---- amplitude, script 18's estimator (published normalisation) ----
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + factor(SM$grp)); colnames(design) <- levels(factor(SM$grp))
E <- voom(dge, design, plot = FALSE)$E
col_of <- function(g, c, r) which(SM$knockout == g & SM$cond == c & SM$replicate == r)
reps <- sort(unique(SM$replicate))
I_rep <- function(c, r) { i <- vapply(c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"),
  function(g) { k <- col_of(g, c, r); if (length(k) != 1L) NA_integer_ else k }, 0L)
  if (anyNA(i)) return(NULL); E[, i[1]] - E[, i[2]] - E[, i[3]] + E[, i[4]] }
rms_true <- function(s, n) sqrt(max(0, mean(s^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE)))
PR <- combn(reps, 2L)
amp <- function(label, f) { per <- rbindlist(lapply(seq_len(ncol(PR)), function(j) {
    a <- f(PR[1, j]); b <- f(PR[2, j]); sig <- (a + b) / 2; nul <- (b - a) / 2
    data.table(rms_sig = sqrt(mean(sig^2)), rms_null = sqrt(mean(nul^2)), rms_tru = rms_true(sig, nul)) }))
  data.table(quantity = label, n_pairs = nrow(per), rms_sig = mean(per$rms_sig), rms_null = mean(per$rms_null),
             rms_tru = mean(per$rms_tru), rms_tru_min = min(per$rms_tru), rms_tru_max = max(per$rms_tru)) }
HEAT <- c("37_10", "37_30", "42_10", "42_30")
AMPH <- rbindlist(lapply(HEAT, function(c_) amp(paste0("I(", c_, ")"), function(r) I_rep(c_, r))))
LA <- fread(file.path(hsd, "p2_interaction_amplitude.tsv"))
stopifnot(max(abs(AMPH$rms_tru - LA[match(sub("^I\\((.*)\\)$", "\\1", AMPH$quantity), condition), rms_tru])) < 1e-10)
say("validation: 4 saved heat-condition amplitudes reproduced")
AMP <- rbind(AMPH, amp("I(25)", function(r) I_rep("25_0", r)),
             rbindlist(lapply(HEAT, function(c_) amp(paste0("I(", c_, ") - I(25)"), function(r) I_rep(c_, r) - I_rep("25_0", r)))))
AMP[, snr := rms_tru / rms_null]
fwrite(SUMB, file.path(hsd, "p2_i25_summary.tsv"), sep = "\t")
fwrite(CLB, file.path(hsd, "p2_i25_classes.tsv"), sep = "\t")
fwrite(AMP, file.path(hsd, "p2_i25_amplitude.tsv"), sep = "\t")
say("=== 72_p2_i25 | %s ===", format(Sys.time()))
say("\ngenome-wide median, model-based:"); print(SUMB[, .(regime = substr(regime, 1, 12), contrast, n, med_I = round(med_I, 4), se = round(se_of_med, 4), z = round(z, 1))])
say("\nclass medians at the anchor:"); print(CLB[order(regime, cls, contrast), .(regime = substr(regime, 1, 12), contrast, cls, n, med = round(med_I, 4), vs_ref = round(vs_ref, 4), z_ref = round(z_vs_ref, 1))])
say("\namplitude (published normalisation):"); print(AMP[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 4) else x)])
