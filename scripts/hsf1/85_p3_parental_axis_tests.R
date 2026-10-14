#!/usr/bin/env Rscript
## 85_p3_parental_axis_tests.R -- P3 parental-axis tests, run once. The tests, their uncertainty and the multiplicity family below were
## fixed before this analysis was run.
##
## Every quadratic statistic here is a ratio of noise-corrected rms divergences. Its uncertainty:
##   SE_rho       : on the log scale, SE^2 = VIF_num SE_num^2 + VIF_den SE_den^2; SE_num and SE_den from separate gene
##                  bootstraps of the numerator and denominator amplitudes; VIF = 1 + (n - 1) rho, rho the set's mean
##                  inter-gene residual correlation of the parental log fold changes (computed below).
##   SE_jk        : leave-one-species-out jackknife of log ratio (4 variants, each removing the 3 pairs with that species).
##                  Leave-one-block-out is not computable: the noise-corrected amplitude needs both replicate blocks
##                  (the null contrast is (b2 - b1)/2), so a single block has no noise floor. No substitute is used.
##   SE_used      = max(SE_rho, SE_jk); z = log(ratio) / SE_used; two-sided p; interval exp(log r +/- 1.96 SE_used).
##   SE_rho2      : the same with VIF = 1 + (n - 1) mean(rho^2), bias-corrected by 1/df; its verdict is reported beside.
## A1 whole-transcriptome amplitude matching (procedure of 30_amplitude_matched_module_test.R, pool = every
##    non-metabolism gene with divergence data);
## A2 uncollapsed module-label permutation at beta = 0.657 (1000 permutations), reported; the module inference
##    uses the SE_rho test instead;
## A3 the measured exponent is primary rather than beta = 0 (the measured exponent was named primary), so the family is
##    5 modules x 3 timepoints x the 6 exponents the paper reports (0, 0.25, 0.5, 0.657, 0.75, 1): 90 tests, Holm;
## A4 (i) leave-one-child-out of GO:0005975 (direct children with >= 1 module gene; genes whose direct annotations in the
##    subtree lie under that child only are removed); (ii) KEGG central carbon pathways as the module definition.
## Cold confirmation: the cold arm's metabolism ratio at 30 min, beta = 0, with the same adjusted interval.
suppressPackageStartupMessages({ library(data.table); library(limma) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd <- file.path(root, "data/02_processed/hsf1"); csd <- file.path(root, "data/02_processed/crossspecies")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
## ---- transplant script 28 up to DV_heat / DV_cold (skips its sink, log lines and writes) ----
EX <- parse(file.path(root, "scripts/hsf1/28_cold_arm_specificity.R"))
for (e in EX) {
  s <- paste(deparse(e), collapse = " ")
  if (grepl("^sink\\(|quit\\(|run_arm|COLD <-|fwrite\\(|^say\\(", s)) next
  eval(e)
  if (exists("TESTMODS") && exists("DV_cold")) break
}
stopifnot(exists("DV_heat"), exists("DV_cold"), exists("MOD"), is.function(ratio_beta), is.function(parental_lfc))
PUB <- fread(file.path(hsd, "s4_module_divergence_beta.tsv"))
chk <- rbindlist(lapply(TP, function(t) data.table(time = t, module = TESTMODS, r = as.numeric(ratio_beta(DV_heat[time == t], TESTMODS, 0)))))
chk <- merge(chk, PUB[abs(beta) < 1e-9], by = c("time", "module")); stopifnot(max(abs(chk$r - chk$ratio)) < 1e-9)
say("fidelity: %d saved beta = 0 module ratios reproduced", nrow(chk))
BETA_FAM <- c(0, 0.25, 0.5, 0.657, 0.75, 1)
## 134 vs 154
say("metabolism: %d orthologous groups in the module table, %d with divergence data in the heat arm (the %d others are outside the parental gene universe)",
    MOD[module == "metabolism", .N], uniqueN(DV_heat[module == "metabolism", hog]), MOD[module == "metabolism", .N] - uniqueN(DV_heat[module == "metabolism", hog]))
## ---- inter-gene correlation on parental log fold changes (species x block x time, residual 12 df) ----
LFCW <- function(cond) { L <- parental_lfc(cond); L[, col := paste(species, block, time, sep = "_")]
  W <- dcast(L, hog ~ col, value.var = "lfc"); M <- as.matrix(W[, -1]); rownames(M) <- W$hog; M[complete.cases(M), ] }
DES_of <- function(M) { g <- factor(sub("_block[12]_", "_", colnames(M))); model.matrix(~ 0 + g) }
RES <- list(heat = LFCW("heat"), cold = LFCW("cold"))
set.seed(341)
rho_of <- function(arm, genes) {
  M <- RES[[arm]]; g <- intersect(genes, rownames(M)); if (length(g) < 5) return(c(rho = NA, rho2 = NA, n = length(g)))
  if (length(g) > 1500) g <- sample(g, 1500)
  des <- DES_of(M); rho <- interGeneCorrelation(M[g, ], des)$correlation
  R <- M[g, ] - M[g, ] %*% des %*% solve(crossprod(des)) %*% t(des); df <- ncol(M) - ncol(des)
  pr <- cbind(sample(g, 20000, TRUE), sample(g, 20000, TRUE)); pr <- pr[pr[, 1] != pr[, 2], , drop = FALSE]
  r <- vapply(seq_len(nrow(pr)), function(k) cor(R[pr[k, 1], ], R[pr[k, 2], ]), 0)
  c(rho = max(0, rho), rho2 = max(0, mean(r^2) - 1 / df), n = length(g))
}
RHO_CACHE <- new.env()
rho_set <- function(arm, genes) { key <- paste(arm, digest_key(genes)); if (!is.null(RHO_CACHE[[key]])) return(RHO_CACHE[[key]])
  v <- rho_of(arm, genes); RHO_CACHE[[key]] <- v; v }
digest_key <- function(g) digest::digest(sort(unique(g)), algo = "md5")   # md5 of the sorted set; a key of size + character-code sum can collide
## ---- the uncertainty engine for "set A against set B" on one arm, one time, one beta ----
NB <- 1000L
amp_set <- function(D, beta) { sc <- pmax(D$amp, 1e-9)^beta; rms_true(D$d_sig / sc, D$d_nul / sc) }
ratio_unc <- function(D, inA, inB, beta, arm, seed = 1) {
  A <- D[inA]; B <- D[inB]; ra <- amp_set(A, beta); rb <- amp_set(B, beta)
  r <- if (rb > 0) ra / rb else NA_real_
  out <- list(ratio = r)
  if (!is.finite(r) || r <= 0) return(c(out, list(se_used = NA, lo = NA, hi = NA, z = NA, p = NA, se_rho = NA, se_jk = NA, se_rho2 = NA,
                                                 lo_boot = NA, hi_boot = NA, lo_rho2 = NA, hi_rho2 = NA, nA = uniqueN(A$hog), nB = uniqueN(B$hog))))
  set.seed(seed)
  bs <- function(X) { hs <- unique(X$hog); hi <- split(seq_len(nrow(X)), X$hog)
    vapply(seq_len(NB), function(k) amp_set(X[unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)], beta), 0) }
  ba <- bs(A); bb <- bs(B)
  la <- log(ba[ba > 0]); lb <- log(bb[bb > 0]); se_a <- sd(la); se_b <- sd(lb)
  lr <- log(ba / bb); lr <- lr[is.finite(lr)]
  pa <- rho_set(arm, unique(A$hog)); pb <- rho_set(arm, unique(B$hog))
  vif <- function(p, n, which) 1 + (n - 1) * p[[which]]
  se_rho <- sqrt(vif(pa, uniqueN(A$hog), "rho") * se_a^2 + vif(pb, uniqueN(B$hog), "rho") * se_b^2)
  se_rho2 <- sqrt(vif(pa, uniqueN(A$hog), "rho2") * se_a^2 + vif(pb, uniqueN(B$hog), "rho2") * se_b^2)
  jk <- vapply(SPECIES, function(s) { keep <- !grepl(s, D$pair, fixed = TRUE)
    a2 <- amp_set(D[inA & keep], beta); b2 <- amp_set(D[inB & keep], beta); if (a2 > 0 && b2 > 0) log(a2 / b2) else NA_real_ }, 0)
  se_jk <- if (all(is.finite(jk))) sqrt(3 / 4 * sum((jk - mean(jk))^2)) else NA_real_
  se_used <- max(se_rho, se_jk, na.rm = TRUE)
  c(out, list(se_used = se_used, lo = exp(log(r) - 1.96 * se_used), hi = exp(log(r) + 1.96 * se_used),
              z = log(r) / se_used, p = 2 * pnorm(-abs(log(r) / se_used)), se_rho = se_rho, se_jk = se_jk, se_rho2 = se_rho2,
              lo_boot = exp(quantile(lr, .025, names = FALSE)), hi_boot = exp(quantile(lr, .975, names = FALSE)),
              lo_rho2 = exp(log(r) - 1.96 * se_rho2), hi_rho2 = exp(log(r) + 1.96 * se_rho2),
              jk_min = exp(min(jk)), jk_max = exp(max(jk)), nA = uniqueN(A$hog), nB = uniqueN(B$hog)))
}
as_row <- function(x) as.data.table(lapply(x, function(v) if (is.null(v)) NA_real_ else v))
## ---- A3 family (90 tests): heat arm, module vs the other four ----
FAM <- rbindlist(lapply(TP, function(t) { D <- DV_heat[time == t]
  rbindlist(lapply(BETA_FAM, function(be) rbindlist(lapply(TESTMODS, function(m) {
    x <- ratio_unc(D, D$module == m, D$module != m, be, "heat", seed = 1000 * t + round(100 * be)); cbind(data.table(time = t, beta = be, module = m), as_row(x)) })))) }))
FAM[, p_holm := p.adjust(p, "holm")]
FAM[, p_holm_rho2 := p.adjust(2 * pnorm(-abs(log(ratio) / se_rho2)), "holm")]
fwrite(FAM, file.path(hsd, "p3_parental_axis_family.tsv"), sep = "\t")
say("A3 family: %d tests (5 modules x 3 timepoints x %d exponents)", nrow(FAM), length(BETA_FAM))
## ---- A2 permutation at beta = 0.657, uncollapsed ----
A2 <- rbindlist(lapply(TP, function(t) { DTt <- DV_heat[time == t]; hmod <- unique(DTt[, .(hog, module)])
  obs <- ratio_beta(DTt, TESTMODS, 0.657); set.seed(657 + t)
  nl <- t(vapply(seq_len(1000), function(k) { hp <- copy(hmod)[, module := sample(module)]
    ratio_beta(merge(DTt[, .(hog, d_sig, d_nul, amp)], hp, by = "hog"), TESTMODS, 0.657) }, numeric(length(TESTMODS))))
  data.table(time = t, module = TESTMODS, ratio = as.numeric(obs),
             p_perm = vapply(seq_along(TESTMODS), function(j) mean(abs(nl[, j] - 1) >= abs(obs[j] - 1), na.rm = TRUE), 0)) }))
A2[, p_perm_holm15 := p.adjust(pmax(p_perm, 1 / 1000), "holm")]
fwrite(A2, file.path(hsd, "p3_parental_axis_A2_permutation_0657.tsv"), sep = "\t")
## ---- A1 whole-transcriptome amplitude matching ----
MOD_saved <- MOD
MOD <- data.table(hog = UNIV)[, module := ifelse(hog %in% MOD_saved[module == "metabolism", hog], "metabolism", "other")]
DV_all <- build_DV(parental_lfc("heat")); MOD <- MOD_saved
TOL <- 0.10; A1 <- list(); A1P <- list()
for (t in TP) { D <- DV_all[time == t]; AG <- unique(D[, .(hog, module, amp)])
  tgt <- AG[module == "metabolism"][order(amp)]; pool <- AG[module != "metabolism"][order(amp)]; used <- rep(FALSE, nrow(pool)); mt <- mc <- integer(0)
  for (i in seq_len(nrow(tgt))) { cand <- which(!used & abs(pool$amp - tgt$amp[i]) <= TOL * tgt$amp[i]); if (!length(cand)) next
    j <- cand[which.min(abs(pool$amp[cand] - tgt$amp[i]))]; used[j] <- TRUE; mt <- c(mt, i); mc <- c(mc, j) }
  MT <- tgt[mt]; MC <- pool[mc]
  x <- ratio_unc(D, D$hog %in% MT$hog, D$hog %in% MC$hog, 0, "heat", seed = 30 + t)
  ## the pair bootstrap of 30_amplitude_matched_module_test.R (seed 30), reported beside
  set.seed(30); bt <- replicate(2000, { k <- sample(seq_len(nrow(MT)), nrow(MT), TRUE); rm2 <- D[hog %in% MT$hog[k]]; rc2 <- D[hog %in% MC$hog[k]]
    v <- rms_true(rc2$d_sig, rc2$d_nul); if (v > 0) rms_true(rm2$d_sig, rm2$d_nul) / v else NA_real_ })
  A1[[length(A1) + 1]] <- cbind(data.table(time = t, n_metab = nrow(tgt), n_matched = nrow(MT), pool = nrow(pool),
    amp_med_metab = median(MT$amp), amp_med_ctrl = median(MC$amp)), as_row(x),
    data.table(lo_pairboot = quantile(bt, .025, na.rm = TRUE), hi_pairboot = quantile(bt, .975, na.rm = TRUE)))
  A1P[[length(A1P) + 1]] <- data.table(time = t, hog_metab = MT$hog, hog_ctrl = MC$hog, amp_metab = MT$amp, amp_ctrl = MC$amp) }
A1 <- rbindlist(A1); fwrite(A1, file.path(hsd, "p3_parental_axis_A1_matched.tsv"), sep = "\t"); fwrite(rbindlist(A1P), file.path(hsd, "p3_parental_axis_A1_pairs.tsv"), sep = "\t")
## ---- cold confirmation ----
COLD <- rbindlist(lapply(c(30L, 60L), function(t) { D <- DV_cold[time == t]
  cbind(data.table(time = t, beta = 0, module = "metabolism"), as_row(ratio_unc(D, D$module == "metabolism", D$module != "metabolism", 0, "cold", seed = 7 + t))) }))
fwrite(COLD, file.path(hsd, "p3_parental_axis_cold.tsv"), sep = "\t")
## ---- A4 (i) leave-one-child-out of GO:0005975 ----
obo <- readLines(file.path(root, "data/03_reference/go-basic.obo")); starts <- grep("^\\[Term\\]", obo)
PAR <- rbindlist(lapply(seq_along(starts), function(k) { b <- obo[starts[k]:(if (k < length(starts)) starts[k + 1] - 1 else length(obo))]
  id <- sub("^id: ", "", grep("^id: ", b, value = TRUE)[1])
  ps <- c(sub("^is_a: (GO:[0-9]+).*", "\\1", grep("^is_a: ", b, value = TRUE)), sub("^relationship: part_of (GO:[0-9]+).*", "\\1", grep("^relationship: part_of ", b, value = TRUE)))
  if (length(ps)) data.table(child = id, parent = ps) else NULL }))
desc_of <- function(root_id) { out <- root_id; frontier <- root_id
  while (length(frontier)) { nx <- setdiff(PAR[parent %in% frontier, child], out); out <- c(out, nx); frontier <- nx }; out }
ROOT <- "GO:0005975"; KIDS <- PAR[parent == ROOT, unique(child)]; SUB <- desc_of(ROOT)
GAF <- fread(cmd = paste("zcat", shQuote(file.path(root, "data/03_reference/sgd.gaf.gz"))), sep = "\t", header = FALSE, skip = "SGD\t", quote = "", fill = TRUE)
GAF <- GAF[V9 == "P" & V5 %in% SUB & !grepl("NOT", V4)]
GAF[, sys := vapply(strsplit(V11, "|", fixed = TRUE), function(x) { y <- grep("^Y[A-P][LR][0-9]{3}[WC](-[A-Z])?$", x, value = TRUE); if (length(y)) y[1] else NA_character_ }, "")]
HB <- fread(file.path(csd, "hog_to_scer.tsv")); h2s <- unique(HB[, .(hog = gene, sys)])
MET <- MOD[module == "metabolism", hog]; MS <- h2s[hog %in% MET]
KID_DESC <- lapply(KIDS, desc_of); names(KID_DESC) <- KIDS
gene_kids <- GAF[sys %in% MS$sys, .(kids = list(unique(unlist(lapply(V5, function(g) KIDS[vapply(KID_DESC, function(d) g %in% d, TRUE)])))),
                                   direct_root = any(V5 == ROOT)), by = sys]
used_kids <- unique(unlist(gene_kids$kids))
say("A4(i): GO:0005975 has %d direct children; %d carry >= 1 module gene", length(KIDS), length(used_kids))
module_ratios <- function(DV, beta, drop = character(0), redefine = NULL) {
  D <- copy(DV); if (length(drop)) D <- D[!(module == "metabolism" & hog %in% drop)]
  if (!is.null(redefine)) { D <- D[module != "metabolism" | hog %in% redefine]
    extra <- setdiff(redefine, D$hog)   # KEGG genes outside every module: brought in only if present in DV_all
    if (length(extra)) D <- rbind(D, DV_all[hog %in% extra][, module := "metabolism"], fill = TRUE) }
  as.list(ratio_beta(D, TESTMODS, beta)) }
A4 <- list()
for (k in used_kids) { only <- gene_kids[vapply(kids, function(x) identical(sort(x), k), TRUE) & !direct_root, sys]
  drop <- unique(MS[sys %in% only, hog])
  for (t in TP) for (be in c(0, 0.657)) A4[[length(A4) + 1]] <- cbind(data.table(check = "leave_one_child_out", variant = k, n_removed = length(drop),
     time = t, beta = be), as.data.table(module_ratios(DV_heat[time == t], be, drop = drop))) }
## ---- A4 (ii) KEGG ----
KG <- unique(unlist(lapply(list.files(file.path(root, "data/03_external/kegg2026"), full.names = TRUE), function(f) sub("^sce:", "", fread(f, header = FALSE)$V2))))
KH <- unique(h2s[sys %in% KG, hog]); in_other <- MOD[hog %in% KH & module != "metabolism"]
say("A4(ii): KEGG central-carbon genes %d; one-to-one HOGs %d; of these in another module %d (kept there, by the exclusivity priority), in the GO module %d",
    length(KG), length(KH), nrow(in_other), sum(KH %in% MET))
KDEF <- setdiff(KH, in_other$hog)
for (t in TP) for (be in c(0, 0.657)) {
  D <- DV_heat[time == t]; D <- D[module != "metabolism" | hog %in% KDEF]
  ex <- setdiff(KDEF, D$hog); if (length(ex)) D <- rbind(D, DV_all[time == t & hog %in% ex][, module := "metabolism"], fill = TRUE)
  A4[[length(A4) + 1]] <- cbind(data.table(check = "KEGG", variant = "KEGG central carbon", n_removed = NA_integer_, time = t, beta = be,
    n_module = uniqueN(D[module == "metabolism", hog])), as.data.table(as.list(ratio_beta(D, TESTMODS, be)))) }
A4 <- rbindlist(A4, fill = TRUE); fwrite(A4, file.path(hsd, "p3_parental_axis_A4.tsv"), sep = "\t")
fwrite(data.table(hog = in_other$hog, module = in_other$module), file.path(hsd, "p3_parental_axis_A4_kegg_overlap.tsv"), sep = "\t")
say("\n=== 85_p3_parental_axis_tests | %s ===", format(Sys.time()))
say("\nA1 whole-transcriptome matching (beta = 0):"); print(A1[, .(time, n_metab, n_matched, pool, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3),
  se_rho = round(se_rho, 3), se_jk = round(se_jk, 3), lo_rho2 = round(lo_rho2, 3), lo_pairboot = round(lo_pairboot, 3), hi_pairboot = round(hi_pairboot, 3))])
say("\ncold confirmation:"); print(COLD[, .(time, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3), se_rho = round(se_rho, 3), se_jk = round(se_jk, 3), lo_rho2 = round(lo_rho2, 3))])
say("\nA3 family, carbohydrate metabolism:"); print(FAM[module == "metabolism", .(time, beta, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3), p = signif(p, 2), p_holm = signif(p_holm, 2), p_holm_rho2 = signif(p_holm_rho2, 2), se_rho = round(se_rho, 3), se_jk = round(se_jk, 3))])
say("\nA3 family, any cell with Holm p < 0.05:"); print(FAM[p_holm < 0.05, .(time, beta, module, ratio = round(ratio, 3), p_holm = signif(p_holm, 2))])
say("\nA2 permutation at beta 0.657:"); print(A2)
say("\nA4:"); print(A4[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
