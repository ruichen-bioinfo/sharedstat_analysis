#!/usr/bin/env Rscript
## 86_p3_parental_axis_scope.R -- scope extension of script 85's correlation-adjusted uncertainty on the parental
## axis, same uncertainty engine as 85:
## (a) Table 2's three-way split of carbohydrate metabolism (Msn2/4-associated / Hsf1 core or joint / neither) against
##     everything else, at beta 0, 0.657, 0.75, 1; (b) the P3 Results' core definitions the unanimous exclusive Hsf1-only definition
##     and the GO proteostasis definition against other induced genes, with the Results' statistic
##     (median over the six species pairs of the per-pair ratio), at beta 0, 0.383, 1. For (b) the numerator and the
##     denominator are bootstrapped separately (core genes resampled with the denominator fixed, and vice versa), each
##     SE is inflated by its set's VIF (the SE_rho form of 85), and the leave-one-species-out jackknife uses the three remaining pairs.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
EX <- parse(file.path(root, "scripts/hsf1/85_p3_parental_axis_tests.R"))
for (e in EX) { s <- paste(deparse(e), collapse = " "); if (grepl("^FAM <- ", s)) break; eval(e) }
NB <- 500L
HB <- fread(file.path(csd, "hog_to_scer.tsv")); h2s <- unique(HB[, .(hog = gene, sys)])
## ---- (a) three-way split ----
STR <- fread(file.path(hsd, "a5_strata_membership_37_30.tsv")); J <- merge(MOD, STR, by = "hog", all.x = TRUE)
msn <- J[module == "metabolism" & !is.na(Msn24_only) & Msn24_only, hog]
hsf <- J[module == "metabolism" & !(hog %in% msn) & ((!is.na(Hsf1_core) & Hsf1_core) | (!is.na(Hsf1_joint) & Hsf1_joint)), hog]
DVs <- copy(DV_heat)[, grp3 := fifelse(module != "metabolism", module, fifelse(hog %in% msn, "metab_Msn24", fifelse(hog %in% hsf, "metab_Hsf1", "metab_neither")))]
G3 <- c("metab_Msn24", "metab_Hsf1", "metab_neither"); LS <- fread(file.path(hsd, "s4_module_stratum_split3_beta.tsv"))
SPL <- rbindlist(lapply(TP, function(t) { D <- DVs[time == t]; rbindlist(lapply(c(0, 0.657, 0.75, 1), function(be) rbindlist(lapply(G3, function(g)
  cbind(data.table(time = t, beta = be, grp = g), as_row(ratio_unc(D, D$grp3 == g, D$grp3 != g, be, "heat", seed = t + 7 * which(G3 == g))))))))}))
chk <- merge(SPL[, .(time, beta, grp, ratio)], LS[, .(time, beta, grp, landed = ratio)], by = c("time", "beta", "grp"))
stopifnot(nrow(chk) >= 24, max(abs(chk$ratio - chk$landed), na.rm = TRUE) < 1e-9)
say("fidelity: %d landed split ratios reproduced", nrow(chk))
fwrite(SPL, file.path(hsd, "p3_parental_axis_split3.tsv"), sep = "\t")
## ---- (b) core definitions of the P3 Results ----
L <- unique(parental_lfc("heat")[, .(hog, species, block, time, lfc)])
W <- dcast(L, hog + species + time ~ block, value.var = "lfc"); setnames(W, c("block1", "block2"), c("b1", "b2"), skip_absent = TRUE)
W <- W[!is.na(b1) & !is.na(b2)][, sig := (b1 + b2) / 2][, nul := (b2 - b1) / 2]
S <- dcast(W, hog + time ~ species, value.var = "sig"); N <- dcast(W, hog + time ~ species, value.var = "nul")
DVU <- rbindlist(lapply(PAIRS, function(p) { a <- substr(p, 1, 2); b <- substr(p, 4, 5)
  data.table(hog = S$hog, time = S$time, pair = p, d_sig = S[[a]] - S[[b]], d_nul = N[[a]] - N[[b]]) }))
DVU <- merge(DVU, W[, .(amp = mean(abs(sig), na.rm = TRUE)), by = .(hog, time)], by = c("hog", "time"))[is.finite(d_sig) & is.finite(d_nul) & is.finite(amp)]
CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv")); UN <- CL[, .(nd = uniqueN(class), cls = class[1]), by = gene][nd == 1L]
GO <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "go_annotations_propagated_universe.tsv.gz"))))
IND <- unique(h2s[sys %in% CL[, .(k = sum(induced), n = .N), by = gene][k == n, gene], hog])
DEFS <- list(D1_unanimous_exclusive = unique(h2s[sys %in% UN[cls == "Hsf1_only", gene], hog]),
             D5_GO_proteostasis = unique(h2s[sys %in% unique(GO[go %in% MODDEF[["proteostasis"]], sys]), hog]))
pp_stat <- function(Dc, Do, beta, pairs = PAIRS, minp = 4L) {
  r <- vapply(pairs, function(p) { a <- Dc[pair == p]; o <- Do[pair == p]; if (length(unique(a$hog)) < 10 || length(unique(o$hog)) < 10) return(NA_real_)
    vo <- amp_set(o, beta); if (vo > 0) amp_set(a, beta) / vo else NA_real_ }, 0)
  if (sum(is.finite(r)) < minp) NA_real_ else median(r, na.rm = TRUE) }
bs_rows <- function(X) { hs <- unique(X$hog); hi <- split(seq_len(nrow(X)), X$hog); X[unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)] }
LR <- fread(file.path(hsd, "s4_core_definition_robustness.tsv"))
CORE <- rbindlist(lapply(names(DEFS), function(k) rbindlist(lapply(TP, function(t) rbindlist(lapply(c(0, 0.383, 1), function(be) {
  D <- DVU[time == t]; Dc <- D[hog %in% DEFS[[k]]]; Do <- D[!(hog %in% DEFS[[k]]) & hog %in% IND]
  r <- pp_stat(Dc, Do, be); set.seed(86 + t)
  bn <- vapply(seq_len(NB), function(i) pp_stat(bs_rows(Dc), Do, be), 0); bd <- vapply(seq_len(NB), function(i) pp_stat(Dc, bs_rows(Do), be), 0)
  ba <- vapply(seq_len(NB), function(i) pp_stat(bs_rows(Dc), bs_rows(Do), be), 0)
  se_n <- sd(log(bn[bn > 0]), na.rm = TRUE); se_d <- sd(log(bd[bd > 0]), na.rm = TRUE)
  pc <- rho_set("heat", unique(Dc$hog)); po <- rho_set("heat", unique(Do$hog)); nc <- uniqueN(Dc$hog); no <- uniqueN(Do$hog)
  se_rho <- sqrt((1 + (nc - 1) * pc[["rho"]]) * se_n^2 + (1 + (no - 1) * po[["rho"]]) * se_d^2)
  se_rho2 <- sqrt((1 + (nc - 1) * pc[["rho2"]]) * se_n^2 + (1 + (no - 1) * po[["rho2"]]) * se_d^2)
  jk <- vapply(SPECIES, function(s) { pr <- PAIRS[!grepl(s, PAIRS, fixed = TRUE)]; v <- pp_stat(Dc, Do, be, pr, 3L); if (is.finite(v) && v > 0) log(v) else NA_real_ }, 0)
  se_jk <- if (all(is.finite(jk))) sqrt(3 / 4 * sum((jk - mean(jk))^2)) else NA_real_
  se <- max(se_rho, se_jk, na.rm = TRUE)
  data.table(definition = k, time = t, beta = be, n_core = nc, ratio = r, lo = exp(log(r) - 1.96 * se), hi = exp(log(r) + 1.96 * se),
             se_rho = se_rho, se_jk = se_jk, se_rho2 = se_rho2, lo_rho2 = exp(log(r) - 1.96 * se_rho2), hi_rho2 = exp(log(r) + 1.96 * se_rho2),
             lo_boot = quantile(ba, .025, na.rm = TRUE), hi_boot = quantile(ba, .975, na.rm = TRUE)) }))))))
chk <- merge(CORE[beta %in% c(0, 1), .(definition, time, beta, ratio)], LR[, .(definition, time, beta, landed = ratio)], by = c("definition", "time", "beta"))
say("fidelity: reported core ratios vs landed, max |diff| %.2e over %d cells", max(abs(chk$ratio - chk$landed)), nrow(chk))
stopifnot(max(abs(chk$ratio - chk$landed)) < 1e-9)
fwrite(CORE, file.path(hsd, "p3_parental_axis_core_definitions.tsv"), sep = "\t")
say("\n(a) three-way split, adjusted intervals:"); print(SPL[, .(time, beta, grp, nA, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3), lo_boot = round(lo_boot, 3), hi_boot = round(hi_boot, 3), lo_jk = round(exp(log(ratio) - 1.96 * se_jk), 3))])
say("\n(b) reported core definitions, adjusted intervals:"); print(CORE[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])
