#!/usr/bin/env Rscript
## 87a_p3_core_adjusted.R -- hybrid (cis) axis, with the correlation-adjusted uncertainty of script 85 and the reading of
## the result both fixed before this analysis was run.
## The six core definitions of the P3 Results (incl.
## reported_core = the equivalence-bound quantity, and D5 = the Gene Ontology core) and the Figure 3 strata (each stratum
## against its denominator in DENOM), at
## beta 0, 0.383444 and 1. Statistic: script 04/45's -- median over the six hybrids of the per-cell ratio of
## noise-corrected cis divergence, per timepoint.
## Primary interval: separate gene bootstraps of numerator and denominator; SE on log scale
##   SE^2 = VIF_num SE_num^2 + VIF_den SE_den^2, VIF = 1 + (n - 1) rho-bar, rho-bar on per-library allele differences.
## No design-level jackknife was specified for this axis (leave-one-block is not computable: the noise floor needs
## both blocks), so the primary interval is the rho-bar form. Sensitivities beside it (not verdicts): leave-one-hybrid-out
## jackknife of the log statistic with a t(5) multiplier; the mean-rho^2 form; the joint gene bootstrap (as published).
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
EX <- parse(file.path(root, "scripts/hsf1/45_core_definition_at_headline_beta.R"))
for (e in EX) { s <- paste(deparse(e), collapse = " "); if (grepl("^RES <- list\\(\\)", s)) break
  if (grepl("^sink\\(|^fwrite\\(|quit\\(", s)) next; eval(e, envir = globalenv()) }
stopifnot(exists("DEFS"), exists("cell_ratio_beta"), exists("S_M"), exists("tp_of"))
say("script 45 validation checks passed (beta = 0 strict table; equivalence-bound result at beta = 0.383444)")
source(file.path(root, "scripts/hsf1/87_allele_rho.R"))
NB <- 1000L; T5 <- qt(0.975, 5)
hyb_of <- sub("_[0-9]+$", "", CELLS); HYBS <- unique(hyb_of); stopifnot(length(HYBS) == 6L)
stat <- function(i_n, i_d, beta, keep = rep(TRUE, length(CELLS))) {
  r <- cell_ratio_beta(i_n, i_d, beta); vapply(TP, function(t) median(r[tp_of == t & keep], na.rm = TRUE), 0) }
engine <- function(label, what, fo, de, beta, seed) {
  o <- stat(fo, de, beta); set.seed(seed)
  bn <- t(replicate(NB, stat(sample(fo, length(fo), TRUE), de, beta)))
  bd <- t(replicate(NB, stat(fo, sample(de, length(de), TRUE), beta)))
  bj <- t(replicate(NB, stat(sample(fo, length(fo), TRUE), sample(de, length(de), TRUE), beta)))
  pn <- allele_rho(universe[fo]); pd <- allele_rho(universe[de])
  jk <- sapply(HYBS, function(h) log(stat(fo, de, beta, hyb_of != h)))            # TP x 6
  rbindlist(lapply(seq_along(TP), function(j) {
    sn <- sd(log(bn[bn[, j] > 0, j]), na.rm = TRUE); sd_ <- sd(log(bd[bd[, j] > 0, j]), na.rm = TRUE)
    se <- sqrt((1 + (length(fo) - 1) * pn[["rho"]]) * sn^2 + (1 + (length(de) - 1) * pd[["rho"]]) * sd_^2)
    se2 <- sqrt((1 + (length(fo) - 1) * pn[["rho2"]]) * sn^2 + (1 + (length(de) - 1) * pd[["rho2"]]) * sd_^2)
    jv <- jk[j, ]; sej <- if (all(is.finite(jv))) sqrt(5 / 6 * sum((jv - mean(jv))^2)) else NA_real_
    r <- o[j]; lr <- log(r)
    data.table(set = label, what = what, time = TP[j], beta = beta, n_num = length(fo), n_den = length(de), ratio = r,
      lo = exp(lr - 1.96 * se), hi = exp(lr + 1.96 * se), se_rho = se, rho_num = pn[["rho"]], rho_den = pd[["rho"]],
      lo_rho2 = exp(lr - 1.96 * se2), hi_rho2 = exp(lr + 1.96 * se2),
      se_jk_hyb = sej, lo_jk_t5 = exp(lr - T5 * sej), hi_jk_t5 = exp(lr + T5 * sej),
      lo_boot = quantile(bj[, j], .025, na.rm = TRUE), hi_boot = quantile(bj[, j], .975, na.rm = TRUE))
  }))
}
OUT <- list()
for (k in names(DEFS)) { fo <- to_idx(DEFS[[k]]); de <- to_idx(setdiff(i_den_base, DEFS[[k]]))
  if (length(fo) < 20L || length(de) < 20L) next
  for (be in BETAS_REPORT) OUT[[length(OUT) + 1]] <- engine(k, "reported core definition vs other induced", fo, de, be, 45L) }
for (s in names(DENOM)) { fo <- idx_of(s); de <- idx_of(DENOM[[s]]); fo <- fo[!is.na(fo)]; de <- de[!is.na(de)]
  for (be in c(0, 0.383444)) OUT[[length(OUT) + 1]] <- engine(s, paste("Figure 3 stratum vs", DENOM[[s]]), fo, de, be, 46L) }
R <- rbindlist(OUT)
## fidelity: point estimates equal the stored tables
LB <- fread(file.path(root, "data/02_processed/hsf1/a5_core_definition_by_beta.tsv"))
ck <- merge(R[what %like% "reported core", .(definition = set, time, beta = round(beta, 4), ratio)], LB[, .(definition, time, beta = round(beta, 4), landed = ratio)], by = c("definition", "time", "beta"))
stopifnot(nrow(ck) >= 45L, max(abs(ck$ratio - ck$landed)) < 1e-9)
say("fidelity: %d landed reported core-definition ratios reproduced", nrow(ck))
fwrite(R, file.path(root, "data/02_processed/hsf1/p3_adjusted_core.tsv"), sep = "\t")
## the planned reading: reported_core, 30 min, widest adjusted upper bound across the three exponents
Z <- R[set == "reported_core" & time == 30]; UB <- max(Z$hi)
band <- if (UB <= 1.25) "INFORMATIVE" else if (UB < 1.5) "LIMITED" else "UNINFORMATIVE"
say("\nINTERPRETATION RULE: Hsf1 core vs other induced, 30 min, widest adjusted upper bound over beta = %.3f -> %s", UB, band)
say("  (sensitivities: widest upper bound, hybrid jackknife t(5) %.3f; mean-rho^2 %.3f; joint gene bootstrap %.3f)",
    max(Z$hi_jk_t5), max(Z$hi_rho2), max(Z$hi_boot))
G <- R[set == "D5_GO_proteostasis"]
say("GO core 'resolvably less divergent' kept only if an adjusted interval lies below 1: cells with hi < 1: %d of %d", G[hi < 1, .N], nrow(G))
print(R[, .(set, time, beta = round(beta, 3), n_num, ratio = round(ratio, 3), lo = round(lo, 3), hi = round(hi, 3), rho_num = round(rho_num, 4),
            jk_t5 = sprintf("[%.3f, %.3f]", lo_jk_t5, hi_jk_t5), rho2 = sprintf("[%.3f, %.3f]", lo_rho2, hi_rho2), boot = sprintf("[%.3f, %.3f]", lo_boot, hi_boot))])
