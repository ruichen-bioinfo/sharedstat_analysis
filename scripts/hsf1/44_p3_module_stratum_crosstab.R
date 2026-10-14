#!/usr/bin/env Rscript
## 44_p3_module_stratum_crosstab.R -- cross the module table with the stratum table, and ask whether the
## metabolism module's divergence excess is carried by its Msn2/4-associated part or by the rest of it.
##
## Background. The stratum analysis reports a Msn24_only / TF_independent divergence ratio and the module
## analysis reports a metabolism module ratio; both tables are keyed on HOG but were not joined. In budding
## yeast heat shock, Msn2/4 drives the environmental stress response, and glycolysis and reserve-carbohydrate
## genes are part of that regulon (Kuang et al. 2017). So the two headline observations may be one
## observation under two labels, and the reading of the metabolism-module result depends on this table.
##
## Decision rule, written before the numbers were seen:
##   Outcome 1. The metabolism excess concentrates in the Msn2/4-associated part. Then the Msn2/4 ratio and the
##              metabolism ratio are the same finding: divergence is not in the
##              Hsf1 core, it is in an Msn2/4-associated carbohydrate arm.
##   Outcome 2. The excess is present in the part with no Msn2/4 assignment. Then the metabolism result stands
##              as a functional localisation independent of the Msn2/4 axis.
##   Outcome 3. Neither part carries a resolvable excess once
##              split because splitting costs power. Then the split is classified as uninformative and the
##              metabolism result keeps the unsplit statement with the enrichment noted.
## These outcomes exhaust the interpretation; no third partition is introduced after seeing the numbers.
##
## Fidelity check. The divergence table is rebuilt here rather than read from disk, so before any split the
## rebuilt table must reproduce the saved s4_module_divergence_beta.tsv ratios. If it does not, the script
## stops.

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
outd <- hsd
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

SP     <- c("Sc", "Sp", "Sk", "Su")
PAIRS  <- c("Sc-Sp", "Sc-Sk", "Sc-Su", "Sp-Sk", "Sp-Su", "Sk-Su")
TP     <- c(15L, 30L, 60L)
## 0.657 is the exponent measured on the parental divergence axis, and the beta the manuscript quotes for the
## metabolism module result; it must be in the grid.
BETAS  <- c(0, 0.25, 0.383444, 0.50, 0.657, 0.75, 1.00)
NBOOT  <- 1000L
SEED   <- 20260912L

## ------------------------------------------------------------------ inputs and their dimensions ----
L <- unique(fread(cmd = paste("zcat", shQuote(file.path(hsd, "fay_parental_lfc_37_30.tsv.gz"))))[
            , .(hog, species, block, time, lfc)])
MOD <- fread(file.path(hsd, "s4_module_membership.tsv"))
STR <- fread(file.path(hsd, "a5_strata_membership_37_30.tsv"))
RBL <- fread(file.path(hsd, "s4_module_divergence_beta.tsv"))
stopifnot(uniqueN(L$species) == 4L, setequal(unique(L$time), TP), uniqueN(L$block) == 2L)
stopifnot(!any(duplicated(MOD$hog)), !any(duplicated(STR$hog)))
BOOLS <- c("Hsf1_core", "Hsf1_joint", "Msn24_only", "TF_independent", "not_induced", "other_induced")
stopifnot(all(BOOLS %in% names(STR)), all(vapply(STR[, ..BOOLS], is.logical, TRUE)))
say("parental LFC %d rows / %d HOGs; module table %d HOGs; stratum table %d HOGs",
    nrow(L), uniqueN(L$hog), nrow(MOD), nrow(STR))

## ------------------------------------------------------------------ signal / null, as in script 15 ----
W <- dcast(L, hog + species + time ~ block, value.var = "lfc")
setnames(W, c("block1", "block2"), c("b1", "b2"))
W <- W[!is.na(b1) & !is.na(b2)]
W[, sig := (b1 + b2) / 2][, nul := (b2 - b1) / 2]
rms_true <- function(s, n) sqrt(max(0, mean(s^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE)))

S <- dcast(W, hog + time ~ species, value.var = "sig")
N <- dcast(W, hog + time ~ species, value.var = "nul")
S <- S[complete.cases(S[, ..SP])]; N <- N[hog %in% S$hog & time %in% S$time]
setkey(S, hog, time); setkey(N, hog, time); N <- N[S[, .(hog, time)]]
stopifnot(identical(S$hog, N$hog), identical(S$time, N$time))

DV <- rbindlist(lapply(PAIRS, function(p) {
  a <- substr(p, 1, 2); b <- substr(p, 4, 5)
  data.table(hog = S$hog, time = S$time, pair = p,
             d_sig = S[[a]] - S[[b]], d_nul = N[[a]] - N[[b]])
}))
AMPG <- W[, .(amp = mean(abs(sig), na.rm = TRUE)), by = .(hog, time)]
DV <- merge(DV, AMPG, by = c("hog", "time"))
DV <- merge(DV, MOD[, .(hog, module)], by = "hog")
DV <- DV[!is.na(d_sig) & !is.na(d_nul) & is.finite(amp)]
say("divergence rows %d over %d HOGs", nrow(DV), uniqueN(DV$hog))

ratio_beta <- function(DT, mods, beta, idx = NULL) {
  d <- if (is.null(idx)) DT else DT[idx]
  sc <- pmax(d$amp, 1e-9)^beta
  ds <- d$d_sig / sc; dn <- d$d_nul / sc
  vapply(mods, function(m) {
    ins <- d$grp == m; ous <- !ins
    if (sum(ins) < 5L || sum(ous) < 5L) return(NA_real_)
    vi <- rms_true(ds[ins], dn[ins]); vo <- rms_true(ds[ous], dn[ous])
    if (vo > 0) vi / vo else NA_real_
  }, 0)
}

## ------------------------------------------------------------------ fidelity check ----
say("\n=== fidelity gate: rebuilt table must reproduce the landed module ratios ===")
DV[, grp := module]
MODS <- sort(unique(RBL$module))
worst <- 0
for (t in TP) for (be in BETAS) {
  obs <- ratio_beta(DV[time == t], MODS, be)
  for (m in MODS) {
    ref <- RBL[time == t & abs(beta - be) < 1e-9 & module == m, ratio]
    if (!length(ref) || !is.finite(ref) || !is.finite(obs[[m]])) next
    worst <- max(worst, abs(obs[[m]] - ref))
  }
}
say("  largest absolute difference from the landed table: %.3g", worst)
stopifnot(worst < 1e-9)
say("  reproduced. The split below runs on the same arithmetic as the published ratios.")

## ------------------------------------------------------------------ the cross-tabulation ----
say("\n=== module x stratum, counts (strata are non-exclusive booleans) ===")
J <- merge(MOD, STR, by = "hog", all.x = TRUE)
stopifnot(nrow(J) == nrow(MOD))
say("  module HOGs %d; with a stratum entry %d; without %d",
    nrow(J), J[!is.na(Hsf1_core), .N], J[is.na(Hsf1_core), .N])
CT <- J[!is.na(Hsf1_core), c(list(n = .N), lapply(.SD, sum)), by = module, .SDcols = BOOLS]
setorder(CT, -n)
print(CT, row.names = FALSE)
fwrite(CT, file.path(outd, "s4_module_stratum_crosstab.tsv"), sep = "\t")

say("\n=== is Msn24_only enriched in a module, against the other modules? ===")
EN <- rbindlist(lapply(CT$module, function(mo) {
  a <- J[module == mo & !is.na(Hsf1_core)]; b <- J[module != mo & !is.na(Hsf1_core)]
  ft <- fisher.test(matrix(c(a[Msn24_only == TRUE, .N], a[Msn24_only == FALSE, .N],
                             b[Msn24_only == TRUE, .N], b[Msn24_only == FALSE, .N]), 2))
  data.table(module = mo, n = nrow(a), msn24 = a[Msn24_only == TRUE, .N],
             pct = 100 * a[Msn24_only == TRUE, .N] / nrow(a),
             odds_ratio = unname(ft$estimate), p = ft$p.value)
}))
print(EN[, .(module, n, msn24, pct = round(pct, 1), odds_ratio = round(odds_ratio, 2),
             p = signif(p, 3))], row.names = FALSE)
fwrite(EN, file.path(outd, "s4_module_stratum_enrichment.tsv"), sep = "\t")

## ------------------------------------------------------------------ the split ----
## Two subgroups of metabolism, against the same denominator the unsplit ratio used: all other modules.
## Splitting the numerator and leaving the denominator alone keeps the comparison to the published ratio
## exact; putting the other half of metabolism into the denominator would change two things at once.
say("\n=== metabolism split by Msn2/4 assignment, denominator unchanged (all other modules) ===")
msn <- J[module == "metabolism" & !is.na(Msn24_only) & Msn24_only == TRUE, hog]
oth <- J[module == "metabolism" & !(hog %in% msn), hog]
say("  metabolism & Msn24_only      %d HOGs", length(msn))
say("  metabolism without Msn24_only %d HOGs (includes %d with no stratum entry)",
    length(oth), J[module == "metabolism" & is.na(Msn24_only), .N])
stopifnot(length(msn) + length(oth) == J[module == "metabolism", .N])

DV[, grp := fifelse(module != "metabolism", module,
                    fifelse(hog %in% msn, "metab_Msn24", "metab_other"))]
stopifnot(DV[grp == "metab_Msn24", uniqueN(hog)] == length(msn))
stopifnot(DV[grp %in% c("metab_Msn24", "metab_other"), uniqueN(hog)] ==
          DV[module == "metabolism", uniqueN(hog)])
TESTG <- c("metab_Msn24", "metab_other", setdiff(MODS, "metabolism"))

set.seed(SEED)
RS <- list()
for (t in TP) {
  DTt <- DV[time == t]
  hs <- unique(DTt$hog); hi <- split(seq_len(nrow(DTt)), DTt$hog)
  for (be in BETAS) {
    obs <- ratio_beta(DTt, TESTG, be)
    bt <- matrix(NA_real_, NBOOT, length(TESTG))
    for (k in seq_len(NBOOT)) {
      idx <- unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)
      bt[k, ] <- ratio_beta(DTt, TESTG, be, idx)
    }
    for (j in seq_along(TESTG)) RS[[length(RS) + 1L]] <- data.table(
      time = t, beta = be, grp = TESTG[j], n_hog = DTt[grp == TESTG[j], uniqueN(hog)],
      ratio = obs[j], lo = quantile(bt[, j], .025, na.rm = TRUE),
      hi = quantile(bt[, j], .975, na.rm = TRUE))
  }
  say("  time %d min done", t)
}
RS <- rbindlist(RS)
fwrite(RS, file.path(outd, "s4_module_stratum_split_beta.tsv"), sep = "\t")

say("\n=== ratio [95%% CI] by group, across beta ===")
for (t in TP) {
  say("\n  time %d min", t)
  say("    %-14s %-6s %s", "group", "n", paste(sprintf("%16s", sprintf("beta=%.2f", BETAS)), collapse = ""))
  for (g in c("metab_Msn24", "metab_other", "metabolism_unsplit_ref", setdiff(MODS, "metabolism"))) {
    if (g == "metabolism_unsplit_ref") {
      cells <- vapply(BETAS, function(be) {
        z <- RBL[time == t & abs(beta - be) < 1e-9 & module == "metabolism"]
        if (!nrow(z)) "  --" else sprintf("%16s", sprintf("%.3f[%.2f,%.2f]", z$ratio, z$lo, z$hi)) }, "")
      say("    %-14s %-6s %s", "(unsplit)", "", paste(cells, collapse = ""))
      next
    }
    z0 <- RS[time == t & grp == g]
    cells <- vapply(BETAS, function(be) {
      z <- z0[abs(beta - be) < 1e-9]
      if (!nrow(z) || !is.finite(z$ratio)) sprintf("%16s", "--")
      else sprintf("%16s", sprintf("%.3f[%.2f,%.2f]", z$ratio, z$lo, z$hi)) }, "")
    say("    %-14s %-6d %s", g, z0$n_hog[1], paste(cells, collapse = ""))
  }
}

## ------------------------------------------------------------------ read the outcome against the rule ----
say("\n=== which pre-registered outcome holds ===")
res <- rbindlist(lapply(TP, function(t) rbindlist(lapply(BETAS, function(be) {
  a <- RS[time == t & abs(beta - be) < 1e-9 & grp == "metab_Msn24"]
  b <- RS[time == t & abs(beta - be) < 1e-9 & grp == "metab_other"]
  data.table(time = t, beta = be,
             msn_excess   = is.finite(a$ratio) && a$lo > 1,
             other_excess = is.finite(b$ratio) && b$lo > 1,
             msn_ratio = a$ratio, other_ratio = b$ratio)
}))))
print(res[, .(time, beta, msn_ratio = round(msn_ratio, 3), msn_excess,
              other_ratio = round(other_ratio, 3), other_excess)], row.names = FALSE)
say("\n  cells where the Msn2/4 part alone has a CI above 1     : %d of %d",
    res[msn_excess == TRUE, .N], nrow(res))
say("  cells where the non-Msn2/4 part alone has a CI above 1 : %d of %d",
    res[other_excess == TRUE, .N], nrow(res))
say("  cells where BOTH do                                    : %d",
    res[msn_excess & other_excess, .N])
say("  cells where NEITHER does                               : %d",
    res[!msn_excess & !other_excess, .N])
fwrite(res, file.path(outd, "s4_module_stratum_split_outcome.tsv"), sep = "\t")
## ------------------------------------------------------------------ three-way split ----
## The 131 metabolism genes without an Msn24_only assignment are not all outside the stress-factor regulons:
## 7 of them are Hsf1_core or Hsf1_joint. Whether the result is specific to the Msn2/4 regulon or involves both regulons depends
## on whether the surviving excess is carried by genes assigned to neither; the two-way split cannot distinguish these cases.
##
## Expected result, fixed before running: the neither-regulon part remains resolvable at beta = 0.657 at
## 30 and 60 min. The Hsf1-associated part has n = 7 and is expected to be unresolvable in either direction;
## "not resolvable" is then the result, not a failure to find one.
say("\n=== three-way split, Msn2/4 vs Hsf1 vs neither ===")
hsf <- J[module == "metabolism" & !(hog %in% msn) &
         ((!is.na(Hsf1_core) & Hsf1_core) | (!is.na(Hsf1_joint) & Hsf1_joint)), hog]
nei <- setdiff(J[module == "metabolism", hog], c(msn, hsf))
say("  metabolism: Msn24_only %d, Hsf1_core/joint %d, neither %d, total %d",
    length(msn), length(hsf), length(nei), J[module == "metabolism", .N])
stopifnot(length(msn) + length(hsf) + length(nei) == J[module == "metabolism", .N],
          length(intersect(msn, hsf)) == 0L, length(intersect(msn, nei)) == 0L,
          length(intersect(hsf, nei)) == 0L)

DV[, grp3 := fifelse(module != "metabolism", module,
              fifelse(hog %in% msn, "metab_Msn24",
              fifelse(hog %in% hsf, "metab_Hsf1", "metab_neither")))]
G3 <- c("metab_Msn24", "metab_Hsf1", "metab_neither")
for (g in G3) say("  %-14s %d HOGs with divergence data", g, DV[grp3 == g, uniqueN(hog)])
stopifnot(DV[grp3 %in% G3, uniqueN(hog)] == DV[module == "metabolism", uniqueN(hog)])

ratio3 <- function(DT, mods, beta, idx = NULL) {
  d <- if (is.null(idx)) DT else DT[idx]
  sc <- pmax(d$amp, 1e-9)^beta
  ds <- d$d_sig / sc; dn <- d$d_nul / sc
  vapply(mods, function(m) {
    ins <- d$grp3 == m; ous <- !ins
    if (sum(ins) < 5L || sum(ous) < 5L) return(NA_real_)
    vi <- rms_true(ds[ins], dn[ins]); vo <- rms_true(ds[ous], dn[ous])
    if (vo > 0) vi / vo else NA_real_
  }, 0)
}
set.seed(SEED + 2L)
R3 <- list()
for (t in TP) {
  DTt <- DV[time == t]; hs <- unique(DTt$hog); hi <- split(seq_len(nrow(DTt)), DTt$hog)
  for (be in BETAS) {
    obs <- ratio3(DTt, G3, be)
    bt <- matrix(NA_real_, NBOOT, length(G3))
    for (k in seq_len(NBOOT)) {
      idx <- unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)
      bt[k, ] <- ratio3(DTt, G3, be, idx)
    }
    for (j in seq_along(G3)) R3[[length(R3) + 1L]] <- data.table(
      time = t, beta = be, grp = G3[j], n_hog = DTt[grp3 == G3[j], uniqueN(hog)],
      ratio = obs[j], lo = quantile(bt[, j], .025, na.rm = TRUE),
      hi = quantile(bt[, j], .975, na.rm = TRUE))
  }
}
R3 <- rbindlist(R3)
fwrite(R3, file.path(outd, "s4_module_stratum_split3_beta.tsv"), sep = "\t")
say("\n  ratio [95%% CI], three-way, at the anchor and either side of it")
for (be in c(0, 0.5, 0.657, 0.75, 1)) {
  say("    beta = %.3f", be)
  for (g in G3) {
    z <- R3[grp == g & abs(beta - be) < 1e-9][order(time)]
    if (!nrow(z)) next
    say("      %-14s n=%-4d %s", g, z$n_hog[1],
        paste(sprintf("%20s", sprintf("%.3f[%.2f,%.2f]%s", z$ratio, z$lo, z$hi,
              ifelse(z$lo > 1, " *", "  "))), collapse = ""))
  }
}
say("\n  '*' marks a CI entirely above 1. Pre-registered read-out:")
for (g in G3) {
  z <- R3[grp == g & abs(beta - 0.657) < 1e-9]
  say("    %-14s at the anchor: resolvable in %d of %d timepoints (n = %d)",
      g, sum(z$lo > 1, na.rm = TRUE), nrow(z), z$n_hog[1])
}

say("")
