#!/usr/bin/env Rscript
## 47_cold_arm_regulon_split.R -- the cold arm at the anchor exponent, split by regulon assignment.
##
## Background. The manuscript says the metabolic excess also appears after a shift to 12 C, which is what makes
## the finding a property of temperature shift rather than of heat. Script 28 runs the cold arm unsplit and
## only on the beta grid {0, 0.25, 0.5, 0.75, 1}, so two things are unknown: whether cold reproduces the
## localisation at the exponent the manuscript anchors on, and whether it reproduces the split, i.e. that the
## surviving excess sits in genes assigned to neither stress-factor regulon. A module-level reproduction with a
## split-level failure would narrow the interpretation.
##
## Decision rule, written before running:
##   Module level: cold metabolism is known to be resolvable at 30 min for beta <= 0.5. If that holds at the
##     anchor too, the module-level reproduction stands.
##   Split level: if the neither-regulon part is resolvable at 30 min, the reproduction holds at split level.
##     If it is not, the result is classified as module-level but not split-level reproduction. Both outcomes use
##     the same interpretation criteria.
##   The Hsf1-associated arm has n = 7 in the heat arm and will be no larger here; it is classified as
##     unresolvable rather than treated as evidence either way.
##
## Transplanted, not re-typed. Script 28 builds both arms through one code path and checks itself against
## the published heat numbers. Its `parental_lfc`, `build_DV`, module definitions and `rms_true` are taken
## unchanged; the only thing added is the regulon grouping and the anchor exponent.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

ANCHOR <- 0.657
BETAS_RUN <- c(0, 0.25, 0.50, ANCHOR, 0.75, 1.00)
NBOOT <- 1000L
SEED  <- 20260912L

SRC <- file.path(root, "scripts/hsf1/28_cold_arm_specificity.R")
EX  <- parse(SRC)
need <- c("parental_lfc", "build_DV", "rms_true", "MOD", "TP", "PAIRS", "SPECIES")
have <- function() all(vapply(need, exists, TRUE, envir = globalenv()))
n_eval <- 0L
for (e in EX) {
  s <- paste(deparse(e), collapse = " ")
  if (grepl("^sink\\(|^fwrite\\(|quit\\(|^writeLines\\(", s)) next
  if (grepl("^say\\(", s)) next
  try(eval(e, envir = globalenv()), silent = TRUE)
  n_eval <- n_eval + 1L
  if (have()) break
}
say("transplant: evaluated %d of %d expressions from %s", n_eval, length(EX), basename(SRC))
stopifnot(have(), is.function(parental_lfc), is.function(build_DV), is.function(rms_true))

DVc <- build_DV(parental_lfc("cold"))
stopifnot("module" %in% names(DVc), all(c("d_sig", "d_nul", "amp", "hog", "time") %in% names(DVc)))
say("cold divergence rows %d over %d HOGs", nrow(DVc), uniqueN(DVc$hog))

## ------------------------------------------------------------------ fidelity check ----
say("\n=== fidelity gate: reproduce the landed cold module ratios ===")
ratio_grp <- function(DT, col, mods, beta, idx = NULL) {
  d <- if (is.null(idx)) DT else DT[idx]
  sc <- pmax(d$amp, 1e-9)^beta
  ds <- d$d_sig / sc; dn <- d$d_nul / sc
  g <- d[[col]]
  vapply(mods, function(m) {
    ins <- g == m; ous <- !ins
    if (sum(ins) < 5L || sum(ous) < 5L) return(NA_real_)
    vi <- rms_true(ds[ins], dn[ins]); vo <- rms_true(ds[ous], dn[ous])
    if (vo > 0) vi / vo else NA_real_
  }, 0)
}
REF <- fread(file.path(hsd, "s4_module_divergence_beta_cold.tsv"))
MODS <- sort(unique(REF$module))
worst <- 0
for (t in TP) for (be in unique(REF$beta)) {
  obs <- ratio_grp(DVc[time == t], "module", MODS, be)
  for (m in MODS) {
    r <- REF[time == t & abs(beta - be) < 1e-9 & module == m, ratio]
    if (!length(r) || !is.finite(r) || !is.finite(obs[[m]])) next
    worst <- max(worst, abs(obs[[m]] - r))
  }
}
say("  worst absolute difference from s4_module_divergence_beta_cold.tsv: %.3g", worst)
stopifnot(worst < 1e-9)
say("  reproduced. The cold arm below is script 28's, not a reimplementation.")

## ------------------------------------------------------------------ the regulon grouping ----
M <- fread(file.path(hsd, "s4_module_membership.tsv"))
S <- fread(file.path(hsd, "a5_strata_membership_37_30.tsv"))
J <- merge(M[module == "metabolism"], S, by = "hog", all.x = TRUE)
msn <- J[!is.na(Msn24_only) & Msn24_only == TRUE, hog]
hsf <- J[!(hog %in% msn) & ((!is.na(Hsf1_core) & Hsf1_core) | (!is.na(Hsf1_joint) & Hsf1_joint)), hog]
nei <- setdiff(J$hog, c(msn, hsf))
stopifnot(length(msn) + length(hsf) + length(nei) == nrow(J))
say("\nmetabolism grouping carried over from the heat arm: Msn24 %d, Hsf1 %d, neither %d",
    length(msn), length(hsf), length(nei))

DVc[, grp3 := fifelse(module != "metabolism", module,
             fifelse(hog %in% msn, "metab_Msn24",
             fifelse(hog %in% hsf, "metab_Hsf1", "metab_neither")))]
G3 <- c("metab_Msn24", "metab_Hsf1", "metab_neither")
for (g in G3) say("  %-14s %d HOGs with cold divergence data", g, DVc[grp3 == g, uniqueN(hog)])
## The cold arm's expression filter is its own, so these counts need not match the heat arm and are asserted
## to be non-empty rather than equal.
stopifnot(all(vapply(G3, function(g) DVc[grp3 == g, uniqueN(hog)] > 0L, TRUE)))

set.seed(SEED)
RES <- list()
for (t in TP) {
  DTt <- DVc[time == t]; hs <- unique(DTt$hog); hi <- split(seq_len(nrow(DTt)), DTt$hog)
  for (be in BETAS_RUN) {
    for (col in c("module", "grp3")) {
      mods <- if (col == "module") "metabolism" else G3
      obs <- ratio_grp(DTt, col, mods, be)
      bt <- matrix(NA_real_, NBOOT, length(mods))
      for (k in seq_len(NBOOT)) {
        idx <- unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)
        bt[k, ] <- ratio_grp(DTt, col, mods, be, idx)
      }
      for (j in seq_along(mods)) RES[[length(RES) + 1L]] <- data.table(
        arm = "cold", level = col, time = t, beta = be, grp = mods[j],
        n_hog = DTt[get(col) == mods[j], uniqueN(hog)], ratio = obs[j],
        lo = quantile(bt[, j], .025, na.rm = TRUE), hi = quantile(bt[, j], .975, na.rm = TRUE))
    }
  }
  say("  cold %d min done", t)
}
R <- rbindlist(RES)
## A cell is only evidence if the point estimate itself exists. The cold arm at 15 min has a near-zero
## denominator -- Fay report an average of 11 cold-responsive genes with an allele-by-time interaction -- so
## its ratios come back NA with intervals like [6.79, 6.79]. Note: a flag that tests only that the lower bound
## is finite and above 1 marks those cells resolvable, so a finite ratio and a finite upper bound are also
## required, and degeneracy is recorded explicitly.
R[, degenerate := !is.finite(ratio) | !is.finite(lo) | !is.finite(hi) | (hi - lo) < 1e-9]
R[, resolvable := !degenerate & lo > 1]
# The column order is fixed with `resolvable` before `degenerate`, the order 48_flag_degenerate_cells.R produces.
# Without it the two trailing columns swap, and a table comparison reports a difference with a maximum numeric
# deviation of zero.
setcolorder(R, c("arm", "level", "time", "beta", "grp", "n_hog", "ratio", "lo", "hi",
                 "resolvable", "degenerate"))
fwrite(R, file.path(hsd, "s4_cold_regulon_split_beta.tsv"), sep = "\t")

say("\n=== cold arm: unsplit metabolism, and the three-way split ===")
for (be in BETAS_RUN) {
  say("  beta = %.3f", be)
  for (g in c("metabolism", G3)) {
    z <- R[grp == g & abs(beta - be) < 1e-9][order(time)]
    if (!nrow(z)) next
    say("    %-14s n=%-4d %s", g, z$n_hog[1],
        paste(sprintf("%20s", sprintf("%.3f[%.2f,%.2f]%s", z$ratio, z$lo, z$hi,
              ifelse(z$resolvable, " *", "  "))), collapse = ""))
  }
}

say("\n=== degeneracy, reported before any read-out ===")
DG <- R[, .(cells = .N, degenerate = sum(degenerate)), by = time]
for (i in seq_len(nrow(DG))) say("  %2d min: %d of %d cells degenerate", DG$time[i], DG$degenerate[i], DG$cells[i])
say("  All 15-minute cold cells are excluded from the read-out below on that basis, which is the same")
say("  exclusion already applied to the cold 15-minute timepoint.")

say("\n=== the pre-registered read-out, degenerate cells excluded ===")
for (g in c("metabolism", G3)) {
  z <- R[grp == g & abs(beta - ANCHOR) < 1e-9]
  if (!nrow(z)) next
  say("  %-14s at beta = %.3f: resolvable in %d of %d timepoints (n = %d)",
      g, ANCHOR, sum(z$resolvable, na.rm = TRUE), nrow(z), z$n_hog[1])
}
z30 <- R[grp == "metab_neither" & time == 30L & abs(beta - ANCHOR) < 1e-9]
say("\n  the deciding cell, neither-regulon at 30 min and the anchor: %.3f [%.2f, %.2f] -> %s",
    z30$ratio, z30$lo, z30$hi,
    if (isTRUE(z30$resolvable)) "reproduced at split level" else
      "module level only; the split does not reproduce in cold")
say("")
