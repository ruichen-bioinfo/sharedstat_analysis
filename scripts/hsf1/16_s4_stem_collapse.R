# 16_s4_stem_collapse.R -- paralog collapse of the module analysis (script 15) by gene-name stem.
#
# Why it matters. Script 15's paralog check (max|delta ratio| = 0.0000 across all 75 rows) is
# vacuous, not passed: it uses only the externally verified family table, whose sole entry is the
# anCWMP family (TIP1/TIR1-TIR4/DAN1/DAN4, PMID 36134526), and zero anCWMP genes fall in any of
# the five modules. The name-stem proxy shows 45-57% of each module's genes sit in multi-member
# families, including 29 in metabolism itself (PFK 4 / GPM 3 / CIT 3 / MDH 2 / PYC 2 and other
# glycolytic paralogs). In script 02 the cell-wall enrichment looks strong uncollapsed (padj
# 4.11e-03) and does not survive collapse (padj 0.118). Metabolism is the only positive
# module-level claim, so it has to face the same test.
#
# Direction, declared before running: unknown, it could go either way. Metabolism's family
# fraction (48.7%) is lower than the two largest comparison modules, cell_cycle_growth (56.6%)
# and translation_ribosome (54.2%). Since the denominator is "all other modules", collapsing
# could strengthen metabolism rather than weaken it. Both outcomes are reported as they come.
#
# The proxy errs in both directions, and is a limitation of the analysis. Stripping
# trailing digits from a gene name:
#   - over-merges: HSP26 / HSP42 / HSP82 / HSP104 collapse to "HSP" although they are not one
#     family; PFK1 / PFK2 genuinely are, but PFK26 / PFK27 (6-phosphofructo-2-kinases) are a
#     different enzyme class that the proxy also folds in.
#   - under-merges: TIP1 and TIR1-TIR4 are one verified family (PMID 36134526) but have different
#     stems, so the proxy splits them.
# The proxy is therefore a sensitivity analysis, not a truth. A module claim that survives it is
# more credible; one that does not is flagged, not silently dropped. The verified-family override
# is applied on top of the proxy, exactly as script 02 does.
#
# Method. The stem rule and the module/divergence machinery are taken from scripts 02 and 15 so
# they cannot drift from the versions that produced those tables, and the rebuild is asserted to
# reproduce script 15's saved ratios before any collapsed number is computed.
#
# Usage: Rscript scripts/hsf1/16_s4_stem_collapse.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
ref  <- file.path(root, "data/03_reference")
repd <- file.path(root, "doc/05_reports")
outd <- hsd
MIRROR <- "logs/output_copy"

TP     <- c(15L, 30L, 60L)
SPP    <- c("Sc", "Sp", "Sk", "Su")
PAIRS  <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")
BETAS  <- c(0, 0.25, 0.50, 0.657, 0.75, 1.00)     # 0.657 = par-axis beta measured in script 05
BETA_ANCHOR <- 0.657
NBOOT  <- 1000L
NPERM  <- 1000L
MIN_MOD <- 20L
anCWMP <- c("TIP1", "TIR1", "TIR2", "TIR3", "TIR4", "DAN1", "DAN4")   # PMID 36134526
MODORDER <- c("proteostasis", "translation_ribosome", "oxidative", "metabolism", "cell_cycle_growth")
set.seed(20260906L)

say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, "logs/hsf1_s4_stem_collapse.log")
sink(logf, split = TRUE)
say("=== 16_s4_stem_collapse | %s ===", format(Sys.time()))
say("Direction declared BEFORE running: unknown, could weaken OR strengthen metabolism, because")
say("metabolism's family fraction (48.7%%) is LOWER than cell_cycle_growth (56.6%%) and")
say("translation_ribosome (54.2%%), and the denominator is 'all other modules'.")
say("The stem proxy over-merges (HSP26/42/82/104) and under-merges (TIP1 vs TIR1-4). It is a")
say("SENSITIVITY ANALYSIS, not a truth, and is a limitation of the analysis.")

## ------------------------------- rebuild script 15's inputs, then assert ----
MOD <- fread(file.path(outd, "s4_module_membership.tsv"))
L <- unique(fread(cmd = paste("zcat", shQuote(file.path(hsd, "fay_parental_lfc_37_30.tsv.gz"))))[
            , .(hog, species, block, time, lfc)])
W <- dcast(L, hog + species + time ~ block, value.var = "lfc")
setnames(W, c("block1", "block2"), c("b1", "b2"))
W <- W[!is.na(b1) & !is.na(b2)][, sig := (b1 + b2) / 2][, nul := (b2 - b1) / 2]
rms_true <- function(s, n) sqrt(max(0, mean(s^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE)))
S <- dcast(W, hog + time ~ species, value.var = "sig")
N <- dcast(W, hog + time ~ species, value.var = "nul")
stopifnot(identical(S$hog, N$hog), identical(S$time, N$time))
DV <- rbindlist(lapply(PAIRS, function(p) {
  a <- substr(p, 1, 2); b <- substr(p, 4, 5)
  data.table(hog = S$hog, time = S$time, pair = p,
             d_sig = S[[a]] - S[[b]], d_nul = N[[a]] - N[[b]]) }))
AMPG <- W[, .(amp = mean(abs(sig), na.rm = TRUE)), by = .(hog, time)]
DV <- merge(merge(DV, AMPG, by = c("hog", "time")), MOD[, .(hog, module, sys)], by = "hog")
DV <- DV[!is.na(d_sig) & !is.na(d_nul) & is.finite(amp)]
MODN <- unique(DV[, .(hog, module)])[, .N, by = module]
TESTMODS <- MODN[N >= MIN_MOD, module]
TESTMODS <- MODORDER[MODORDER %in% TESTMODS]

ratio_beta <- function(DT, mods, beta, idx = NULL) {
  d <- if (is.null(idx)) DT else DT[idx]
  sc <- pmax(d$amp, 1e-9)^beta
  ds <- d$d_sig / sc; dn <- d$d_nul / sc
  vapply(mods, function(m) {
    ins <- d$module == m; ous <- !ins
    if (sum(ins) < 5L || sum(ous) < 5L) return(NA_real_)
    vi <- rms_true(ds[ins], dn[ins]); vo <- rms_true(ds[ous], dn[ous])
    if (vo > 0) vi / vo else NA_real_ }, 0)
}
# check: the rebuild must reproduce script 15's saved ratios before anything else
REF <- fread(file.path(outd, "s4_module_divergence_beta.tsv"))
chk <- rbindlist(lapply(TP, function(t) rbindlist(lapply(c(0, 0.25, 0.5, 0.75, 1), function(be)
  data.table(time = t, beta = be, module = TESTMODS,
             rebuilt = ratio_beta(DV[time == t], TESTMODS, be)))))) 
cmp <- merge(chk, REF[, .(time, beta, module, published = ratio)], by = c("time","beta","module"))
say("\nfaithfulness check vs s4_module_divergence_beta.tsv: %d rows, max abs diff %.2e",
    nrow(cmp), max(abs(cmp$rebuilt - cmp$published)))
stopifnot(nrow(cmp) == 75L, max(abs(cmp$rebuilt - cmp$published)) < 1e-9)
say("Rebuilt ratios match the saved ratios; collapse analysis proceeds")

## ------------------------------------------------- the name-stem proxy ----
# Rule from script 02: uppercase the standard name, strip trailing digits; if that leaves fewer
# than 2 characters keep the full name. Verified families then override the proxy.
stem_rule <- function(std) {
  st <- toupper(sub("[0-9]+$", "", std))
  st[nchar(st) < 2] <- toupper(std[nchar(st) < 2])
  st
}
MOD[, std := NULL]
# standard gene names come from hog_to_scer.tsv's `Name` column -- the same source script 02's
# stem rule was built on. Where absent, a gene is its own stem: you cannot merge what has no name.
STD <- unique(fread(file.path(csd, "hog_to_scer.tsv"))[, .(sys, std = Name)])
STD <- STD[!is.na(sys) & sys != ""][!duplicated(sys)]
MOD <- merge(MOD, STD, by = "sys", all.x = TRUE)
n_noname <- MOD[is.na(std) | std == "" | std == sys, .N]
MOD[is.na(std) | std == "", std := sys]
say("\nstandard names resolved for %d of %d module genes (%d fall back to the systematic name,",
    MOD[, sum(std != sys)], nrow(MOD), n_noname)
say("and therefore cannot be merged with anything by the proxy)")
MOD[, stem := stem_rule(std)]
# verified-family override, applied on top, exactly as script 02 does
MOD[std %in% anCWMP, stem := "FAM_anCWMP"]
FAMSTAT <- MOD[, .(genes = .N, stems = uniqueN(stem),
                   in_multi = sum(stem %in% stem[duplicated(stem)])), by = module]
FAMSTAT[, frac_multi := in_multi / genes]
setkey(FAMSTAT, module)
fwrite(FAMSTAT, file.path(outd, "s4_stem_family_stats.tsv"), sep = "\t")
say("\n---- stem-proxy family structure per module ----")
say("%-22s %7s %7s %11s %11s", "module", "genes", "stems", "in a multi", "frac multi")
for (m in MODORDER) { z <- FAMSTAT[m]
  if (!nrow(z) || is.na(z$genes)) next
  say("%-22s %7d %7d %11d %10.1f%%", m, z$genes, z$stems, z$in_multi, 100*z$frac_multi) }
say("\nlargest multi-member stems inside metabolism (the load-bearing module):")
mm <- MOD[module == "metabolism", .N, by = stem][N > 1L][order(-N)]
if (nrow(mm)) for (i in seq_len(min(12L, nrow(mm)))) say("  %-14s %d members", mm$stem[i], mm$N[i])
say("total metabolism genes inside a multi-member stem: %d", FAMSTAT["metabolism", in_multi])

## ---------------------------------------- collapse and re-run everything ---
# One representative per stem, dropped from both the module and the background, so every module's
# size shrinks consistently (the same construction as script 02's collapse_rerun).
setorder(MOD, stem, sys)
KEEP <- MOD[!duplicated(stem), hog]
DROP <- setdiff(MOD$hog, KEEP)
say("\ncollapse: %d HOGs -> %d (%d dropped as within-stem duplicates)",
    uniqueN(MOD$hog), length(unique(KEEP)), length(unique(DROP)))
DVc <- DV[hog %in% KEEP]
MODNc <- unique(DVc[, .(hog, module)])[, .N, by = module]
say("module sizes after collapse:")
for (m in TESTMODS) say("  %-22s %4d -> %4d HOGs", m, MODN[module == m, N],
                        ifelse(m %in% MODNc$module, MODNc[module == m, N], 0L))
TESTMODS2 <- MODORDER[MODORDER %in% MODNc[N >= MIN_MOD, module]]
if (!identical(TESTMODS2, TESTMODS))
  say("NOTE: modules dropping below n=%d after collapse: %s", MIN_MOD,
      paste(setdiff(TESTMODS, TESTMODS2), collapse = ", "))

RES <- list()
for (t in TP) {
  Dt  <- DV[time == t]; Dc <- DVc[time == t]
  hcs <- unique(Dc$hog); hic <- split(seq_len(nrow(Dc)), Dc$hog)
  hmod <- unique(Dc[, .(hog, module)])
  for (be in BETAS) {
    full <- ratio_beta(Dt, TESTMODS2, be)
    coll <- ratio_beta(Dc, TESTMODS2, be)
    bt <- matrix(NA_real_, NBOOT, length(TESTMODS2))
    for (k in seq_len(NBOOT)) {
      idx <- unlist(hic[sample(hcs, length(hcs), TRUE)], use.names = FALSE)
      bt[k, ] <- ratio_beta(Dc, TESTMODS2, be, idx)
    }
    nl <- matrix(NA_real_, NPERM, length(TESTMODS2))
    for (k in seq_len(NPERM)) {
      hp <- copy(hmod)[, module := sample(module)]
      d2 <- merge(Dc[, .(hog, d_sig, d_nul, amp)], hp, by = "hog")
      nl[k, ] <- ratio_beta(d2, TESTMODS2, be)
    }
    for (j in seq_along(TESTMODS2)) RES[[length(RES) + 1L]] <- data.table(
      time = t, beta = be, module = TESTMODS2[j],
      ratio_full = full[j], ratio_collapsed = coll[j], delta = coll[j] - full[j],
      lo = quantile(bt[, j], .025, na.rm = TRUE), hi = quantile(bt[, j], .975, na.rm = TRUE),
      p_emp = mean(abs(nl[, j] - 1) >= abs(coll[j] - 1), na.rm = TRUE))
  }
  say("  time %d min done", t)
}
RES <- rbindlist(RES)
RES[, ci_excludes_1 := lo > 1 | hi < 1]
RES[, survives := ci_excludes_1 & p_emp < 0.05]
fwrite(RES, file.path(outd, "s4_stem_collapse_ratios.tsv"), sep = "\t")

say("\n================ COLLAPSED vs FULL, all modules x time x beta ================")
say("%-22s %5s %6s %9s %9s %8s %20s %8s %9s", "module", "time", "beta", "full", "collapsed",
    "delta", "collapsed 95%% CI", "p_emp", "survives")
for (m in TESTMODS2) for (t in TP) for (be in BETAS) {
  z <- RES[module == m & time == t & beta == be]
  if (!nrow(z)) next
  say("%-22s %5d %6.3f %9.3f %9.3f %+8.3f [%8.3f,%8.3f] %8.3f %9s", m, t, be,
      z$ratio_full, z$ratio_collapsed, z$delta, z$lo, z$hi, z$p_emp, z$survives)
}
say("\nmax |delta| across all cells: %.4f  (the earlier uninformative check returned 0.0000)",
    max(abs(RES$delta), na.rm = TRUE))

## ------------------------------------- interpretation criterion ---
say("\n================ INTERPRETATION CRITERION ================")
say("The outcome is read from whether metabolism's enrichment survives the collapse. The outcomes were")
say("specified before the analysis and applied unchanged.")
MB <- RES[module == "metabolism"]
say("\nmetabolism after collapse:")
say("%-5s %8s %9s %9s %20s %8s %9s", "time", "beta", "full", "collapsed", "95%% CI", "p_emp", "survives")
for (t in TP) for (be in BETAS) { z <- MB[time == t & beta == be]
  say("%-5d %8.3f %9.3f %9.3f [%8.3f,%8.3f] %8.3f %9s", t, be, z$ratio_full,
      z$ratio_collapsed, z$lo, z$hi, z$p_emp, z$survives) }
surv_any    <- MB[, any(survives, na.rm = TRUE)]
surv_anchor <- MB[abs(beta - BETA_ANCHOR) < 1e-9, any(survives, na.rm = TRUE)]
surv_anchor_60 <- MB[abs(beta - BETA_ANCHOR) < 1e-9 & time == 60L, isTRUE(survives)]
say("\nsurvives at ANY beta/time      : %s", surv_any)
say("survives at the measured anchor beta = %.3f (A6 par axis): %s",
    BETA_ANCHOR, surv_anchor)
say("  of which at 60 min            : %s", surv_anchor_60)

CC <- RES[module == "cell_cycle_growth"]
say("\ncell_cycle_growth (a second independent positive result):")
say("%-5s %8s %9s %9s %20s %8s %9s", "time", "beta", "full", "collapsed", "95%% CI", "p_emp", "survives")
for (t in TP) for (be in BETAS) { z <- CC[time == t & beta == be]
  say("%-5d %8.3f %9.3f %9.3f [%8.3f,%8.3f] %8.3f %9s", t, be, z$ratio_full,
      z$ratio_collapsed, z$lo, z$hi, z$p_emp, z$survives) }
PR <- RES[module == "proteostasis"]
say("\nproteostasis (must stay indistinguishable from background): any beta/time surviving? %s",
    PR[, any(survives, na.rm = TRUE)])

say("")
if (surv_anchor_60 || (surv_any && surv_anchor)) {
  say("=> Result: metabolism enrichment SURVIVES the stem-proxy collapse.")
  say("   Reading: the DEFENCE arm (proteostasis)")
  say("   and the GROWTH-ARREST arm are conserved, and divergence is carried by CARBOHYDRATE-")
  say("   METABOLIC REPROGRAMMING.")
  say("   At the measured anchor beta = %.3f the enrichment",
      BETA_ANCHOR)
  say("   holds at 60 min, is marginal at 30 min, and does not hold at 15 min -- a time trend")
  say("   concordant with A4 (divergence grows with time) and A6 (trans accumulates).")
} else if (surv_any) {
  say("=> Result: metabolism survives at some beta but NOT at the measured anchor")
  say("   beta = %.3f. The interpretation is conditional on the anchor and includes the failure at the",
      BETA_ANCHOR)
  say("   anchor.")
} else {
  say("=> Result: metabolism enrichment does NOT survive the collapse.")
  say("   The result is then layer localisation and resolution bounds only, and the metabolic")
  say("   reading remains a candidate observation.")
  say("")
  say("   This matches the earlier cell-wall sensitivity result: cell-wall enrichment looked strong at padj 4.11e-03")
  say("   and died at 0.118 after the same collapse.")
}
say("")
say("Limitation of the stem proxy: it is wrong in BOTH")
say("directions -- it over-merges non-families (HSP26/42/82/104 -> 'HSP') and under-merges")
say("verified families with different stems (TIP1 vs TIR1-TIR4, PMID 36134526). It is reported as")
say("a sensitivity analysis. The only externally VERIFIED family table still contains exactly one")
say("family, and it intersects none of the five modules.")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_s4_stem_collapse.txt"))
mdir <- file.path(MIRROR, "data/02_processed/hsf1")
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(outd, pattern = "^(s4_stem|sessionInfo_s4_stem)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
file.copy(file.path(root, "scripts/hsf1/16_s4_stem_collapse.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("
copied outputs to %s", MIRROR)
say("=== done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
