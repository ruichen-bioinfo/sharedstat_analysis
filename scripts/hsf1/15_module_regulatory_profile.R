# 15_module_regulatory_profile.R -- which cellular programs diverged? Decision rule fixed before
# the analysis was run.
#
# The question. Which cellular programs diverged, not merely which genes. This is the only
# cell-biological endpoint in the analysis, linking transcriptional response to
# cellular process.
#
# Prediction, stated before the numbers. The load-driven candidate explanation says divergence is
# driven by proteome-level thermal stability differences. If so:
#   (i) divergence should concentrate in the chaperone/proteostasis module rather than be
#       uniform across modules, and
#   (ii) the proteostasis module's response amplitude should order the species the way
#        thermotolerance does (Sc, Sp thermophilic; Sk, Su cryophilic; Sc vs Su has the
#        PMID 40476605 anchor of 8 C thermotolerance difference and 1.6 C mean melting-point
#        shift; only the abstract was used).
# Outcome that can fail: if module-level divergence is uniform across modules, the load-driven
# hypothesis is not supported. This outcome and the supporting outcome are interpreted under
# the same criteria. Neither outcome demonstrates that upstream regulation evolved or that
# Model A holds.
#
# Amplitude confounding is at its worst here and is handled throughout. Modules differ
# enormously in response amplitude (chaperones are the most strongly induced genes in the cell),
# and rms-type divergence statistics scale with amplitude. The same mechanism turns script 04's
# apparent Msn2/4 positive (1.27) into ~1.0 once corrected. So there is no uncorrected version of
# this analysis: every module ratio is reported across the whole beta in [0, 1] range, where the
# divergence axis is divided by amp^beta before the statistic is formed. A ratio that survives
# only at beta = 0 is classified as an amplitude artefact.
#
# Module definitions are fixed in advance from external sources. Terms come from go-basic.obo and
# the SGD GAF, both sha256-recorded in manifests/EXTERNAL_INPUTS.tsv, via the propagated
# annotation table of script 02. No module is invented and no term is chosen after seeing a
# result. GO:0051082 (unfolded protein binding) is not used because it is obsolete in the
# ontology release used.
#   proteostasis         GO:0006457 protein folding            + GO:0009408 response to heat
#   translation_ribosome GO:0002181 cytoplasmic translation    + GO:0042254 ribosome biogenesis
#   oxidative            GO:0006979 response to oxidative stress
#   metabolism           GO:0005975 carbohydrate metabolic proc + GO:0006096 glycolytic process
#   cell_cycle_growth    GO:0007049 cell cycle
# Modules are made mutually exclusive in the fixed priority order above, and exclusivity is
# asserted. Reason: comparing a subset against a superset dilutes automatically, so the
# denominator must be mutually exclusive and checked. Overlap counts are reported so the cost of
# the exclusivity rule is visible.
#
# Paralog collapse. The one verified expanded family in this dataset is the anCWMP family
# TIP1/TIR1-TIR4/DAN1/DAN4 (PMID 36134526). Every module ratio is re-run on the collapsed
# background and reported side by side, because a category-level cross-species claim that does
# not survive collapse is not supported.
#
# n = 4 species. Species-level statements are descriptive only; no significance is claimed at
# the species or pair level. All testing is at the gene level (n ~ 4000).
#
# Usage: Rscript scripts/hsf1/15_module_regulatory_profile.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
repd <- file.path(root, "doc/05_reports")
outd <- hsd
MIRROR <- "logs/output_copy"

TP     <- c(15L, 30L, 60L)
SPP    <- c("Sc", "Sp", "Sk", "Su")
THERMAL <- c(Sc = "thermophilic", Sp = "thermophilic", Sk = "cryophilic", Su = "cryophilic")
PAIRS  <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")
BETAS  <- c(0, 0.25, 0.50, 0.75, 1.00)
NBOOT  <- 1000L
NPERM  <- 1000L
MIN_MOD <- 20L        # a module below this is reported but not tested
anCWMP <- c("TIP1", "TIR1", "TIR2", "TIR3", "TIR4", "DAN1", "DAN4")   # PMID 36134526
set.seed(20260906L)

MODDEF <- list(
  proteostasis         = c("GO:0006457", "GO:0009408"),
  translation_ribosome = c("GO:0002181", "GO:0042254"),
  oxidative            = c("GO:0006979"),
  metabolism           = c("GO:0005975", "GO:0006096"),
  cell_cycle_growth    = c("GO:0007049"))
MODORDER <- names(MODDEF)      # exclusivity priority, fixed in advance

say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, "logs/hsf1_module_profile.log")
sink(logf, split = TRUE)
say("=== 15_module_regulatory_profile | %s | pre-registered ===",
    format(Sys.time()))
say("Pre-registered prediction: load-driven => divergence CONCENTRATED in proteostasis and the")
say("proteostasis amplitude ordering tracks thermotolerance. Alternative outcome: UNIFORM across")
say("modules => load-driven NOT supported. Neither outcome implies Model A.")
say("Amplitude confounding is handled by evaluating every ratio across beta in [0,1].")

## ------------------------------------------------------ modules ------------
GO <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "go_annotations_propagated_universe.tsv.gz"))))
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
L  <- unique(fread(cmd = paste("zcat", shQuote(file.path(hsd, "fay_parental_lfc_37_30.tsv.gz"))))[
             , .(hog, species, block, time, lfc)])
say("\nparental LFC: %d rows, %d HOGs, %d species x %d blocks x %d timepoints",
    nrow(L), uniqueN(L$hog), uniqueN(L$species), uniqueN(L$block), uniqueN(L$time))
h2s <- unique(HB[, .(hog = gene, sys)])
say("HOG -> Sc systematic bridge: %d pairs covering %d HOGs", nrow(h2s), uniqueN(h2s$hog))

say("\n---- module assignment (locked GO terms, pre-registered exclusivity priority) ----")
raw <- rbindlist(lapply(MODORDER, function(m)
  data.table(module = m, sys = unique(GO[go %in% MODDEF[[m]], sys]))))
say("raw (non-exclusive) module sizes in Sc systematic names:")
for (m in MODORDER) say("  %-22s %5d genes  (%s)", m, raw[module == m, .N],
                        paste(MODDEF[[m]], collapse = " + "))
ovl <- raw[, .N, by = sys][N > 1L, .N]
say("genes in more than one module before exclusivity: %d", ovl)
raw[, prio := match(module, MODORDER)]
setorder(raw, sys, prio)
EX <- raw[, .SD[1L], by = sys][, .(sys, module)]
stopifnot(!any(duplicated(EX$sys)))          # the denominator must be exclusive
say("after exclusivity (priority %s):", paste(MODORDER, collapse = " > "))
for (m in MODORDER) say("  %-22s %5d genes", m, EX[module == m, .N])
say("exclusivity asserted: no gene appears in two modules.")

MOD <- merge(h2s, EX, by = "sys", allow.cartesian = TRUE)
MOD <- unique(MOD[, .(hog, module, sys)])
MOD <- MOD[!duplicated(hog)]                 # one module per HOG
say("\nmapped onto HOGs present in the Fay universe:")
MODN <- MOD[hog %in% unique(L$hog), .N, by = module]
for (m in MODORDER) say("  %-22s %5d HOGs %s", m,
  ifelse(m %in% MODN$module, MODN[module == m, N], 0L),
  ifelse((if (m %in% MODN$module) MODN[module == m, N] else 0L) < MIN_MOD, "<- below test threshold", ""))
fwrite(MOD, file.path(outd, "s4_module_membership.tsv"), sep = "\t")

## --------------------------- signal / null from the two replicate blocks ----
# Same estimator as script 03: two complete blocks give signal = mean and null = half-difference.
# The null carries no heat signal by construction, so rms_true = sqrt(max(0, E[sig^2]-E[null^2]))
# removes the species-specific noise floor, which matters because the four libraries differ
# about 30-fold in depth.
W <- dcast(L, hog + species + time ~ block, value.var = "lfc")
setnames(W, c("block1", "block2"), c("b1", "b2"))
W <- W[!is.na(b1) & !is.na(b2)]
W[, sig := (b1 + b2) / 2][, nul := (b2 - b1) / 2]
say("\nsignal/null built for %d HOG x species x time cells", nrow(W))
rms_true <- function(s, n) sqrt(max(0, mean(s^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE)))

## ------------------- (A) amplitude profile: module x species x timepoint ----
say("\n================ (A) MODULE AMPLITUDE PROFILE (noise-corrected) ================")
WM <- merge(W, MOD[, .(hog, module)], by = "hog")
AMP <- WM[, .(n = .N, rms_sig = sqrt(mean(sig^2)), rms_null = sqrt(mean(nul^2)),
              rms_true = rms_true(sig, nul)), by = .(module, species, time)]
AMP[, snr := rms_true / rms_null]
fwrite(AMP, file.path(outd, "s4_module_amplitude.tsv"), sep = "\t")
for (t in TP) {
  say("\ntime %d min -- rms_true (noise-corrected response amplitude):", t)
  say("%-22s %9s %9s %9s %9s", "module", SPP[1], SPP[2], SPP[3], SPP[4])
  for (m in MODORDER) {
    z <- AMP[module == m & time == t]
    if (!nrow(z)) next
    v <- setNames(z$rms_true, z$species)
    say("%-22s %9.3f %9.3f %9.3f %9.3f", m, v[SPP[1]], v[SPP[2]], v[SPP[3]], v[SPP[4]])
  }
}
say("\nmodules differ in amplitude by a large factor, which is exactly why every divergence")
say("ratio below is reported across beta in [0,1] rather than at beta = 0 alone.")

## --- (A2) check (ii), fixed in advance: does proteostasis amplitude track thermotolerance? ---
say("\n---- pre-registered check (ii): species ordering of proteostasis amplitude ----")
say("DESCRIPTIVE ONLY, n = 4 species; no significance is claimed.")
ORD <- rbindlist(lapply(TP, function(t) {
  z <- AMP[module == "proteostasis" & time == t]
  o <- z[order(-rms_true)]
  data.table(time = t, ordering = paste(o$species, collapse = " > "),
             thermal = paste(THERMAL[o$species], collapse = " > "),
             thermo_first_two = all(THERMAL[o$species[1:2]] == "thermophilic"),
             cryo_first_two   = all(THERMAL[o$species[1:2]] == "cryophilic"))
}))
fwrite(ORD, file.path(outd, "s4_thermal_ordering.tsv"), sep = "\t")
for (i in seq_len(nrow(ORD))) say("  %2d min: %-22s (%s)", ORD$time[i], ORD$ordering[i], ORD$thermal[i])
say("prediction (ii) would require the two thermophiles (Sc, Sp) to sit on one end at every")
say("timepoint. thermophiles-first at 15/30/60: %s / %s / %s ; cryophiles-first: %s / %s / %s",
    ORD$thermo_first_two[1], ORD$thermo_first_two[2], ORD$thermo_first_two[3],
    ORD$cryo_first_two[1], ORD$cryo_first_two[2], ORD$cryo_first_two[3])

## ---------- (B) divergence per module, beta-calibrated, gene-level test -----
say("\n================ (B) MODULE-LEVEL DIVERGENCE, beta-CALIBRATED ================")
say("Divergence per gene x species pair x time: d_sig = sig_i - sig_j, d_null = nul_i - nul_j,")
say("using the same two-block construction, so the noise floor of BOTH species is subtracted.")
S <- dcast(W, hog + time ~ species, value.var = "sig")
N <- dcast(W, hog + time ~ species, value.var = "nul")
stopifnot(identical(S$hog, N$hog), identical(S$time, N$time))
DV <- rbindlist(lapply(PAIRS, function(p) {
  a <- substr(p, 1, 2); b <- substr(p, 4, 5)
  data.table(hog = S$hog, time = S$time, pair = p,
             d_sig = S[[a]] - S[[b]], d_nul = N[[a]] - N[[b]])
}))
# amplitude per gene x time = mean |sig| across the four species (species-symmetric, so it
# cannot itself encode which pair is being compared)
AMPG <- W[, .(amp = mean(abs(sig), na.rm = TRUE)), by = .(hog, time)]
DV <- merge(DV, AMPG, by = c("hog", "time"))
DV <- merge(DV, MOD[, .(hog, module)], by = "hog")
DV <- DV[!is.na(d_sig) & !is.na(d_nul) & is.finite(amp)]
say("divergence rows (HOG x pair x time, module-assigned): %d over %d HOGs",
    nrow(DV), uniqueN(DV$hog))

# ratio_beta(module) = rms_true(module | d/amp^beta) / rms_true(all other modules | d/amp^beta)
# The denominator is the union of the other exclusive modules, never the superset.
ratio_beta <- function(DT, mods, beta, idx = NULL) {
  d <- DT
  if (!is.null(idx)) d <- DT[idx]
  sc <- pmax(d$amp, 1e-9)^beta
  ds <- d$d_sig / sc; dn <- d$d_nul / sc
  out <- vapply(mods, function(m) {
    ins <- d$module == m; ous <- !ins
    if (sum(ins) < 5L || sum(ous) < 5L) return(NA_real_)
    vi <- rms_true(ds[ins], dn[ins]); vo <- rms_true(ds[ous], dn[ous])
    if (vo > 0) vi / vo else NA_real_
  }, 0)
  out
}
TESTMODS <- MODN[N >= MIN_MOD, module]
say("modules with >= %d HOGs, entering the test: %s", MIN_MOD, paste(TESTMODS, collapse = ", "))

RB <- list()
for (t in TP) {
  DTt <- DV[time == t]
  for (be in BETAS) {
    obs <- ratio_beta(DTt, TESTMODS, be)
    # gene-level bootstrap: resample HOGs (not rows), so all pairs of a gene move together
    hs <- unique(DTt$hog); hi <- split(seq_len(nrow(DTt)), DTt$hog)
    bt <- matrix(NA_real_, NBOOT, length(TESTMODS))
    for (k in seq_len(NBOOT)) {
      idx <- unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)
      bt[k, ] <- ratio_beta(DTt, TESTMODS, be, idx)
    }
    for (j in seq_along(TESTMODS)) RB[[length(RB) + 1L]] <- data.table(
      time = t, beta = be, module = TESTMODS[j], ratio = obs[j],
      lo = quantile(bt[, j], .025, na.rm = TRUE), hi = quantile(bt[, j], .975, na.rm = TRUE))
  }
  say("  time %d min done", t)
}
RB <- rbindlist(RB); fwrite(RB, file.path(outd, "s4_module_divergence_beta.tsv"), sep = "\t")
say("\nratio = rms_true(module) / rms_true(all other modules), divergence divided by amp^beta")
for (t in TP) {
  say("\ntime %d min:", t)
  say("%-22s %s", "module", paste(sprintf("%18s", sprintf("beta=%.2f", BETAS)), collapse = ""))
  for (m in TESTMODS) {
    cells <- vapply(BETAS, function(be) {
      z <- RB[time == t & beta == be & module == m]
      sprintf("%6.3f[%.2f,%.2f]", z$ratio, z$lo, z$hi) }, "")
    say("%-22s %s", m, paste(sprintf("%18s", cells), collapse = ""))
  }
}
say("\nA ratio that excludes 1 only at beta = 0 and collapses toward 1 as beta rises IS an")
say("amplitude artefact -- the A5 Msn2/4 signature. A ratio stable across the whole beta range is")
say("not explained by amplitude.")

## ------- (C) the control that can fail: is divergence uniform across modules? ---
say("\n================ (C) FAILABLE CONTROL: UNIFORMITY ACROSS MODULES ================")
say("if module-level divergence is uniform, the load-driven hypothesis is NOT")
say("supported; both outcomes are evaluated under the same criteria.")
say("Null = permute module labels among genes at FIXED module sizes (gene-level, so it is a")
say("real can-fail control and not the observed data relabelled as a null).")
UNI <- list()
for (t in TP) {
  DTt <- DV[time == t]
  hmod <- unique(DTt[, .(hog, module)])
  for (be in BETAS) {
    obs <- ratio_beta(DTt, TESTMODS, be)
    obs_spread <- max(obs, na.rm = TRUE) - min(obs, na.rm = TRUE)
    nl <- matrix(NA_real_, NPERM, length(TESTMODS)); nsp <- numeric(NPERM)
    for (k in seq_len(NPERM)) {
      hp <- copy(hmod); hp[, module := sample(module)]
      d2 <- merge(DTt[, .(hog, d_sig, d_nul, amp)], hp, by = "hog")
      v <- ratio_beta(d2, TESTMODS, be)
      nl[k, ] <- v; nsp[k] <- max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
    }
    for (j in seq_along(TESTMODS)) {
      p <- mean(abs(nl[, j] - 1) >= abs(obs[j] - 1), na.rm = TRUE)
      UNI[[length(UNI) + 1L]] <- data.table(time = t, beta = be, module = TESTMODS[j],
        ratio = obs[j], p_emp = p, null_med = median(nl[, j], na.rm = TRUE))
    }
    UNI[[length(UNI) + 1L]] <- data.table(time = t, beta = be, module = "__SPREAD__",
      ratio = obs_spread, p_emp = mean(nsp >= obs_spread, na.rm = TRUE),
      null_med = median(nsp, na.rm = TRUE))
  }
  say("  time %d min done", t)
}
UNI <- rbindlist(UNI); fwrite(UNI, file.path(outd, "s4_uniformity_test.tsv"), sep = "\t")
say("\nmodule-label permutation, p_emp (two-sided about 1; __SPREAD__ is max-min across modules)")
for (t in TP) {
  say("\ntime %d min:", t)
  say("%-22s %s", "module", paste(sprintf("%14s", sprintf("beta=%.2f", BETAS)), collapse = ""))
  for (m in c(TESTMODS, "__SPREAD__")) {
    cells <- vapply(BETAS, function(be) {
      z <- UNI[time == t & beta == be & module == m]; sprintf("%.3f(p=%.3f)", z$ratio, z$p_emp) }, "")
    say("%-22s %s", m, paste(sprintf("%14s", cells), collapse = ""))
  }
}
SPREAD_SIG <- UNI[module == "__SPREAD__" & beta >= 0.25, any(p_emp < 0.05)]
PROT_SIG   <- UNI[module == "proteostasis" & beta >= 0.25, any(p_emp < 0.05)]
say("\nacross beta >= 0.25 (i.e. with at least partial amplitude correction):")
say("  any timepoint where the across-module SPREAD exceeds its permutation null: %s", SPREAD_SIG)
say("  any timepoint where PROTEOSTASIS itself differs from the other modules:    %s", PROT_SIG)

## ------------------- (D) paralog-collapse sensitivity --------
say("\n================ (D) PARALOG-COLLAPSE SENSITIVITY ================")
say("The one verified expanded family here is anCWMP TIP1/TIR1-TIR4/DAN1/DAN4 (PMID 36134526).")
say("A category-level cross-species result that does not survive collapse is not supported.")
drop_hogs <- MOD[sys %in% anCWMP[-1L], hog]
say("anCWMP HOGs present: %d ; dropped to a single representative: %d",
    MOD[sys %in% anCWMP, .N], length(drop_hogs))
COL <- rbindlist(lapply(TP, function(t) {
  d <- DV[time == t & !(hog %in% drop_hogs)]
  rbindlist(lapply(BETAS, function(be) {
    v <- ratio_beta(d, TESTMODS, be)
    data.table(time = t, beta = be, module = TESTMODS, ratio_collapsed = v) }))
}))
COL <- merge(COL, RB[, .(time, beta, module, ratio_full = ratio)],
             by = c("time", "beta", "module"))
COL[, delta := ratio_collapsed - ratio_full]
fwrite(COL, file.path(outd, "s4_paralog_collapse.tsv"), sep = "\t")
say("max |ratio change| after collapse, across all module x time x beta: %.4f",
    max(abs(COL$delta), na.rm = TRUE))
say("=> %s", ifelse(max(abs(COL$delta), na.rm = TRUE) < 0.02,
  "collapse changes nothing material; every module conclusion survives it",
  "collapse moves at least one ratio; the affected module must be FLAGGED, not reported plainly"))

## ------------------------------------------------- what this analysis concludes ---
say("\n================ WHAT CAN BE CONCLUDED ================")
prot_stable <- RB[module == "proteostasis"][, .(any_excl = any(lo > 1 | hi < 1)), by = beta]
say("proteostasis ratio excludes 1 at each beta: %s",
    paste(sprintf("%.2f:%s", prot_stable$beta, prot_stable$any_excl), collapse = "  "))
say("")
if (!SPREAD_SIG) {
  say("VERDICT -- THE PRE-REGISTERED FAILABLE OUTCOME FIRED:")
  say("module-level divergence is UNIFORM ACROSS MODULES once amplitude is even partially")
  say("corrected. The across-module spread does not exceed its gene-level module-label")
  say("permutation null at any timepoint for beta >= 0.25.")
  say("=> THE LOAD-DRIVEN HYPOTHESIS IS NOT SUPPORTED by module-level divergence.")
  say("This is a real, pre-registered negative, not an absence of analysis: the control could")
  say("have fired the other way and did not.")
} else {
  say("VERDICT: the across-module spread exceeds its permutation null with amplitude correction")
  say("in place, so divergence is NOT uniform across cellular programs. Read which module carries")
  say("it from the table above, and only at beta values where its CI excludes 1.")
}
say("")
say("The analysis does not support:")
say(" - an inference that upstream regulation evolved or that Model A holds")
say(" - species- or hybrid-level significance (n = 4 / n = 6)")
say(" - a biological interpretation of a ratio that excludes 1 only at beta = 0")
say("")
say("Notes:")
say(" - censoring of the detection bound. Pooled over 52 substrates, the number of")
say("   substrates NOT detected anywhere up to the top of the tested grid (c = 0.80) is ZERO at")
say("   all three timepoints. The 90th-percentile bound is therefore exactly 0.80 and is NOT")
say("   right-censored; the interval is c ~ 0.02 - 0.80.")
say(" - MID2. Mid2 identities (Sp 90.2 / Sk 80.0 / Su 75.3) lack the background needed")
say("   for comparison with Hsf1's (Sp 91.0 / Sk 81.4 / Su 78.5) as evidence that Mid2 is unusually")
say("   diverged: that needs a genome-wide pairwise-identity background distribution, which was")
say("   not computed. A two-point comparison without a background is inconclusive.")
say(" - For PMID 41390490 (Mid2/Ssa1-T492) and PMID 40476605 (Sc/Su thermostability,")
say("   8 C, 85%%, 1.6 C), only the abstracts were consulted.")
say("")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_s4_module.txt"))
mdir <- file.path(MIRROR, "data/02_processed/hsf1")
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(outd, pattern = "^(s4_|sessionInfo_s4)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
dir.create(file.path(MIRROR, "scripts/hsf1"), showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(root, "scripts/hsf1/15_module_regulatory_profile.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("
copied outputs to %s", MIRROR)
say("=== done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
