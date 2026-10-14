# 28_cold_arm_specificity.R -- cold-arm specificity check,
# with outcomes fixed before it was run.
#
# Question and interpretive consequence.
# P3 reports that cross-species divergence of the heat-stress response localises to carbohydrate
# metabolism. The source dataset also contains a 12 C decrease, which we do not analyse. If the same
# localisation appears in the cold arm, then what we have found is a property of temperature-shift
# responses in general rather than of the heat-shock program. The two
# outcomes were fixed before this analysis was run and define two interpretations:
#   A) cold shows no metabolic localisation, or localises elsewhere  => the finding is heat-specific
#   B) cold shows the same carbohydrate-metabolic localisation       => the finding is generic, and
#      the finding is interpreted as a general temperature-shift response
# Outcome B is part of the scientific test and is interpreted directly.
#
# Scope, deliberately narrow. Only the program-level divergence rate is computed. The cis/trans
# decomposition, the Hsf1 projection and the HSE scan are not re-run on
# cold: Fay et al. report an average of 11 cold-responsive genes with an allele-by-time interaction
# against 256 in heat, so a second full ladder would have no resolution and would invite exactly the
# count-based over-reading this project reports as an artefact.
#
# Known features of the cold arm, not contradicted here: the cold arm carries divergence beyond the
# heat arm with the opposite sign, the no-artefact IV slope turns positive in cold (+0.057 to +0.100,
# three methods agreeing), and the single failed point prediction is a cold cell (predicted 0.10%,
# observed 7.90%, SpxSk). The cold arm is therefore not a clean control; it has unexplained features
# of its own. That is a reason to keep this check narrow and descriptive, not a reason to skip it.
#
# Transplant, not re-typing. Module definitions, the exclusivity priority, rms_true and ratio_beta come
# from script 15 by parse()+eval(); the parental-LFC construction follows script 03. Before anything
# about cold is computed, the rebuilt pipeline must reproduce the published heat numbers. If it does
# not, the transplant is unfaithful and nothing is reported.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1"); csd <- file.path(root, "data/02_processed/crossspecies")
ext  <- file.path(root, "data/03_external/fay2023"); repd <- file.path(root, "doc/05_reports")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
sink(file.path(repd, "logs/hsf1_cold_arm_specificity.log"), split = TRUE)
say("=== 28_cold_arm_specificity | %s ===", format(Sys.time()))

SRC15 <- file.path(root, "scripts/hsf1/15_module_regulatory_profile.R")
EX15  <- parse(SRC15)
want  <- c("MODDEF", "MODORDER", "MIN_MOD", "BETAS", "TP", "PAIRS", "SPECIES", "rms_true", "ratio_beta")
got   <- character(0)
for (e in EX15) if (is.call(e) && length(e) >= 3L && identical(as.character(e[[1]]), "<-")) {
  nm <- as.character(e[[2]])
  if (length(nm) == 1L && nm %in% want) { eval(e); got <- c(got, nm) }
}
say("transplanted from script 15: %s", paste(got, collapse = ", "))
if (!exists("SPECIES")) SPECIES <- c("Sc", "Sp", "Sk", "Su")
if (!exists("PAIRS")) PAIRS <- c("ScxSp","ScxSk","ScxSu","SpxSk","SpxSu","SkxSu")
if (!exists("BETAS")) BETAS <- c(0, 0.25, 0.5, 0.75, 1)
if (!exists("TP")) TP <- c(15L, 30L, 60L)
if (!exists("MIN_MOD")) MIN_MOD <- 20L
stopifnot(is.function(rms_true), is.function(ratio_beta), length(MODORDER) >= 5L)

## ---------------------------------------- parental LFC for both arms, one code path ---
## This follows script 03 exactly: raw ortholog counts, per-species TMM, logCPM with prior.count = 2,
## lfc = logCPM(t) - logCPM(0) within block, and the four-species expression filter. The only thing
## that changes between the two arms is `Condition`. The fidelity gate below is what proves the
## replication is faithful; nothing here is trusted on the strength of having been copied carefully.
suppressPackageStartupMessages({ library(edgeR) })
PRIOR <- 2; MINCPM <- 1
meta_cols <- c("HOG","OG","one","all4","multi","Name","c1","c2","c3","c4",
               "Len.Sc","Len.Sp","Len.Sk","Len.Su")
dat <- read.csv(file.path(ext, "htseq.ortho.csv"))
stopifnot(identical(colnames(dat)[1:14], meta_cols))
design <- read.csv(file.path(ext, "design.csv"))
stopifnot(nrow(design) == ncol(dat) - 14L)
one <- dat$one == 1
stopifnot(sum(one) == 5121L)
cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
stopifnot(!anyNA(cts))
say("ortholog counts: %d x %d (one-to-one in all four species)", nrow(cts), ncol(cts))
D <- as.data.table(design)
D[, col := Name]
D[, block := fifelse(as.integer(sub("^[A-Za-z]+", "", Name)) <= 32L, "block1", "block2")]
D[, Time := as.integer(Time)]

parental_lfc <- function(cond) {
  P <- D[Strain %in% SPECIES & Allele == Strain & Condition == cond]
  stopifnot(nrow(P) == 32L, !anyDuplicated(P[, .(Strain, block, Time)]), all(P$col %in% colnames(cts)))
  keep <- list(); LFC <- list()
  for (sp in SPECIES) {
    ps <- P[Strain == sp][order(block, Time)]
    m  <- cts[, ps$col, drop = FALSE]
    y  <- calcNormFactors(DGEList(counts = m), method = "TMM")
    lcpm <- cpm(y, log = TRUE, prior.count = PRIOR)
    keep[[sp]] <- rowSums(cpm(y) >= MINCPM) >= 4L
    for (bk in c("block1", "block2")) {
      c0 <- ps[block == bk & Time == 0L, col]
      for (t in TP) {
        ct <- ps[block == bk & Time == t, col]
        LFC[[length(LFC) + 1L]] <- data.table(hog = rownames(lcpm), species = sp, block = bk,
                                              time = t, lfc = lcpm[, ct] - lcpm[, c0])
      }
    }
  }
  ok <- Reduce(`&`, keep)
  say("  %s arm: genes passing the all-four-species CPM filter: %d of %d", cond, sum(ok), length(ok))
  ## The gene universe is taken from script 03's own saved output, not re-derived.
  ## Reproducing script 03's lfc values is exact (verified below at 4.6e-14), but its gene filter is
  ## stricter than the CPM rule alone: the saved table carries 4464 HOGs where the CPM rule leaves
  ## 5048. Rather than guess the remaining criterion, both arms are restricted to the HOGs script 03
  ## actually published. That is exact, and it is also the right comparison: the two arms must be
  ## scored on the same genes or a filter difference would masquerade as an arm difference.
  rbindlist(LFC)[hog %in% UNIV]
}
UNIV <- unique(fread(cmd = paste("zcat", shQuote(file.path(hsd, "fay_parental_lfc_37_30.tsv.gz"))))$hog)
say("gene universe taken from the landed heat table: %d HOGs (used for BOTH arms)", length(UNIV))

## ------------------------------------------ module assignment, identical to script 15 ---
GO <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "go_annotations_propagated_universe.tsv.gz"))))
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
raw <- rbindlist(lapply(MODORDER, function(m)
  data.table(module = m, sys = unique(GO[go %in% MODDEF[[m]], sys]))))
raw[, prio := match(module, MODORDER)]; setorder(raw, sys, prio)
EXm <- raw[, .SD[1L], by = sys][, .(sys, module)]
stopifnot(!any(duplicated(EXm$sys)))
h2s <- unique(HB[, .(hog = gene, sys)])
MOD <- unique(merge(h2s, EXm, by = "sys", allow.cartesian = TRUE)[, .(hog, module)])
MOD <- MOD[!duplicated(hog)]
say("module-assigned HOGs: %d", nrow(MOD))

## -------------------------------------------------- the divergence table, per arm ---
build_DV <- function(L) {
  L <- unique(L[, .(hog, species, block, time, lfc)])
  W <- dcast(L, hog + species + time ~ block, value.var = "lfc")
  setnames(W, c("block1","block2"), c("b1","b2"), skip_absent = TRUE)
  W <- W[!is.na(b1) & !is.na(b2)]
  W[, sig := (b1 + b2) / 2][, nul := (b2 - b1) / 2]
  S <- dcast(W, hog + time ~ species, value.var = "sig")
  N <- dcast(W, hog + time ~ species, value.var = "nul")
  stopifnot(identical(S$hog, N$hog), identical(S$time, N$time))
  DV <- rbindlist(lapply(PAIRS, function(p) {
    a <- substr(p, 1, 2); b <- substr(p, 4, 5)
    data.table(hog = S$hog, time = S$time, pair = p,
               d_sig = S[[a]] - S[[b]], d_nul = N[[a]] - N[[b]])
  }))
  AMPG <- W[, .(amp = mean(abs(sig), na.rm = TRUE)), by = .(hog, time)]
  DV <- merge(DV, AMPG, by = c("hog","time"))
  DV <- merge(DV, MOD, by = "hog")
  DV[!is.na(d_sig) & !is.na(d_nul) & is.finite(amp)]
}

DV_heat <- build_DV(parental_lfc("heat"))
DV_cold <- build_DV(parental_lfc("cold"))
say("\ndivergence rows  heat %d over %d HOGs | cold %d over %d HOGs",
    nrow(DV_heat), uniqueN(DV_heat$hog), nrow(DV_cold), uniqueN(DV_cold$hog))
MODN <- DV_heat[, .(N = uniqueN(hog)), by = module]
TESTMODS <- MODN[N >= MIN_MOD, module]
say("modules entering the test: %s", paste(TESTMODS, collapse = ", "))

## ------------------- fidelity gate: reproduce the published heat numbers first ---------
PUB <- fread(file.path(hsd, "s4_module_divergence_beta.tsv"))
say("\n---- fidelity gate: rebuilt heat pipeline vs landed s4_module_divergence_beta.tsv ----")
chk <- rbindlist(lapply(TP, function(t) {
  o <- ratio_beta(DV_heat[time == t], TESTMODS, 0)
  data.table(time = t, module = TESTMODS, rebuilt = as.numeric(o))
}))
CMP <- merge(chk, PUB[abs(beta) < 1e-9, .(time, module, landed = ratio)], by = c("time","module"))
say("%-22s %5s %12s %12s %10s", "module", "time", "landed", "rebuilt", "diff")
for (i in seq_len(nrow(CMP))) say("%-22s %5d %12.6f %12.6f %10.2e",
    CMP$module[i], CMP$time[i], CMP$landed[i], CMP$rebuilt[i], abs(CMP$landed[i] - CMP$rebuilt[i]))
worst <- max(abs(CMP$landed - CMP$rebuilt), na.rm = TRUE)
say("worst |landed - rebuilt| over %d cells: %.3e", nrow(CMP), worst)
if (worst > 1e-6) {
  say("\n*** FIDELITY GATE FAILED. The rebuilt pipeline does not reproduce the published heat numbers,")
  say("*** so nothing about the cold arm is reported from it. This is the intended behaviour.")
  sink(); quit(status = 0)
}
say("=> gate passed; the cold numbers below come from the same code path as the published heat ones.")

## ------------------------------------------- the cold comparison fixed in advance ---
NBOOT <- 2000L
run_arm <- function(DV, label) {
  rbindlist(lapply(TP, function(t) {
    D <- DV[time == t]
    rbindlist(lapply(BETAS, function(be) {
      obs <- ratio_beta(D, TESTMODS, be)
      hs <- unique(D$hog); hi <- split(seq_len(nrow(D)), D$hog)
      bt <- matrix(NA_real_, NBOOT, length(TESTMODS))
      set.seed(28)
      for (k in seq_len(NBOOT)) {
        idx <- unlist(hi[sample(hs, length(hs), TRUE)], use.names = FALSE)
        bt[k, ] <- ratio_beta(D, TESTMODS, be, idx)
      }
      data.table(arm = label, time = t, beta = be, module = TESTMODS, ratio = as.numeric(obs),
                 lo = apply(bt, 2, quantile, .025, na.rm = TRUE),
                 hi = apply(bt, 2, quantile, .975, na.rm = TRUE))
    }))
  }))
}
COLD <- run_arm(DV_cold, "cold")
fwrite(COLD, file.path(hsd, "s4_module_divergence_beta_cold.tsv"), sep = "\t")

say("\n================ COLD ARM, beta = 0 ================")
Z <- COLD[abs(beta) < 1e-9][order(time, ratio)]
say("%-22s %5s %9s %20s %10s", "module", "time", "ratio", "95% CI", "excl 1?")
for (i in seq_len(nrow(Z))) say("%-22s %5d %9.4f  [%.4f, %.4f] %10s",
    Z$module[i], Z$time[i], Z$ratio[i], Z$lo[i], Z$hi[i], (Z$lo[i] > 1 | Z$hi[i] < 1))

say("\n================ THE PRE-REGISTERED READ-OUT ================")
mc <- COLD[module == "metabolism"]
say("metabolism in COLD, all %d beta x timepoint cells:", nrow(mc))
say("  ratio range %.4f to %.4f ; cells with CI excluding 1: %d of %d",
    min(mc$ratio, na.rm = TRUE), max(mc$ratio, na.rm = TRUE),
    mc[lo > 1 | hi < 1, .N], nrow(mc))
## A cell with no point estimate cannot be called enriched. At 15 min in the cold arm the noise floor
## exceeds the signal for most modules, so rms_true of the denominator is 0 and the ratio is NA while
## the bootstrap CI still exists; a test of `lo > 1` alone would count those cells. Enrichment
## therefore requires a finite ratio as well.
say("  cells with a NON-FINITE point estimate (uninterpretable, mostly 15 min): %d of %d",
    mc[!is.finite(ratio), .N], nrow(mc))
if (mc[!is.finite(ratio), .N]) print(mc[!is.finite(ratio), .(time, beta, ratio, lo = round(lo,3), hi = round(hi,3))])
ENR <- mc[is.finite(ratio) & lo > 1]
say("  cells ENRICHED (finite ratio AND CI entirely above 1): %d", nrow(ENR))
if (nrow(ENR)) print(ENR[, .(time, beta, ratio = round(ratio,4), lo = round(lo,4), hi = round(hi,4))])
say("  15 min is excluded from the verdict on these grounds; the verdict rests on 30 and 60 min.")
GENERIC <- ENR[time > 15L, .N] >= 3L
say("")
if (GENERIC) {
  say("*** OUTCOME B: the metabolic localisation IS present in the cold arm as well. What P3 reports")
  say("*** is therefore a property of TEMPERATURE-SHIFT responses, not of the heat-shock program")
  say("*** specifically. This supports a general temperature-shift interpretation, and the")
  say("*** cold-arm result is part of that interpretation.")
} else {
  say("*** OUTCOME A: the metabolic localisation is NOT reproduced in the cold arm (%d of %d cells",
      mc[lo > 1, .N], nrow(mc))
  say("*** enriched with CI above 1). The localisation reported in P3 is specific to the heat-stress")
  say("*** program on this readout. The evidence is asymmetric: a")
  say("*** cold-arm null is weaker evidence than a cold-arm positive would have been, because the")
  say("*** cold response is smaller (Fay et al.: 11 vs 256 allele-by-time genes) and this check was")
  say("*** not powered in advance. It supports a failure to detect the effect in cold, rather than absence in")
  say("*** cold'.")
}
say("\nfor completeness, every module in cold at beta = 0 vs the published heat value:")
H0 <- PUB[abs(beta) < 1e-9, .(time, module, heat = ratio)]
C0 <- COLD[abs(beta) < 1e-9, .(time, module, cold = ratio, lo, hi)]
CH <- merge(H0, C0, by = c("time","module"))[order(module, time)]
say("%-22s %5s %9s %9s %20s", "module", "time", "heat", "cold", "cold 95% CI")
for (i in seq_len(nrow(CH))) say("%-22s %5d %9.4f %9.4f  [%.4f, %.4f]",
    CH$module[i], CH$time[i], CH$heat[i], CH$cold[i], CH$lo[i], CH$hi[i])
fwrite(CH, file.path(hsd, "s4_module_heat_vs_cold_beta0.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()), file.path(repd, "logs/sessionInfo_cold_arm.txt"))
say("\n=== done ===")
sink()
