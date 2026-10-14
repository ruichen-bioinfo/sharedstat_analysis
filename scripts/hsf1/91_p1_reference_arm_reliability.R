#!/usr/bin/env Rscript
## 91_p1_reference_arm_reliability.R -- reference-arm reliability, not the DoD-contrast reliability. The measurement
## and its reporting were fixed before this analysis was run.
##
## The manuscript identity is gamma_hat_naive = gamma * h, with h = var(R)/var(R-hat), R the
## reference arm's true heat response (first-named allele). Script 36 measures the
## reliability of the DoD contrast e (its own header says so). This script measures h on the
## reference arm, using the TMM, logCPM, and detection filter from script 04, alongside
## the DoD reliability with explicit labels. Outputs are written to separate files.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
source(file.path(root, "pkg/sharedstat/R/reliability.R"))

## output helper (as in script 89)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

## ---- constants from 04_divergence_by_dependency.R ----
ext  <- file.path(root, "data/03_external/fay2023")
csd  <- file.path(root, "data/02_processed/crossspecies")
hsd  <- file.path(root, "data/02_processed/hsf1")
ANCHOR  <- "37_30"
HYB     <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")
TP      <- c(15L, 30L, 60L)
PRIOR   <- 2
MINCPM  <- 1

## ---- inputs shared with script 04 ----
dat <- read.csv(file.path(ext, "htseq.ortho.csv"))
design <- read.csv(file.path(ext, "design.csv"))
one <- dat$one == 1
stopifnot(sum(one) == 5121L)
cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
stopifnot(!anyNA(cts))

D <- as.data.table(design); D[, col := Name]
D[, num := as.integer(sub("^[A-Za-z]+", "", Name))]
D[, block := fifelse(num <= 32L, "block1", "block2")]
H <- D[Strain %in% HYB & Condition == "heat"]; H[, Time := as.integer(Time)]
stopifnot(nrow(H) == 96L, !anyDuplicated(H[, .(Strain, Allele, block, Time)]))
POL <- data.table(hybrid = HYB, a1 = substr(HYB, 1, 2), a2 = substr(HYB, 4, 5))

## ---- gene strata shared with script 04 ----
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
CS <- fread(file.path(hsd, sprintf("hsf1_core_set_%s.tsv", ANCHOR)))
CL <- fread(file.path(hsd, sprintf("hsf1_epistasis_classes_%s.tsv", ANCHOR)))
UC <- dcast(CL[, .(gene, method, class)], gene ~ method, value.var = "class")
UC[, unanimous := fifelse(DESeq2 == edgeR & edgeR == limma, DESeq2, NA_character_)]
core  <- CS[set == "Hsf1_core_unanimous",  gene]
joint <- CS[set == "Hsf1_joint_candidate", gene]
map_hog <- function(g) { h <- HB[Name %in% g | sys %in% g, unique(gene)]; h[h %in% rownames(cts)] }
STRATA <- list(
  Hsf1_core      = core,
  Hsf1_joint     = joint,
  Msn24_only     = UC[unanimous == "Msn24_only",    gene],
  TF_independent = UC[unanimous == "TF_independent", gene],
  other_induced  = UC[unanimous %in% c("Msn24_only", "TF_independent"), gene],
  not_induced    = UC[unanimous == "not_induced",    gene]
)
HOGS <- lapply(STRATA, map_hog)
## Note: a gene may belong to several strata (other_induced is the union of Msn24_only and
## TF_independent, exactly as 04:183,304). Assigning each gene to the first matching stratum
## would leave other_induced empty (18 NA cells) and make the denominator 90 instead of the
## 108, so each stratum is queried independently against HOGS, as 04 does.
stopifnot(setequal(HOGS$other_induced, union(HOGS$Msn24_only, HOGS$TF_independent)))

## ---- per-hybrid logCPM + reference-arm two-block responses (04's block, arm form) ----
ARM <- list()   # per gene: reference-arm response in each block, per hybrid x time
DOD <- list()   # per gene: the DoD sig/nul, to reproduce script 36's quantity beside it
for (g in HYB) {
  hs <- H[Strain == g][order(Allele, block, Time)]
  a1 <- POL[hybrid == g, a1]; a2 <- POL[hybrid == g, a2]
  m  <- cts[, hs$col, drop = FALSE]
  y  <- calcNormFactors(DGEList(counts = m), method = "TMM")   # verbatim
  lcpm <- cpm(y, log = TRUE, prior.count = PRIOR)              # verbatim
  cpm_raw <- cpm(y)
  k1 <- rowSums(cpm_raw[, hs[Allele == a1, col], drop = FALSE] >= MINCPM) >= 4L
  k2 <- rowSums(cpm_raw[, hs[Allele == a2, col], drop = FALSE] >= MINCPM) >= 4L
  keep <- k1 & k2
  gv <- function(al, bk, t) { cc <- hs[Allele == al & block == bk & Time == t, col]
                              stopifnot(length(cc) == 1L); lcpm[, cc] }
  say("  %-6s a1=%s detected both alleles: %d", g, a1, sum(keep))
  for (t in TP) {
    rA1 <- gv(a1, "block1", t) - gv(a1, "block1", 0L)   # reference arm, block 1
    rA2 <- gv(a1, "block2", t) - gv(a1, "block2", 0L)   # reference arm, block 2
    rB1 <- gv(a2, "block1", t) - gv(a2, "block1", 0L)
    rB2 <- gv(a2, "block2", t) - gv(a2, "block2", 0L)
    hog <- rownames(lcpm)
    ARM[[length(ARM)+1L]] <- data.table(hog = hog, hybrid = g, time = t,
                                         rA1 = rA1, rA2 = rA2, keep = keep)
    DOD[[length(DOD)+1L]] <- data.table(hog = hog, hybrid = g, time = t,
                                         d1 = rB1 - rA1, d2 = rB2 - rA2, keep = keep)
  }
}
ARM <- rbindlist(ARM); DOD <- rbindlist(DOD)
## ARM/DOD are keyed by gene only; stratum membership is applied per-stratum below (a gene
## can be in several strata, e.g. Msn24_only and other_induced), matching 04.
setkey(ARM, hybrid, time); setkey(DOD, hybrid, time)

## ---- reliability per stratum x hybrid x timepoint ----
## reference arm: pass the two block responses to reliability_from_blocks (center=TRUE);
## it forms sig=(rA1+rA2)/2, nul=(rA2-rA1)/2 internally, so h_arm = 1 - ms_nul/ms_sig.
rel <- function(dt, b1col, b2col) {
  x <- dt[keep == TRUE & is.finite(get(b1col)) & is.finite(get(b2col))]
  if (nrow(x) < 3L) return(data.table(h = NA_real_, snr = NA_real_, n = nrow(x)))
  r <- reliability_from_blocks(x[[b1col]], x[[b2col]], center = TRUE)
  data.table(h = r$h, snr = r$snr, n = r$n_pairs)
}
grid <- CJ(stratum = names(HOGS), hybrid = HYB, time = TP, sorted = FALSE)
res <- grid[, {
  gg <- HOGS[[.BY$stratum]]
  a <- ARM[.(.BY$hybrid, .BY$time)][hog %in% gg]
  d <- DOD[.(.BY$hybrid, .BY$time)][hog %in% gg]
  ra <- rel(a, "rA1", "rA2")
  # DoD reliability: its two "blocks" are d1, d2, so reliability_from_blocks(d1,d2) reproduces script 36
  rd <- if (nrow(d)) { dk <- d[keep == TRUE & is.finite(d1) & is.finite(d2)]
                       if (nrow(dk) >= 3L) { rr <- reliability_from_blocks(dk$d1, dk$d2, center = TRUE)
                                             data.table(h = rr$h, n = rr$n_pairs) }
                       else data.table(h = NA_real_, n = nrow(dk)) } else data.table(h = NA_real_, n = 0L)
  .(h_arm = ra$h, snr_arm = ra$snr, n_arm = ra$n, h_dod = rd$h, n_dod = rd$n)
}, by = .(stratum, hybrid, time)]

fwrite(res, file.path(hsd, "p1_reference_arm_reliability.tsv"), sep = "\t")
say("\nwrote p1_reference_arm_reliability.tsv (%d rows)", nrow(res))

## ---- reversal-regime recount: h_arm < 1/gamma, point estimate ----
E <- res[is.finite(h_arm) & h_arm > 0]
say("\n=== reference-arm reliability h_arm: distribution over estimable cells (n=%d) ===", nrow(E))
say("min %.4f  median %.4f  max %.4f", min(E$h_arm), median(E$h_arm), max(E$h_arm))
say("\n=== reversal-regime counts (POINT ESTIMATE), h_arm < 1/gamma ===")
say("%-8s %10s %14s %14s", "gamma", "1/gamma", "cells below", "of estimable")
for (g in c(1.1, 1.2, 1.3, 1.5, 2.0))
  say("%-8.1f %10.4f %14d %14d", g, 1/g, sum(E$h_arm < 1/g), nrow(E))

## ---- side-by-side: DoD-contrast reliability (script 36's quantity), explicitly labelled ----
Ed <- res[is.finite(h_dod) & h_dod > 0]
say("\n=== FOR COMPARISON: DoD-contrast reliability (script 36's quantity; NOT read as reversal) ===")
say("min %.4f  median %.4f  max %.4f  | cells h<1/1.1: %d of %d",
    min(Ed$h_dod), median(Ed$h_dod), max(Ed$h_dod), sum(Ed$h_dod < 1/1.1), nrow(Ed))

## ---- sanity anchor, a diagnostic only (it does not affect the estimates) ----
say("\n=== diagnostic anchor (reported either way) ===")
say("Under IV validity Table 1 implies reference-arm h ~ 0.76-0.94 per hybrid.")
say("median h_arm by hybrid:")
for (g in HYB) { x <- E[hybrid == g]; if (nrow(x)) say("  %-6s %.4f (n=%d)", g, median(x$h_arm), nrow(x)) }
say("=== done ===")
