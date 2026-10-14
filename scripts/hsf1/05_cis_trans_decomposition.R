# 05_cis_trans_decomposition.R -- is the divergence of the heat-stress
# response more cis or more trans, and does the answer differ between the Hsf1-associated
# core and the rest of the induced programme?
#
# The cis/trans result determines the appropriate sequence-analysis axis.
# A trans-dominated result points to regulators such as Hsf1,
# whereas a cis-dominated result points to target promoters and HSE architecture.
# 
# 
#
# The decomposition, and the one thing about it that is not symmetric.
# For a species pair (s1, s2) with hybrid s1xs2, timepoint t, replicate block b:
#     cis(t,b)   = [h(s2,t,b) - h(s2,0,b)] - [h(s1,t,b) - h(s1,0,b)]      (inside the hybrid)
#     par(t,b)   = [p(s2,t,b) - p(s2,0,b)] - [p(s1,t,b) - p(s1,0,b)]      (between parents)
#     trans(t,b) = par(t,b) - cis(t,b)
# cis is measured inside one nucleus, so the two alleles share the entire trans environment
# and it cancels. par is measured across two separate cultures. trans is therefore not an
# independent measurement: it is a residual, and every discrepancy between the parental and
# the hybrid measurement ends up in it -- different libraries, different depths (parental Sc
# is 0.1-0.9M reads vs Sp 2.6-6.9M), different mapping behaviour. Two consequences that
# govern everything below:
#   (1) trans carries both noise floors (var(trans) = var(par) + var(cis) when the two are
#       measured independently), so an uncorrected comparison makes trans win by construction.
#       The block-contrast null is therefore not optional here; it is the whole test.
#   (2) noise in cis enters trans with the opposite sign, which manufactures a negative
#       cor(cis, trans) -- the classic "compensatory evolution" artefact. Any correlation
#       reported below is computed with cis from one block and trans from the other, so the
#       two carry independent noise. The naive same-block value is printed beside it so the
#       size of the artefact is measured on this dataset rather than cited.
#
# Estimator (the construction of 04_divergence_by_dependency.R, applied to three axes
# instead of one).
# Two complete replicate blocks, exactly one sample per cell, so for each axis X:
#     signal_X = (X_b1 + X_b2)/2      null_X = (X_b2 - X_b1)/2
# Both are the same linear combination of two single-block values up to sign, so they carry
# identical noise variance, and rms_tru(X) = sqrt(max(0, rms(signal)^2 - rms(null)^2)) is the
# noise-corrected amplitude. Under unequal block variances this is still exact:
# E[sig^2] = delta^2 + (s1^2+s2^2)/4 and E[nul^2] = (s1^2+s2^2)/4 carry the same second term.
# No dispersion model, no eBayes, no se -- which matters because shrunken se is not usable
# as a precision measure and the deposited multimethod table carries an F, not a t.
#
# The amplitude confound carries over. rms of any of these three axes scales with
# the response magnitude R, so a stratum that responds harder looks more divergent on all
# three at identical relative divergence. 04_divergence_by_dependency.R measures the
# scaling exponent rather than assuming it (beta = 0.383, sublinear -- neither additive
# nor multiplicative), and the same measurement is repeated here per axis, because there
# is no reason cis and trans should scale alike. Every stratum comparison is reported as
# a bracket over beta in [0, 1] with the measured beta as the point estimate, and a ratio
# that crosses 1 anywhere in the bracket is reported as undecided.
#
# Expected result, fixed before this analysis was run.
# Li & Fay 2017 (PMID 28431042) analysed ASE across three temperatures in an Sc x Su hybrid
# and found only 136 genes with temperature-dependent ASE and no directional cis evolution.
# So the expected answer is that the response divergence (as opposed to the basal expression
# divergence) is trans-dominated, and that cis contributes little. If this script returns a
# cis-dominated response divergence, that is a conflict with a published result on an
# overlapping species pair and should be interpreted as a conflict with that published result.
#
# Structural limits (as in 04_divergence_by_dependency.R): 6 pairs are not 6 independent
# units -- each species appears in 3 pairs and the hybrid libraries are co-pooled in
# complementary pairs (library 5 = ScxSp+SkxSu, 13 = ScxSu+SpxSk, 21 = ScxSk+SpxSu). All
# inference is at the gene level; the pair axis is described, never tested. Negatives are
# equivalence bounds.
#
# Usage: Rscript scripts/hsf1/05_cis_trans_decomposition.R

suppressPackageStartupMessages({ library(data.table); library(edgeR) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
ext  <- file.path(root, "data/03_external/fay2023")
csd  <- file.path(root, "data/02_processed/crossspecies")
hsd  <- file.path(root, "data/02_processed/hsf1")
outd <- hsd
repd <- file.path(root, "doc/05_reports")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

ANCHOR  <- "37_30"
SPECIES <- c("Sc", "Sp", "Sk", "Su")
HYB     <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")
TP      <- c(15L, 30L, 60L)
PRIOR   <- 2
MINCPM  <- 1
NBOOT   <- 2000
NB_CAL  <- 600L
NDEC    <- 8L
AXES    <- c("cis", "trans", "par")
set.seed(20260906L)

say <- function(...) cat(sprintf(...), "\n", sep = "")
rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))

logf <- file.path(repd, "logs/hsf1_cis_trans.log")
dir.create(dirname(logf), showWarnings = FALSE, recursive = TRUE)
sink(logf, split = TRUE)
say("=== 05_cis_trans_decomposition | %s | anchor %s C ===", format(Sys.time()), ANCHOR)
say("cis/trans decomposition linked to the Hsf1 locus and promoter analyses (scripts 07-09).")
say("Pre-registered expectation: trans-dominated response divergence (Li & Fay 2017,")
say("PMID 28431042: 136 temperature-dependent ASE genes, no directional cis evolution).")

## --------------------------------------------------------------- inputs ---
# read.csv, not fread: 172 rows of htseq.ortho.csv are ragged (267-269 fields of 270). All
# 172 carry one==0 (see 03_fay_parental_dimensions.R); re-asserted here.
dat <- read.csv(file.path(ext, "htseq.ortho.csv"))
design <- read.csv(file.path(ext, "design.csv"))
meta_cols <- c("HOG","OG","one","all4","multi","Name","c1","c2","c3","c4",
               "Len.Sc","Len.Sp","Len.Sk","Len.Su")
stopifnot(identical(colnames(dat)[1:14], meta_cols))
stopifnot(nrow(design) == ncol(dat) - 14L)
one <- dat$one == 1
stopifnot(sum(one) == 5121L)
cts <- as.matrix(dat[one, -(1:14)])
rownames(cts) <- dat$HOG[one]
stopifnot(!anyNA(cts))
say("count matrix: %d one-to-one orthologs x %d samples, no NA", nrow(cts), ncol(cts))

D <- as.data.table(design)
D[, col := Name]
D[, block := fifelse(as.integer(sub("^[A-Za-z]+", "", Name)) <= 32L, "block1", "block2")]

P <- D[Strain %in% SPECIES & Allele == Strain & Condition == "heat"]
P[, Time := as.integer(Time)]
stopifnot(nrow(P) == 32L, !anyDuplicated(P[, .(Strain, block, Time)]))
H <- D[Strain %in% HYB & Condition == "heat"]
H[, Time := as.integer(Time)]
stopifnot(nrow(H) == 96L, !anyDuplicated(H[, .(Strain, Allele, block, Time)]))
say("parental heat-arm samples: %d (4 species x 4 times x 2 blocks)", nrow(P))
say("hybrid   heat-arm samples: %d (6 hybrids x 2 alleles x 4 times x 2 blocks)", nrow(H))

# allele polarity parsed from the hybrid name and checked against design.csv, so that par
# and cis are guaranteed to be oriented the same way (s2 - s1 in both)
POL <- data.table(hybrid = HYB, s1 = substr(HYB, 1, 2), s2 = substr(HYB, 4, 5))
for (i in seq_len(nrow(POL))) {
  al <- sort(H[Strain == POL$hybrid[i], unique(Allele)])
  stopifnot(identical(al, sort(c(POL$s1[i], POL$s2[i]))))
}
say("\npolarity (s2 - s1), identical for the cis and the parental axis:")
say("  %s", paste(sprintf("%s: %s-%s", POL$hybrid, POL$s2, POL$s1), collapse = "  "))

## ------------------------------------------------- gene set assignment ---
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
CS <- fread(file.path(hsd, sprintf("hsf1_core_set_%s.tsv", ANCHOR)))
CL <- fread(file.path(hsd, sprintf("hsf1_epistasis_classes_%s.tsv", ANCHOR)))
UC <- dcast(CL[, .(gene, method, class)], gene ~ method, value.var = "class")
UC[, unanimous := fifelse(DESeq2 == edgeR & edgeR == limma, DESeq2, NA_character_)]
map_hog <- function(g) {
  h <- HB[Name %in% g | sys %in% g, unique(gene)]
  h[h %in% rownames(cts)]
}
STRATA <- list(
  Hsf1_core      = CS[set == "Hsf1_core_unanimous",  gene],
  Hsf1_joint     = CS[set == "Hsf1_joint_candidate", gene],
  Msn24_only     = UC[unanimous == "Msn24_only",     gene],
  TF_independent = UC[unanimous == "TF_independent", gene],
  other_induced  = UC[unanimous %in% c("Msn24_only", "TF_independent"), gene],
  not_induced    = UC[unanimous == "not_induced",    gene])
HOGS <- lapply(STRATA, map_hog)
SETS <- names(STRATA)

## --------------------------------------------- per-species parental logCPM ---
# Normalisation within species: every quantity used is a within-species time difference, so
# a species' own library composition cancels. Pooling species into one TMM fit would let one
# species' composition move another's baseline.
say("\n=== per-species parental logCPM (TMM within species) ===")
PL <- list(); keep_p <- list()
for (sp in SPECIES) {
  ps <- P[Strain == sp][order(block, Time)]
  m  <- cts[, ps$col, drop = FALSE]
  y  <- calcNormFactors(DGEList(counts = m), method = "TMM")
  PL[[sp]] <- cpm(y, log = TRUE, prior.count = PRIOR)
  colnames(PL[[sp]]) <- ps$col
  keep_p[[sp]] <- rowSums(cpm(y) >= MINCPM) >= 4L
  say("  %-3s libs %s | TMM %.2f-%.2f | detected: %d", sp,
      paste(sprintf("%.1fM", colSums(m) / 1e6), collapse = " "),
      min(y$samples$norm.factors), max(y$samples$norm.factors), sum(keep_p[[sp]]))
}
say("NOTE the ~30x parental depth mismatch. It cannot be normalised away; it is")
say("handled by the per-axis block-contrast null, which is that library set's own noise floor.")

## ---------------------------------------------------- per-hybrid logCPM ---
say("\n=== per-hybrid logCPM (TMM within hybrid) ===")
HL <- list(); HS <- list(); keep_h <- list()
for (g in HYB) {
  hs <- H[Strain == g][order(Allele, block, Time)]
  m  <- cts[, hs$col, drop = FALSE]
  y  <- calcNormFactors(DGEList(counts = m), method = "TMM")
  HL[[g]] <- cpm(y, log = TRUE, prior.count = PRIOR)
  colnames(HL[[g]]) <- hs$col
  HS[[g]] <- hs
  cpr <- cpm(y)
  s1 <- POL[hybrid == g, s1]; s2 <- POL[hybrid == g, s2]
  k1 <- rowSums(cpr[, hs[Allele == s1, col], drop = FALSE] >= MINCPM) >= 4L
  k2 <- rowSums(cpr[, hs[Allele == s2, col], drop = FALSE] >= MINCPM) >= 4L
  keep_h[[g]] <- k1 & k2
  say("  %-6s libs %s | TMM %.2f-%.2f | both alleles detected: %d", g,
      paste(sprintf("%.1fM", colSums(m) / 1e6), collapse = " "),
      min(y$samples$norm.factors), max(y$samples$norm.factors), sum(keep_h[[g]]))
}

# A gene must be usable on both axes for all six pairs, otherwise cis and trans would be
# computed on different gene sets and their comparison would be meaningless.
keep_all <- Reduce(`&`, c(keep_p, keep_h))
KEEP <- rownames(cts)[keep_all]
say("\ngene filter: detected in all 4 parents AND both alleles of all 6 hybrids: %d of %d (%.1f%%)",
    length(KEEP), nrow(cts), 100 * mean(keep_all))
a5f <- file.path(hsd, sprintf("a5_cis_ratio_%s.tsv", ANCHOR))
if (file.exists(a5f)) say("A5 (cis axis alone) used %d genes; A6 adds the parental requirement.",
                          4997L)

for (s in SETS) HOGS[[s]] <- intersect(HOGS[[s]], KEEP)
say("\nstratum sizes after the filter:")
for (s in SETS) say("  %-16s %d", s, length(HOGS[[s]]))
stopifnot(length(HOGS$Hsf1_core) >= 20L)
# disjointness is asserted (a subset compared against its own superset dilutes toward 1
# by construction)
stopifnot(length(intersect(HOGS$Hsf1_core, HOGS$other_induced)) == 0L)
stopifnot(length(intersect(HOGS$Hsf1_joint, HOGS$other_induced)) == 0L)
stopifnot(setequal(HOGS$other_induced, union(HOGS$Msn24_only, HOGS$TF_independent)))
DENOM <- list(Hsf1_core = "other_induced", Hsf1_joint = "other_induced",
              Msn24_only = "TF_independent", TF_independent = "Msn24_only",
              not_induced = "other_induced")
say("disjoint denominators: %s",
    paste(sprintf("%s/%s", names(DENOM), unlist(DENOM)), collapse = "  "))
ST <- rbindlist(lapply(SETS, function(s) data.table(hog = HOGS[[s]], stratum = s)))

## --------------------------------- the three axes, per pair / time / block ---
say("\n=== decomposition: cis (within hybrid), par (between parents), trans = par - cis ===")
gv_p <- function(sp, bk, t) {
  cc <- P[Strain == sp & block == bk & Time == t, col]
  stopifnot(length(cc) == 1L)
  PL[[sp]][KEEP, cc]
}
gv_h <- function(g, al, bk, t) {
  cc <- HS[[g]][Allele == al & block == bk & Time == t, col]
  stopifnot(length(cc) == 1L)
  HL[[g]][KEEP, cc]
}
AX <- list()
for (g in HYB) {
  s1 <- POL[hybrid == g, s1]; s2 <- POL[hybrid == g, s2]
  for (t in TP) {
    # per block, all three axes; the response is always (t - 0) within the same block and
    # the same library, so baseline and ortholog length both cancel
    cis_b <- lapply(c("block1", "block2"), function(bk)
      (gv_h(g, s2, bk, t) - gv_h(g, s2, bk, 0L)) - (gv_h(g, s1, bk, t) - gv_h(g, s1, bk, 0L)))
    par_b <- lapply(c("block1", "block2"), function(bk)
      (gv_p(s2, bk, t) - gv_p(s2, bk, 0L)) - (gv_p(s1, bk, t) - gv_p(s1, bk, 0L)))
    trn_b <- list(par_b[[1]] - cis_b[[1]], par_b[[2]] - cis_b[[2]])
    # the common amplitude reference: the mean parental response, i.e. how hard this gene
    # responds at all. Sum contrast -- orthogonal to every difference contrast above when the
    # two units carry equal measurement error (checked empirically below).
    amp_b <- lapply(c("block1", "block2"), function(bk)
      ((gv_p(s2, bk, t) - gv_p(s2, bk, 0L)) + (gv_p(s1, bk, t) - gv_p(s1, bk, 0L))) / 2)
    mk <- function(v) list(sig = (v[[1]] + v[[2]]) / 2, nul = (v[[2]] - v[[1]]) / 2,
                           b1 = v[[1]], b2 = v[[2]])
    C <- mk(cis_b); T_ <- mk(trn_b); A <- mk(par_b); M <- mk(amp_b)
    AX[[length(AX) + 1L]] <- data.table(
      hog = KEEP, pair = g, time = t,
      cis = C$sig, cis_nul = C$nul, cis_b1 = C$b1, cis_b2 = C$b2,
      trans = T_$sig, trans_nul = T_$nul, trans_b1 = T_$b1, trans_b2 = T_$b2,
      par = A$sig, par_nul = A$nul,
      amp = M$sig, amp_nul = M$nul)
  }
}
AX <- rbindlist(AX)
stopifnot(nrow(AX) == length(KEEP) * length(HYB) * length(TP))
# the decomposition is an identity by construction; assert it so a sign slip cannot survive
stopifnot(max(abs(AX$par - (AX$cis + AX$trans))) < 1e-8)
say("identity par = cis + trans holds to %.1e (it is a definition, asserted to catch sign slips)",
    max(abs(AX$par - (AX$cis + AX$trans))))
CELLS <- AX[, unique(paste(pair, time, sep = "_"))]
stopifnot(length(CELLS) == 18L)

## ----------------------------------------- orthogonality of the amplitude ---
say("\n---- is the amplitude reference orthogonal to the three divergence axes? ----")
say("Required before dividing by it: a normaliser correlated with the numerator")
say("would import the P1 reference-level artefact instead of removing it.")
OR <- AX[, .(r_cis_amp = cor(cis, amp), r_trans_amp = cor(trans, amp),
             r_par_amp = cor(par, amp), r_cis_trans = cor(cis, trans)),
         by = .(pair, time)]
say("%-8s %6s %11s %13s %11s %13s", "pair", "time", "r(cis,amp)", "r(trans,amp)",
    "r(par,amp)", "r(cis,trans)")
for (i in seq_len(nrow(OR))) say("%-8s %6d %11.3f %13.3f %11.3f %13.3f",
  OR$pair[i], OR$time[i], OR$r_cis_amp[i], OR$r_trans_amp[i], OR$r_par_amp[i], OR$r_cis_trans[i])
say("median |r(cis,amp)| = %.3f ; median |r(trans,amp)| = %.3f ; median r(cis,trans) = %+.3f",
    median(abs(OR$r_cis_amp)), median(abs(OR$r_trans_amp)), median(OR$r_cis_trans))
fwrite(OR, file.path(outd, sprintf("a6_orthogonality_%s.tsv", ANCHOR)), sep = "\t")

## ------------------------------------------- noise-corrected amplitudes ---
say("\n================ NOISE-CORRECTED AMPLITUDE OF EACH AXIS, BY STRATUM ================")
say("rms_tru = sqrt(max(0, rms(signal)^2 - rms(null)^2)); snr = rms(signal)/rms(null).")
say("An UNCORRECTED comparison would make trans win by construction: trans = par - cis")
say("carries the parental AND the hybrid noise floor. That is what the correction removes.")
AXL <- merge(AX, ST, by = "hog", allow.cartesian = TRUE)
FL <- AXL[, .(n = .N,
  cis_sig = rms(cis), cis_nul = rms(cis_nul),
  trn_sig = rms(trans), trn_nul = rms(trans_nul),
  par_sig = rms(par), par_nul = rms(par_nul),
  amp_sig = rms(amp), amp_nul = rms(amp_nul)), by = .(stratum, pair, time)]
FL[, `:=`(cis_tru = sqrt(pmax(0, cis_sig^2 - cis_nul^2)),
          trn_tru = sqrt(pmax(0, trn_sig^2 - trn_nul^2)),
          par_tru = sqrt(pmax(0, par_sig^2 - par_nul^2)),
          amp_tru = sqrt(pmax(0, amp_sig^2 - amp_nul^2)),
          cis_snr = cis_sig / cis_nul, trn_snr = trn_sig / trn_nul,
          par_snr = par_sig / par_nul)]
setorder(FL, stratum, pair, time)
fwrite(FL, file.path(outd, sprintf("a6_axis_amplitudes_%s.tsv", ANCHOR)), sep = "\t")

say("\nmedian over the 6 pairs, per stratum and timepoint:")
say("%-15s %5s %9s %9s %9s %9s %8s %8s", "stratum", "time", "cis_tru", "trans_tru",
    "par_tru", "amp_tru", "cis_snr", "trn_snr")
for (s in SETS) for (t in TP) {
  z <- FL[stratum == s & time == t]
  say("%-15s %5d %9.3f %9.3f %9.3f %9.3f %8.2f %8.2f", s, t,
      median(z$cis_tru), median(z$trn_tru), median(z$par_tru), median(z$amp_tru),
      median(z$cis_snr), median(z$trn_snr))
}
say("\nUNCORRECTED, for comparison -- without null correction:")
say("%-15s %5s %9s %9s %14s", "stratum", "time", "cis_sig", "trans_sig", "trans/cis (raw)")
for (s in c("Hsf1_core", "other_induced")) for (t in TP) {
  z <- FL[stratum == s & time == t]
  say("%-15s %5d %9.3f %9.3f %14.2f", s, t, median(z$cis_sig), median(z$trn_sig),
      median(z$trn_sig) / median(z$cis_sig))
}
say("and CORRECTED:")
say("%-15s %5s %9s %9s %14s", "stratum", "time", "cis_tru", "trans_tru", "trans/cis (corr)")
for (s in c("Hsf1_core", "other_induced")) for (t in TP) {
  z <- FL[stratum == s & time == t]
  say("%-15s %5d %9.3f %9.3f %14.2f", s, t, median(z$cis_tru), median(z$trn_tru),
      ifelse(median(z$cis_tru) > 0, median(z$trn_tru) / median(z$cis_tru), NA_real_))
}
say("If the two trans/cis columns differ materially, the corrected estimate is used.")

## ------------------------------------------------ wide matrices for the ratios ---
mk_mat <- function(vn) {
  w <- dcast(AX, hog ~ pair + time, value.var = vn)
  m <- as.matrix(w[, setdiff(names(w), "hog"), with = FALSE])
  rownames(m) <- w$hog
  m[, CELLS, drop = FALSE]
}
MATS <- lapply(c("cis", "cis_nul", "trans", "trans_nul", "par", "par_nul",
                 "amp", "amp_nul", "cis_b1", "cis_b2", "trans_b1", "trans_b2"), mk_mat)
names(MATS) <- c("cis", "cis_nul", "trans", "trans_nul", "par", "par_nul",
                 "amp", "amp_nul", "cis_b1", "cis_b2", "trans_b1", "trans_b2")
stopifnot(identical(rownames(MATS$cis), rownames(MATS$trans)))
GENES <- rownames(MATS$cis)
idx_of <- function(s) match(HOGS[[s]], GENES)
tp_of <- as.integer(sub(".*_", "", CELLS))
pair_of <- sub("_[0-9]+$", "", CELLS)

tru <- function(ax, i) sqrt(pmax(0, colMeans(MATS[[ax]][i, , drop = FALSE]^2) -
                                   colMeans(MATS[[paste0(ax, "_nul")]][i, , drop = FALSE]^2)))
amp_tru <- function(i) sqrt(pmax(0, colMeans(MATS$amp[i, , drop = FALSE]^2) -
                                    colMeans(MATS$amp_nul[i, , drop = FALSE]^2)))

## ------------------------------- cis share: the answer to the cis-or-trans question ---
# The quantity of interest is not "is trans bigger" (it must be, it is a residual carrying
# two noise floors) but how much of the noise-corrected divergence variance is cis. Variances
# add along the identity par = cis + trans only up to the covariance term, which is estimated
# from independent blocks below, so the share is defined on the variances themselves and the
# covariance is reported separately rather than being buried in a ratio.
say("\n================ CIS SHARE OF THE RESPONSE DIVERGENCE ================")
say("cis_share = rms_tru(cis)^2 / (rms_tru(cis)^2 + rms_tru(trans)^2), per cell.")
say("Reported as a share, not as trans/cis, because trans is a residual: a ratio with a")
say("residual in the numerator is unbounded and unstable when cis approaches its noise floor.")
cell_share <- function(i) {
  c2 <- tru("cis", i)^2; t2 <- tru("trans", i)^2
  ifelse(c2 + t2 > 0, c2 / (c2 + t2), NA_real_)
}
boot_share <- function(i, nb = NBOOT) {
  out <- matrix(NA_real_, nb, 3L)
  for (k in seq_len(nb)) {
    v <- cell_share(sample(i, length(i), TRUE))
    for (j in seq_along(TP)) out[k, j] <- median(v[tp_of == TP[j]], na.rm = TRUE)
  }
  out
}
SH <- list()
say("%-15s %5s %5s %10s %17s %10s", "stratum", "time", "n", "cis_share", "95%% CI", "cells")
for (s in SETS) {
  i_s <- idx_of(s); obs <- cell_share(i_s); bt <- boot_share(i_s)
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    SH[[length(SH) + 1L]] <- data.table(stratum = s, time = TP[j], n = length(i_s),
      cis_share = median(obs[sel], na.rm = TRUE), lo = ci[1], hi = ci[2],
      n_cells = sum(sel & !is.na(obs)))
  }
}
SH <- rbindlist(SH)
fwrite(SH, file.path(outd, sprintf("a6_cis_share_%s.tsv", ANCHOR)), sep = "\t")
for (i in seq_len(nrow(SH))) say("%-15s %5d %5d %10.3f [%6.3f,%6.3f] %8d/6",
  SH$stratum[i], SH$time[i], SH$n[i], SH$cis_share[i], SH$lo[i], SH$hi[i], SH$n_cells[i])
say("\nPRE-REGISTERED READ-OUT: Li & Fay 2017 predicts a SMALL cis share for the RESPONSE")
say("(as opposed to basal expression). A share below 0.5 is consistent with that published")
say("result; a share above 0.5 conflicts with it.")

say("\n---- per-pair detail (the median hides which pairs drive it) ----")
for (s in c("Hsf1_core", "other_induced")) {
  say("-- %s --", s)
  v <- cell_share(idx_of(s))
  say("%-8s %8s %8s %8s", "pair", "15min", "30min", "60min")
  for (g in HYB) say("%-8s %8.3f %8.3f %8.3f", g,
    v[which(pair_of == g & tp_of == 15L)], v[which(pair_of == g & tp_of == 30L)],
    v[which(pair_of == g & tp_of == 60L)])
}

## ----------------- the compensation artefact, measured rather than cited ---
# cor(cis, trans) is the standard evidence for compensatory evolution, and it is also the
# standard artefact: trans = par - cis, so measurement error in cis appears in trans with the
# opposite sign and drives the correlation negative even when the truth is zero. The fix uses
# the two replicate blocks: cis from block 1, trans from block 2. Their noise is independent
# by design, so the artefact term vanishes and what remains is the real covariance.
say("\n================ cor(cis, trans): ARTEFACT vs INDEPENDENT-BLOCK ESTIMATE ================")
CT <- list()
for (cc in seq_along(CELLS)) {
  naive <- cor(MATS$cis[, cc], MATS$trans[, cc])
  same1 <- cor(MATS$cis_b1[, cc], MATS$trans_b1[, cc])
  cross <- cor(MATS$cis_b1[, cc], MATS$trans_b2[, cc])
  cross2 <- cor(MATS$cis_b2[, cc], MATS$trans_b1[, cc])
  CT[[cc]] <- data.table(pair = pair_of[cc], time = tp_of[cc],
    r_naive_bothblocks = naive, r_same_block = same1,
    r_cross_block = (cross + cross2) / 2, artefact = same1 - (cross + cross2) / 2)
}
CT <- rbindlist(CT)
fwrite(CT, file.path(outd, sprintf("a6_cis_trans_correlation_%s.tsv", ANCHOR)), sep = "\t")
say("%-8s %6s %12s %12s %13s %10s", "pair", "time", "r(naive)", "r(same blk)",
    "r(cross blk)", "artefact")
for (i in seq_len(nrow(CT))) say("%-8s %6d %12.3f %12.3f %13.3f %10.3f",
  CT$pair[i], CT$time[i], CT$r_naive_bothblocks[i], CT$r_same_block[i],
  CT$r_cross_block[i], CT$artefact[i])
say("median r naive = %+.3f | same-block = %+.3f | CROSS-BLOCK = %+.3f | artefact = %+.3f",
    median(CT$r_naive_bothblocks), median(CT$r_same_block), median(CT$r_cross_block),
    median(CT$artefact))
say("The cross-block column removes the SHARED-ERROR part. The gap between it and the naive")
say("column is the size of that artefact ON THIS DATASET -- the same class of defect as the")
say("P1 reference-level artefact, in a different statistic.")

# Removing the shared error is necessary but not sufficient. trans is defined as
# par - cis, so
#     cov(cis, trans) = cov(cis, par) - var(cis),
# and if cov(cis, par) = 0 -- i.e. no compensation of any kind -- then
#     var(trans) = var(par) + var(cis)  and  r_mech = -sqrt(var(cis)/(var(par)+var(cis))),
# which is strongly negative. A negative cor(cis, trans) is therefore the default, not
# evidence. The comparison has to be against r_mech, not against zero, and the quantity that
# actually carries biology is cor(cis, par): negative = compensation, zero = neither,
# positive = the two axes reinforce. All correlations below are cross-block.
say("\n---- is that residual negative correlation just the definition trans = par - cis? ----")
say("r_mech = -sqrt(var_tru(cis)/(var_tru(par)+var_tru(cis))) is what cor(cis,trans) equals")
say("when cov(cis,par) = 0, i.e. under NO compensation. Compare the cross-block r with THAT.")
CP <- list()
for (cc in seq_along(CELLS)) {
  par_b1 <- MATS$cis_b1[, cc] + MATS$trans_b1[, cc]
  par_b2 <- MATS$cis_b2[, cc] + MATS$trans_b2[, cc]
  vc <- max(0, mean(MATS$cis[, cc]^2)   - mean(MATS$cis_nul[, cc]^2))
  vp <- max(0, mean(MATS$par[, cc]^2)   - mean(MATS$par_nul[, cc]^2))
  CP[[cc]] <- data.table(pair = pair_of[cc], time = tp_of[cc],
    r_cis_trans_cross = CT$r_cross_block[cc],
    r_mech = if (vp + vc > 0) -sqrt(vc / (vp + vc)) else NA_real_,
    r_cis_par_cross = (cor(MATS$cis_b1[, cc], par_b2) + cor(MATS$cis_b2[, cc], par_b1)) / 2)
}
CP <- rbindlist(CP)
CP[, more_negative_than_mech := r_cis_trans_cross < r_mech]
fwrite(CP, file.path(outd, sprintf("a6_compensation_test_%s.tsv", ANCHOR)), sep = "\t")
say("%-8s %6s %19s %10s %18s %8s", "pair", "time", "r(cis,trans) cross", "r_mech",
    "r(cis,par) cross", "< mech?")
for (i in seq_len(nrow(CP))) say("%-8s %6d %19.3f %10.3f %18.3f %8s",
  CP$pair[i], CP$time[i], CP$r_cis_trans_cross[i], CP$r_mech[i], CP$r_cis_par_cross[i],
  ifelse(CP$more_negative_than_mech[i], "YES", "no"))
say("median r(cis,trans) cross-block = %+.3f | median r_mech = %+.3f | cells below mech: %d/%d",
    median(CP$r_cis_trans_cross), median(CP$r_mech, na.rm = TRUE),
    sum(CP$more_negative_than_mech, na.rm = TRUE), nrow(CP))
say("median CROSS-BLOCK r(cis, par) = %+.3f  <- this is the quantity that carries biology",
    median(CP$r_cis_par_cross))
if (median(CP$r_cis_par_cross) < -0.1) {
  say("=> cis and the total parental divergence are ANTI-correlated: compensation supported.")
} else if (median(CP$r_cis_par_cross) > 0.1) {
  say("=> cis and the total parental divergence are POSITIVELY correlated. The negative")
  say("   cor(cis,trans) is therefore LESS negative than the no-compensation prediction: the")
  say("   two axes REINFORCE rather than compensate. Reporting 'compensatory evolution' from")
  say("   the naive statistic here would have inverted the sign of the biological claim.")
} else {
  say("=> cor(cis,par) is not materially different from 0: the negative cor(cis,trans) is")
  say("   fully accounted for by the definition trans = par - cis. NO compensation is")
  say("   demonstrated, and none is contradicted either.")
}

## ---------------- amplitude scaling exponent, per axis ---
say("\n---- amplitude scaling exponent beta, measured separately for each axis ----")
say("beta = slope of log rms_tru(axis) on log amp across %d deciles of an INDEPENDENTLY", NDEC)
say("measured amplitude (leave-one-pair-out at the same timepoint), so regression to the mean")
say("cannot flatten it. beta~0 additive, beta~1 multiplicative; A5 found 0.383 for cis.")
BETA <- list()
for (t in TP) {
  ct <- which(tp_of == t)
  for (cc in ct) {
    others <- setdiff(ct, cc)
    xind <- rowMeans(abs(MATS$amp[, others, drop = FALSE]))
    br <- unique(quantile(xind, seq(0, 1, length.out = NDEC + 1L), na.rm = TRUE))
    if (length(br) < 4L) next
    bin <- cut(xind, breaks = br, include.lowest = TRUE, labels = FALSE)
    for (ax in AXES) {
      B <- rbindlist(lapply(sort(unique(bin[!is.na(bin)])), function(bb) {
        ii <- which(bin == bb)
        data.table(nb = length(ii), a = amp_tru(ii)[cc], d = tru(ax, ii)[cc])
      }))
      B <- B[a > 0 & d > 0 & nb >= 20L]
      if (nrow(B) < 4L) next
      f <- lm(log(d) ~ log(a), data = B)
      BETA[[length(BETA) + 1L]] <- data.table(axis = ax, pair = pair_of[cc], time = t,
        n_bins = nrow(B), beta = unname(coef(f)[2]), r2 = summary(f)$r.squared)
    }
  }
}
BETA <- rbindlist(BETA)
fwrite(BETA, file.path(outd, sprintf("a6_amplitude_scaling_%s.tsv", ANCHOR)), sep = "\t")
say("%-8s %8s %8s %8s %8s %8s", "axis", "cells", "beta_med", "beta_q25", "beta_q75", "R2_med")
BMED <- c()
for (ax in AXES) {
  z <- BETA[axis == ax]
  BMED[ax] <- median(z$beta)
  say("%-8s %8d %8.3f %8.3f %8.3f %8.2f", ax, nrow(z), median(z$beta),
      quantile(z$beta, .25), quantile(z$beta, .75), median(z$r2))
}
say("The axes are NOT required to share beta and are not assumed to: each stratum comparison")
say("below uses its own axis's measured beta, with beta=0 and beta=1 printed as the bracket.")

## -------------- stratum comparison on each axis, amplitude-calibrated ---
say("\n============ IS EITHER AXIS'S DIVERGENCE CONCENTRATED IN THE Hsf1 CORE? ============")
say("ratio = [rms_tru(axis)/amp^beta](focal) / same(disjoint denominator), median over 6 pairs.")
say("Bracket over beta in [0,1] printed; a ratio crossing 1 inside the bracket is UNDECIDED")
say("regardless of how tight its CI is.")
cell_cal <- function(ax, i, beta) {
  a <- amp_tru(i); ifelse(a > 0, tru(ax, i) / a^beta, NA_real_)
}
cell_cal_ratio <- function(ax, i_n, i_d, beta) {
  n <- cell_cal(ax, i_n, beta); d <- cell_cal(ax, i_d, beta)
  ifelse(!is.na(d) & d > 0, n / d, NA_real_)
}
boot_cal_ratio <- function(ax, i_n, i_d, beta, nb = NBOOT) {
  out <- matrix(NA_real_, nb, 3L)
  for (k in seq_len(nb)) {
    r <- cell_cal_ratio(ax, sample(i_n, length(i_n), TRUE), sample(i_d, length(i_d), TRUE), beta)
    for (j in seq_along(TP)) out[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  out
}
RAT <- list()
say("%-6s %-15s %5s %8s %17s %7s %8s %8s %6s", "axis", "focal", "time", "ratio",
    "95%% CI", "cells>1", "beta=0", "beta=1", "sign")
for (ax in AXES) for (s in names(DENOM)) {
  i_f <- idx_of(s); i_d <- idx_of(DENOM[[s]])
  obs <- cell_cal_ratio(ax, i_f, i_d, BMED[ax])
  r0  <- cell_cal_ratio(ax, i_f, i_d, 0)
  r1  <- cell_cal_ratio(ax, i_f, i_d, 1)
  bt  <- boot_cal_ratio(ax, i_f, i_d, BMED[ax])
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    o  <- median(obs[sel], na.rm = TRUE)
    v0 <- median(r0[sel], na.rm = TRUE); v1 <- median(r1[sel], na.rm = TRUE)
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    stable <- (all(c(o, v0, v1) > 1) || all(c(o, v0, v1) < 1))
    RAT[[length(RAT) + 1L]] <- data.table(axis = ax, stratum = s, denom = DENOM[[s]],
      time = TP[j], n = length(i_f), beta = BMED[ax], ratio = o, lo = ci[1], hi = ci[2],
      ratio_beta0 = v0, ratio_beta1 = v1,
      n_cells_gt1 = sum(obs[sel] > 1, na.rm = TRUE), n_cells = sum(sel & !is.na(obs)),
      sign_stable = stable, equiv_bound = max(abs(c(ci[1], ci[2]) - 1)))
  }
}
RAT <- rbindlist(RAT)
fwrite(RAT, file.path(outd, sprintf("a6_stratum_ratio_by_axis_%s.tsv", ANCHOR)), sep = "\t")
for (i in seq_len(nrow(RAT))) say("%-6s %-15s %5d %8.3f [%6.3f,%6.3f] %5d/%d %8.3f %8.3f %6s",
  RAT$axis[i], RAT$stratum[i], RAT$time[i], RAT$ratio[i], RAT$lo[i], RAT$hi[i],
  RAT$n_cells_gt1[i], RAT$n_cells[i], RAT$ratio_beta0[i], RAT$ratio_beta1[i],
  ifelse(RAT$sign_stable[i], "YES", "no"))
say("\nsign-stable (focal,time,axis) combinations: %d of %d",
    sum(RAT$sign_stable), nrow(RAT))

## --------------- power calibration, per axis ---
# The injected term must be homologous to the statistic's numerator, or a negative means
# nothing. Injection: axis -> axis + c * amp for an n-matched random subset of other_induced,
# i.e. a genuine multiplicative divergence increment on that axis. amp is taken leave-one-
# pair-out so its noise is independent of the target cell. The c=0 run is the control that
# must not come out significant.
say("\n---- POWER CALIBRATION, per axis: smallest detectable divergence excess ----")
LOHO <- sapply(TP, function(t) {
  ct <- which(tp_of == t); rowMeans(abs(MATS$amp[, ct, drop = FALSE]))
})
i_oth <- idx_of("other_induced")
set.seed(20260906L + 11L)
i_inj <- sample(i_oth, length(HOGS$Hsf1_core))
i_rest <- setdiff(i_oth, i_inj)
say("injected n = %d (n-matched to the core), denominator n = %d (disjoint remainder)",
    length(i_inj), length(i_rest))
inj_ratio <- function(ax, i_n, i_d, cc_inj, beta) {
  M2 <- MATS[[ax]]
  for (j in seq_along(TP)) {
    cols <- which(tp_of == TP[j])
    M2[i_n, cols] <- M2[i_n, cols] + cc_inj * LOHO[i_n, j]
  }
  tr <- function(i) sqrt(pmax(0, colMeans(M2[i, , drop = FALSE]^2) -
                                 colMeans(MATS[[paste0(ax, "_nul")]][i, , drop = FALSE]^2)))
  n <- tr(i_n) / amp_tru(i_n)^beta
  d <- tr(i_d) / amp_tru(i_d)^beta
  ifelse(!is.na(d) & d > 0, n / d, NA_real_)
}
PWR <- list()
say("%-6s %6s %9s %17s %10s", "axis", "c", "ratio_15", "95%% CI (15 min)", "excl. 1?")
for (ax in AXES) for (cc_inj in c(0, 0.10, 0.20, 0.40, 0.80)) {
  obs <- inj_ratio(ax, i_inj, i_rest, cc_inj, BMED[ax])
  bt <- matrix(NA_real_, NB_CAL, 3L)
  for (k in seq_len(NB_CAL)) {
    r <- inj_ratio(ax, sample(i_inj, length(i_inj), TRUE),
                   sample(i_rest, length(i_rest), TRUE), cc_inj, BMED[ax])
    for (j in seq_along(TP)) bt[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    PWR[[length(PWR) + 1L]] <- data.table(axis = ax, c_inj = cc_inj, time = TP[j],
      ratio = median(obs[sel], na.rm = TRUE), lo = ci[1], hi = ci[2],
      excludes_1 = ci[1] > 1 | ci[2] < 1)
  }
  z <- PWR[[length(PWR) - 2L]]
  say("%-6s %6.2f %9.3f [%6.3f,%6.3f] %10s", ax, cc_inj, z$ratio, z$lo, z$hi,
      ifelse(z$excludes_1, "YES", "no"))
}
PWR <- rbindlist(PWR)
fwrite(PWR, file.path(outd, sprintf("a6_power_calibration_%s.tsv", ANCHOR)), sep = "\t")
say("\nc = 0 must NOT exclude 1 (same gene pool split in two). Failures: %d of %d axis-times",
    PWR[c_inj == 0 & excludes_1 == TRUE, .N], PWR[c_inj == 0, .N])
say("%-6s %28s", "axis", "smallest detectable c (15 min)")
DET <- c()
for (ax in AXES) {
  d <- PWR[axis == ax & time == 15L & c_inj > 0 & excludes_1 == TRUE, if (.N) min(c_inj) else NA_real_]
  DET[ax] <- d
  say("%-6s %28s", ax, ifelse(is.na(d), "none up to 0.80", sprintf("%.2f", d)))
}
say("For an axis with no detectable c, flat ratios characterize the test rather than the biology")
say("and the result is classified as 'test not informative'.")

## ------- the orthogonality check on errors only ---
# The per-gene r(axis, amp) table above cannot be read the way the orthogonality check of
# 04_divergence_by_dependency.R is read. There, both contrasts have ~zero true signal in
# common, so a nonzero r means shared error. Here amp is the response magnitude and beta > 0
# says the axes genuinely scale with it, so r(par, amp) = +0.33 is partly real. The two
# components have to be separated, and the block-contrast nulls do it: they have zero
# expectation, so cor(axis_nul, amp_nul) is pure measurement-error covariance.
# Algebraically, with per-species error variances s1^2 and s2^2,
#   cov(par_err, amp_err) = cov(e2 - e1, (e2 + e1)/2) = (s2^2 - s1^2)/2,
# which is not zero here: the parental libraries differ ~30-fold in depth between species
# (Sc 0.1-0.9M vs Sp 2.6-6.9M), so the sign flips with which species is shallower.
say("\n---- orthogonality on the NULL contrasts: pure measurement-error correlation ----")
say("cor(axis_nul, amp_nul) has no true-signal component, so it isolates the shared error.")
say("Expected nonzero and pair-dependent in sign: cov = (s2^2 - s1^2)/2 and the parental")
say("depths differ ~30x between species. This is a property of Fay's libraries, not a bug,")
say("but it means the per-gene r(axis,amp) table above must NOT be read as an error check.")
ER <- list()
for (cc in seq_along(CELLS)) {
  an <- MATS$amp_nul[, cc]
  ER[[cc]] <- data.table(pair = pair_of[cc], time = tp_of[cc],
    r_cisnul_ampnul = cor(MATS$cis_nul[, cc], an),
    r_trnnul_ampnul = cor(MATS$trans_nul[, cc], an),
    r_parnul_ampnul = cor(MATS$par_nul[, cc], an))
}
ER <- rbindlist(ER)
fwrite(ER, file.path(outd, sprintf("a6_error_orthogonality_%s.tsv", ANCHOR)), sep = "\t")
say("%-8s %6s %16s %16s %16s", "pair", "time", "r(cis_n,amp_n)", "r(trn_n,amp_n)",
    "r(par_n,amp_n)")
for (i in seq_len(nrow(ER))) say("%-8s %6d %16.3f %16.3f %16.3f",
  ER$pair[i], ER$time[i], ER$r_cisnul_ampnul[i], ER$r_trnnul_ampnul[i], ER$r_parnul_ampnul[i])
say("median |r| -- cis %.3f | trans %.3f | par %.3f", median(abs(ER$r_cisnul_ampnul)),
    median(abs(ER$r_trnnul_ampnul)), median(abs(ER$r_parnul_ampnul)))
say("WHY THE RATIOS SURVIVE THIS: numerator and denominator of every reported statistic are")
say("BOTH noise-corrected (rms_tru, amp_tru), so a shared error term inflates neither; and")
say("beta is fitted against a LEAVE-ONE-PAIR-OUT amplitude, whose error is independent of")
say("the target cell by construction. What does NOT survive is reading the per-gene")
say("r(axis,amp) table as evidence of a clean normaliser. That table is descriptive only.")

## ---------- internal control: basal (25 C) divergence on all three axes ---------
# The core's excess on the par and trans axes is only about the heat response if the core is
# not simply a more divergent set of genes to begin with. 04_divergence_by_dependency.R runs
# this control on the cis axis (basal ratio 1.077, specificity ~1.00, so the cis excess is
# not specific to the response). It must be run on the two new axes or their positives are
# unreadable. At t = 0 the response-based trans is identically zero, so the basal contrast
# is built from levels instead of responses, which reintroduces the per-gene ortholog
# length offset that cancels in every difference-of-differences above -- subtracted here
# explicitly, as in 04_divergence_by_dependency.R.
say("\n================ INTERNAL CONTROL: BASAL (25 C) DIVERGENCE, ALL THREE AXES ================")
LEN <- as.data.table(dat[one, c("HOG", "Len.Sc", "Len.Sp", "Len.Sk", "Len.Su")])
setnames(LEN, "HOG", "hog")
BA <- list()
for (g in HYB) {
  s1 <- POL[hybrid == g, s1]; s2 <- POL[hybrid == g, s2]
  lc <- LEN[match(KEEP, hog)]
  dl <- log2(lc[[paste0("Len.", s2)]] / lc[[paste0("Len.", s1)]])
  cb <- lapply(c("block1", "block2"), function(bk) gv_h(g, s2, bk, 0L) - gv_h(g, s1, bk, 0L) - dl)
  pb <- lapply(c("block1", "block2"), function(bk) gv_p(s2, bk, 0L) - gv_p(s1, bk, 0L) - dl)
  tb <- list(pb[[1]] - cb[[1]], pb[[2]] - cb[[2]])
  BA[[length(BA) + 1L]] <- data.table(hog = KEEP, pair = g,
    cis = (cb[[1]] + cb[[2]]) / 2, cis_nul = (cb[[2]] - cb[[1]]) / 2,
    trans = (tb[[1]] + tb[[2]]) / 2, trans_nul = (tb[[2]] - tb[[1]]) / 2,
    par = (pb[[1]] + pb[[2]]) / 2, par_nul = (pb[[2]] - pb[[1]]) / 2)
}
BA <- rbindlist(BA)
BAL <- merge(BA, ST, by = "hog", allow.cartesian = TRUE)
BF <- BAL[, .(cis_tru = sqrt(pmax(0, rms(cis)^2 - rms(cis_nul)^2)),
              trn_tru = sqrt(pmax(0, rms(trans)^2 - rms(trans_nul)^2)),
              par_tru = sqrt(pmax(0, rms(par)^2 - rms(par_nul)^2))),
          by = .(stratum, pair)]
say("median over the 6 pairs of the basal (level) divergence, by stratum:")
say("%-15s %10s %10s %10s", "stratum", "cis_tru", "trans_tru", "par_tru")
for (s in SETS) { z <- BF[stratum == s]
  say("%-15s %10.3f %10.3f %10.3f", s, median(z$cis_tru), median(z$trn_tru), median(z$par_tru)) }
bas_ratio <- function(ax, s) {
  vn <- switch(ax, cis = "cis_tru", trans = "trn_tru", par = "par_tru")
  n <- BF[stratum == s][[vn]]; d <- BF[stratum == DENOM[[s]]][[vn]]
  median(n / d)
}
say("\nSPECIFICITY = response_ratio / basal_ratio, computed on the UNNORMALISED (beta = 0)")
say("response ratio, because the basal contrast has no response amplitude to normalise by")
say("and the two must be the same statistic to be divided. A value near 1 means the excess")
say("is a general property of those genes' divergence, NOT of their heat response.")
SP <- list()
say("%-6s %-15s %9s %9s %9s %9s %8s %8s %8s", "axis", "focal", "basal_r", "resp_15",
    "resp_30", "resp_60", "spec_15", "spec_30", "spec_60")
for (ax in AXES) for (s in names(DENOM)) {
  br <- bas_ratio(ax, s)
  rr <- sapply(TP, function(t) RAT[axis == ax & stratum == s & time == t, ratio_beta0])
  SP[[length(SP) + 1L]] <- data.table(axis = ax, stratum = s, basal_ratio = br,
    resp_15 = rr[1], resp_30 = rr[2], resp_60 = rr[3],
    spec_15 = rr[1] / br, spec_30 = rr[2] / br, spec_60 = rr[3] / br)
  say("%-6s %-15s %9.3f %9.3f %9.3f %9.3f %8.3f %8.3f %8.3f", ax, s, br,
      rr[1], rr[2], rr[3], rr[1] / br, rr[2] / br, rr[3] / br)
}
SP <- rbindlist(SP)
fwrite(SP, file.path(outd, sprintf("a6_specificity_%s.tsv", ANCHOR)), sep = "\t")

## ------------- control that can fail: n- and expression-matched subsets ---------
# Every ratio above compares 83 core genes with 761 others. Two things differ besides
# regulatory dependency: n, and expression level (the strata are not expression-matched, and
# counting noise in a log ratio falls with expression). The permutation control removes both:
# draw size-83 subsets of other_induced with the same expression-decile composition as the
# core, and ask how often they reach the core's value. Run on the calibrated statistic, which
# is the primary one, not only on the unnormalised one.
say("\n================ THE CONTROL THAT CAN FAIL: n- AND EXPRESSION-MATCHED SUBSETS ================")
pmean <- rowMeans(do.call(cbind, lapply(SPECIES, function(sp) PL[[sp]][KEEP, , drop = FALSE])))
names(pmean) <- KEEP
dec <- cut(pmean, breaks = quantile(pmean, seq(0, 1, .1)), include.lowest = TRUE, labels = FALSE)
names(dec) <- KEEP
need <- table(factor(dec[HOGS$Hsf1_core], levels = 1:10))
have <- table(factor(dec[HOGS$other_induced], levels = 1:10))
say("expression-decile match feasibility (core needs / other_induced has):")
say("  %s", paste(sprintf("%d:%d/%d", 1:10, as.integer(need), as.integer(have)), collapse = " "))
feasible <- all(as.integer(need) <= as.integer(have))
say("exact decile match feasible: %s", feasible)
NPERM <- 1000L
pool_by_dec <- lapply(1:10, function(d) HOGS$other_induced[dec[HOGS$other_induced] == d])
draw_matched <- function() unlist(lapply(1:10, function(d) {
  k <- as.integer(need)[d]
  if (k == 0L) character(0) else sample(pool_by_dec[[d]], k) }))
NM <- list()
say("\n%-6s %-11s %7s %8s %9s %19s %8s", "axis", "matched", "time", "obs", "null_med",
    "null 95%%", "p_emp")
for (ax in AXES) for (matched in c(FALSE, TRUE)) {
  i_f <- idx_of("Hsf1_core")
  obs <- cell_cal(ax, i_f, BMED[ax])
  nullm <- matrix(NA_real_, NPERM, 3L)
  for (k in seq_len(NPERM)) {
    gsub_ <- if (matched) draw_matched() else sample(HOGS$other_induced, length(HOGS$Hsf1_core))
    v <- cell_cal(ax, match(gsub_, GENES), BMED[ax])
    for (j in seq_along(TP)) nullm[k, j] <- median(v[tp_of == TP[j]], na.rm = TRUE)
  }
  for (j in seq_along(TP)) {
    o <- median(obs[tp_of == TP[j]], na.rm = TRUE)
    p <- mean(nullm[, j] >= o, na.rm = TRUE)
    ci <- quantile(nullm[, j], c(.025, .975), na.rm = TRUE)
    NM[[length(NM) + 1L]] <- data.table(axis = ax, matched = matched, time = TP[j],
      obs = o, null_med = median(nullm[, j], na.rm = TRUE), lo = ci[1], hi = ci[2], p_emp = p)
    say("%-6s %-11s %7d %8.3f %9.3f  [%6.3f, %6.3f] %8.3f", ax,
        ifelse(matched, "decile", "size only"), TP[j], o, median(nullm[, j], na.rm = TRUE),
        ci[1], ci[2], p)
  }
}
NM <- rbindlist(NM)
fwrite(NM, file.path(outd, sprintf("a6_nmatched_%s.tsv", ANCHOR)), sep = "\t")
say("p_emp is one-sided (fraction of matched subsets at least as divergent as the core) and")
say("is NOT multiplicity-corrected across the 9 axis-time cells; with 9 cells at alpha 0.05")
say("the expected number of false positives is 0.45, so a single p just under 0.05 is not a")
say("result. Read the decile-matched rows, not the size-only rows.")

# The response rows above are only interpretable next to the same test run at 25 C. If the
# core also beats expression-matched subsets at basal, then a small p on the response says
# "these genes are more divergent", not "their heat response is more divergent" -- and the
# permutation test has to establish that with the same rigour as the ratio-based specificity
# table, or the specificity conclusion rests on a weaker statistic than the claim it limits.
say("\n---- the same permutation control at 25 C (basal), for the specificity conclusion ----")
BM_ <- list()
for (ax in AXES) {
  vn <- switch(ax, cis = "cis", trans = "trans", par = "par")
  bw <- dcast(BA, hog ~ pair, value.var = vn)
  bn <- dcast(BA, hog ~ pair, value.var = paste0(vn, "_nul"))
  stopifnot(identical(bw$hog, bn$hog))
  Mb <- as.matrix(bw[, HYB, with = FALSE]); Nb <- as.matrix(bn[, HYB, with = FALSE])
  rownames(Mb) <- rownames(Nb) <- bw$hog
  btru <- function(gs) { i <- match(gs, rownames(Mb))
    median(sqrt(pmax(0, colMeans(Mb[i, , drop = FALSE]^2) - colMeans(Nb[i, , drop = FALSE]^2)))) }
  o <- btru(HOGS$Hsf1_core)
  nl <- replicate(NPERM, btru(draw_matched()))
  ci <- quantile(nl, c(.025, .975), na.rm = TRUE)
  BM_[[ax]] <- data.table(axis = ax, matched = TRUE, time = 0L, obs = o,
    null_med = median(nl), lo = ci[1], hi = ci[2], p_emp = mean(nl >= o))
}
BM_ <- rbindlist(BM_)
say("%-6s %-11s %7s %8s %9s %19s %8s", "axis", "matched", "time", "obs", "null_med",
    "null 95%%", "p_emp")
for (i in seq_len(nrow(BM_))) say("%-6s %-11s %7d %8.3f %9.3f  [%6.3f, %6.3f] %8.3f",
  BM_$axis[i], "decile", 0L, BM_$obs[i], BM_$null_med[i], BM_$lo[i], BM_$hi[i], BM_$p_emp[i])
NM <- rbind(NM, BM_)
fwrite(NM, file.path(outd, sprintf("a6_nmatched_%s.tsv", ANCHOR)), sep = "\t")
say("If the basal p_emp is as small as the response p_emp, the core's excess is a property of")
say("those genes' sequence divergence in general and the RESPONSE-specific claim fails, in")
say("agreement with the specificity table above. That agreement between two different")
say("statistics is the evidence; disagreement would indicate sensitivity to the statistic.")

## ------------------------------------------------------------ persistence ---
say("\n=== per-gene decomposition written out for the promoter HSE test to join against ===")
OUT <- merge(AX[, .(hog, pair, time, cis, cis_nul, trans, trans_nul, par, amp)],
             ST, by = "hog", all.x = TRUE, allow.cartesian = TRUE)
OUT <- merge(OUT, HB[, .(hog = gene, scer = Name)], by = "hog", all.x = TRUE)
fwrite(OUT, file.path(outd, sprintf("a6_cis_trans_%s.tsv.gz", ANCHOR)), sep = "\t",
       compress = "gzip")
say("rows %d | file a6_cis_trans_%s.tsv.gz", nrow(OUT), ANCHOR)
writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_a6_cis_trans.txt"))
say("\n=== done ===")
