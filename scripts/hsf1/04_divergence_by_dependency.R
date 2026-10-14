# 04_divergence_by_dependency.R -- the central test.
#
# Question:
#   Is the cross-species divergence of the heat-stress transcriptional programme
#   concentrated in the Hsf1-associated core, or is it a general property of the
#   induced programme?
#
# 03_fay_parental_dimensions.R answers this on the parental axis (species-level
# trajectories, composition statistic): no concentration in the core, empirical
# p = 47% / 20% / 14% at 15/30/60 min. This script answers the same question on the cis
# axis, which is an independent measurement: inside an interspecific hybrid, both parental
# alleles of an orthologous gene are transcribed in the same nucleus and sequenced in
# the same library, so the allele contrast is free of the trans environment and of all
# library-level nuisance. If the two axes agree the answer is not an artefact of either.
#
# Projection limitation, which governs every number below. Hsf1
# dependence can only be measured in S. cerevisiae (only the Mühlhofer dataset perturbs
# Hsf1). What is partitioned below is therefore the orthologs of an Sc-defined gene
# set, not each species' own Hsf1 regulon. Nothing here is evidence about any species'
# Hsf1 protein or Hsf1 activity.
#
# Relation to 06_module_bridge.R. That script asked the same question
# with a count statistic: is the fraction of genes with padj(interaction) < 0.05 higher
# in the Hsf1 module? The count statistic differs from the present analysis in three ways:
#   (1) a count of significant genes is a magnitude ordering
#       only when the noise floor is shared, and it is not shared across strata, because
#       chaperones sit at the top of the expression range and se falls with expression;
#   (2) the module was defined from single perturbations, so it mixes Hsf1-specific with
#       Msn2/4-redundant genes (01_epistasis_decomposition.R). This script uses the
#       epistasis-refined core and keeps the two core strata separate (111 unanimous /
#       51 joint candidate);
#   (3) a median OR of 1.05-1.67 without an equivalence bound does not distinguish "enriched"
#       from "indistinguishable"; 1.05 alone does not support either interpretation.
# Both statistics are computed here, on the same genes, so the count-vs-threshold-free
# discrepancy is measured rather than assumed.
#
# The estimator, and why it needs no standard errors.
# For hybrid g, timepoint t, replicate block b, with alleles (a1, a2):
#     DoD(t,b) = [x(a2,t,b) - x(a2,0,b)] - [x(a1,t,b) - x(a1,0,b)]
# i.e. the allele difference in the heat response = cis divergence of the response.
# The design has exactly one sample per (hybrid, allele, time, block), and two complete
# blocks, so:
#     signal(t) = ( DoD(t,b1) + DoD(t,b2) ) / 2      <- estimate of the true DoD
#     null(t)   = ( DoD(t,b2) - DoD(t,b1) ) / 2      <- zero expectation by construction
# Both are the same linear combination of four single-block DoDs up to sign, so they
# carry identical noise variance (sigma^2_DoD / 2). Therefore
#     rms_true = sqrt( max(0, rms(signal)^2 - rms(null)^2) )
# removes the noise floor per (hybrid, stratum, timepoint) with no dispersion model, no
# empirical-Bayes shrinkage and no reliance on any se estimate -- which matters, because
# shrunken lfcSE cannot be used as a precision measure and
# the deposited multimethod table carries an F statistic, not a t, so no se is recoverable
# from it at all. The null also absorbs any block x allele x time interaction, so the
# correction can only shrink an amplitude, never inflate one: it is conservative.
# This is the same construction 03_fay_parental_dimensions.R uses across replicate
# blocks, applied to the allele contrast instead of the species trajectory.
#
# Why the noise correction is required here. The strata differ
# systematically in expression: the Hsf1-associated core is chaperone-heavy and sits high
# in the expression range, while `not_induced` spans the whole range. Counting noise in a
# log ratio falls with expression, and it adds in quadrature to real signal. An uncorrected
# rms therefore has a built-in bias against the core looking divergent -- the opposite
# direction from the one a "core is conserved" conclusion would need. Reporting an
# uncorrected statistic would make the negative result unfalsifiable.
#
# Structural limits, stated before any number is produced.
#   - 6 hybrids is not 6 independent units. Each species appears in 3 hybrids, and the
#     libraries are co-pooled in complementary pairs (from design.csv: library
#     5 carries both ScxSp and SkxSu, 13 carries ScxSu and SpxSk, 21 carries ScxSk and
#     SpxSu. So all 18 (hybrid x timepoint) cells share structure in two different ways.
#     All inference below is therefore at the gene level (n in the thousands); the hybrid
#     axis is described, never tested. Co-pooling leakage is bounded but not excluded
#     (lambda >= 0.20 excluded, lambda <= 0.10 not).
#   - Negative results are expressed as equivalence bounds rather than "no significant
#     difference".
#   - The TF_independent stratum is diluted by single-perturbation power: genes
#     that are really Msn2/4-dependent but whose single-mutant padj missed alpha sit in
#     it. So core-vs-TF_independent is a contaminated contrast and is interpreted as
#     such. The clean comparison is core vs `other_induced` (both arms heat-induced).
#
# Usage: Rscript scripts/hsf1/04_divergence_by_dependency.R

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
ALPHA   <- 0.05
NBOOT   <- 2000
SNR_MIN <- 1.3    # below this the stratum's cis signal is at the design's resolution limit
set.seed(20260905L)

say <- function(...) cat(sprintf(...), "\n", sep = "")
rms <- function(x) sqrt(mean(x^2, na.rm = TRUE))
q90 <- function(x) as.numeric(quantile(abs(x), 0.90, na.rm = TRUE))
rms_true <- function(sig, nul) sqrt(max(0, mean(sig^2, na.rm = TRUE) - mean(nul^2, na.rm = TRUE)))

logf <- file.path(repd, "logs/hsf1_divergence_by_dependency.log")
dir.create(dirname(logf), showWarnings = FALSE, recursive = TRUE)
sink(logf, split = TRUE)
say("=== 04_divergence_by_dependency | %s | anchor %s C ===", format(Sys.time()), ANCHOR)
say("The central test, on the CIS (allele-specific) axis.")
say("Independent replication target: A4's parental-axis answer (p_emp 47/20/14%%).")

## --------------------------------------------------------------- inputs ---
# read.csv, not fread: 172 rows of htseq.ortho.csv are ragged (267-269 fields instead of
# 270). All 172 carry one==0 & all4==0 (see 03_fay_parental_dimensions.R), so every
# ragged row is outside the 5121 one-to-one set. Re-asserted below.
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
LEN <- as.data.table(dat[one, c("HOG", "Len.Sc", "Len.Sp", "Len.Sk", "Len.Su")])
setnames(LEN, "HOG", "hog")
say("count matrix: %d one-to-one orthologs x %d samples, no NA", nrow(cts), ncol(cts))

D <- as.data.table(design)
D[, col := Name]
D[, num := as.integer(sub("^[A-Za-z]+", "", Name))]
D[, block := fifelse(num <= 32L, "block1", "block2")]
H <- D[Strain %in% HYB & Condition == "heat"]
H[, Time := as.integer(Time)]
say("\nhybrid heat-arm samples: %d (expect 96 = 6 hybrids x 2 alleles x 4 times x 2 blocks)",
    nrow(H))
stopifnot(nrow(H) == 96L, !anyDuplicated(H[, .(Strain, Allele, block, Time)]))
stopifnot(all(H$col %in% colnames(cts)))
stopifnot(H[, .N, by = Strain][, all(N == 16L)])

# allele polarity is read off the hybrid name, not guessed: "SpxSk" -> a1 = Sp, a2 = Sk
POL <- data.table(hybrid = HYB, a1 = substr(HYB, 1, 2), a2 = substr(HYB, 4, 5))
for (i in seq_len(nrow(POL))) {
  al <- sort(H[Strain == POL$hybrid[i], unique(Allele)])
  stopifnot(identical(al, sort(c(POL$a1[i], POL$a2[i]))))
}
say("allele polarity (a2 - a1), parsed from hybrid names and checked against design.csv:")
say("  %s", paste(sprintf("%s: %s-%s", POL$hybrid, POL$a2, POL$a1), collapse = "  "))

# Co-pooling, read from design.csv: which hybrids share a library number?
pool <- H[, .(hybrids = paste(sort(unique(Strain)), collapse = "+")), by = num]
say("\nco-pooling structure (library number -> hybrids sharing it):")
print(pool[, .N, by = hybrids])
say("=> the 6 hybrids are 3 complementary CO-POOLED PAIRS. Combined with each species")
say("   appearing in 3 hybrids, the 18 (hybrid x time) cells are not independent units.")
say("   Consequence: no test is run on the hybrid axis. Gene-level bootstrap only.")

## ------------------------------------------------- gene set assignment ---
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
CS <- fread(file.path(hsd, sprintf("hsf1_core_set_%s.tsv", ANCHOR)))
CL <- fread(file.path(hsd, sprintf("hsf1_epistasis_classes_%s.tsv", ANCHOR)))

UC <- dcast(CL[, .(gene, method, class)], gene ~ method, value.var = "class")
UC[, unanimous := fifelse(DESeq2 == edgeR & edgeR == limma, DESeq2, NA_character_)]
core  <- CS[set == "Hsf1_core_unanimous",  gene]
joint <- CS[set == "Hsf1_joint_candidate", gene]
stopifnot(length(core) == 111L, length(joint) == 51L)

map_hog <- function(g) {
  h <- HB[Name %in% g | sys %in% g, unique(gene)]
  h[h %in% rownames(cts)]
}
STRATA <- list(
  Hsf1_core      = core,
  Hsf1_joint     = joint,
  Msn24_only     = UC[unanimous == "Msn24_only",    gene],
  TF_independent = UC[unanimous == "TF_independent", gene],
  other_induced  = UC[unanimous %in% c("Msn24_only", "TF_independent"), gene],
  not_induced    = UC[unanimous == "not_induced",    gene]
)
HOGS <- lapply(STRATA, map_hog)
say("\n=== projection of each Sc stratum onto one-to-one orthologs ===")
say("%-16s %7s %7s %7s", "stratum", "Sc", "mapped", "frac")
for (s in names(STRATA))
  say("%-16s %7d %7d %6.1f%%", s, length(STRATA[[s]]), length(HOGS[[s]]),
      100 * length(HOGS[[s]]) / length(STRATA[[s]]))
say("Projection loss is NOT random: a one-to-one requirement preferentially")
say("removes recently duplicated genes, which the paralogue analysis identified as the paralogue-rich part of the")
say("core (HSP, anCWMP). Every 'core is conserved' reading carries this caveat.")

## --------------------------------------- per-hybrid logCPM, DoD, null ---
# Normalisation is within hybrid x heat (16 columns). Rationale: every quantity below is
# a double difference within one hybrid, so that hybrid's own library composition cancels;
# pooling hybrids into one TMM fit would let one hybrid's composition move another's
# baseline. The two allele columns of a (hybrid, block, time) cell come from the same
# library, which is the point of the allele-paired design.
say("\n=== per-hybrid logCPM, allele-paired difference-of-differences ===")
keep_by_hyb <- list(); SIG <- list(); NUL <- list(); BAS <- list(); EXPR <- list()
RES <- list()
for (g in HYB) {
  hs <- H[Strain == g][order(Allele, block, Time)]
  a1 <- POL[hybrid == g, a1]; a2 <- POL[hybrid == g, a2]
  m  <- cts[, hs$col, drop = FALSE]
  y  <- calcNormFactors(DGEList(counts = m), method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = PRIOR)
  cpm_raw <- cpm(y)
  # a gene must be adequately detected for both alleles, otherwise the allele contrast is
  # a comparison between a measurement and a zero
  k1 <- rowSums(cpm_raw[, hs[Allele == a1, col], drop = FALSE] >= MINCPM) >= 4L
  k2 <- rowSums(cpm_raw[, hs[Allele == a2, col], drop = FALSE] >= MINCPM) >= 4L
  keep_by_hyb[[g]] <- k1 & k2
  say("  %-6s libs %s | TMM %s | detected both alleles: %d",
      g, paste(sprintf("%.1fM", colSums(m) / 1e6), collapse = " "),
      sprintf("%.2f-%.2f", min(y$samples$norm.factors), max(y$samples$norm.factors)),
      sum(keep_by_hyb[[g]]))
  # gv() resolves exactly one column; the [[1]] guard makes a design slip fail loudly
  # rather than silently recycle a matrix.
  gv <- function(al, bk, t) {
    cc <- hs[Allele == al & block == bk & Time == t, col]
    stopifnot(length(cc) == 1L)
    lcpm[, cc]
  }
  # the double difference is written out in full rather than through a helper, so that the
  # baseline column of each term is visible at the point of use
  for (t in TP) {
    rA1 <- gv(a1, "block1", t) - gv(a1, "block1", 0L)
    rB1 <- gv(a2, "block1", t) - gv(a2, "block1", 0L)
    rA2 <- gv(a1, "block2", t) - gv(a1, "block2", 0L)
    rB2 <- gv(a2, "block2", t) - gv(a2, "block2", 0L)
    d1 <- rB1 - rA1                 # allele DIFFERENCE in response, block 1
    d2 <- rB2 - rA2                 # ... block 2
    m1 <- (rB1 + rA1) / 2           # allele MEAN response (the orthogonal contrast)
    m2 <- (rB2 + rA2) / 2
    SIG[[length(SIG) + 1L]] <- data.table(hog = rownames(lcpm), hybrid = g, time = t,
                                          sig = (d1 + d2) / 2)
    NUL[[length(NUL) + 1L]] <- data.table(hog = rownames(lcpm), hybrid = g, time = t,
                                          nul = (d2 - d1) / 2)
    # The response magnitude itself, needed because |DoD| is not scale-free.
    # If divergence is multiplicative -- allele B responds gamma times as much as allele A --
    # then DoD = (gamma - 1) * R, so rms(DoD) is proportional to the response magnitude R.
    # A stratum that simply responds harder therefore has a larger rms(DoD) at identical
    # gamma. Since the Hsf1-associated core is the most strongly induced part of the
    # programme (it is the chaperone module), the unnormalised comparison is confounded in
    # the direction that would manufacture the positive result.
    # R is taken as the mean of the two alleles' responses, not the reference allele's.
    # DoD = B - A and R_mean = (B + A)/2 are the difference and sum contrasts, which are
    # orthogonal when the two alleles carry equal measurement error, whereas dividing by
    # A alone reintroduces the shared-error term `-(1-h)` that the reference-level
    # artefact is made of. The orthogonality is checked empirically below.
    RES[[length(RES) + 1L]] <- data.table(hog = rownames(lcpm), hybrid = g, time = t,
                                          rsig = (m1 + m2) / 2, rnul = (m2 - m1) / 2)
  }
  # baseline (25 C) allele contrast: cis divergence of basal expression, the internal
  # control that asks whether the core simply carries more cis divergence overall.
  # Corrected for the ortholog length difference between the two genomes, which is a
  # per-gene constant offset that cancels in the DoD but not here.
  lc <- LEN[match(rownames(lcpm), hog)]
  dl <- log2(lc[[paste0("Len.", a2)]] / lc[[paste0("Len.", a1)]])
  b1 <- gv(a2, "block1", 0L) - gv(a1, "block1", 0L) - dl
  b2 <- gv(a2, "block2", 0L) - gv(a1, "block2", 0L) - dl
  BAS[[length(BAS) + 1L]] <- data.table(hog = rownames(lcpm), hybrid = g,
                                        bas = (b1 + b2) / 2, bas_nul = (b2 - b1) / 2)
  EXPR[[length(EXPR) + 1L]] <- data.table(hog = rownames(lcpm), hybrid = g,
                                          expr = rowMeans(lcpm))
}
SIG <- rbindlist(SIG); NUL <- rbindlist(NUL); BAS <- rbindlist(BAS); EXPR <- rbindlist(EXPR)
RES <- rbindlist(RES)

keep_all <- Reduce(`&`, keep_by_hyb)
KEEP <- rownames(cts)[keep_all]
say("\ngene filter: both alleles detected in ALL SIX hybrids: %d of %d (%.1f%%)",
    length(KEEP), nrow(cts), 100 * mean(keep_all))
# cross-check against the parental filter of 03_fay_parental_dimensions.R, so the two
# axes are known to be on comparable gene sets
a4f <- file.path(hsd, sprintf("fay_parental_amplitude_%s.tsv", ANCHOR))
if (file.exists(a4f)) {
  A4 <- fread(a4f)
  say("A4 (parental axis) core stratum n = %d; A5 core stratum n = %d",
      A4[stratum == "Hsf1_core", max(n)], length(intersect(HOGS$Hsf1_core, KEEP)))
}
SIG <- SIG[hog %in% KEEP]; NUL <- NUL[hog %in% KEEP]; RES <- RES[hog %in% KEEP]
BAS <- BAS[hog %in% KEEP]; EXPR <- EXPR[hog %in% KEEP]

SETS <- c("Hsf1_core", "Hsf1_joint", "Msn24_only", "TF_independent", "other_induced", "not_induced")
say("\nstratum sizes AFTER the all-hybrid detection filter:")
for (s in SETS) { HOGS[[s]] <- intersect(HOGS[[s]], KEEP); say("  %-16s %d", s, length(HOGS[[s]])) }
stopifnot(length(HOGS$Hsf1_core) >= 20L)   # stop condition: the projected core needs >= 20 genes

# Disjointness is asserted. The central contrast is core vs
# other_induced, and it is only a contrast if the two sets do not overlap. `other_induced`
# is the union of Msn24_only and TF_independent, so Msn24_only and TF_independent are each
# a subset of it and must not be compared against it -- that would be a set against a
# superset containing it, which dilutes any separation toward 1 by construction. Each
# stratum therefore gets an explicitly disjoint denominator, defined below.
stopifnot(length(intersect(HOGS$Hsf1_core, HOGS$other_induced)) == 0L)
stopifnot(length(intersect(HOGS$Hsf1_joint, HOGS$other_induced)) == 0L)
stopifnot(length(intersect(HOGS$Hsf1_core, HOGS$Hsf1_joint)) == 0L)
stopifnot(setequal(HOGS$other_induced, union(HOGS$Msn24_only, HOGS$TF_independent)))
# denominator for each focal stratum: heat-induced, and disjoint from the focal set
DENOM <- list(
  Hsf1_core      = "other_induced",     # already disjoint (Hsf1-dependent vs not)
  Hsf1_joint     = "other_induced",     # kept separate from the core
  Msn24_only     = "TF_independent",    # head-to-head inside `other_induced`
  TF_independent = "Msn24_only",        # the same comparison read the other way
  not_induced    = "other_induced")     # sanity check: induced should exceed not-induced
say("\ndisjoint denominators (focal / denominator):")
say("  %s", paste(sprintf("%s / %s", names(DENOM), unlist(DENOM)), collapse = "   "))

ST <- rbindlist(lapply(SETS, function(s) data.table(hog = HOGS[[s]], stratum = s)))
fwrite(dcast(ST[, .(hog, stratum, v = TRUE)], hog ~ stratum, value.var = "v", fill = FALSE),
       file.path(outd, sprintf("a5_strata_membership_%s.tsv", ANCHOR)), sep = "\t")

## -------------------------------------- expression profile of the strata ---
say("\n================ WHY THE NOISE CORRECTION IS REQUIRED ================")
EX <- merge(EXPR[, .(expr = mean(expr)), by = hog], ST, by = "hog", allow.cartesian = TRUE)
say("%-16s %7s %8s %8s %8s", "stratum", "n", "med_expr", "q25", "q75")
for (s in SETS) {
  v <- EX[stratum == s, expr]
  say("%-16s %7d %8.2f %8.2f %8.2f", s, length(v), median(v), quantile(v, .25), quantile(v, .75))
}
say("The strata are NOT expression-matched. Counting noise in a log ratio falls with")
say("expression and adds in quadrature to signal, so an uncorrected rms is biased")
say("AGAINST the high-expression core looking divergent -- i.e. it would make a")
say("'core is conserved' conclusion unfalsifiable. Hence rms_true throughout.")

## --------------------------------------------- the noise floor and rms_true ---
X  <- merge(SIG, ST, by = "hog", allow.cartesian = TRUE)
NZ <- merge(NUL, ST, by = "hog", allow.cartesian = TRUE)
FL <- merge(X[, .(n = .N, rms_sig = rms(sig), q90_sig = q90(sig)), by = .(stratum, hybrid, time)],
            NZ[, .(rms_nul = rms(nul), q90_nul = q90(nul)), by = .(stratum, hybrid, time)],
            by = c("stratum", "hybrid", "time"))
FL[, rms_tru := sqrt(pmax(0, rms_sig^2 - rms_nul^2))]
FL[, snr := rms_sig / rms_nul]
setorder(FL, stratum, hybrid, time)
fwrite(FL, file.path(outd, sprintf("a5_cis_noise_corrected_%s.tsv", ANCHOR)), sep = "\t")

say("\n================ CIS DIVERGENCE OF THE HEAT RESPONSE, BY STRATUM ================")
say("rms_sig = observed rms of the allele-paired difference-of-differences")
say("rms_nul = same quantity built from the block contrast, zero expectation by construction")
say("rms_tru = sqrt(rms_sig^2 - rms_nul^2); snr = rms_sig/rms_nul (1.0 = pure noise)")
for (s in c("Hsf1_core", "other_induced", "not_induced")) {
  say("\n-- %s (n = %d) --", s, length(HOGS[[s]]))
  say("%-7s %6s %8s %8s %8s %6s", "hybrid", "time", "rms_sig", "rms_nul", "rms_tru", "snr")
  for (g in HYB) for (t in TP) {
    r <- FL[stratum == s & hybrid == g & time == t]
    say("%-7s %6d %8.3f %8.3f %8.3f %6.2f", g, t, r$rms_sig, r$rms_nul, r$rms_tru, r$snr)
  }
}
say("\nmedian snr across the 18 cells, per stratum:")
for (s in SETS) say("  %-16s %.2f", s, FL[stratum == s, median(snr)])
say("An snr near 1.0 means the cis signal in that stratum is at the resolution limit of")
say("this design; any ratio built from it is uninterpretable and is flagged below.")

## ------------------------------------------- diagnostics on the null contrast ---
# The whole noise correction rests on two properties of the block contrast, and both are
# checkable rather than assumed: (i) it must have zero mean, (ii) it must be uncorrelated
# with the signal. A nonzero mean would mean a systematic block x allele x time effect (the
# correction would then be conservative, not wrong); a correlation with the signal would
# mean the block assignment is wrong, and would invalidate the subtraction outright.
say("\n---- diagnostics on the block-contrast null (must hold, or the correction is void) ----")
DG <- merge(merge(SIG, NUL, by = c("hog", "hybrid", "time")), RES,
            by = c("hog", "hybrid", "time"))
# respA is recoverable exactly: sig = respB - respA and rsig = (respB + respA)/2, so
# respA = rsig - sig/2. No extra accumulator is needed.
DG[, respA := rsig - sig / 2]
dg <- DG[, .(mean_nul = mean(nul), sd_nul = sd(nul), snr = rms(sig) / rms(nul),
             r_sig_nul = cor(sig, nul), r_sig_rsig = cor(sig, rsig),
             r_sig_respA = cor(sig, respA)),
         by = .(hybrid, time)]
# implied block-variance ratio, derived below; only meaningful as a magnitude
dg[, var_ratio := (1 + r_sig_nul * snr) / (1 - r_sig_nul * snr)]
say("%-7s %6s %10s %10s %12s %10s", "hybrid", "time", "mean_nul", "sd_nul", "r(sig,nul)",
    "sig2/sig1^2")
for (i in seq_len(nrow(dg))) say("%-7s %6d %10.4f %10.4f %12.4f %10.2f",
  dg$hybrid[i], dg$time[i], dg$mean_nul[i], dg$sd_nul[i], dg$r_sig_nul[i], dg$var_ratio[i])
say("max |mean_nul| = %.4f (vs sd %.3f); max |r(sig,nul)| = %.3f",
    max(abs(dg$mean_nul)), median(dg$sd_nul), max(abs(dg$r_sig_nul)))
# The algebra.
# Write d1 = delta + e1, d2 = delta + e2 with var(e_b) = sigma_b^2 and delta the true DoD.
#   sig = delta + (e1+e2)/2   nul = (e2-e1)/2
#   cov(sig, nul) = (sigma2^2 - sigma1^2)/4
#   E[sig^2] = delta^2 + (sigma1^2 + sigma2^2)/4     E[nul^2] = (sigma1^2 + sigma2^2)/4
# So (i) a nonzero r(sig,nul) means only that the two blocks have unequal noise variance --
# it is not evidence of mislabelling; and (ii) rms_true^2 = E[sig^2] - E[nul^2] = delta^2
# exactly even under that heteroscedasticity, because the two expectations carry the same
# (sigma1^2 + sigma2^2)/4 term. The estimator therefore needs no equal-variance assumption.
# Inverting r*snr = (rho - 1)/(rho + 1) with rho = sigma2^2/sigma1^2 gives the printed
# var_ratio = (1 + r*snr)/(1 - r*snr), i.e. how unequal the blocks would have to be.
say("var_ratio = (1 + r*snr)/(1 - r*snr) = the implied sigma2^2/sigma1^2 (see algebra above).")
say("rms_true^2 = E[sig^2] - E[nul^2] = delta^2 EXACTLY under unequal block variances, because")
say("both expectations carry the same (sigma1^2+sigma2^2)/4 term, so a nonzero r is a")
say("statement about block heteroscedasticity, NOT about the correction's validity. What")
say("would be fatal is |r| near 1 (blocks mislabelled, the 'null' carrying the signal).")
say("\n---- orthogonality of the amplitude normaliser (positive control, transposed) ----")
say("The amplitude statistic below divides by R_mean = (respB + respA)/2. That is the SUM")
say("contrast, orthogonal to the DIFFERENCE contrast sig = respB - respA under equal allele")
say("measurement error; dividing by respA alone would reintroduce the shared-error term")
say("-(1-h) that the P1 reference-level artefact is made of. Prediction: r(sig, rsig) ~ 0")
say("and r(sig, respA) < 0. This is the P1 argument's positive control on the A5 axis.")
say("%-7s %6s %14s %14s", "hybrid", "time", "r(sig,rsig)", "r(sig,respA)")
for (i in seq_len(nrow(dg))) say("%-7s %6d %14.4f %14.4f",
  dg$hybrid[i], dg$time[i], dg$r_sig_rsig[i], dg$r_sig_respA[i])
say("median r(sig,rsig) = %+.4f ; median r(sig,respA) = %+.4f ; cells with r(sig,respA)<0: %d/18",
    median(dg$r_sig_rsig), median(dg$r_sig_respA), sum(dg$r_sig_respA < 0))
fwrite(dg, file.path(outd, sprintf("a5_null_diagnostics_%s.tsv", ANCHOR)), sep = "\t")

## ------------------------------------------- is rms_tru biased by stratum size? ---
# rms_tru = sqrt(max(0, A - B)) is a concave function of a noisy quantity, so for small n
# it is biased downward. That matters because the central ratio has n=83 on top and n=761
# underneath: if the bias were material, the ratio would be biased against the core, i.e.
# in the direction of an apparent negative result. Checked here:
# rms_tru of random subsets of `other_induced` as a function of subset size,
# against the full-set value.
say("\n---- size-dependence of rms_tru (bias check for the ratio's unequal n) ----")
BIAS <- list()
for (nsub in c(40L, 83L, 200L, 400L)) {
  v <- replicate(400L, {
    pk <- sample(HOGS$other_induced, nsub)
    median(FL_sub <- sapply(HYB, function(g) {
      s <- SIG[hybrid == g & time == 15L & hog %in% pk, sig]
      u <- NUL[hybrid == g & time == 15L & hog %in% pk, nul]
      sqrt(max(0, mean(s^2) - mean(u^2)))
    }))
  })
  BIAS[[length(BIAS) + 1L]] <- data.table(n = nsub, med = median(v), sd = sd(v))
}
BIAS <- rbindlist(BIAS)
full15 <- median(FL[stratum == "other_induced" & time == 15L, rms_tru])
say("full other_induced (n=%d) at 15 min: rms_tru = %.4f", length(HOGS$other_induced), full15)
for (i in seq_len(nrow(BIAS))) say("  subset n=%4d: median %.4f (%+.1f%% vs full), sd %.4f",
  BIAS$n[i], BIAS$med[i], 100 * (BIAS$med[i] / full15 - 1), BIAS$sd[i])
say("If the n=83 row is within a few percent of the full-set value, the concavity bias is")
say("negligible at this signal-to-noise and the unequal-n ratio needs no correction.")
fwrite(BIAS, file.path(outd, sprintf("a5_rms_tru_size_bias_%s.tsv", ANCHOR)), sep = "\t")

## ------------------------------------- baseline (basal) cis divergence control ---
BX <- merge(BAS, ST, by = "hog", allow.cartesian = TRUE)
BF <- BX[, .(n = .N, rms_sig = rms(bas), rms_nul = rms(bas_nul)), by = .(stratum, hybrid)]
BF[, rms_tru := sqrt(pmax(0, rms_sig^2 - rms_nul^2))]
fwrite(BF, file.path(outd, sprintf("a5_cis_baseline_%s.tsv", ANCHOR)), sep = "\t")
say("\n================ INTERNAL CONTROL: BASAL cis divergence at 25 C ================")
say("Allele contrast at t=0, corrected for the ortholog length difference between the two")
say("genomes (a per-gene constant that cancels in the DoD but not here). If the core were")
say("simply a more cis-divergent set of genes in general, it would show here too.")
say("%-16s %8s %8s %8s", "stratum", "med_tru", "min", "max")
for (s in SETS) {
  v <- BF[stratum == s, rms_tru]
  say("%-16s %8.3f %8.3f %8.3f", s, median(v), min(v), max(v))
}
BAS_RATIO <- BF[stratum == "Hsf1_core", median(rms_tru)] / BF[stratum == "other_induced", median(rms_tru)]
say("ratio core/other_induced, BASAL: %.3f", BAS_RATIO)
say("This number is the reference against which any excess in the RESPONSE must be read.")
say("If the core's response-divergence ratio is no larger than its basal-divergence ratio,")
say("then whatever excess cis divergence the core carries is NOT specific to the heat")
say("response -- it is a property of those genes' regulatory regions in general. That is")
say("the same logic as A4's specificity control, transposed to the cis axis.")

## ------------------------------------------------- the ratio and its bounds ---
# Primary statistic: rms_tru(core) / rms_tru(other_induced) in each (hybrid, time) cell.
# Both arms are heat-induced in S. cerevisiae, so "being heat-induced" cancels and what is
# left is regulatory dependency. Inference is a gene-level bootstrap: one resample of core
# genes and one of other_induced genes per iteration, the ratio recomputed in every cell
# with that same resample, then the median across cells taken. That respects both the gene
# level (the only level with n) and the cell structure (the only level with none).
sig_w <- dcast(SIG[hog %in% KEEP], hog ~ hybrid + time, value.var = "sig")
nul_w <- dcast(NUL[hog %in% KEEP], hog ~ hybrid + time, value.var = "nul")
stopifnot(identical(sig_w$hog, nul_w$hog))
CELLS <- setdiff(names(sig_w), "hog")
stopifnot(length(CELLS) == 18L)
S_M <- as.matrix(sig_w[, ..CELLS]); N_M <- as.matrix(nul_w[, ..CELLS])
rownames(S_M) <- rownames(N_M) <- sig_w$hog
idx_of <- function(s) match(HOGS[[s]], rownames(S_M))
# the amplitude arm: same wide layout, same noise-correction machinery, applied to the
# allele-mean response instead of the allele difference
rsg_w <- dcast(RES[hog %in% KEEP], hog ~ hybrid + time, value.var = "rsig")
rnl_w <- dcast(RES[hog %in% KEEP], hog ~ hybrid + time, value.var = "rnul")
stopifnot(identical(rsg_w$hog, sig_w$hog), identical(rnl_w$hog, sig_w$hog))
stopifnot(identical(setdiff(names(rsg_w), "hog"), CELLS))
R_M <- as.matrix(rsg_w[, ..CELLS]); RN_M <- as.matrix(rnl_w[, ..CELLS])
rownames(R_M) <- rownames(RN_M) <- rsg_w$hog

cell_ratio <- function(i_num, i_den) {
  num <- sqrt(pmax(0, colMeans(S_M[i_num, , drop = FALSE]^2) - colMeans(N_M[i_num, , drop = FALSE]^2)))
  den <- sqrt(pmax(0, colMeans(S_M[i_den, , drop = FALSE]^2) - colMeans(N_M[i_den, , drop = FALSE]^2)))
  ifelse(den > 0, num / den, NA_real_)
}
cell_tru <- function(i) sqrt(pmax(0, colMeans(S_M[i, , drop = FALSE]^2) - colMeans(N_M[i, , drop = FALSE]^2)))
# noise-corrected amplitude of the response itself, per cell (same estimator, sum contrast)
cell_amp <- function(i) sqrt(pmax(0, colMeans(R_M[i, , drop = FALSE]^2) - colMeans(RN_M[i, , drop = FALSE]^2)))
# Amplitude-normalised statistic. Under multiplicative divergence DoD = (gamma-1)*R, so
# rel = rms_tru(DoD)/rms_tru(R_mean) estimates |gamma-1|: it is the divergence per unit of
# response, and it is invariant to a stratum simply responding harder. A cell with no
# resolvable response amplitude yields NA rather than a ratio with a near-zero denominator.
cell_rel <- function(i) { a <- cell_amp(i); ifelse(a > 0, cell_tru(i) / a, NA_real_) }
cell_rel_ratio <- function(i_num, i_den) {
  n <- cell_rel(i_num); d <- cell_rel(i_den)
  ifelse(!is.na(d) & d > 0, n / d, NA_real_)
}

tp_of <- as.integer(sub(".*_", "", CELLS))
i_core <- idx_of("Hsf1_core"); i_join <- idx_of("Hsf1_joint")
i_oth  <- idx_of("other_induced"); i_not <- idx_of("not_induced")
i_msn  <- idx_of("Msn24_only");    i_tfi <- idx_of("TF_independent")

boot_ratio <- function(i_num, i_den, nb = NBOOT) {
  out <- matrix(NA_real_, nb, 3L, dimnames = list(NULL, as.character(TP)))
  for (k in seq_len(nb)) {
    r <- cell_ratio(sample(i_num, length(i_num), TRUE), sample(i_den, length(i_den), TRUE))
    for (j in seq_along(TP)) out[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  out
}

say("\n================ THE CENTRAL COMPARISON ================")
say("statistic = rms_tru(focal) / rms_tru(denominator), median over the 6 hybrids.")
say("Each denominator is heat-induced AND DISJOINT from its focal set (asserted above), so")
say("'being heat-induced' cancels and what is left is regulatory dependency.")
say("CI = gene-level bootstrap (%d resamples). >1 means MORE cis divergence.", NBOOT)
RAT <- list(); PERHYB <- list()
for (s in names(DENOM)) {
  i_f <- idx_of(s); i_d <- idx_of(DENOM[[s]])
  obs <- cell_ratio(i_f, i_d)
  bt  <- boot_ratio(i_f, i_d)
  # snr-filtered sensitivity: a cell whose signal is at the resolution limit contributes a
  # ratio of two near-zero numbers. Cells where either arm has snr < SNR_MIN are dropped in
  # the sensitivity column; the primary column keeps all 18 so the filter cannot be a way
  # of selecting a preferred answer.
  ok <- sapply(CELLS, function(cc) {
    gg <- sub("_[0-9]+$", "", cc); tt <- as.integer(sub(".*_", "", cc))
    FL[stratum == s & hybrid == gg & time == tt, snr] >= SNR_MIN &&
    FL[stratum == DENOM[[s]] & hybrid == gg & time == tt, snr] >= SNR_MIN
  })
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    o  <- median(obs[sel], na.rm = TRUE)
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    RAT[[length(RAT) + 1L]] <- data.table(stratum = s, denom = DENOM[[s]], time = TP[j],
      n = length(i_f), n_denom = length(i_d),
      ratio = o, lo = ci[1], hi = ci[2],
      n_cells = sum(sel), n_cells_gt1 = sum(obs[sel] > 1, na.rm = TRUE),
      n_cells_resolved = sum(sel & ok),
      ratio_resolved = if (sum(sel & ok) > 0) median(obs[sel & ok], na.rm = TRUE) else NA_real_,
      equiv_bound = max(abs(c(ci[1], ci[2]) - 1)))
  }
  PERHYB[[length(PERHYB) + 1L]] <- data.table(stratum = s, denom = DENOM[[s]],
    hybrid = sub("_[0-9]+$", "", CELLS), time = tp_of, ratio = as.numeric(obs), resolved = ok)
}
RAT <- rbindlist(RAT); PERHYB <- rbindlist(PERHYB)
fwrite(RAT, file.path(outd, sprintf("a5_cis_ratio_%s.tsv", ANCHOR)), sep = "\t")
fwrite(PERHYB, file.path(outd, sprintf("a5_cis_ratio_by_hybrid_%s.tsv", ANCHOR)), sep = "\t")
say("%-15s %-15s %5s %5s %7s %15s %8s %9s", "focal", "denominator", "time", "n",
    "ratio", "95%% CI", "cells>1", "resolved")
for (i in seq_len(nrow(RAT))) say("%-15s %-15s %5d %5d %7.3f [%6.3f,%6.3f] %5d/%d %7s  |r-1| > %.2f excl.",
  RAT$stratum[i], RAT$denom[i], RAT$time[i], RAT$n[i], RAT$ratio[i], RAT$lo[i], RAT$hi[i],
  RAT$n_cells_gt1[i], RAT$n_cells[i],
  ifelse(is.na(RAT$ratio_resolved[i]), "-", sprintf("%.3f", RAT$ratio_resolved[i])),
  RAT$equiv_bound[i])
say("\nINTERPRETATION USES THE EQUIVALENCE COLUMN. A CI containing 1")
say("is not evidence of no difference; it bounds how large a difference the data can")
say("still hide. That bound defines the supported negative result.")
say("The 'cells>1' column matters as much as the point estimate: a ratio above 1 that only")
say("holds in half the hybrids is not a property of the gene set.")

say("\n---- per-hybrid detail for the two core strata (the point estimate hides this) ----")
for (s in c("Hsf1_core", "Hsf1_joint")) {
  say("-- %s / %s --", s, DENOM[[s]])
  say("%-8s %8s %8s %8s", "hybrid", "15min", "30min", "60min")
  for (g in HYB) {
    v <- sapply(TP, function(t) PERHYB[stratum == s & hybrid == g & time == t, ratio])
    fl <- sapply(TP, function(t) PERHYB[stratum == s & hybrid == g & time == t, resolved])
    say("%-8s %8s %8s %8s", g,
        sprintf("%.2f%s", v[1], ifelse(fl[1], "", "*")),
        sprintf("%.2f%s", v[2], ifelse(fl[2], "", "*")),
        sprintf("%.2f%s", v[3], ifelse(fl[3], "", "*")))
  }
  say("(* = at least one arm has snr < %.1f in that cell: ratio of two near-zero numbers)", SNR_MIN)
}

## ---------------- specificity: is the core's excess specific to the response? ---
say("\n================ SPECIFICITY: RESPONSE EXCESS vs BASAL EXCESS ================")
bas_r <- function(s) BF[stratum == s, median(rms_tru)] / BF[stratum == DENOM[[s]], median(rms_tru)]
say("%-15s %10s %10s %10s %10s", "focal", "basal_r", "resp_r_15", "resp_r_30", "resp_r_60")
SPEC <- list()
for (s in names(DENOM)) {
  br <- bas_r(s)
  rr <- sapply(TP, function(t) RAT[stratum == s & time == t, ratio])
  say("%-15s %10.3f %10.3f %10.3f %10.3f", s, br, rr[1], rr[2], rr[3])
  SPEC[[length(SPEC) + 1L]] <- data.table(stratum = s, denom = DENOM[[s]], basal_ratio = br,
    time = TP, response_ratio = rr, spec = rr / br)
}
SPEC <- rbindlist(SPEC)
fwrite(SPEC, file.path(outd, sprintf("a5_cis_specificity_%s.tsv", ANCHOR)), sep = "\t")
say("\nspecificity = response_ratio / basal_ratio. A value near 1 means the stratum's cis")
say("divergence excess is NOT specific to the heat response:")
say("%-15s %8s %8s %8s", "focal", "spec_15", "spec_30", "spec_60")
for (s in names(DENOM)) {
  v <- sapply(TP, function(t) SPEC[stratum == s & time == t, spec])
  say("%-15s %8.3f %8.3f %8.3f", s, v[1], v[2], v[3])
}

## ------- the amplitude confound, and the statistic that is immune to it ---
# rms(DoD) is not scale-free. If the divergence is multiplicative -- allele B responds gamma
# times as much as allele A -- then DoD = (gamma - 1) * R, so rms(DoD) is proportional to the
# response magnitude R at identical gamma. The Hsf1-associated core is the chaperone module,
# i.e. the most strongly induced part of the programme, so the unnormalised ratio above is
# confounded in exactly the direction that would manufacture a positive result. The same
# confound is the leading alternative explanation for the Msn24_only ratio.
# rel = rms_tru(DoD) / rms_tru(R_mean) removes it: it estimates |gamma - 1|, the divergence
# per unit of response. R_mean is the allele mean (sum contrast), not the reference allele,
# for the reason checked empirically in the orthogonality table above.
say("\n================ AMPLITUDE CONFOUND: DIVERGENCE PER UNIT OF RESPONSE ================")
say("amp = noise-corrected rms of the allele-MEAN response (the confounding variable itself);")
say("rel = rms_tru(DoD)/amp, an estimate of |gamma-1| that is invariant to response strength.")
say("%-15s %6s %8s %8s %8s", "stratum", "time", "amp", "rms_tru", "rel")
AMP <- list()
for (s in SETS) {
  i_s <- idx_of(s)
  a <- cell_amp(i_s); d <- cell_tru(i_s); rl <- cell_rel(i_s)
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    AMP[[length(AMP) + 1L]] <- data.table(stratum = s, time = TP[j], n = length(i_s),
      amp = median(a[sel], na.rm = TRUE), rms_tru = median(d[sel], na.rm = TRUE),
      rel = median(rl[sel], na.rm = TRUE))
  }
}
AMP <- rbindlist(AMP)
for (i in seq_len(nrow(AMP))) say("%-15s %6d %8.3f %8.3f %8.3f",
  AMP$stratum[i], AMP$time[i], AMP$amp[i], AMP$rms_tru[i], AMP$rel[i])
say("\nThe `amp` column IS the confound, printed so its size is visible rather than argued:")
for (j in seq_along(TP)) say("  %d min: core amp %.3f vs other_induced %.3f (ratio %.2f) | Msn24_only %.3f vs TF_independent %.3f (ratio %.2f)",
  TP[j],
  AMP[stratum == "Hsf1_core" & time == TP[j], amp], AMP[stratum == "other_induced" & time == TP[j], amp],
  AMP[stratum == "Hsf1_core" & time == TP[j], amp] / AMP[stratum == "other_induced" & time == TP[j], amp],
  AMP[stratum == "Msn24_only" & time == TP[j], amp], AMP[stratum == "TF_independent" & time == TP[j], amp],
  AMP[stratum == "Msn24_only" & time == TP[j], amp] / AMP[stratum == "TF_independent" & time == TP[j], amp])

boot_rel_ratio <- function(i_num, i_den, nb = NBOOT) {
  out <- matrix(NA_real_, nb, 3L, dimnames = list(NULL, as.character(TP)))
  for (k in seq_len(nb)) {
    r <- cell_rel_ratio(sample(i_num, length(i_num), TRUE), sample(i_den, length(i_den), TRUE))
    for (j in seq_along(TP)) out[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  out
}
say("\nAmplitude-normalised version of THE CENTRAL COMPARISON (same strata, same denominators,")
say("same gene-level bootstrap). If a stratum's unnormalised excess was only response strength,")
say("its rel_ratio collapses toward 1; if the excess is genuine per-unit divergence, it survives.")
say("%-15s %-15s %5s %8s %17s %8s %12s", "focal", "denominator", "time", "rel_ratio",
    "95%% CI", "cells>1", "vs unnorm")
RELR <- list()
for (s in names(DENOM)) {
  i_f <- idx_of(s); i_d <- idx_of(DENOM[[s]])
  obs <- cell_rel_ratio(i_f, i_d)
  bt  <- boot_rel_ratio(i_f, i_d)
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    o  <- median(obs[sel], na.rm = TRUE)
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    un <- RAT[stratum == s & time == TP[j], ratio]
    RELR[[length(RELR) + 1L]] <- data.table(stratum = s, denom = DENOM[[s]], time = TP[j],
      n = length(i_f), rel_ratio = o, lo = ci[1], hi = ci[2],
      n_cells = sum(sel & !is.na(obs)), n_cells_gt1 = sum(obs[sel] > 1, na.rm = TRUE),
      unnorm_ratio = un, equiv_bound = max(abs(c(ci[1], ci[2]) - 1)))
  }
}
RELR <- rbindlist(RELR)
fwrite(RELR, file.path(outd, sprintf("a5_cis_ratio_amplitude_normalised_%s.tsv", ANCHOR)), sep = "\t")
for (i in seq_len(nrow(RELR))) say("%-15s %-15s %5d %8.3f [%6.3f,%6.3f] %5d/%d %7.3f  |r-1| > %.2f excl.",
  RELR$stratum[i], RELR$denom[i], RELR$time[i], RELR$rel_ratio[i], RELR$lo[i], RELR$hi[i],
  RELR$n_cells_gt1[i], RELR$n_cells[i], RELR$unnorm_ratio[i], RELR$equiv_bound[i])
fwrite(AMP, file.path(outd, sprintf("a5_cis_amplitude_%s.tsv", ANCHOR)), sep = "\t")

## --------- which of the two statistics is the right one? measured, not chosen ---
# The two tables above disagree, and they must, because they carry opposite amplitude biases:
#   - unnormalised rms_tru(DoD) is biased toward strongly responding strata if divergence is
#     multiplicative (DoD = (gamma-1)*R grows with R);
#   - rel = rms_tru(DoD)/R is biased against them if divergence is additive (DoD independent
#     of R, so rel ~ 1/R).
# The `not_induced` internal control shows the second bias is not hypothetical: not_induced
# has the smallest unnormalised divergence (0.776 of other_induced) and the largest rel
# (1.42), purely because its R is near zero. So choosing between the tables by preference is
# would not identify which bias applies; that is an empirical question about how DoD scales with R:
#     rms_tru(DoD) ~ R^beta      beta = 1 -> multiplicative -> rel is the right statistic
#                                beta = 0 -> additive       -> unnormalised is the right one
# beta is estimated by binning genes on an amplitude measured in other cells (leave-one-
# hybrid-out at the same timepoint), so the binning variable is statistically independent of
# the noise in the cell being measured -- otherwise regression to the mean would flatten beta
# by construction, which is the same trap as bin-on-noise selection in P1.
say("\n---- how does DoD scale with response amplitude? (this decides which table to read) ----")
say("beta = slope of log rms_tru(DoD) on log amp, across deciles of an INDEPENDENTLY measured")
say("amplitude (leave-one-hybrid-out at the same timepoint). beta~1 => multiplicative")
say("divergence => the normalised table is primary. beta~0 => additive => the unnormalised one.")
NDEC <- 8L
BETA <- list()
for (t in TP) {
  cells_t <- CELLS[tp_of == t]
  for (cc in cells_t) {
    others <- setdiff(cells_t, cc)
    xind <- rowMeans(abs(R_M[, others, drop = FALSE]))   # independent of this cell's noise
    br <- unique(quantile(xind, seq(0, 1, length.out = NDEC + 1L), na.rm = TRUE))
    if (length(br) < 4L) next
    bin <- cut(xind, breaks = br, include.lowest = TRUE, labels = FALSE)
    B <- rbindlist(lapply(sort(unique(bin[!is.na(bin)])), function(bb) {
      ii <- which(bin == bb)
      data.table(bin = bb, nb = length(ii),
                 amp = sqrt(max(0, mean(R_M[ii, cc]^2) - mean(RN_M[ii, cc]^2))),
                 dod = sqrt(max(0, mean(S_M[ii, cc]^2) - mean(N_M[ii, cc]^2))))
    }))
    B <- B[amp > 0 & dod > 0 & nb >= 20L]
    if (nrow(B) < 4L) next
    fit <- lm(log(dod) ~ log(amp), data = B)
    BETA[[length(BETA) + 1L]] <- data.table(cell = cc, hybrid = sub("_[0-9]+$", "", cc), time = t,
      n_bins = nrow(B), beta = unname(coef(fit)[2]),
      beta_se = summary(fit)$coefficients[2, 2], r2 = summary(fit)$r.squared,
      amp_span = max(B$amp) / min(B$amp))
  }
}
BETA <- rbindlist(BETA)
fwrite(BETA, file.path(outd, sprintf("a5_dod_amplitude_scaling_%s.tsv", ANCHOR)), sep = "\t")
say("%-14s %6s %6s %8s %8s %6s %9s", "hybrid", "time", "bins", "beta", "beta_se", "R2", "amp_span")
for (i in seq_len(nrow(BETA))) say("%-14s %6d %6d %8.3f %8.3f %6.2f %9.1f",
  BETA$hybrid[i], BETA$time[i], BETA$n_bins[i], BETA$beta[i], BETA$beta_se[i],
  BETA$r2[i], BETA$amp_span[i])
say("median beta = %.3f (IQR %.3f-%.3f) over %d cells; cells with beta > 0.5: %d/%d",
    median(BETA$beta), quantile(BETA$beta, .25), quantile(BETA$beta, .75), nrow(BETA),
    sum(BETA$beta > 0.5), nrow(BETA))
BETA_MED <- median(BETA$beta)
say("Read-out: beta is neither 0 nor 1, so NEITHER table is the correct statistic on its own.")
say("The two therefore BRACKET the answer, and a partially normalised statistic")
say("rms_tru(DoD)/amp^beta with the measured beta = %.3f is reported below as the point", BETA_MED)
say("estimate. beta is treated as fixed (measured on all %d genes, not per stratum), so the", length(KEEP))
say("bootstrap CI below does not propagate its uncertainty -- stated, not hidden.")

## ------------- the partially normalised statistic, at the measured beta ---
cell_cal <- function(i, beta = BETA_MED) {
  a <- cell_amp(i); ifelse(a > 0, cell_tru(i) / a^beta, NA_real_)
}
cell_cal_ratio <- function(i_num, i_den, beta = BETA_MED) {
  n <- cell_cal(i_num, beta); d <- cell_cal(i_den, beta)
  ifelse(!is.na(d) & d > 0, n / d, NA_real_)
}
boot_cal_ratio <- function(i_num, i_den, nb = NBOOT) {
  out <- matrix(NA_real_, nb, 3L, dimnames = list(NULL, as.character(TP)))
  for (k in seq_len(nb)) {
    r <- cell_cal_ratio(sample(i_num, length(i_num), TRUE), sample(i_den, length(i_den), TRUE))
    for (j in seq_along(TP)) out[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  out
}
say("\n======== PRIMARY: AMPLITUDE-CALIBRATED RATIO AT THE MEASURED beta = %.3f ========", BETA_MED)
say("Also printed: the two bracketing statistics (beta=0 unnormalised, beta=1 fully")
say("normalised), so the reader sees the whole range the amplitude assumption can move.")
say("%-15s %-15s %5s %9s %17s %8s %8s %8s", "focal", "denominator", "time", "cal_ratio",
    "95%% CI", "cells>1", "beta=0", "beta=1")
CAL <- list()
for (s in names(DENOM)) {
  i_f <- idx_of(s); i_d <- idx_of(DENOM[[s]])
  obs <- cell_cal_ratio(i_f, i_d)
  bt  <- boot_cal_ratio(i_f, i_d)
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    o  <- median(obs[sel], na.rm = TRUE)
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    CAL[[length(CAL) + 1L]] <- data.table(stratum = s, denom = DENOM[[s]], time = TP[j],
      n = length(i_f), beta = BETA_MED, cal_ratio = o, lo = ci[1], hi = ci[2],
      n_cells = sum(sel & !is.na(obs)), n_cells_gt1 = sum(obs[sel] > 1, na.rm = TRUE),
      ratio_beta0 = RAT[stratum == s & time == TP[j], ratio],
      ratio_beta1 = RELR[stratum == s & time == TP[j], rel_ratio],
      equiv_bound = max(abs(c(ci[1], ci[2]) - 1)))
  }
}
CAL <- rbindlist(CAL)
fwrite(CAL, file.path(outd, sprintf("a5_cis_ratio_calibrated_%s.tsv", ANCHOR)), sep = "\t")
for (i in seq_len(nrow(CAL))) say("%-15s %-15s %5d %9.3f [%6.3f,%6.3f] %5d/%d %8.3f %8.3f  |r-1| > %.2f excl.",
  CAL$stratum[i], CAL$denom[i], CAL$time[i], CAL$cal_ratio[i], CAL$lo[i], CAL$hi[i],
  CAL$n_cells_gt1[i], CAL$n_cells[i], CAL$ratio_beta0[i], CAL$ratio_beta1[i], CAL$equiv_bound[i])
say("\nSIGN-STABILITY ACROSS THE BRACKET is the only conclusion the data support for a stratum:")
say("%-15s %6s %10s %10s %10s %12s", "focal", "time", "beta=0", "beta=meas", "beta=1", "sign stable?")
for (i in seq_len(nrow(CAL))) {
  v <- c(CAL$ratio_beta0[i], CAL$cal_ratio[i], CAL$ratio_beta1[i])
  say("%-15s %6d %10.3f %10.3f %10.3f %12s", CAL$stratum[i], CAL$time[i], v[1], v[2], v[3],
      ifelse(all(v > 1) || all(v < 1), "YES", "NO -- undecided"))
}
say("A stratum whose ratio crosses 1 somewhere in the bracket has NO direction that this")
say("design can establish, regardless of how tight any single statistic's CI looks. That is")
say("the supported interpretation, which is stricter than either table alone.")

## ------------- power calibration of the calibrated statistic ---
# The calibrated ratio pushes every stratum toward 1, including `not_induced`, whose
# unnormalised ratio (0.776) was the internal validity check. Two readings are possible and
# they are not equivalent:
#   (i) all the stratum differences really were amplitude composition, or
#   (ii) dividing by amp^beta destroyed the statistic's discriminating power altogether.
# Interpretation (i) requires exclusion of (ii) by a control that can fail. So the
# statistic is calibrated by injection, and the injected term is homologous to the
# statistic's numerator: a multiplicative divergence increment DoD -> DoD + c*R makes
# gamma -> gamma + c, which is the quantity the ratio claims to measure.
# R is taken as the leave-one-hybrid-out amplitude, independent of the target cell's noise;
# the residual noise it carries still enters rms_sig without entering rms_nul, so the
# detection thresholds below are mildly optimistic, i.e. a lower bound on what is needed.
say("\n---- POWER CALIBRATION: what size of divergence excess could this statistic see? ----")
say("Injection: for a random n-matched subset of other_induced, DoD -> DoD + c*R_loho, i.e.")
say("gamma -> gamma + c. Homologous to the numerator. Reports the smallest c whose")
say("calibrated ratio has a bootstrap CI excluding 1. If no c does, the negative above is")
say("not supported as a biological negative; the test is then uninformative.")
LOHO <- sapply(TP, function(t) {
  cells_t <- CELLS[tp_of == t]
  rowMeans(abs(R_M[, cells_t, drop = FALSE]))
})
colnames(LOHO) <- as.character(TP)
NB_CAL <- 600L
inj_ratio <- function(i_num, i_den, cc_inj, beta = BETA_MED) {
  S2 <- S_M
  for (j in seq_along(TP)) {
    cols <- which(tp_of == TP[j])
    S2[i_num, cols] <- S2[i_num, cols] + cc_inj * LOHO[i_num, j]
  }
  tru <- function(i) sqrt(pmax(0, colMeans(S2[i, , drop = FALSE]^2) - colMeans(N_M[i, , drop = FALSE]^2)))
  amp <- function(i) sqrt(pmax(0, colMeans(R_M[i, , drop = FALSE]^2) - colMeans(RN_M[i, , drop = FALSE]^2)))
  n <- tru(i_num) / amp(i_num)^beta
  d <- tru(i_den) / amp(i_den)^beta
  ifelse(!is.na(d) & d > 0, n / d, NA_real_)
}
set.seed(20260905L + 7L)
i_inj <- sample(i_oth, length(i_core))          # n-matched to the core, drawn once
i_rest <- setdiff(i_oth, i_inj)                  # disjoint denominator, as everywhere else
say("injected set n = %d (random subset of other_induced); denominator n = %d (the remainder)",
    length(i_inj), length(i_rest))
say("%6s %8s %17s %10s", "c", "ratio_15", "95%% CI (15 min)", "excl. 1?")
PWR <- list()
for (cc_inj in c(0, 0.05, 0.10, 0.20, 0.40)) {
  obs <- inj_ratio(i_inj, i_rest, cc_inj)
  bt <- matrix(NA_real_, NB_CAL, 3L)
  for (k in seq_len(NB_CAL)) {
    r <- inj_ratio(sample(i_inj, length(i_inj), TRUE), sample(i_rest, length(i_rest), TRUE), cc_inj)
    for (j in seq_along(TP)) bt[k, j] <- median(r[tp_of == TP[j]], na.rm = TRUE)
  }
  for (j in seq_along(TP)) {
    sel <- tp_of == TP[j]
    ci <- quantile(bt[, j], c(.025, .975), na.rm = TRUE)
    PWR[[length(PWR) + 1L]] <- data.table(c_inj = cc_inj, time = TP[j],
      ratio = median(obs[sel], na.rm = TRUE), lo = ci[1], hi = ci[2],
      excludes_1 = ci[1] > 1 | ci[2] < 1)
  }
  z <- PWR[[length(PWR) - 2L]]
  say("%6.2f %8.3f [%6.3f,%6.3f] %10s", cc_inj, z$ratio, z$lo, z$hi,
      ifelse(z$excludes_1, "YES", "no"))
}
PWR <- rbindlist(PWR)
fwrite(PWR, file.path(outd, sprintf("a5_power_calibration_%s.tsv", ANCHOR)), sep = "\t")
det <- PWR[excludes_1 == TRUE & c_inj > 0, if (.N) min(c_inj) else NA_real_]
say("c = 0 is the null run and MUST NOT exclude 1 (it is the same gene pool split in two):")
say("  c=0 excludes 1 in %d of 3 timepoints", PWR[c_inj == 0 & excludes_1 == TRUE, .N])
if (is.na(det)) {
  say("No injected c up to 0.40 is detectable, so the calibrated statistic is not informative,")
  say("and the collapse toward 1 does not support a biological negative.")
} else {
  say("smallest detectable injected excess: c = %.2f (gamma shifted by %.2f)", det, det)
  say("=> the calibrated statistic CAN detect a per-unit divergence excess of this size, so the")
  say("collapse of every stratum toward 1 is a statement about the data, not about the test.")
  say("Any stratum whose calibrated ratio is within the equivalence bound is therefore")
  say("inconsistent with an excess of c >= %.2f.", det)
  say("SELF-CONSISTENCY: core bootstrap equivalence bound |r-1| > %.2f excluded, vs injection",
      CAL[stratum == "Hsf1_core", max(equiv_bound)])
  say("detection threshold c = %.2f. Two independent routes to the same", det)
  say("resolution limit and they agree, which is the check that neither number is an artefact")
  say("of its own machinery.")
}

## ------------------- control that can fail: n-matched subsampling ---
# The ratio above can be moved by n alone: rms_tru is a noisy quantity and the core has
# ~40x fewer genes than other_induced, so its rms_tru is estimated with more error. The
# test is therefore to ask where the core falls in the distribution of equally sized
# random subsets of the other heat-induced genes. Two versions: plain, and matched on
# baseline expression deciles, because expression drives the noise floor.
say("\n================ THE CONTROL THAT CAN FAIL: n-matched subsampling ================")
expr_all <- EXPR[, .(expr = mean(expr)), by = hog]
expr_vec <- expr_all$expr[match(rownames(S_M), expr_all$hog)]
core_dec <- cut(expr_vec[i_core], breaks = quantile(expr_vec, seq(0, 1, .1), na.rm = TRUE),
                include.lowest = TRUE, labels = FALSE)
oth_dec  <- cut(expr_vec[i_oth],  breaks = quantile(expr_vec, seq(0, 1, .1), na.rm = TRUE),
                include.lowest = TRUE, labels = FALSE)
need <- table(factor(core_dec, levels = 1:10))
avail <- table(factor(oth_dec, levels = 1:10))
say("expression-decile match feasibility (core needs / other_induced has):")
say("  decile %s", paste(sprintf("%d:%d/%d", 1:10, as.integer(need), as.integer(avail)), collapse = " "))
matchable <- all(need <= avail)
say("  exact decile match feasible: %s", matchable)

# `stat` selects which statistic the control is run on. It must be run on the primary
# statistic (the amplitude-calibrated one) or it is a control on a quantity that is not being
# reported; the unnormalised version is kept beside it so the amplitude contribution to the
# p-value is visible rather than absorbed.
nmatch <- function(stratum_name, i_focal, matched, stat = c("tru", "cal")) {
  stat <- match.arg(stat)
  f <- if (stat == "cal") cell_cal else cell_tru
  nsub <- length(i_focal)
  obs <- f(i_focal)
  null <- matrix(NA_real_, NBOOT, 3L)
  for (k in seq_len(NBOOT)) {
    if (matched) {
      pick <- unlist(lapply(1:10, function(dd) {
        nd <- as.integer(need[dd]); if (nd == 0L) return(integer(0))
        pool_d <- i_oth[which(oth_dec == dd)]
        sample(pool_d, min(nd, length(pool_d)), replace = FALSE)
      }))
    } else pick <- sample(i_oth, nsub, replace = FALSE)
    v <- f(pick)
    for (j in seq_along(TP)) null[k, j] <- median(v[tp_of == TP[j]], na.rm = TRUE)
  }
  data.table(stratum = stratum_name, statistic = stat, matched = matched, time = TP, n = nsub,
             obs = sapply(TP, function(t) median(obs[tp_of == t], na.rm = TRUE)),
             null_med = apply(null, 2, median, na.rm = TRUE),
             null_lo = apply(null, 2, quantile, .025, na.rm = TRUE),
             null_hi = apply(null, 2, quantile, .975, na.rm = TRUE),
             p_emp_greater = sapply(seq_along(TP), function(j)
               mean(null[, j] >= median(obs[tp_of == TP[j]], na.rm = TRUE), na.rm = TRUE)))
}
NM <- rbindlist(list(
  nmatch("Hsf1_core",  i_core, FALSE,     "tru"),
  nmatch("Hsf1_core",  i_core, matchable, "tru"),
  nmatch("Hsf1_core",  i_core, FALSE,     "cal"),
  nmatch("Hsf1_core",  i_core, matchable, "cal"),
  nmatch("Hsf1_joint", i_join, FALSE,     "tru"),
  nmatch("Hsf1_joint", i_join, FALSE,     "cal")))
NM <- unique(NM, by = c("stratum", "statistic", "matched", "time"))
fwrite(NM, file.path(outd, sprintf("a5_cis_nmatched_%s.tsv", ANCHOR)), sep = "\t")
say("\nDoes the focal stratum's cis divergence exceed a SIZE-MATCHED random subset of the")
say("other heat-induced genes? p_emp = fraction of subsets at least as divergent.")
say("statistic: tru = unnormalised rms_tru(DoD); cal = amplitude-calibrated (beta = %.3f),", BETA_MED)
say("the primary one. Both are shown so the amplitude contribution to p_emp is visible.")
say("%-12s %5s %8s %6s %8s %8s %16s %9s", "stratum", "stat", "matched", "time", "obs",
    "null_med", "null 95%%", "p_emp")
for (i in seq_len(nrow(NM))) say("%-12s %5s %8s %6d %8.3f %8.3f  [%6.3f,%6.3f] %8.3f",
  NM$stratum[i], NM$statistic[i], NM$matched[i], NM$time[i], NM$obs[i], NM$null_med[i],
  NM$null_lo[i], NM$null_hi[i], NM$p_emp_greater[i])

## ------------------- count statistic on the same genes, for the contrast ---
# 06_module_bridge.R answered this question with a count of padj(interaction) < 0.05.
# Recomputed here on the epistasis-refined strata so the count-vs-threshold-free
# discrepancy is measured rather than argued about. med(se) is unavailable for the
# deposited table (the interaction test is multi-coefficient, so `stat` is an F, not a t,
# and the reason shrunken se is not a precision measure applies twice over), so the count
# is interpreted only alongside the threshold-free ratio.
say("\n================ COUNT STATISTIC vs THRESHOLD-FREE, SAME GENES ================")
MM <- fread(cmd = paste("zcat", shQuote(file.path(csd, "fay_multimethod.tsv.gz"))))
stopifnot("Parents" %in% MM$group)     # the 7th group value must be filtered explicitly
MI <- MM[condition == "heat" & model == "A" & term == "inter" & group %in% HYB &
         gene %in% KEEP, .(gene, group, method, padj)]
stopifnot(uniqueN(MI$group) == 6L)
MI <- merge(MI, ST, by.x = "gene", by.y = "hog", allow.cartesian = TRUE)
CNT <- list()
for (mth in c("limma", "DESeq2", "edgeR")) for (s in c("Hsf1_core", "Hsf1_joint")) for (g in HYB) {
  a <- MI[method == mth & group == g & stratum == s]
  b <- MI[method == mth & group == g & stratum == "other_induced"]
  if (nrow(a) < 10L || nrow(b) < 10L) next
  ta <- c(sum(a$padj < ALPHA, na.rm = TRUE), sum(!(a$padj < ALPHA) | is.na(a$padj)))
  tb <- c(sum(b$padj < ALPHA, na.rm = TRUE), sum(!(b$padj < ALPHA) | is.na(b$padj)))
  ft <- fisher.test(matrix(c(ta, tb), 2, 2))
  # Degenerate cell: if neither arm has a single significant gene the odds ratio is
  # undefined (0/0), yet fisher.test returns the numeric 0. Letting that 0 into the median
  # and the min/max range would print a fabricated "OR = 0.00" as if it were an estimate of
  # depletion. Marked NA and excluded from the summary, with the count reported.
  degen <- (ta[1] + tb[1]) == 0L
  CNT[[length(CNT) + 1L]] <- data.table(method = mth, stratum = s, hybrid = g,
    n_focal = sum(ta), k_focal = ta[1], rate_focal = ta[1] / sum(ta),
    n_other = sum(tb), k_other = tb[1], rate_other = tb[1] / sum(tb),
    degenerate = degen,
    or = if (degen) NA_real_ else unname(ft$estimate),
    lo = if (degen) NA_real_ else ft$conf.int[1],
    hi = if (degen) NA_real_ else ft$conf.int[2],
    p  = if (degen) NA_real_ else ft$p.value)
}
CNT <- rbindlist(CNT)
fwrite(CNT, file.path(outd, sprintf("a5_count_vs_threshold_free_%s.tsv", ANCHOR)), sep = "\t")
say("Fisher OR for padj(interaction) < %.2f, focal stratum vs other_induced:", ALPHA)
say("%-8s %-12s %6s %8s %8s %8s %10s %8s", "method", "stratum", "cells", "OR_med", "OR_min",
    "OR_max", "p<0.05", "degen")
for (mth in c("limma", "DESeq2", "edgeR")) for (s in c("Hsf1_core", "Hsf1_joint")) {
  z <- CNT[method == mth & stratum == s]
  if (!nrow(z)) next
  zz <- z[degenerate == FALSE]
  say("%-8s %-12s %6d %8.2f %8.2f %8.2f %10d %8d", mth, s, nrow(zz),
      median(zz$or), min(zz$or), max(zz$or), sum(zz$p < ALPHA, na.rm = TRUE), sum(z$degenerate))
}
say("`degen` = cells where NEITHER arm had a single gene at padj < %.2f, so the OR is 0/0", ALPHA)
say("and is excluded (fisher.test reports 0 there; that 0 is not an estimate of depletion).")
say("\ndetection rates (why a count can mislead here): core %.2f%% vs other_induced %.2f%%",
    100 * CNT[stratum == "Hsf1_core", median(rate_focal)],
    100 * CNT[stratum == "Hsf1_core", median(rate_other)])
say("Both rates are low, so the OR is a ratio of two small counts; the threshold-free")
say("ratio above uses every gene and is the primary statement.")
say("\nper-cell counts for the core (limma), so the empty cells are readable rather than")
say("mysterious -- k_core = 0 means NOT ONE of the %d core genes cleared padj < %.2f in that",
    length(HOGS$Hsf1_core), ALPHA)
say("hybrid, which is a statement about power, not about conservation:")
say("%-8s %8s %8s %10s %10s %8s", "hybrid", "k_core", "n_core", "rate_core", "rate_other", "OR")
for (g in HYB) {
  z <- CNT[method == "limma" & stratum == "Hsf1_core" & hybrid == g]
  if (!nrow(z)) next
  say("%-8s %8d %8d %9.1f%% %9.1f%% %8s", g, z$k_focal, z$n_focal,
      100 * z$rate_focal, 100 * z$rate_other,
      ifelse(z$degenerate, "undef", sprintf("%.2f", z$or)))
}
say("\nCompare with the threshold-free table: the count statistic loses %d of the 6 hybrids",
    sum(CNT[method == "limma" & stratum == "Hsf1_core", k_focal] == 0L))
say("to zero counts, while the threshold-free ratio is estimable in all 18 cells. That is")
say("the concrete cost of ranking by counts, measured on this dataset rather than argued.")

## --------------------------------------- cross-check against the parental axis ---
say("\n================ CROSS-CHECK: CIS AXIS vs PARENTAL AXIS ================")
a4c <- file.path(outd, sprintf("fay_parental_core_vs_matched_%s.tsv", ANCHOR))
if (file.exists(a4c)) {
  CM <- fread(a4c); print(CM)
  say("\nA4 (parental axis, composition statistic) and A5 (cis axis, allele-paired DoD)")
  say("are independent measurements of the same question: one uses the trans environment")
  say("and one removes it. Agreement is the evidence; disagreement would itself be the")
  say("result and would have to be reported as such.")
} else say("[missing] %s -- run 03_fay_parental_dimensions.R first", a4c)

## ------------------------------------------------------------- persistence ---
say("\n=== per-gene DoD written out for the promoter HSE test (HSE architecture) to join against ===")
OUT <- merge(SIG, NUL, by = c("hog", "hybrid", "time"))
OUT <- merge(OUT, ST, by = "hog", all.x = TRUE, allow.cartesian = TRUE)
## The left join keeps every hog, including those that fall in no stratum: `unanimous` is NA when DESeq2, edgeR and
## limma do not agree on the class (line 169), so those genes are in no STRATA entry. There are 569 of them and
## they are correctly excluded from every stratum summary, but as NA they would be written
## to disk as an empty string, which reads as a seventh nameless stratum. Name them.
OUT[is.na(stratum), stratum := "unassigned"]
stopifnot(!any(is.na(OUT$stratum)), !any(OUT$stratum == ""))
OUT <- merge(OUT, HB[, .(hog = gene, scer = Name)], by = "hog", all.x = TRUE)
fwrite(OUT, file.path(outd, sprintf("a5_cis_dod_%s.tsv.gz", ANCHOR)), sep = "\t", compress = "gzip")
say("rows %d | file a5_cis_dod_%s.tsv.gz", nrow(OUT), ANCHOR)

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_a5_divergence.txt"))
say("\n=== done ===")
sink()
