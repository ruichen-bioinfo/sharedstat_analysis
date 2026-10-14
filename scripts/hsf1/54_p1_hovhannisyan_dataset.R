#!/usr/bin/env Rscript
## 54_p1_hovhannisyan_dataset.R -- the reliability h and the shared-reference artefact, measured on a third
## dataset.
##
## Background. P1's argument is that the reliability of a shared reference arm can be measured from the
## replicate structure of the design being analysed, and that at the reliabilities real designs have, the naive
## estimator inverts the sign of small divergences. The first two datasets show this. Hovhannisyan et al. 2020
## adds a third from another laboratory, another species pair, three biological replicates instead of two, and
## the opposite temperature direction. Nothing about the argument is specific to heat or to the Fay design, and
## a third dataset is the cheapest way to show that.
##
## What h is here. h is the reliability of the quantity that appears on both sides of the shared-reference
## regression: one allele's own response. In this design the response is log2 of the 12 C expression over the
## 30 C expression within one allele, and three biological replicates give three estimates of it, so h follows
## from the same signal-and-null construction the package uses, generalised to three replicates exactly as in
## the cold-arm decomposition:
##     sig = (d1+d2+d3)/3 ; n1 = (d2-d1)/sqrt(6) ; n2 = (2*d3-d1-d2)/sqrt(18)
## all with variance sigma^2/3, sig carrying the signal and n1, n2 carrying none. Then
##     h = 1 - mean(c(mean(n1^2), mean(n2^2))) / mean(sig^2), each contrast centred over features first,
## which reduces to `reliability_from_blocks(..., center = TRUE)` when there are two replicates. FIDELITY GATE: on two of the three
## replicates the function in the package must return the same number as this code, so the generalisation is
## checked against the published implementation rather than asserted to match it.
##
## The artefact. The naive shared-reference estimator regresses a contrast on the arm it contains. Here the
## contrast is the allelic difference in response inside the hybrid, cis = r(Su allele) - r(Sc allele), and the
## axis is the Sc allele's own response r(Sc). Since cis contains -r(Sc), the slope of cis on r(Sc) is biased
## downward by (1-h) even when the two alleles diverge by a constant factor. The instrumented estimator uses a
## quantity correlated with the reference's signal but not its error: here, the Sc parent's response, measured
## in different libraries from the hybrid.
##
## Decision rule, written before the numbers were seen:
##   Expected: h in the same range as the Fay design (median about 0.7, some cells much lower); the naive slope
##     below the instrumented one by roughly (1-h) times the gain; the instrumented estimate closer to the
##     parental gain than the naive one.
##   Failable: h near 1 in every arm, which would make the artefact negligible here and would mean the two-block
##     Fay reliabilities are unusually low rather than typical. That outcome would limit how
##     far P1's numbers generalise.

suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hov  <- file.path(root, "data/02_processed/hov2020")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
MINCPM <- 1; MINLIB <- 4L; PRIOR <- 2; NBOOT <- 1000L; SEED <- 20260913L

D <- fread(file.path(hov, "design.tsv"))
LENT <- fread(file.path(hov, "cds_lengths.tsv"))
KEEP <- LENT[len_frac < 0.05, sys]
say("genes with alleles of comparable CDS length: %d of %d", length(KEEP), nrow(LENT))

read_counts <- function(al) {
  q <- fread(file.path(hov, "quant", al, "quant.sf"))
  q[, c("sys", "allele") := tstrsplit(Name, "|", fixed = TRUE)]
  q[sys %in% KEEP, .(sys, allele, n = NumReads)]
}
LIB <- list()
for (i in seq_len(nrow(D))) {
  al <- D$alias[i]; g <- D$geno[i]; q <- read_counts(al)
  if (g == "SCSU") {
    LIB[[paste0(al, "|Sc")]] <- q[allele == "Sc", .(sys, n)]
    LIB[[paste0(al, "|Su")]] <- q[allele == "Su", .(sys, n)]
  } else {
    own <- if (g == "SC") "Sc" else "Su"
    LIB[[paste0(al, "|", own)]] <- q[allele == own, .(sys, n)]
  }
}
## Counts are TMM-normalised and analysed as log-CPM, as for the Fay hybrids (P1 Methods). Each allele of a
## hybrid library is a separate column, so the normalisation factors are computed over all allele columns together.
CM <- dcast(rbindlist(lapply(names(LIB), function(k) LIB[[k]][, .(sys, key = k, n)])), sys ~ key, value.var = "n", fill = 0)
CMm <- as.matrix(CM[, -1]); rownames(CMm) <- CM$sys
DG <- calcNormFactors(DGEList(CMm), method = "TMM")
CPMm <- cpm(DG, log = FALSE); LCPMm <- cpm(DG, log = TRUE, prior.count = PRIOR)
CPM <- merge(melt(as.data.table(CPMm, keep.rownames = "sys"), id.vars = "sys", variable.name = "key", value.name = "cpm"),
             melt(as.data.table(LCPMm, keep.rownames = "sys"), id.vars = "sys", variable.name = "key", value.name = "lcpm"),
             by = c("sys", "key"))
CPM[, key := as.character(key)]
CPM[, c("alias", "allele") := tstrsplit(key, "|", fixed = TRUE)]
CPM[, c("geno", "tempc", "rep") := tstrsplit(alias, "_", fixed = TRUE)]
CPM[, temp := as.integer(sub("C$", "", tempc))]
ok_h <- CPM[geno == "SCSU", .(k = sum(cpm >= MINCPM)), by = sys][k >= MINLIB, sys]
ok_p <- CPM[geno != "SCSU", .(k = sum(cpm >= MINCPM)), by = sys][k >= MINLIB, sys]
UNIV <- intersect(intersect(ok_h, ok_p), KEEP)
say("expression filter: %d genes", length(UNIV))
stopifnot(length(UNIV) > 2000L)

## ------------------------------------------------------------------ per-arm responses ----
## An "arm" is one allele in one genotype: the Sc parent, the Su parent, and each of the two alleles inside the
## hybrid. Temperature is consumed by the response, so it does not index an arm.
W <- dcast(CPM[sys %in% UNIV], sys + geno + allele ~ temp + rep, value.var = "lcpm")
ARMS <- unique(W[, .(geno, allele)])
say("\narms: %s", paste(sprintf("%s/%s", ARMS$geno, ARMS$allele), collapse = ", "))

resp3 <- function(geno_, allele_) {
  z <- W[geno == geno_ & allele == allele_]
  d <- data.table(sys = z$sys)
  for (i in 1:3) d[[paste0("d", i)]] <- z[[sprintf("12_rep%d", i)]] - z[[sprintf("30_rep%d", i)]]
  d[complete.cases(d)]
}
## The three-replicate reliability comes from the package, not from a formula written here, so the
## third-dataset calculation is performed by the package. The check below runs
## in that direction: the package function must reproduce the inline formula, not the other way round.
source(file.path(root, "pkg/sharedstat/R/reliability_replicates.R"))
source(file.path(root, "pkg/sharedstat/R/reliability.R"))
h_from3 <- function(d) {
  # centred over features, as in P1 (Methods): h is a ratio of variances
  r <- reliability_from_replicates(as.matrix(d[, .(d1, d2, d3)]), center = TRUE)
  list(h = r$h, snr = r$snr, n = r$n_pairs, rms_sig = r$rms_signal, rms_nul = r$rms_null,
       h_pairwise = r$h_pairwise)
}
h_from3_inline <- function(d) {          # the inline formula, kept only for the fidelity check
  sig <- (d$d1 + d$d2 + d$d3) / 3
  n1  <- (d$d2 - d$d1) / sqrt(6)
  n2  <- (2 * d$d3 - d$d1 - d$d2) / sqrt(18)
  sig <- sig - mean(sig); n1 <- n1 - mean(n1); n2 <- n2 - mean(n2)
  1 - mean(c(mean(n1^2), mean(n2^2))) / mean(sig^2)
}

## ------------------------------------------------------------------ implementation validation ----
say("\n=== fidelity gate: on two replicates this code must equal the package function ===")
rel2 <- function(b1, b2) reliability_from_blocks(b1, b2, center = TRUE)$h
d_test <- resp3("SC", "Sc")
## 1. two columns of the package function must equal the published two-block function
h_two_rep <- reliability_from_replicates(as.matrix(d_test[, .(d1, d2)]), center = TRUE)$h
h_two_blk <- reliability_from_blocks(d_test$d1, d_test$d2, center = TRUE)$h
say("  two replicates: reliability_from_replicates %.12f, reliability_from_blocks %.12f, difference %.2e",
    h_two_rep, h_two_blk, abs(h_two_rep - h_two_blk))
stopifnot(abs(h_two_rep - h_two_blk) < 1e-12)
## 2. three columns of the package function must equal the inline formula
h_three_pkg <- h_from3(d_test)$h
h_three_inl <- h_from3_inline(d_test)
say("  three replicates: package %.12f, the formula this script carried %.12f, difference %.2e",
    h_three_pkg, h_three_inl, abs(h_three_pkg - h_three_inl))
stopifnot(abs(h_three_pkg - h_three_inl) < 1e-12)
say("  both implementation checks passed: the package reduces to the two-block form and reproduces the")
say("  three-replicate values this analysis reports.")

## ------------------------------------------------------------------ h per arm ----
say("\n=== reliability of each arm's response, three replicates ===")
set.seed(SEED)
H <- rbindlist(lapply(seq_len(nrow(ARMS)), function(i) {
  g <- ARMS$geno[i]; a <- ARMS$allele[i]
  d <- resp3(g, a); r <- h_from3(d)
  hb <- replicate(NBOOT, { s <- sample(nrow(d), nrow(d), TRUE); h_from3(d[s])$h })
  ## all three two-replicate subsets, so the three-replicate value can be compared with the pairs
  pairs <- combn(3, 2)
  h2 <- apply(pairs, 2, function(p) rel2(d[[paste0("d", p[1])]], d[[paste0("d", p[2])]]))
  data.table(arm = sprintf("%s/%s", g, a), n = r$n, h = r$h, snr = r$snr,
             lo = quantile(hb, .025), hi = quantile(hb, .975),
             h_pairs = paste(sprintf("%.3f", h2), collapse = "; "),
             inv_gain_threshold = 1 / r$h)
}))
print(H[, .(arm, n, h = round(h, 4), snr = round(snr, 2), lo = round(lo, 3), hi = round(hi, 3),
            h_pairs, threshold_gain = round(inv_gain_threshold, 3))], row.names = FALSE)
fwrite(H, file.path(hsd, "p1x_hov2020_reliability.tsv"), sep = "\t")

say("\n  the reversal threshold is a true gain of 1/h; a true gain below it is reported as a loss")
for (i in seq_len(nrow(H))) say("    %-12s h = %.3f  =>  gains up to %.2f-fold invert", H$arm[i], H$h[i],
                               H$inv_gain_threshold[i])

## ------------------------------------------------------------------ the artefact ----
## The naive estimator: regress the allelic response difference inside the hybrid on the Sc allele's own
## response. The instrument: the Sc parent's response, from different libraries.
say("\n=== the shared-reference artefact on this dataset ===")
hSc <- resp3("SCSU", "Sc"); hSu <- resp3("SCSU", "Su"); pSc <- resp3("SC", "Sc")
M <- merge(merge(hSc[, .(sys, hsc = (d1 + d2 + d3) / 3)],
                 hSu[, .(sys, hsu = (d1 + d2 + d3) / 3)], by = "sys"),
           pSc[, .(sys, psc = (d1 + d2 + d3) / 3)], by = "sys")
say("  genes with all three arms: %d", nrow(M))
M[, cis := hsu - hsc]
naive <- 1 + cov(M$cis, M$hsc) / var(M$hsc)
iv    <- 1 + cov(M$cis, M$psc) / cov(M$hsc, M$psc)
hh    <- H[arm == "SCSU/Sc", h]
say("  naive gain, regressing the allelic difference on the reference arm itself : %.4f", naive)
say("  instrumented gain, using the Sc parent's response as the instrument        : %.4f", iv)
say("  the closed form predicts naive = gain x h; with h = %.3f that is %.4f from the", hh, iv * hh)
say("  instrumented gain, against %.4f observed", naive)
set.seed(SEED + 1L)
bt <- replicate(NBOOT, {
  s <- sample(nrow(M), nrow(M), TRUE); m <- M[s]
  c(1 + cov(m$cis, m$hsc) / var(m$hsc), 1 + cov(m$cis, m$psc) / cov(m$hsc, m$psc))
})
say("  naive        %.4f [%.4f, %.4f]", naive, quantile(bt[1, ], .025), quantile(bt[1, ], .975))
say("  instrumented %.4f [%.4f, %.4f]", iv, quantile(bt[2, ], .025), quantile(bt[2, ], .975))
A <- data.table(dataset = "Hovhannisyan 2020", n = nrow(M), h_reference_arm = hh,
                naive_gain = naive, naive_lo = quantile(bt[1, ], .025), naive_hi = quantile(bt[1, ], .975),
                iv_gain = iv, iv_lo = quantile(bt[2, ], .025), iv_hi = quantile(bt[2, ], .975),
                predicted_naive = iv * hh)
A[, artefact_share := (iv - naive) / abs(iv - 1)]
fwrite(A, file.path(hsd, "p1x_hov2020_artefact.tsv"), sep = "\t")

say("\n=== results for the pre-specified criteria ===")
say("  h across the four arms: %.3f to %.3f, median %.3f", min(H$h), max(H$h), median(H$h))
say("  is h near 1 in every arm (the failable outcome): %s", all(H$h > 0.95))
say("  naive below instrumented: %s (%.4f against %.4f)", naive < iv, naive, iv)
say("  closed-form prediction within %.4f of the observed naive value", abs(iv * hh - naive))
FAY <- fread(file.path(hsd, "p1_measured_reliability.tsv"))
say("\n  for comparison, the Fay design: %d cells, h from %.3f to %.3f, median %.3f",
    nrow(FAY), min(FAY$h), max(FAY$h), median(FAY$h))
say("")

## ------------------------------------------------------------------ Table 1 row, same formula ----
## The artefact percentage is 100 * |artefact| / |slope| with slope = naive gain - 1, taken from
## the formula implemented in script 41. Note: dividing by the gain instead of the slope gives
## 4.3% instead of 28.4%, since the gain is near 1 and the slope is not.
say("\n=== artefact percentage for this dataset, using script 41's formula ===")
t41 <- parse(file.path(root, "scripts/hsf1/41_table1_per_hybrid.R"))
pct_formula <- "100 * abs(artefact) / abs(slope)"
src41 <- paste(readLines(file.path(root, "scripts/hsf1/41_table1_per_hybrid.R")), collapse = "\n")
stopifnot(grepl("artefact_pct := 100 * abs(artefact) / abs(slope)", src41, fixed = TRUE))
stopifnot(grepl("naive := 1 + slope", src41, fixed = TRUE))
say("  formula from script 41: artefact_pct = %s, naive = 1 + slope", pct_formula)
slope_naive <- naive - 1
artefact_h  <- iv - naive
pct_h       <- 100 * abs(artefact_h) / abs(slope_naive)
say("  naive gain %.4f  =>  slope %.4f", naive, slope_naive)
say("  corrected gain %.4f  =>  artefact %.4f", iv, artefact_h)
say("  artefact as %% of naive slope: %.1f", pct_h)
say("  (dividing by the gain instead of the slope would have given %.1f -- a different denominator)",
    100 * abs(artefact_h) / abs(naive))
T1 <- data.table(hybrid = "S. cerevisiae x S. uvarum (Hovhannisyan)", naive = naive, gain_iv = iv,
                 slope = slope_naive, artefact = artefact_h, artefact_pct = pct_h,
                 h_reference_arm = hh, n_genes = nrow(M), n_rep = 3L)
fwrite(T1, file.path(hsd, "p1x_hov2020_table1_row.tsv"), sep = "\t")
say("  wrote p1x_hov2020_table1_row.tsv")

## Range check against the six Fay hybrids.
F1 <- fread(file.path(root, "data/02_processed/crossspecies/gain_iv.tsv"))[model == "A"]
FT <- F1[, .(slope = median(slope_ols), gain_iv = median(gain_iv), artefact = median(artefact)), by = group]
FT[, pct := 100 * abs(artefact) / abs(slope)]
say("\n  Fay hybrids: artefact %.1f%% to %.1f%% of the naive slope", min(FT$pct), max(FT$pct))
say("  this dataset: %.1f%%  =>  inside that range: %s", pct_h, pct_h >= min(FT$pct) && pct_h <= max(FT$pct))
say("  Fay reliabilities: %d cells, median %.3f, from %.3f to %.3f", nrow(FAY), median(FAY$h), min(FAY$h), max(FAY$h))
say("  this dataset: 4 arms, median %.3f, from %.3f to %.3f", median(H$h), min(H$h), max(H$h))
say("  every arm at or below the inversion threshold for a 10%% gain (1/1.1 = %.6f): %s",
    1/1.1, all(H$h <= 1/1.1))
say("  highest arm sits %.5f below that threshold", 1/1.1 - max(H$h))
