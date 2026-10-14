#!/usr/bin/env Rscript
## 53_p3_r3_per_hybrid.R -- does the P3 cis/trans time trend hold in each hybrid separately?
##
## Background. In the cold-arm decomposition (script 52), leave-one-replicate-out shows that a gene-level
## bootstrap interval says nothing about replicate variance: the cis share moves from 0.337 to 0.555 across
## three subsets while the bootstrap interval is [0.405, 0.442]. The P3 time-trend intervals are gene-level
## bootstraps too, and the Fay design has two replicate blocks, which cannot be left out one at a time. The Fay
## data do have six pairwise hybrids, which are independent biological units. If the time trend holds in each
## of them separately, that is replicate-level evidence of the kind the cold arm cannot provide.
##
## Decision rule, fixed before running:
##   (1) the cis share falls with time in at least 5 of the 6 hybrids;
##   (2) trans exceeds cis at 60 min in at least 5 of the 6.
##   If both hold, the P3 time-trend result is reported with the count. If either fails, its statement
##   is changed to the count that actually holds. No hybrid is excluded, and the criterion is not adjusted
##   after the counts are seen.
##
## The share is the variance share of script 05, rms_tru(cis)^2 / (rms_tru(cis)^2 + rms_tru(trans)^2), so
## that the numbers here are the same quantity as the published a6_cis_share table rather than an
## amplitude share that happens to share its name.
## This reads the saved axis table and computes nothing new: cis_tru and trn_tru per stratum, hybrid and
## timepoint are already there. The share is their ratio.

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

A <- fread(file.path(hsd, "a6_axis_amplitudes_37_30.tsv"))
stopifnot(uniqueN(A$pair) == 6L, setequal(unique(A$time), c(15L, 30L, 60L)),
          all(c("cis_tru", "trn_tru") %in% names(A)))
say("axis table: %d rows, %d hybrids, %d strata, %d timepoints",
    nrow(A), uniqueN(A$pair), uniqueN(A$stratum), uniqueN(A$time))

## The table has no whole-transcriptome row, so the all-gene quantity per hybrid is recombined from the strata
## weighted by their n: rms_true is a root-mean-square, and squares add in proportion to n.
POOL <- A[, .(n = sum(n),
              cis_tru = sqrt(sum(n * cis_tru^2) / sum(n)),
              trn_tru = sqrt(sum(n * trn_tru^2) / sum(n))), by = .(pair, time)]
POOL[, cis_share := cis_tru^2 / (cis_tru^2 + trn_tru^2)]
say("\n=== cis share per hybrid and timepoint, strata pooled by gene count ===")
W <- dcast(POOL, pair ~ time, value.var = "cis_share")
setnames(W, c("15", "30", "60"), c("t15", "t30", "t60"))
W[, falls_monotonically := t15 > t30 & t30 > t60]
W[, falls_15_to_60 := t15 > t60]
print(W[, .(pair, t15 = round(t15, 3), t30 = round(t30, 3), t60 = round(t60, 3),
            falls_monotonically, falls_15_to_60)], row.names = FALSE)

say("\n=== trans against cis at 60 min, per hybrid ===")
S60 <- POOL[time == 60L, .(pair, cis_tru = round(cis_tru, 3), trn_tru = round(trn_tru, 3),
                           trans_exceeds_cis = trn_tru > cis_tru)]
print(S60, row.names = FALSE)

say("\n=== the pre-registered read-out ===")
say("  (1) cis share falls from 15 to 60 min in %d of 6 hybrids (monotone in all three points: %d of 6)",
    W[falls_15_to_60 == TRUE, .N], W[falls_monotonically == TRUE, .N])
say("      criterion was at least 5 of 6: %s", W[falls_15_to_60 == TRUE, .N] >= 5L)
say("  (2) trans exceeds cis at 60 min in %d of 6 hybrids", S60[trans_exceeds_cis == TRUE, .N])
say("      criterion was at least 5 of 6: %s", S60[trans_exceeds_cis == TRUE, .N] >= 5L)

## The same question inside each stratum, since a trend that holds genome-wide but not in the Hsf1-dependent
## genes would matter for how the time-trend result and the Hsf1-stratum result are joined.
say("\n=== the same two counts within each stratum ===")
ST <- A[, .(cis_share = cis_tru^2 / (cis_tru^2 + trn_tru^2), cis_tru, trn_tru), by = .(stratum, pair, time)]
SUM <- ST[, {
  w <- dcast(.SD, pair ~ time, value.var = "cis_share")
  setnames(w, c("15", "30", "60"), c("t15", "t30", "t60"))
  s60 <- .SD[time == 60L]
  .(n_falls = w[t15 > t60, .N], n_hyb = nrow(w),
    n_trans_gt_cis_60 = s60[trn_tru > cis_tru, .N])
}, by = stratum]
print(SUM, row.names = FALSE)

fwrite(W, file.path(hsd, "p3_r3_per_hybrid_share.tsv"), sep = "\t")
fwrite(POOL, file.path(hsd, "p3_r3_per_hybrid_axes.tsv"), sep = "\t")
fwrite(SUM, file.path(hsd, "p3_r3_per_hybrid_by_stratum.tsv"), sep = "\t")
say("")
