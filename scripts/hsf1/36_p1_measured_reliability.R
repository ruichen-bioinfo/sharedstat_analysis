#!/usr/bin/env Rscript
## 36_p1_measured_reliability.R -- measure h in the real data, so the bias formula stops being hypothetical.
##
## Purpose. The paper derives gamma_hat = gamma*h and shows the direction inverts once h < 1/gamma. Both
## statements are about a quantity, h, that has to be measured in the data the paper analyses. Without that
## measurement the reader is told a hazard exists but not whether their own design is in the regime where
## it bites. With it, the threshold becomes checkable.
##
## How h is identified, with no extra assumptions. The design has two replicate blocks, so for a
## quantity measured as b1 = R + u1 and b2 = R + u2 with independent errors of equal variance,
##      sig = (b1 + b2)/2 = R + (u1 + u2)/2      var(sig) = var(R) + sigma^2/2
##      nul = (b2 - b1)/2 =      (u2 - u1)/2      var(nul) =            sigma^2/2
## and therefore
##      h = 1 - mean(nul^2)/mean(sig^2) = 1 - 1/snr^2,   snr = rms(sig)/rms(nul)  (uncentred mean squares).
## Nothing here is fitted. The null contrast is the same one the amplitude estimator already subtracts, so
## the reliability is read off quantities the pipeline computes anyway.
##
## What this h is and is not. It is the reliability of the allele-paired difference-of-differences, the
## contrast this design actually forms. It is not a reliability of "the data" in general, and it is not
## transferable to another design (P1 reports the centred reference-arm h of scripts 91/93). Reported per stratum, hybrid and timepoint, because a single number
## would hide the spread, and the spread is the point: an artefact whose size depends on h is not constant
## across a design whose h is not constant.

suppressPackageStartupMessages(library(data.table))
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")

## ---- route 1: from the saved rms table (h = 1 - 1/snr^2) --------------------------------------------
FL <- fread(file.path(hsd, "a5_cis_noise_corrected_37_30.tsv"))
FL[, h_from_snr := 1 - 1 / snr^2]

## ---- route 2: recomputed from the per-gene contrasts, independently of the saved rms values ----------
D <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a5_cis_dod_37_30.tsv.gz"))))
H <- D[, {
  keep <- is.finite(sig) & is.finite(nul)          # paired, as in rms_true
  vs <- mean(sig[keep]^2); vn <- mean(nul[keep]^2)
  .(n = sum(keep), h_recomputed = 1 - vn / vs)
}, by = .(stratum, hybrid, time)]

M <- merge(FL[, .(stratum, hybrid, time, snr, rms_sig, rms_nul, rms_tru, h_from_snr)],
           H, by = c("stratum", "hybrid", "time"))
say("=== validation check: the two routes to h must agree ===")
say("max |h_from_snr - h_recomputed| over %d cells = %.3e", nrow(M),
    max(abs(M$h_from_snr - M$h_recomputed)))
stopifnot(max(abs(M$h_from_snr - M$h_recomputed)) < 1e-6)
say("agreement is exact to numerical precision, so h is a property of the measured contrasts and not of")
say("either route's arithmetic.\n")

setnames(M, "h_recomputed", "h")

## A reliability is a variance ratio and cannot be negative. h = 1 - 1/snr^2 goes negative exactly when the
## null contrast exceeds the signal contrast, i.e. when the estimated true variance is at or below zero --
## the same condition under which the amplitude estimator clamps at zero. Those cells are not reliability
## estimates and are reported separately rather than folded into a median, for the same reason a clamped
## amplitude must be reported as "not distinguishable from noise" rather than as a small positive number.
M[, estimable := h > 0]
say("=== cells where h is not estimable (null contrast exceeds signal, so var(R) <= 0) ===")
say("%d of %d cells (%.1f%%)%s", sum(!M$estimable), nrow(M), 100 * mean(!M$estimable),
    if (sum(!M$estimable)) sprintf(", by stratum: %s", paste(sprintf("%s=%d",
      M[estimable == FALSE, .N, by = stratum]$stratum, M[estimable == FALSE, .N, by = stratum]$N),
      collapse = ", ")) else "")
say("These are excluded from every summary below. Excluding them is conservative in the direction that")
say("matters: they are the LEAST reliable cells, so dropping them raises the reported median h and")
say("understates how often the sign can invert.\n")
E <- M[estimable == TRUE]

say("=== measured reliability of the allele-paired contrast, by stratum (estimable cells only) ===")
say("%-16s %5s %8s %8s %8s %8s %10s", "stratum", "cells", "min h", "median", "max h", "median snr", "h<1/1.3?")
for (s in sort(unique(E$stratum))) {
  x <- E[stratum == s]
  say("%-16s %5d %8.4f %8.4f %8.4f %8.2f %10s", s, nrow(x), min(x$h), median(x$h), max(x$h),
      median(x$snr), sprintf("%d/%d", sum(x$h < 1/1.3), nrow(x)))
}

say("\n=== where these data sit relative to the sign-inversion threshold h < 1/gamma ===")
say("%-8s %10s %14s %14s %16s", "gamma", "1/gamma", "cells below", "of total", "median gamma*h")
for (g in c(1.1, 1.2, 1.3, 1.5, 2.0)) {
  below <- sum(E$h < 1/g)
  say("%-8.1f %10.4f %14d %14d %16.4f", g, 1/g, below, nrow(E), g * median(E$h))
}
say("\nReading. A cell with h below 1/gamma is one where a true gain of gamma would be reported as a loss.")
say("The fraction of cells in that state is a property of this design, and it is not small.")

say("\n=== the implied naive bias, cell by cell, at three illustrative gains ===")
for (g in c(1.1, 1.3, 1.5)) {
  b <- g * E$h - g
  say("gamma = %.1f : bias median %+.4f, range %+.4f to %+.4f, sign inverted in %d of %d cells",
      g, median(b), min(b), max(b), sum(g * E$h < 1), nrow(E))
}

fwrite(M[order(stratum, hybrid, time)], file.path(hsd, "p1_measured_reliability.tsv"), sep = "\t")
say("\nwrote %s (%d rows)", file.path(hsd, "p1_measured_reliability.tsv"), nrow(M))
