#!/usr/bin/env Rscript
## 48_flag_degenerate_cells.R -- add a degeneracy flag to the saved module-ratio tables.
##
## Background. A ratio of noise-corrected quantities is undefined when its denominator is at or near zero, and the cold
## arm at 15 minutes is exactly that case: Fay et al. report a mean of eleven cold-responsive genes with an
## allele-by-time interaction, so the noise floor reaches the signal and rms_true collapses. The tables record
## what the arithmetic returned, which is correct, but what it returned includes cells a reader can mistake for
## a result:
##   translation_ribosome, cold, 15 min: ratio 0 with interval [0, 0] at four exponents, readable as
##     "resolvably below 1" when it is in fact an empty denominator;
##   cell_cycle_growth, cold, 15 min, beta 0.75: no point estimate but a lower bound of 1.31, readable as
##     "resolvably above 1";
##   proteostasis, cold, 15 min, beta 1: interval [1.212, 1.212], zero width.
## Note: a bootstrap read-out that tests only that the lower bound is finite and above 1 marks such cells
## resolvable. The flag marks them explicitly for downstream scripts and for readers of the supplementary tables.
##
## Scope. It adds one boolean column. It does not recompute, round or alter any ratio
## or interval bound, and it asserts that afterwards.
##
## Criterion: a cell is degenerate if its point estimate, lower bound or upper bound is not finite, or if the
## interval has zero width. Cells that are merely wide are not degenerate: the cold 15-minute oxidative cells
## at beta 0 and 0.25, 3.413 [1.64, 12.81] and 5.142 [1.78, 15.74], are real estimates and are reported.

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

TARGETS <- c("s4_module_divergence_beta.tsv", "s4_module_divergence_beta_cold.tsv",
             "s4_module_stratum_split_beta.tsv", "s4_module_stratum_split3_beta.tsv",
             "s4_cold_regulon_split_beta.tsv", "s4_stem_collapse_ratios.tsv",
             "s4_metabolism_amplitude_matched.tsv")
RATIO_COL <- c("ratio", "ratio_full", "cal_ratio")

total <- 0L
for (f in TARGETS) {
  p <- file.path(hsd, f)
  if (!file.exists(p)) { say("  (absent) %s", f); next }
  X <- fread(p)
  rc <- intersect(RATIO_COL, names(X))
  if (!length(rc) || !all(c("lo", "hi") %in% names(X))) { say("  (no ratio/lo/hi) %s", f); next }
  r <- rc[1]
  before <- copy(X[, c(r, "lo", "hi"), with = FALSE])
  X[, degenerate := !is.finite(get(r)) | !is.finite(lo) | !is.finite(hi) | (hi - lo) < 1e-9]
  after <- X[, c(r, "lo", "hi"), with = FALSE]
  stopifnot(identical(before, after))            # the flag must not touch a number
  n <- sum(X$degenerate); total <- total + n
  fwrite(X, p, sep = "\t")
  say("  %-40s %4d cells, %3d flagged degenerate", f, nrow(X), n)
}
say("\n  total flagged: %d", total)
say("  Ratios and interval bounds were asserted unchanged in every table written.")

## Read the flagged tables back and report what a reader would now see, per timepoint, for the cold arm.
C <- fread(file.path(hsd, "s4_module_divergence_beta_cold.tsv"))
say("\n  cold arm, degenerate cells by timepoint:")
print(C[, .(cells = .N, degenerate = sum(degenerate)), by = time][order(time)], row.names = FALSE)
say("\n  cold 15 min cells that are wide but NOT degenerate, and therefore reportable:")
print(C[time == 15L & !degenerate, .(module, beta, ratio = round(ratio, 3),
                                     lo = round(lo, 2), hi = round(hi, 2))], row.names = FALSE)
