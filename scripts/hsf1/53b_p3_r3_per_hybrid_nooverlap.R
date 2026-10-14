#!/usr/bin/env Rscript
## 53b_p3_r3_per_hybrid_nooverlap.R -- the genome-wide per-hybrid pooling of 53_p3_r3_per_hybrid.R without the
## overlapping stratum.
##
## Pooling every stratum row of a6_axis_amplitudes_37_30.tsv by gene count. One of those strata, other_induced
## (761 genes), is the union of TF_independent (536) and Msn24_only (225), so those 761 genes enter the original pool twice
## and the pooled n is 5184 instead of 4423. This script repeats the pooling with other_induced left out, reproduces the original script's
## saved values first (fidelity check), and writes the corrected values to new tables beside the original tables, which are not
## overwritten. It also writes a comparison for every pooled value.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
A <- fread(file.path(hsd, "a6_axis_amplitudes_37_30.tsv"))
pool <- function(d) {
  P <- d[, .(n = sum(n), cis_tru = sqrt(sum(n * cis_tru^2) / sum(n)), trn_tru = sqrt(sum(n * trn_tru^2) / sum(n))),
         by = .(pair, time)]
  P[, cis_share := cis_tru^2 / (cis_tru^2 + trn_tru^2)][]
}
OLD <- pool(A); NEW <- pool(A[stratum != "other_induced"])
L <- fread(file.path(hsd, "p3_r3_per_hybrid_axes.tsv"))
chk <- merge(L, OLD, by = c("pair", "time"), suffixes = c("_landed", "_re"))
stopifnot(nrow(chk) == 18L, max(abs(chk$cis_share_landed - chk$cis_share_re)) < 1e-12)
say("fidelity: the 18 original pooled rows reproduced exactly")
stopifnot(all(NEW$n == 4423L))
CMP <- merge(OLD, NEW, by = c("pair", "time"), suffixes = c("_old", "_new"))
fwrite(NEW, file.path(hsd, "p3_r3_per_hybrid_axes_nooverlap.tsv"), sep = "\t")
fwrite(CMP, file.path(hsd, "p3_r3_per_hybrid_nooverlap_comparison.tsv"), sep = "\t")
W <- dcast(NEW, pair ~ time, value.var = "cis_share"); setnames(W, c("15", "30", "60"), c("t15", "t30", "t60"))
W[, falls_monotonically := t15 > t30 & t30 > t60][, falls_15_to_60 := t15 > t60]
fwrite(W, file.path(hsd, "p3_r3_per_hybrid_share_nooverlap.tsv"), sep = "\t")
s60 <- NEW[time == 60L]
say("without the overlap: cis share falls 15 -> 60 min in %d of 6, monotonically in %d; trans > cis at 60 min in %d of 6",
    W[falls_15_to_60 == TRUE, .N], W[falls_monotonically == TRUE, .N], s60[trn_tru > cis_tru, .N])
say("largest change in a pooled cis share: %.3f", max(abs(CMP$cis_share_new - CMP$cis_share_old)))
