# 13_resolution_shard.R -- one shard of the resolution calibration. The equivalence bound is a
# conservative quantile of the per-substrate minimum detectable coupling, not its median.
#
# Across randomized substrates, the median is 15 min 0.30 / 30 min 0.35 / 60 min 0.30, while
# the per-substrate values range 0.10-0.80. An equivalence bound means "an effect below this value
# cannot be detected by this design". A median says only that half the noise realisations could
# detect it, so it overstates sensitivity and understates the resolution floor, the one direction
# an equivalence claim must never err in. The bound is therefore reported as min / median / 90th
# percentile / max, with the 90th percentile as the bound and the median as description.
#
# A 90th percentile estimated from only 10 values is essentially the 9th order statistic and
# is very coarse. This script is shardable: it takes a seed offset and a substrate count from the
# environment, so several shards run in parallel on separate cores and are pooled afterwards
# in a pooling analysis.
#
# The promoters are not rebuilt. The matrices are reconstructed from what script 10 wrote to disk
# (site_turnover.tsv.gz + a6_cis_trans_37_30.tsv.gz), which skips the ~15 min of promoter
# extraction and site pairing per shard. That shortcut is safe only if the rebuild is identical,
# so the script asserts that the reconstructed matrices reproduce script 10's observed
# discriminant to 1e-9 before it calibrates anything, and stops otherwise.
#
# The estimator, the substrate generator and the signal-free assertion are parsed out of
# script 10 rather than copied, as in scripts 11 and 12, so they cannot drift.
#
# Usage: SHARD=1 N_SUB_SHARD=10 Rscript scripts/hsf1/13_resolution_shard.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
repd <- file.path(root, "doc/05_reports")
SRC10 <- file.path(root, "scripts/hsf1/10_promoter_hse_architecture_v2.R")

SHARD       <- as.integer(Sys.getenv("SHARD", "1"))
N_SUB_SHARD <- as.integer(Sys.getenv("N_SUB_SHARD", "10"))
NPERM_SWEEP <- 300L
ALPHA_INFER <- 0.05 / 3
TP          <- c(15L, 30L, 60L)
C_GRID      <- c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.80)   # identical to script 10
HYB <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")
# script 10's observed discriminant, hard-coded as the faithfulness target
OBS_DISC_REF <- c(0.167376135861948, 0.00257636505985581, 0.0641650146394714)

say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, sprintf("logs/hsf1_resolution_shard%d.log", SHARD))
sink(logf, split = TRUE)
say("=== 13_resolution_shard | shard %d | %d substrates | %s ===",
    SHARD, N_SUB_SHARD, format(Sys.time()))
say("the equivalence bound is a CONSERVATIVE quantile (90th), not the median.")

## ------------------------- transplant the estimator out of script 10 --------
ex <- parse(SRC10)
want <- c("perm_G", "arm_cells", "med_by_t", "assert_signal_free", "perm_p", "inject_on")
got <- character(0)
for (e in ex) {
  if (is.call(e) && length(e) >= 3L && as.character(e[[1]]) %in% c("<-", "=") &&
      is.name(e[[2]]) && as.character(e[[2]]) %in% want) {
    eval(e, envir = globalenv()); got <- c(got, as.character(e[[2]]))
  }
}
say("transplanted from script 10: %s", paste(got, collapse = ", "))
missing <- setdiff(want, got)
if (length(missing)) stop(sprintf("script 10 no longer defines: %s", paste(missing, collapse = ", ")))

## ------------------------------- rebuild the matrices from what 10 wrote ----
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
A6 <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
PD <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "b3v2_site_turnover.tsv.gz"))))
sysmap <- unique(HB[gene %in% unique(A6$hog), .(hog = gene, sys)])
A6s <- unique(merge(A6[, .(hog, pair, time, cis, cis_nul, trans, trans_nul, amp)], sysmap, by = "hog"))
D <- merge(A6s, PD, by = c("sys", "pair"))
D <- D[!is.na(cis) & !is.na(cis_nul) & !is.na(trans) & !is.na(trans_nul)]
D[, cell := paste(pair, time, sep = "_")]
mk <- function(vn) {
  w <- dcast(D, sys ~ cell, value.var = vn, fun.aggregate = function(z) z[1])
  m <- as.matrix(w[, setdiff(names(w), "sys"), with = FALSE]); rownames(m) <- w$sys; m
}
CIS <- mk("cis"); CISN <- mk("cis_nul"); TRN <- mk("trans"); TRNN <- mk("trans_nul")
AMPA <- abs(mk("amp"))
CELLS <- colnames(CIS); GEN_ <- rownames(CIS)
pair_of <- sub("_[0-9]+$", "", CELLS); tp_of <- as.integer(sub(".*_", "", CELLS))
PAIRS <- unique(pair_of)
grp_mat <- function(vn) {
  w <- dcast(unique(D[, c("sys", "pair", vn), with = FALSE]), sys ~ pair, value.var = vn,
             fun.aggregate = function(z) z[1])
  m <- as.matrix(w[, PAIRS, with = FALSE]); rownames(m) <- w$sys; m[GEN_, , drop = FALSE]
}
TO <- grp_mat("turnover"); TO[is.na(TO)] <- 0L
say("rebuilt: %d genes x %d cells", length(GEN_), length(CELLS))

## ---- check: the rebuild must reproduce script 10's statistic ----
obs_disc <- med_by_t(arm_cells(CIS, CISN, TO)) - med_by_t(arm_cells(TRN, TRNN, TO))
say("observed discriminant rebuilt: %+.12f / %+.12f / %+.12f", obs_disc[1], obs_disc[2], obs_disc[3])
say("script 10 reference          : %+.12f / %+.12f / %+.12f",
    OBS_DISC_REF[1], OBS_DISC_REF[2], OBS_DISC_REF[3])
stopifnot(all(abs(obs_disc - OBS_DISC_REF) < 1e-9))
say("REBUILD FAITHFUL (max abs diff %.2e) -- calibration may proceed",
    max(abs(obs_disc - OBS_DISC_REF)))

## ------------------------------------------------- the calibration itself ---
set.seed(90610L + 1000L * SHARD)
OUT <- list()
for (sb in seq_len(N_SUB_SHARD)) {
  sid <- sprintf("s%d_%d", SHARD, sb)
  TO_S <- perm_G(TO)
  assert_signal_free(TO_S, TO, sid)                      # the substrate must differ from the observed labels
  stopifnot(all(colSums(TO_S == 1L) == colSums(TO == 1L)))
  for (cc in C_GRID) {
    Xi <- inject_on(cc, TO_S)
    od <- med_by_t(arm_cells(Xi, CISN, TO_S)) - med_by_t(arm_cells(TRN, TRNN, TO_S))
    nl <- matrix(NA_real_, NPERM_SWEEP, length(TP))
    for (k in seq_len(NPERM_SWEEP)) {
      Gp <- perm_G(TO_S)
      nl[k, ] <- med_by_t(arm_cells(Xi, CISN, Gp)) - med_by_t(arm_cells(TRN, TRNN, Gp))
    }
    for (j in seq_along(TP)) {
      p <- perm_p(nl[, j], od[j])
      OUT[[length(OUT) + 1L]] <- data.table(shard = SHARD, substrate = sid, c_inj = cc,
        time = TP[j], delta = od[j], p_emp = p, detected = p < ALPHA_INFER)
    }
  }
  say("  substrate %s done (%d/%d)", sid, sb, N_SUB_SHARD)
}
OUT <- rbindlist(OUT)
of <- file.path(hsd, sprintf("b3v2_resolution_shard%d.tsv", SHARD))
fwrite(OUT, of, sep = "\t")
say("\nwrote %s (%d rows)", basename(of), nrow(OUT))
say("c=0 detections in this shard (must be ~alpha): %d of %d",
    OUT[c_inj == 0 & detected, .N], OUT[c_inj == 0, .N])
writeLines(capture.output(sessionInfo()),
           file.path(hsd, sprintf("sessionInfo_resolution_shard%d.txt", SHARD)))
say("=== done ===")
sink()
