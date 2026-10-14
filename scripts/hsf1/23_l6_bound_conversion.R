# Bound conversion -- convert the detectable coupling of the HSE-architecture test into a
# bound. This is a unit conversion, not a new analysis: the estimator, the substrates and the
# calibrated thresholds are those already on disk; only the units change.
#
# Purpose. The raw resolution range c ~ 0.02 to 0.80 is design-specific, so it is
# not directly transferable to another study. The same
# information stated as a bound is usable by the next study:
#   "with the power available here, the contribution of HSE architecture divergence to
#    allele-specific expression divergence is constrained below X."
#
# Conversion. The calibration injected a coupling of strength c as
#       cis  ->  cis + c * 1{turnover} * |amp|
# so c is a coefficient on the amplitude-scaled turnover indicator and has no external meaning.
# Two conversions give it one, and both are read off the same estimator the arms use:
#
#  (X1) Observable ratio scale, reported as primary because it is what other studies publish:
#       cis_ratio(c) = rms_tru(cis | turnover) / rms_tru(cis | no turnover) under injection c.
#       X1 = cis_ratio(c*) - 1 at the conservative detectable coupling c* (the 90th percentile of the
#       per-substrate minimum, i.e. 0.80 / 0.80 / 0.40 at 15 / 30 / 60 min). Read as: "genes whose
#       HSE architecture differs between species cannot have more than X1 higher allele-specific
#       divergence than genes whose architecture is shared, or we would have detected it."
#
#  (X2) Variance-share scale, reported as secondary because it is the quantity a modeller wants:
#       the injection adds variance c^2 * mean(amp^2) among the affected genes, so
#       X2 = c*^2 * mean(amp^2 | turnover) / rms_tru(cis | turnover)^2, capped at 1.
#       Read as: "HSE architecture divergence can account for at most X2 of the allele-specific
#       divergence variance of the genes it affects."
#
# Scope. The conversion does not create evidence below the detection limit. The bound is an
# upper limit determined by power rather than an estimate of the contribution.
# It constrains contributions below X but does not estimate the contribution or establish that HSE architecture does not
# explain divergence.
#
# Usage: Rscript scripts/hsf1/23_l6_bound_conversion.R

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
hsd <- file.path(root, "data/02_processed/hsf1"); csd <- file.path(root, "data/02_processed/crossspecies")
repd <- file.path(root, "doc/05_reports"); MIRROR <- "logs/output_copy"
SRC10 <- file.path(root, "scripts/hsf1/10_promoter_hse_architecture_v2.R")
TP <- c(15L, 30L, 60L); PAIRS <- c("ScxSp","ScxSk","ScxSu","SpxSk","SpxSu","SkxSu")
C_GRID <- c(0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.60, 0.80)
OBS_DISC_REF <- c(0.167376135861948, 0.00257636505985581, 0.0641650146394714)
set.seed(20260906L)
say <- function(...) cat(sprintf(...), "\n", sep = "")
logf <- file.path(repd, "logs/hsf1_l6_bound.log"); sink(logf, split = TRUE)
say("=== 23_l6_bound_conversion | %s ===", format(Sys.time()))
say("A CONVERSION of the calibrated coupling into a field-usable bound. No new analysis.")

## ---------------------------- reproduce and verify the estimator ---
ex <- parse(SRC10); want <- c("perm_G","arm_cells","med_by_t","inject_on")
for (e in ex) if (is.call(e) && length(e) >= 3L && as.character(e[[1]]) %in% c("<-","=") &&
                  is.name(e[[2]]) && as.character(e[[2]]) %in% want) eval(e, envir = globalenv())
stopifnot(all(want %in% ls(globalenv())))
HB <- fread(file.path(csd, "hog_to_scer.tsv"))
A6 <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
PD <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "b3v2_site_turnover.tsv.gz"))))
sysmap <- unique(HB[gene %in% unique(A6$hog), .(hog = gene, sys)])
D <- merge(unique(merge(A6[, .(hog, pair, time, cis, cis_nul, trans, trans_nul, amp)], sysmap, by = "hog")),
           PD, by = c("sys","pair"))
D <- D[!is.na(cis) & !is.na(cis_nul) & !is.na(trans) & !is.na(trans_nul)][, cell := paste(pair, time, sep = "_")]
mk <- function(vn) { w <- dcast(D, sys ~ cell, value.var = vn, fun.aggregate = function(z) z[1])
  m <- as.matrix(w[, setdiff(names(w), "sys"), with = FALSE]); rownames(m) <- w$sys; m }
CIS <- mk("cis"); CISN <- mk("cis_nul"); TRN <- mk("trans"); TRNN <- mk("trans_nul"); AMPA <- abs(mk("amp"))
CELLS <- colnames(CIS); GEN_ <- rownames(CIS)
pair_of <- sub("_[0-9]+$", "", CELLS); tp_of <- as.integer(sub(".*_", "", CELLS))
PAIRS <- unique(pair_of)
grp_mat <- function(vn) { w <- dcast(unique(D[, c("sys","pair",vn), with = FALSE]), sys ~ pair,
  value.var = vn, fun.aggregate = function(z) z[1])
  m <- as.matrix(w[, PAIRS, with = FALSE]); rownames(m) <- w$sys; m[GEN_, , drop = FALSE] }
TO <- grp_mat("turnover"); TO[is.na(TO)] <- 0L
obs <- med_by_t(arm_cells(CIS, CISN, TO)) - med_by_t(arm_cells(TRN, TRNN, TO))
say("\nrebuild faithfulness vs script 10: max abs diff %.2e", max(abs(obs - OBS_DISC_REF)))
stopifnot(max(abs(obs - OBS_DISC_REF)) < 1e-9)
say("ESTIMATOR REPRODUCED")

## --------------------------- cis_ratio as a function of the injected coupling --
RES <- rbindlist(lapply(C_GRID, function(cc) {
  Xi <- inject_on(cc, TO)
  data.table(c_inj = cc, time = TP, cis_ratio = med_by_t(arm_cells(Xi, CISN, TO)))
}))
fwrite(RES, file.path(hsd, "b3v2_bound_cis_ratio_vs_c.tsv"), sep = "\t")
say("\ncis_ratio = rms_tru(cis | turnover) / rms_tru(cis | no turnover), under injection c")
say("%-8s %12s %12s %12s", "c", "15 min", "30 min", "60 min")
for (cc in C_GRID) { z <- RES[c_inj == cc]; setkey(z, time)
  say("%-8.2f %12.4f %12.4f %12.4f", cc, z[.(15L), cis_ratio], z[.(30L), cis_ratio], z[.(60L), cis_ratio]) }

## ------------------------------------------ the two bounds at the calibrated c* -
RS <- fread(file.path(hsd, "b3v2_power_resolution.tsv"))
setkey(RS, time)
# variance share needs the per-timepoint amplitude of the affected genes and the cis rms there
amp2_to <- vapply(TP, function(t) {
  cs <- which(tp_of == t)
  v <- unlist(lapply(cs, function(c_) { j <- match(pair_of[c_], PAIRS)
    a <- AMPA[, c_]; a[TO[, j] == 1L] }))
  mean(v^2, na.rm = TRUE) }, 0)
cis_to <- vapply(TP, function(t) {
  cs <- which(tp_of == t)
  v <- unlist(lapply(cs, function(c_) { j <- match(pair_of[c_], PAIRS)
    CIS[TO[, j] == 1L, c_] }))
  n <- unlist(lapply(cs, function(c_) { j <- match(pair_of[c_], PAIRS)
    CISN[TO[, j] == 1L, c_] }))
  sqrt(max(0, mean(v^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE))) }, 0)
BND <- rbindlist(lapply(seq_along(TP), function(i) {
  t <- TP[i]; cst <- RS[.(t), p90]
  # cis_ratio at c*, interpolated on the computed grid
  z <- RES[time == t][order(c_inj)]
  r_at <- approx(z$c_inj, z$cis_ratio, xout = cst, rule = 2)$y
  x2 <- min(1, cst^2 * amp2_to[i] / max(cis_to[i]^2, 1e-12))
  data.table(time = t, c_star = cst, cis_ratio_at_cstar = r_at,
             X1_max_excess_ratio = r_at - 1,
             mean_amp2_turnover = amp2_to[i], rms_tru_cis_turnover = cis_to[i],
             X2_max_variance_share = x2)
}))
fwrite(BND, file.path(hsd, "b3v2_l6_bound.tsv"), sep = "\t")
say("\n================ THE BOUND, IN TWO USABLE UNITS ================")
say("%-6s %8s %16s %20s %24s", "time", "c*", "cis_ratio at c*", "X1 = max excess", "X2 = max variance share")
for (i in seq_len(nrow(BND))) say("%-6d %8.2f %16.4f %19.1f%% %23.1f%%",
  BND$time[i], BND$c_star[i], BND$cis_ratio_at_cstar[i], 100*BND$X1_max_excess_ratio[i],
  100*BND$X2_max_variance_share[i])
say("")
say("c* is the 90th percentile of the per-substrate minimum detectable coupling across 52")
say("randomised substrates, using the permutation criterion at the inference's own alpha = 0.05/3.")
say("")
say("What the bound supports:")
say("  'With the power available in this design, divergence of classic or gapped heat-shock-element")
say("   architecture cannot raise the allele-specific expression divergence of the genes it affects")
say("   by more than %.0f%% at 60 min (%.0f%% at 15 and 30 min), and cannot account for more than",
    100*BND[time==60L, X1_max_excess_ratio], 100*BND[time==15L, X1_max_excess_ratio])
say("   %.0f%% of their allele-specific divergence variance. Larger contributions are excluded;",
    100*BND[time==60L, X2_max_variance_share])
say("   smaller ones are not addressed by these data.'")
say("")
say("This is an upper limit determined by power rather than an estimate.")
say("It does not estimate the contribution of HSE architecture or exclude smaller contributions.")
say(". The observed 15 min discriminant is NOT evidence for a contribution: implied")
say("coupling ~0.278 is 35%% of c*, its gene-level bootstrap CI contains 0, and it sits in the")
say("three strain-mismatched pairs.")
writeLines(capture.output(sessionInfo()), file.path(hsd, "sessionInfo_l6_bound.txt"))
for (f in list.files(hsd, pattern = "^(b3v2_bound|b3v2_l6_bound|sessionInfo_l6)", full.names = TRUE))
  file.copy(f, file.path(MIRROR, "data/02_processed/hsf1", basename(f)), overwrite = TRUE)
file.copy(file.path(root, "scripts/hsf1/23_l6_bound_conversion.R"), file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("
copied. === done ===")
sink(); file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
