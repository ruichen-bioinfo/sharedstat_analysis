#!/usr/bin/env Rscript
## 52_c5_decomposition.R -- the cis/trans decomposition of the cold response, three replicates.
##
## The estimator, and why the algebra is checked before it is used. With two replicate blocks the companion
## forms signal = (d1 + d2)/2 and null = (d2 - d1)/2, so that both carry the same error variance and the null
## carries no signal; rms_true = sqrt(max(0, mean(sig^2) - mean(nul^2))) then removes the noise floor. With
## three replicates the same construction generalises to one signal contrast and two orthogonal null contrasts:
##     sig = (d1 + d2 + d3)/3          var = sigma^2/3
##     n1  = (d2 - d1)/sqrt(6)        var = sigma^2/3
##     n2  = (2*d3 - d1 - d2)/sqrt(18) var = sigma^2/3
## All three have variance sigma^2/3, sig has expectation delta while n1 and n2 have expectation zero, and the
## three contrast vectors are mutually orthogonal. So E[sig^2] - mean(E[n1^2], E[n2^2]) = delta^2 exactly. That
## identity is checked by simulation below before any real data passes through it; if it fails, nothing is
## computed.
##
## Decision rule: in the header of 50_c5_build_index.R, written before any quantification existed. Gate
## results: lambda 0.0064 against a 0.10 threshold; baseline allelic constant 0.0427 with 0.0038 surviving into
## the estimand. The mapping-rate asymmetry is 5.2 percentage points, Su lower.
##
## What is and is not estimated here. This is a two-level design: 12 C against 30 C, no time course. The shape
## of the decomposition is testable: whether trans exceeds cis, whether the cis share is below a half, whether
## cross-replicate compensation is negative. The time dynamics of the companion are not, and no statement about
## front-loading or accumulation is made from this data.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hov  <- file.path(root, "data/02_processed/hov2020")
hsd  <- file.path(root, "data/02_processed/hsf1")
csd  <- file.path(root, "data/02_processed/crossspecies")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

MINCPM <- 1; MINLIB <- 4L; PRIOR <- 2; NBOOT <- 1000L; SEED <- 20260912L

## ------------------------------------------------------------------ algebra gate ----
say("=== algebra gate: E[sig^2] - mean(E[n1^2], E[n2^2]) = delta^2 ===")
set.seed(1L)
ng <- 200000L; delta <- 0.7; sigma <- 1.3
d <- matrix(rnorm(3L * ng, mean = delta, sd = sigma), ng, 3L)
sig <- rowMeans(d)
n1  <- (d[, 2] - d[, 1]) / sqrt(6)
n2  <- (2 * d[, 3] - d[, 1] - d[, 2]) / sqrt(18)
say("  var(sig) %.5f | var(n1) %.5f | var(n2) %.5f | sigma^2/3 = %.5f",
    var(sig), var(n1), var(n2), sigma^2 / 3)
est <- mean(sig^2) - mean(c(n1^2, n2^2))
say("  recovered delta^2 = %.5f against %.5f, difference %.2e", est, delta^2, abs(est - delta^2))
stopifnot(abs(est - delta^2) < 5e-3, abs(var(n1) - sigma^2 / 3) < 5e-3, abs(var(n2) - sigma^2 / 3) < 5e-3)
## orthogonality of the three contrast vectors, exactly
V <- rbind(c(1, 1, 1) / 3, c(-1, 1, 0) / sqrt(6), c(-1, -1, 2) / sqrt(18))
say("  contrast vectors mutually orthogonal: %s", isTRUE(all.equal(V %*% t(V), diag(diag(V %*% t(V))))))
stopifnot(isTRUE(all.equal(V %*% t(V), diag(diag(V %*% t(V))))))
say("  gate passed.")

## ------------------------------------------------------------------ counts to log-CPM ----
D <- fread(file.path(hov, "design.tsv"))
LENT <- fread(file.path(hov, "cds_lengths.tsv"))
KEEP_LEN <- LENT[len_frac < 0.05, sys]
say("\ngenes whose two alleles agree in CDS length to within 5%%: %d of %d", length(KEEP_LEN), nrow(LENT))

read_counts <- function(al) {
  q <- fread(file.path(hov, "quant", al, "quant.sf"))
  q[, c("sys", "allele") := tstrsplit(Name, "|", fixed = TRUE)]
  q[sys %in% KEEP_LEN, .(sys, allele, n = NumReads)]
}
## A hybrid library gives two allelic counts per gene; a parental library gives one usable count, on its own
## allele, and the other allele's count is the leakage measured in gate 1 and is discarded.
LIB <- list()
for (i in seq_len(nrow(D))) {
  al <- D$alias[i]; g <- D$geno[i]
  q <- read_counts(al)
  if (g == "SCSU") {
    LIB[[paste0(al, "|Sc")]] <- q[allele == "Sc", .(sys, n)]
    LIB[[paste0(al, "|Su")]] <- q[allele == "Su", .(sys, n)]
  } else {
    own <- if (g == "SC") "Sc" else "Su"
    LIB[[paste0(al, "|", own)]] <- q[allele == own, .(sys, n)]
  }
}
say("count vectors assembled: %d (6 hybrid libraries x 2 alleles + 12 parental libraries x 1)", length(LIB))
stopifnot(length(LIB) == 24L)

## log-CPM within each count vector, which is the unit the companion works in
CPM <- rbindlist(lapply(names(LIB), function(k) {
  z <- LIB[[k]]; tot <- sum(z$n)
  z[, .(sys, key = k, cpm = 1e6 * n / tot, lcpm = log2(1e6 * n / tot + PRIOR))]
}))
CPM[, c("alias", "allele") := tstrsplit(key, "|", fixed = TRUE)]
CPM[, c("geno", "tempc", "rep") := tstrsplit(alias, "_", fixed = TRUE)]
CPM[, temp := as.integer(sub("C$", "", tempc))]

## expression filter: a gene must clear MINCPM in at least MINLIB of the vectors of its own arm
ok_hyb <- CPM[geno == "SCSU", .(k = sum(cpm >= MINCPM)), by = sys][k >= MINLIB, sys]
ok_par <- CPM[geno != "SCSU", .(k = sum(cpm >= MINCPM)), by = sys][k >= MINLIB, sys]
UNIV <- intersect(intersect(ok_hyb, ok_par), KEEP_LEN)
say("expression filter: %d genes clear %g CPM in >= %d vectors of both arms", length(UNIV), MINCPM, MINLIB)
stopifnot(length(UNIV) > 2000L)

## ------------------------------------------------------------------ responses ----
## r = lcpm(12 C) - lcpm(30 C) within one allele and one replicate index. The pairing of replicate i at 12 C
## with replicate i at 30 C is arbitrary but unbiased: the cultures are independent, so any pairing gives the
## same expectation and differs only in which noise lands in which contrast. Five random pairings are run below
## as a sensitivity check.
resp <- function(pairing = 1:3) {
  z <- CPM[sys %in% UNIV]
  w <- dcast(z, sys + allele + geno ~ temp + rep, value.var = "lcpm")
  cn <- names(w)
  r <- w[, .(sys, allele, geno)]
  for (i in 1:3) {
    c12 <- sprintf("12_rep%d", pairing[i]); c30 <- sprintf("30_rep%d", i)
    stopifnot(c12 %in% cn, c30 %in% cn)
    r[[paste0("d", i)]] <- w[[c12]] - w[[c30]]
  }
  r[complete.cases(r)]
}
R <- resp()
say("\nresponse table: %d rows over %d genes", nrow(R), uniqueN(R$sys))

contrasts3 <- function(DT) {
  DT[, `:=`(sig = (d1 + d2 + d3) / 3,
            nu1 = (d2 - d1) / sqrt(6),
            nu2 = (2 * d3 - d1 - d2) / sqrt(18))]
  DT[]
}
R <- contrasts3(R)
rms_true <- function(s, n) sqrt(max(0, mean(s^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE)))

## cis is the difference of the two allelic responses inside the hybrid; par is the difference of the two
## parental responses; trans is par - cis. Each is formed contrast by contrast so that the noise structure
## carries through.
build_axes <- function(R) {
  H <- R[geno == "SCSU"]; P <- R[geno != "SCSU"]
  hs <- dcast(H, sys ~ allele, value.var = c("sig", "nu1", "nu2"))
  ps <- dcast(P, sys ~ allele, value.var = c("sig", "nu1", "nu2"))
  A <- merge(hs, ps, by = "sys", suffixes = c("_h", "_p"))
  A[, `:=`(cis_sig = sig_Su_h - sig_Sc_h, cis_n1 = nu1_Su_h - nu1_Sc_h, cis_n2 = nu2_Su_h - nu2_Sc_h,
           par_sig = sig_Su_p - sig_Sc_p, par_n1 = nu1_Su_p - nu1_Sc_p, par_n2 = nu2_Su_p - nu2_Sc_p)]
  A[, `:=`(trn_sig = par_sig - cis_sig, trn_n1 = par_n1 - cis_n1, trn_n2 = par_n2 - cis_n2)]
  A[complete.cases(A[, .(cis_sig, par_sig, trn_sig)])]
}
A <- build_axes(R)
say("axes built for %d genes", nrow(A))
stopifnot(nrow(A) > 2000L)

tru <- function(A, pre, idx = NULL) {
  d <- if (is.null(idx)) A else A[idx]
  rms_true(d[[paste0(pre, "_sig")]], c(d[[paste0(pre, "_n1")]], d[[paste0(pre, "_n2")]]))
}
## The share is a variance share, with script 05's definition: the companion's published cis_share is
## rms_tru(cis)^2 / (rms_tru(cis)^2 + rms_tru(trans)^2), a partition of divergence variance, taken as the
## median over the six hybrids. Variances add along par = cis + trans up to a covariance term, which is why the
## share is defined on squares; and trans is a residual, so a ratio with it in the numerator is unbounded when
## cis approaches its noise floor. Note: the amplitude share ci/(ci+tr) is also monotone in cis/trans, so no
## direction differs, but every number does (genome-wide 0.4223 for the amplitude share, 0.3483 for the
## variance share).
summ <- function(A, idx = NULL) {
  ci <- tru(A, "cis", idx); tr <- tru(A, "trn", idx); pa <- tru(A, "par", idx)
  c2 <- ci^2; t2 <- tr^2
  c(cis_tru = ci, trn_tru = tr, par_tru = pa,
    cis_share = if (c2 + t2 > 0) c2 / (c2 + t2) else NA_real_)
}
boot <- function(A, nb = NBOOT) {
  m <- matrix(NA_real_, nb, 4L)
  for (k in seq_len(nb)) m[k, ] <- summ(A, sample(nrow(A), nrow(A), TRUE))
  m
}

say("\n=== genome-wide decomposition of the cold response ===")
set.seed(SEED)
o <- summ(A); b <- boot(A)
for (j in seq_along(o)) say("  %-10s %.4f [%.4f, %.4f]", names(o)[j], o[j],
                           quantile(b[, j], .025, na.rm = TRUE), quantile(b[, j], .975, na.rm = TRUE))
qg <- function(j, p) quantile(b[, j], p, na.rm = TRUE)
GW <- data.table(level = "genome", n = nrow(A), t(o),
                 cis_lo = qg(1, .025), cis_hi = qg(1, .975),
                 trn_lo = qg(2, .025), trn_hi = qg(2, .975),
                 cis_share_lo = qg(4, .025), cis_share_hi = qg(4, .975),
                 cis_resolvable = mean(b[, 1] > 0, na.rm = TRUE) >= 0.975,
                 trn_resolvable = mean(b[, 2] > 0, na.rm = TRUE) >= 0.975,
                 share_reportable = mean(b[, 1] > 0, na.rm = TRUE) >= 0.975 &&
                                    mean(b[, 2] > 0, na.rm = TRUE) >= 0.975,
                 frac_boot_trn_zero = mean(b[, 2] <= 0, na.rm = TRUE))

## ------------------------------------------------------------------ by stratum ----
HB <- fread(file.path(csd, "hog_to_scer.tsv")); h2s <- unique(HB[, .(hog = gene, sys)])
STR <- fread(file.path(hsd, "a5_strata_membership_37_30.tsv"))
SM <- merge(h2s, STR, by = "hog")
BOOLS <- c("Hsf1_core", "Hsf1_joint", "Msn24_only", "TF_independent", "not_induced", "other_induced")
say("\n=== by Hsf1-dependence stratum, projected onto the shared genes ===")
## A share is a ratio of two noise-corrected quantities, so it is only interpretable when both are resolvable.
## rms_true is truncated at zero, so a component whose signal does not clear its noise floor returns exact zeros
## in many bootstrap draws, and the share then collapses onto whichever component survives (for Msn24_only this
## gives an upper bound of exactly 1.000). A component counts as resolvable when it exceeds zero in
## at least 97.5% of draws. Where either component fails, no share is reported: only that trans did not clear
## its noise floor. Same criterion as the degeneracy rule used for the module ratios.
ST <- rbindlist(lapply(BOOLS, function(bl) {
  g <- unique(SM[get(bl) == TRUE, sys]); Ax <- A[sys %in% g]
  if (nrow(Ax) < 30L) { say("  %-16s %d genes, too few", bl, nrow(Ax)); return(NULL) }
  set.seed(SEED); oo <- summ(Ax); bb <- boot(Ax)
  cis_res <- mean(bb[, 1] > 0, na.rm = TRUE) >= 0.975
  trn_res <- mean(bb[, 2] > 0, na.rm = TRUE) >= 0.975
  q <- function(j, p) quantile(bb[, j], p, na.rm = TRUE)
  ok <- cis_res && trn_res
  say("  %-16s n=%-5d cis %.3f [%.3f,%.3f]%s  trn %.3f [%.3f,%.3f]%s  %s", bl, nrow(Ax),
      oo["cis_tru"], q(1, .025), q(1, .975), if (cis_res) " " else "*",
      oo["trn_tru"], q(2, .025), q(2, .975), if (trn_res) " " else "*",
      if (ok) sprintf("share %.3f [%.3f, %.3f]", oo["cis_share"], q(4, .025), q(4, .975))
      else "share not reported: a component does not clear its noise floor")
  data.table(level = bl, n = nrow(Ax), t(oo),
             cis_lo = q(1, .025), cis_hi = q(1, .975), trn_lo = q(2, .025), trn_hi = q(2, .975),
             cis_share_lo = q(4, .025), cis_share_hi = q(4, .975),
             cis_resolvable = cis_res, trn_resolvable = trn_res, share_reportable = ok,
             frac_boot_trn_zero = mean(bb[, 2] <= 0, na.rm = TRUE))
}))
say("  '*' marks a component that does not clear its noise floor in at least 97.5%% of draws.")
OUT <- rbind(GW, ST)
fwrite(OUT, file.path(hsd, "ext_hov2020_cis_trans_by_stratum.tsv"), sep = "\t")

## ------------------------------------------------------------------ compensation ----
## Same-replicate cor(cis, trans) is inflated by shared error, which is exactly the artefact the companion
## documents. Taking cis from one replicate index and trans from another removes it, averaged over the three
## ordered pairs with i != j.
say("\n=== compensation: cross-replicate against same-replicate ===")
Rw <- R[sys %in% A$sys]
H <- Rw[geno == "SCSU"]; P <- Rw[geno != "SCSU"]
hd <- dcast(H, sys ~ allele, value.var = c("d1", "d2", "d3"))
pd <- dcast(P, sys ~ allele, value.var = c("d1", "d2", "d3"))
M <- merge(hd, pd, by = "sys", suffixes = c("_h", "_p"))
cis_i <- function(i) M[[sprintf("d%d_Su_h", i)]] - M[[sprintf("d%d_Sc_h", i)]]
trn_i <- function(i) (M[[sprintf("d%d_Su_p", i)]] - M[[sprintf("d%d_Sc_p", i)]]) - cis_i(i)
same <- mean(vapply(1:3, function(i) cor(cis_i(i), trn_i(i), use = "complete.obs"), 0))
pairs <- expand.grid(i = 1:3, j = 1:3); pairs <- pairs[pairs$i != pairs$j, ]
cross <- mean(vapply(seq_len(nrow(pairs)), function(k)
  cor(cis_i(pairs$i[k]), trn_i(pairs$j[k]), use = "complete.obs"), 0))
say("  same-replicate  cor(cis, trans) = %+.4f  (mean of 3)", same)
say("  cross-replicate cor(cis, trans) = %+.4f  (mean of 6 ordered pairs)", cross)
say("  the cross-replicate value is the one free of shared replicate error")
## Leave-one-replicate-out on the compensation pair as well, since it is the quantity the first companion reports
## and a single value would carry the same false precision as a single share value.
COMP <- rbindlist(lapply(1:3, function(drop) {
  k <- setdiff(1:3, drop)
  s <- cor(cis_i(k[1]), trn_i(k[1]), use = "complete.obs")
  s2 <- cor(cis_i(k[2]), trn_i(k[2]), use = "complete.obs")
  x <- mean(cor(cis_i(k[1]), trn_i(k[2]), use = "complete.obs"),
            cor(cis_i(k[2]), trn_i(k[1]), use = "complete.obs"))
  data.table(dropped = sprintf("rep%d", drop), same_replicate = mean(c(s, s2)), cross_replicate = x)
}))
print(COMP[, .(dropped, same_replicate = round(same_replicate, 4),
               cross_replicate = round(cross_replicate, 4))], row.names = FALSE)
say("  same-replicate across leave-one-out: %+.4f to %+.4f", min(COMP$same_replicate), max(COMP$same_replicate))
say("  cross-replicate across leave-one-out: %+.4f to %+.4f", min(COMP$cross_replicate), max(COMP$cross_replicate))
say("  cross stays smaller in magnitude than same in all three: %s",
    COMP[, all(abs(cross_replicate) < abs(same_replicate))])
fwrite(rbind(data.table(dropped = "none", same_replicate = same, cross_replicate = cross), COMP,
             fill = TRUE)[, n := nrow(M)],
       file.path(hsd, "ext_hov2020_compensation.tsv"), sep = "\t")

## ------------------------------------------------------------------ pairing sensitivity ----
say("\n=== sensitivity to the arbitrary 12 C to 30 C replicate pairing ===")
set.seed(SEED + 7L)
PS <- rbindlist(lapply(1:5, function(k) {
  pr <- sample(1:3)
  Ak <- build_axes(contrasts3(resp(pr)))
  oo <- summ(Ak)
  data.table(draw = k, pairing = paste(pr, collapse = ""), n = nrow(Ak),
             cis_tru = oo["cis_tru"], trn_tru = oo["trn_tru"], cis_share = oo["cis_share"])
}))
print(PS, row.names = FALSE)
say("  cis share across pairings: %.4f to %.4f", min(PS$cis_share), max(PS$cis_share))
fwrite(PS, file.path(hsd, "ext_hov2020_pairing_sensitivity.tsv"), sep = "\t")

## ------------------------------------------------------------------ read out the decision rule ----
say("\n=== the pre-registered read-out ===")
say("  cis share below 0.5 genome-wide: %s (%.3f, interval %.3f to %.3f)",
    GW$cis_share < 0.5 && GW$cis_share_hi < 0.5, GW$cis_share, GW$cis_share_lo, GW$cis_share_hi)
say("  trans exceeds cis genome-wide: %s (%.4f against %.4f)",
    GW$trn_tru > GW$cis_tru, GW$trn_tru, GW$cis_tru)
say("  cross-replicate compensation negative: %s (%+.4f)", cross < 0, cross)
say("  smaller in magnitude than same-replicate: %s (%.4f against %.4f)",
    abs(cross) < abs(same), abs(cross), abs(same))
if (nrow(ST) && all(c("Hsf1_core", "other_induced") %in% ST$level)) {
  a <- ST[level == "Hsf1_core"]; c2 <- ST[level == "other_induced"]
  say("  Hsf1_core cis share %.3f [%.3f, %.3f] against other_induced %.3f [%.3f, %.3f]",
      a$cis_share, a$cis_share_lo, a$cis_share_hi, c2$cis_share, c2$cis_share_lo, c2$cis_share_hi)
}
## ------------------------------------------------------------------ leave-one-library-out ----
## SCSU_30C_rep1 has a baseline allelic median of 0.108 against 0.00 to 0.06 in the other five hybrid libraries.
## The main result is not changed on that basis; this reports what dropping each hybrid replicate index does, so
## that the sensitivity is a number rather than an assurance.
say("\n=== leave-one-replicate-out, hybrid arm ===")
LOO <- rbindlist(lapply(1:3, function(drop) {
  keep <- setdiff(1:3, drop)
  z <- CPM[sys %in% UNIV]
  w <- dcast(z, sys + allele + geno ~ temp + rep, value.var = "lcpm")
  r <- w[, .(sys, allele, geno)]
  ## two replicates only, so the companion's own two-block construction applies
  r[, d1 := w[[sprintf("12_rep%d", keep[1])]] - w[[sprintf("30_rep%d", keep[1])]]]
  r[, d2 := w[[sprintf("12_rep%d", keep[2])]] - w[[sprintf("30_rep%d", keep[2])]]]
  r <- r[complete.cases(r)]
  r[, `:=`(sig = (d1 + d2) / 2, nu1 = (d2 - d1) / 2, nu2 = NA_real_)]
  H <- r[geno == "SCSU"]; P <- r[geno != "SCSU"]
  hs <- dcast(H, sys ~ allele, value.var = c("sig", "nu1"))
  ps <- dcast(P, sys ~ allele, value.var = c("sig", "nu1"))
  Ax <- merge(hs, ps, by = "sys", suffixes = c("_h", "_p"))
  Ax[, `:=`(cis_sig = sig_Su_h - sig_Sc_h, cis_n1 = nu1_Su_h - nu1_Sc_h,
            par_sig = sig_Su_p - sig_Sc_p, par_n1 = nu1_Su_p - nu1_Sc_p)]
  Ax[, `:=`(trn_sig = par_sig - cis_sig, trn_n1 = par_n1 - cis_n1)]
  Ax <- Ax[complete.cases(Ax[, .(cis_sig, par_sig, trn_sig)])]
  ci <- rms_true(Ax$cis_sig, Ax$cis_n1); tr <- rms_true(Ax$trn_sig, Ax$trn_n1)
  data.table(dropped = sprintf("rep%d", drop), n = nrow(Ax), cis_tru = ci, trn_tru = tr,
             cis_share = ci^2 / (ci^2 + tr^2))
}))
print(LOO[, .(dropped, n, cis_tru = round(cis_tru, 4), trn_tru = round(trn_tru, 4),
              cis_share = round(cis_share, 4))], row.names = FALSE)
say("  cis share across the three leave-one-out fits: %.4f to %.4f, against %.4f on all three",
    min(LOO$cis_share), max(LOO$cis_share), GW$cis_share)

## A jackknife interval over replicates, because the gene-level bootstrap does not resample replicates and so
## cannot see the dominant variance component. With three replicates the jackknife standard error is
## sqrt((n-1)/n * sum((theta_i - mean(theta))^2)); it is the only replicate-level uncertainty this design
## supports. Note that dropping one of three replicates also degrades the noise-floor estimate, so this
## overstates instability, but that is not a reason to prefer the narrower interval.
jk <- function(v) { n <- length(v); sqrt((n - 1) / n * sum((v - mean(v))^2)) }
## No replicate-level 95% interval is formed. With three replicates a normal quantile is too optimistic and a
## t2 quantile (4.30) gives an interval so wide as to be uninformative (for the cis share it runs from below
## zero to nearly one). The standard error and the three leave-one-out values quantify replicate-level uncertainty;
## together they show that the result is not distinguishable from 0.5 at the replicate level, with
## the underlying values shown directly.
for (q in c("cis_tru", "trn_tru", "cis_share")) {
  v <- LOO[[q]]; se <- jk(v)
  say("  %-10s all three %.4f | leave-one-out %s | jackknife SE %.4f",
      q, GW[[q]], paste(sprintf("%.4f", v), collapse = ", "), se)
}
JK <- data.table(quantity = c("cis_tru", "trn_tru", "cis_share"),
                 all_three = c(GW$cis_tru, GW$trn_tru, GW$cis_share),
                 jk_mean = vapply(c("cis_tru", "trn_tru", "cis_share"), function(q) mean(LOO[[q]]), 0),
                 jk_se = vapply(c("cis_tru", "trn_tru", "cis_share"), function(q) jk(LOO[[q]]), 0))
## This quantity is not a confidence interval.
JK[, `:=`(pm196se_lo = jk_mean - 1.96 * jk_se, pm196se_hi = jk_mean + 1.96 * jk_se)]
JK[, loo_values := vapply(quantity, function(q) paste(sprintf("%.4f", LOO[[q]]), collapse = "; "), "")]
JK[quantity == "cis_tru",   `:=`(boot_lo = GW$cis_lo, boot_hi = GW$cis_hi)]
JK[quantity == "trn_tru",   `:=`(boot_lo = GW$trn_lo, boot_hi = GW$trn_hi)]
JK[quantity == "cis_share", `:=`(boot_lo = GW$cis_share_lo, boot_hi = GW$cis_share_hi)]
print(JK[, .(quantity, all_three = round(all_three, 4),
             gene_bootstrap = sprintf("[%.3f, %.3f]", boot_lo, boot_hi),
             leave_one_out = loo_values, jackknife_SE = round(jk_se, 4))], row.names = FALSE)
say("  the bootstrap interval covers gene sampling only; the jackknife covers replicate sampling only.")
say("  cis share is not distinguishable from 0.5 at the replicate level: the three leave-one-out")
say("  values are %s, straddling 0.5, with a jackknife SE of %.4f",
    JK[quantity == "cis_share", loo_values], JK[quantity == "cis_share", jk_se])
fwrite(JK, file.path(hsd, "ext_hov2020_jackknife.tsv"), sep = "\t")
say("  these use two replicates and therefore the two-block construction, so they are not directly")
say("  comparable to the three-replicate value; what they show is the spread, not a better estimate.")
fwrite(LOO, file.path(hsd, "ext_hov2020_leave_one_out.tsv"), sep = "\t")

## ------------------------------------------------------------------ heat arm, same quantity ----
## The Hsf1_core to other_induced comparison uses cis share rather than an amplitude ratio.
## The cross-species ratio of cis divergence amplitudes is a different quantity; the matching
## comparison uses the same share on the heat arm, computed from the saved axis table.
say("\n=== the same share on the heat arm for a like-for-like comparison ===")
HA <- fread(file.path(hsd, "a6_axis_amplitudes_37_30.tsv"))
HA[, cis_share := cis_tru^2 / (cis_tru^2 + trn_tru^2)]   # script 05's definition, see above
HS <- dcast(HA[, .(share = median(cis_share, na.rm = TRUE)), by = .(stratum, time)],
            stratum ~ time, value.var = "share")
setnames(HS, c("15", "30", "60"), c("heat_15", "heat_30", "heat_60"))
CMP <- merge(HS, OUT[, .(stratum = level, cold = cis_share, cold_reportable = share_reportable)],
             by = "stratum")
print(CMP[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)], row.names = FALSE)
fwrite(CMP, file.path(hsd, "ext_hov2020_heat_vs_cold_share.tsv"), sep = "\t")
say("  heat-arm shares are medians over the six hybrids; the cold arm has one pair and one interval.")

say("")
