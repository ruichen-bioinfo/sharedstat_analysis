#!/usr/bin/env Rscript
## 64_b5_centred_reliability_with_ci.R -- centred reliability estimates with bootstrap intervals.
##
## Recomputes every reported reliability with center = TRUE, attaches a 95% bootstrap interval from resampling
## GENES (not replicates: the features are the sampling unit, the replicate structure is the design), and reports
## how many cells have an upper bound still below 1/1.1 -- the cell count that survives the interval, which is a
## stronger statement than a point estimate below the threshold.
##
## Output columns: h_centred, h_uncentred, lo, hi, upper_below_1_1.
##
## Two guards. The uncentred values of script 36 must be reproduced by center = FALSE before anything else is read off the
## table, and the bootstrap must be seeded, or the interval is not the same interval tomorrow.

suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
for (f in c("reliability.R", "reliability_replicates.R")) source(file.path(root, "pkg/sharedstat/R", f))
hsd <- file.path(root, "data/02_processed/hsf1")
ext <- file.path(root, "data/03_external")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
SEED <- 20260915L
B <- 2000L
THRESH <- 1 / 1.1
OUT <- list()

## Bootstrap over features. The estimator is a ratio of mean squares over genes, so resampling genes with
## replacement is the interval that answers "how much of this reliability is the particular gene set".
boot_ci <- function(make_h, n, B = 2000L, seed = SEED) {
  set.seed(seed)
  v <- vapply(seq_len(B), function(b) make_h(sample.int(n, n, replace = TRUE)), 0)
  v <- v[is.finite(v)]
  if (!length(v)) return(c(NA_real_, NA_real_))
  unname(quantile(v, c(0.025, 0.975)))
}

add <- function(dataset, arm, k, n, hT, hF, ci) {
  OUT[[length(OUT) + 1L]] <<- data.table(
    dataset = dataset, arm = arm, k = k, n_genes = n,
    h_centred = hT, h_uncentred = hF, lo = ci[1], hi = ci[2],
    upper_below_1_1 = ci[2] < THRESH)
}

## ------------------------------------------------------------------ 1. Fay, 108 cells, two blocks ----
say("=== Fay 2023: 108 cells, two replicate blocks ===")
D <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a5_cis_dod_37_30.tsv.gz"))))
P <- fread(file.path(hsd, "p1_measured_reliability.tsv"))
## Eligibility is a named list of the six consensus S. cerevisiae Hsf1-dependence classes, not "whatever label is
## not blank". Note: a test of `stratum != ""` admits genes without a consensus class once they carry the literal
## label "unassigned", and the primary analysis then grows from 108 cells to 126. A scientific exclusion must not
## depend on how a group happens to be spelled.
ELIGIBLE <- c("Hsf1_core", "Hsf1_joint", "Msn24_only", "TF_independent", "not_induced", "other_induced")
## Any label that is not one of the six is the remainder: genes the three calling methods do not agree on. It is
## measured and reported as a sensitivity stratum, never folded into the primary set and never dropped in silence.
D[, eligible := stratum %in% ELIGIBLE]
say("  strata present in the contrast table: %s", paste(sort(unique(D$stratum)), collapse = ", "))
n_rem <- uniqueN(D[eligible == FALSE, hog])
say("  %d of the labels are eligible; the remainder covers %d orthologue groups in %d rows",
    length(ELIGIBLE), n_rem, nrow(D[eligible == FALSE]))
stopifnot(length(setdiff(unique(D$stratum), c(ELIGIBLE, "", NA_character_, "unassigned"))) == 0L)
cells <- unique(D[eligible == TRUE, .(stratum, hybrid, time)])
say("  %d primary cells from %d eligible strata", nrow(cells), uniqueN(cells$stratum))
stopifnot(nrow(cells) == 108L, uniqueN(cells$stratum) == 6L)
for (i in seq_len(nrow(cells))) {
  cl <- cells[i]
  z <- D[stratum == cl$stratum & hybrid == cl$hybrid & time == cl$time]
  k <- is.finite(z$sig) & is.finite(z$nul)
  b1 <- z$sig[k] - z$nul[k]; b2 <- z$sig[k] + z$nul[k]
  hT <- reliability_from_blocks(b1, b2, center = TRUE)$h
  hF <- reliability_from_blocks(b1, b2, center = FALSE)$h
  ci <- boot_ci(function(ix) reliability_from_blocks(b1[ix], b2[ix], center = TRUE)$h, length(b1), B,
                SEED + i)
  add("Fay 2023", sprintf("%s|%s|%s", cl$stratum, cl$hybrid, cl$time), 2L, length(b1), hT, hF, ci)
}
## The remainder, measured the same way through the same function, written to its own file. Reporting it is what
## makes the primary exclusion checkable: a reader can see whether dropping it flattered the result.
SENS <- list()
rem_cells <- unique(D[eligible == FALSE, .(stratum, hybrid, time)])
if (nrow(rem_cells) > 0L) {
  for (i in seq_len(nrow(rem_cells))) {
    cl <- rem_cells[i]
    z <- D[stratum == cl$stratum & hybrid == cl$hybrid & time == cl$time]
    k <- is.finite(z$sig) & is.finite(z$nul)
    b1 <- z$sig[k] - z$nul[k]; b2 <- z$sig[k] + z$nul[k]
    hT <- reliability_from_blocks(b1, b2, center = TRUE)$h
    hF <- reliability_from_blocks(b1, b2, center = FALSE)$h
    ci <- boot_ci(function(ix) reliability_from_blocks(b1[ix], b2[ix], center = TRUE)$h, length(b1), B,
                  SEED + 1000L + i)
    SENS[[length(SENS) + 1L]] <- data.table(
      stratum = "remainder_no_consensus_class", hybrid = cl$hybrid, time = cl$time,
      n_genes = length(b1), h_centred = hT, h_uncentred = hF, lo = ci[1], hi = ci[2],
      upper_below_1_1 = ci[2] < THRESH)
  }
  SE <- rbindlist(SENS)
  fwrite(SE, file.path(hsd, "p1x_b5_sensitivity_remainder.tsv"), sep = "\t")
  say("  sensitivity stratum: %d cells, median h %.3f, %d of %d with the interval below 1/1.1",
      nrow(SE), median(SE$h_centred), sum(SE$upper_below_1_1), nrow(SE))
  say("  landed p1x_b5_sensitivity_remainder.tsv")
}

FY <- rbindlist(OUT)
S <- merge(P[, .(key = sprintf("%s|%s|%s", stratum, hybrid, time), h_pub = h)],
           FY[, .(key = arm, h_uncentred, h_centred, lo, hi, upper_below_1_1)], by = "key")
stopifnot(nrow(S) == 108L)

## Label invariance. The primary cell count must not depend on how the remainder is spelled. Relabelling it in each
## of the ways below must leave the primary set untouched; a different count means eligibility is being decided by
## a string test somewhere.
for (lab in list("", NA_character_, "unassigned", "no_consensus")) {
  Dx <- copy(D)
  Dx[!stratum %in% ELIGIBLE, stratum := lab]
  n_lab <- nrow(unique(Dx[stratum %in% ELIGIBLE, .(stratum, hybrid, time)]))
  stopifnot(n_lab == 108L)
}
say("  label invariance: the primary set is 108 cells whether the remainder is blank, NA, 'unassigned' or otherwise")
worst <- max(abs(S$h_pub - S$h_uncentred))
say("  fidelity gate: max |script 36 (uncentred) - center=FALSE| = %.3e  %s", worst, if (worst < 1e-12) "PASS" else "FAIL")
stopifnot(worst < 1e-12)
say("  h_centred: min %.4f  median %.4f  max %.4f   (uncentred median %.4f)",
    min(S$h_centred), median(S$h_centred), max(S$h_centred), median(S$h_pub))
say("  point estimate below 1/1.1: %d of 108", sum(S$h_centred < THRESH))
say("  UPPER 95%% bound also below 1/1.1: %d of 108", sum(S$upper_below_1_1))
say("  widest interval: %.4f (%s)", max(S$hi - S$lo), S[which.max(hi - lo), key])

## ------------------------------------------------------------------ 2. Keyport Kik, three species ----
say("\n=== Keyport Kik 2024: three species, two replicate blocks ===")
kk <- file.path(ext, "keyportkik2024/GSE234499_counts.tsv.gz")
if (file.exists(kk)) {
  X <- fread(cmd = paste("zcat", shQuote(kk)))[!grepl("^__", ORF)]
  SM <- data.table(sample = paste0("SK", 1:12), cond = rep(c("ctrl", "ctrl", "shock", "shock"), 3),
                   br = rep(c(1L, 2L, 1L, 2L), 3))
  X <- merge(X, SM, by = "sample")
  X[, cpm := 1e6 * counts / sum(counts), by = sample][, lcpm := log2(cpm + 2)]
  for (sp in c("Skud", "Scer", "Kmarx")) {
    Z <- X[species == sp]
    keep <- Z[, .(ok = all(cpm >= 1)), by = ORF][ok == TRUE, ORF]
    W <- dcast(Z[ORF %in% keep], ORF ~ cond + br, value.var = "lcpm")
    r1 <- W$shock_1 - W$ctrl_1; r2 <- W$shock_2 - W$ctrl_2
    hT <- reliability_from_blocks(r1, r2, center = TRUE)$h
    hF <- reliability_from_blocks(r1, r2, center = FALSE)$h
    ci <- boot_ci(function(ix) reliability_from_blocks(r1[ix], r2[ix], center = TRUE)$h, length(r1), B)
    add("Keyport Kik 2024", sp, 2L, length(r1), hT, hF, ci)
    say("  %-6s n %5d  h %.4f (uncentred %.4f)  CI [%.4f, %.4f]  upper < 1/1.1: %s",
        sp, length(r1), hT, hF, ci[1], ci[2], ci[2] < THRESH)
  }
} else say("  counts absent; skipped")

## ------------------------------------------------------------------ 3. GSE310700, two genotypes ----
say("\n=== GSE310700: two genotypes, three replicates ===")
g7 <- file.path(ext, "gse310700/GSE310700_raw_counts.txt.gz")
if (file.exists(g7)) {
  X <- fread(cmd = paste("zcat", shQuote(g7))); setnames(X, 1, "gene")
  M <- data.table(col = names(X)[-1])
  M[, genotype := ifelse(grepl("^SHCHIP", col), "STUB1 knockdown", "control hairpin")]
  M[, temp := as.integer(sub("^SH(CTL|CHIP)([0-9]+)-.*$", "\\2", col))]
  M[, rep  := as.integer(sub(".*-", "", col))]
  stopifnot(all(is.finite(M$temp)), all(is.finite(M$rep)))
  L <- melt(X, id.vars = "gene", variable.name = "col", value.name = "v")
  L[, col := as.character(col)]
  L <- merge(L, M, by = "col")
  L[, cpm := 1e6 * v / sum(v), by = col][, lcpm := log2(cpm + 2)]
  keep <- L[, .(k = sum(cpm >= 1)), by = gene][k == 12L, gene]
  W <- dcast(L[gene %in% keep], gene ~ genotype + temp + rep, value.var = "lcpm")
  for (g in c("control hairpin", "STUB1 knockdown")) {
    m <- cbind(W[[paste0(g, "_42_1")]] - W[[paste0(g, "_37_1")]],
               W[[paste0(g, "_42_2")]] - W[[paste0(g, "_37_2")]],
               W[[paste0(g, "_42_3")]] - W[[paste0(g, "_37_3")]])
    m <- m[complete.cases(m), , drop = FALSE]
    hT <- suppressWarnings(reliability_from_replicates(m, center = TRUE)$h)
    hF <- suppressWarnings(reliability_from_replicates(m, center = FALSE)$h)
    ci <- boot_ci(function(ix) suppressWarnings(reliability_from_replicates(m[ix, , drop = FALSE],
                                                                          center = TRUE)$h), nrow(m), B)
    add("GSE310700", g, 3L, nrow(m), hT, hF, ci)
    say("  %-16s n %5d  h %.4f (uncentred %.4f)  CI [%.4f, %.4f]  upper < 1/1.1: %s",
        g, nrow(m), hT, hF, ci[1], ci[2], ci[2] < THRESH)
  }
} else say("  counts absent; skipped")

## ------------------------------------------------------------------ 4. Hovhannisyan, if quantified ----
say("\n=== Hovhannisyan 2020: four arms, three replicates ===")
qd <- file.path(ext, "hov2020/quant")
if (dir.exists(qd) && length(list.files(qd)) > 0) {
  say("  quantification present: the four arms are computed by 54_p1_hovhannisyan_dataset.R, not duplicated here")
  say("")
  say("")
} else {
  say("  quantification directory absent or empty: the four arms are NOT in this table.")
  say("  They are computed by 54_p1_hovhannisyan_dataset.R.")
}

## ------------------------------------------------------------------ write ----
R <- rbindlist(OUT)
setcolorder(R, c("dataset", "arm", "k", "n_genes", "h_centred", "h_uncentred", "lo", "hi", "upper_below_1_1"))
fwrite(R, file.path(hsd, "p1x_b5_reliability_centred_ci.tsv"), sep = "\t")
say("\n=== summary ===")
say("  %d arms, %d datasets, B = %d, seed = %d", nrow(R), uniqueN(R$dataset), B, SEED)
say("  %-18s %5s %8s %8s %8s %8s", "dataset", "arms", "min h", "median", "max h", "upper<thr")
for (d in unique(R$dataset)) {
  z <- R[dataset == d]
  say("  %-18s %5d %8.4f %8.4f %8.4f %8d", d, nrow(z), min(z$h_centred), median(z$h_centred),
      max(z$h_centred), sum(z$upper_below_1_1))
}
say("  landed: p1x_b5_reliability_centred_ci.tsv")
