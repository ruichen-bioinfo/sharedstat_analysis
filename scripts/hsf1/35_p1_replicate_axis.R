#!/usr/bin/env Rscript
## 35_p1_replicate_axis.R -- the analytic replicate-axis result for P1.
## 
##
## A replicate grid varying n_rep over {2, 3, 6, 12} is circular when the
## per-observation error sd set as
##       sE <- sqrt(sR^2 * (1 - h) / h * n_rep)   and then used   sE / sqrt(n_rep)
## and the n_rep cancels exactly: the effective error sd is identical at every replicate count. That is a
## parameterization makes h and n_rep appear separately interpretable. However,
## since the bias is a function of h alone, holding h fixed along the replicate axis guarantees a constant
## bias, so unchanged bias cannot show that averaging more replicates does not help; the grid
## reflects its parameterization.
##
## What holds instead. The naive estimator has a closed form,
## gamma_hat = gamma * h, so the bias gamma*(h - 1) depends on the reference arm's reliability and on
## nothing else. The two sample sizes an experiment can increase therefore behave differently:
##   features  -- h is a variance ratio, not a sampling error, so the estimate converges to gamma*h.
##                Bias is O(1). More sequencing cannot help at all.
##   replicates-- averaging n_rep observations divides the measurement variance by n_rep, so
##                h = 1/(1 + tau/n_rep) and the bias is O(1/n_rep). Replication helps, through h only,
##                and never reaches zero because h < 1 for every finite n_rep.
## This statement is exact rather than simulated and retains the practical implication,
## because linear improvement against a bias that an instrument removes outright is a poor trade.

suppressPackageStartupMessages(library(data.table))
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
outd <- file.path(root, "data/02_processed/hsf1"); dir.create(outd, showWarnings = FALSE, recursive = TRUE)

## the closed form, and a simulation that is free to disagree with it
h_of  <- function(tau, n_rep) 1 / (1 + tau / n_rep)
naive <- function(gamma, h) gamma * h

sim_cell <- function(gamma, tau, n_rep, n_gene, cor_S = 0.8, nsim = 200L, seed = 1L) {
  set.seed(seed)
  o <- vapply(seq_len(nsim), function(s) {
    vu <- tau / n_rep                       # per-observation variance tau, averaged over n_rep
    R  <- stats::rnorm(n_gene)
    Rh <- R + stats::rnorm(n_gene, 0, sqrt(vu))
    e  <- gamma * R + stats::rnorm(n_gene, 0, sqrt(vu)) - Rh
    S  <- cor_S * R + stats::rnorm(n_gene, 0, sqrt(1 - cor_S^2))
    c(stats::cov(e, Rh) / stats::var(Rh) + 1,
      stats::cov(e, S)  / stats::cov(Rh, S) + 1,
      stats::var(R) / stats::var(Rh))
  }, numeric(3))
  data.table(gamma = gamma, tau = tau, n_rep = n_rep, n_gene = n_gene,
             h_emp = mean(o[3, ]), h_pred = h_of(tau, n_rep),
             naive_emp = mean(o[1, ]), naive_pred = naive(gamma, h_of(tau, n_rep)),
             iv_emp = mean(o[2, ]),
             sign_wrong = mean((gamma > 1 & o[1, ] < 1) | (gamma < 1 & o[1, ] > 1)))
}

say("================ AXIS A: FEATURES. Does the bias fall as more genes are measured? ================")
A <- rbindlist(lapply(c(500L, 2000L, 8000L, 32000L, 128000L),
                      function(ng) sim_cell(1.3, tau = 1, n_rep = 3L, n_gene = ng, seed = ng)))
say("%9s %8s %11s %11s %10s %11s", "n_gene", "h", "naive", "predicted", "IV", "bias")
for (i in seq_len(nrow(A))) say("%9d %8.4f %11.4f %11.4f %10.4f %+11.4f",
  A$n_gene[i], A$h_emp[i], A$naive_emp[i], A$naive_pred[i], A$iv_emp[i], A$naive_emp[i] - A$gamma[i])
say("range of the bias over a 256-fold change in features: %.4f", diff(range(A$naive_emp - A$gamma)))
stopifnot(diff(range(A$naive_emp - A$gamma)) < 0.02)     # O(1): it does not move

say("\n================ AXIS B: REPLICATES. Does the bias fall as replicates are added? ================")
B <- rbindlist(lapply(c(2L, 3L, 6L, 12L, 48L, 192L),
                      function(nr) sim_cell(1.3, tau = 1, n_rep = nr, n_gene = 4000L, seed = 100L + nr)))
say("%7s %8s %9s %11s %11s %10s %11s", "n_rep", "h", "h_pred", "naive", "predicted", "IV", "bias")
for (i in seq_len(nrow(B))) say("%7d %8.4f %9.4f %11.4f %11.4f %10.4f %+11.4f",
  B$n_rep[i], B$h_emp[i], B$h_pred[i], B$naive_emp[i], B$naive_pred[i], B$iv_emp[i], B$naive_emp[i] - B$gamma[i])
bias <- abs(B$naive_emp - B$gamma)
say("bias falls monotonically: %s ; ratio first/last = %.1f-fold",
    all(diff(bias) < 0), bias[1] / bias[length(bias)])
stopifnot(all(diff(bias) < 0), bias[1] / bias[length(bias)] > 10)   # O(1/n_rep): it does move

say("\n=== the closed form reproduces the simulation ===")
say("max |naive_emp - gamma*h_pred| over both axes: %.5f",
    max(abs(c(A$naive_emp - A$naive_pred, B$naive_emp - B$naive_pred))))
stopifnot(max(abs(c(A$naive_emp - A$naive_pred, B$naive_emp - B$naive_pred))) < 0.01)

say("\n=== the constant-h construction, shown arithmetically ===")
for (nr in c(2L, 12L)) {
  h <- 0.5; sE <- sqrt(1 * (1 - h) / h * nr)
  say("  n_rep = %2d: sE = %.6f, and sE/sqrt(n_rep) = %.6f", nr, sE, sE / sqrt(nr))
}
say("  identical at every n_rep, so h was constant along that axis by construction")

OUT <- rbind(A[, axis := "features"], B[, axis := "replicates"], fill = TRUE)
fwrite(OUT, file.path(outd, "p1_replicate_axis.tsv"), sep = "\t")
say("\nwrote %s (%d rows)", file.path(outd, "p1_replicate_axis.tsv"), nrow(OUT))
