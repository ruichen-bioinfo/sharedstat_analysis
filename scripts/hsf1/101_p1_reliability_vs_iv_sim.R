#!/usr/bin/env Rscript
## 101_p1_reliability_vs_iv_sim.R -- simulation comparing the naive, reliability-corrected (naive/h) and instrumented
## gain estimators at true gain 1.1 and 1.3 and reliability 0.7 and 0.9, in five scenarios: clean; gain heterogeneity
## (sd 0.2) independent, correlated 0.3 with R, or correlated 0.3 with the instrument; and leakage of 0.3 of R's
## block-1 error into the instrument. n = 2000 genes, 300 replicates per cell; bias, sd, rmse and percent of estimates
## below 1. The design was fixed before this analysis was run.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
for (f in c("reliability.R", "estimators.R")) source(file.path(root, "pkg/sharedstat/R", f))
set.seed(20261002L); n <- 2000L; NR <- 300L
SC <- list(clean = c(sd = 0, rR = 0, rS = 0, k = 0), `het-indep` = c(0.2, 0, 0, 0), `het-corR` = c(0.2, 0.3, 0, 0),
           `het-corS` = c(0.2, 0, 0.3, 0), leak = c(0, 0, 0, 0.3))
one <- function(g, h, p) {
  su <- sqrt(2 * (1 - h) / h)              # var(R_hat) = 1 + su^2/2, h = 1/(1 + su^2/2)
  R <- rnorm(n); u1 <- rnorm(n, 0, su); u2 <- rnorm(n, 0, su); w <- rnorm(n)
  zR <- (R - mean(R)) / sd(R); zw <- (w - mean(w)) / sd(w)
  d <- p[1] * (p[2] * zR + p[3] * zw + sqrt(1 - p[2]^2 - p[3]^2) * rnorm(n))
  A <- g * R + d + rnorm(n, 0, su / sqrt(2)); Rh <- (R + u1 + R + u2) / 2; e <- A - Rh
  S <- R + w + p[4] * u1
  nv <- 1 + cov(e, Rh) / var(Rh); hh <- reliability_from_blocks(R + u1, R + u2, center = TRUE)$h
  c(naive = nv, reliability_corrected = nv / hh, instrumented = iv_gain(e, Rh, S)$gamma, h_hat = hh) }
OUT <- rbindlist(lapply(c(1.1, 1.3), function(g) rbindlist(lapply(c(0.7, 0.9), function(h) rbindlist(lapply(names(SC), function(s) {
  X <- t(replicate(NR, one(g, h, unname(SC[[s]]))))
  rbindlist(lapply(c("naive", "reliability_corrected", "instrumented"), function(est) { x <- X[, est]
    data.table(gamma = g, h = h, scenario = s, estimator = est, bias = mean(x) - g, sd = sd(x),
               rmse = sqrt(mean((x - g)^2)), pct_wrong_side = 100 * mean(x < 1), mean_h_hat = mean(X[, "h_hat"])) })) }))))))
fwrite(OUT, file.path(root, "data/02_processed/hsf1/p1_reliability_vs_iv_sim.tsv"), sep = "\t")
print(OUT[, .(gamma, h, scenario, estimator, bias = round(bias, 4), sd = round(sd, 4), rmse = round(rmse, 4), wrong = round(pct_wrong_side, 1))], nrows = 200)
