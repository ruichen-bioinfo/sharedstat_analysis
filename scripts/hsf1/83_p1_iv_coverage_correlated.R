#!/usr/bin/env Rscript
## 83_p1_iv_coverage_correlated.R -- coverage of the IV bootstrap interval when gene errors are correlated.
## Same 96-cell grid, same seeds and same estimator as 75_p1_iv_coverage.R; the only change is the error model. Each of
## the three error vectors (reference arm eR, response arm eA, instrument noise) gets its own one-factor structure,
##   e_i = s * ( sqrt(w) * lambda_i * F + sqrt(1 - w) * z_i ),  lambda_i ~ N(1, 1), F ~ N(0, 1), z_i ~ N(0, 1),
## rescaled so its variance is the independent model's. Cov(e_i, e_j) = w E[lambda_i lambda_j] = w and
## Var(e_i) = 2w + 1 - w = 1 + w, so the mean pairwise correlation is w / (1 + w); w = rho / (1 - rho). Factor loadings vary by
## gene, so the correlation is not a common shift that centring removes. Targets: rho = 0.017, the largest stratum
## median of the allele-contrast residual correlation measured in P1's data (82_p1_h_correlation.R), and rho = 0.05 as
## a stress level. Usage: Rscript 83_p1_iv_coverage_correlated.R <rho>
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath(".")); hsd <- file.path(root, "data/02_processed/hsf1")
RHO <- as.numeric(commandArgs(trailingOnly = TRUE)[1]); stopifnot(is.finite(RHO), RHO > 0, RHO < 0.5)
W <- RHO / (1 - RHO)
source(file.path(root, "pkg/sharedstat/R/estimators.R"))
src <- parse(file.path(root, "scripts/hsf1/31_p1_design_grid.R"))
keep <- vapply(src, function(e) grepl("^(N_GENE|NSIM|GRID) <-", paste(deparse(e), collapse = " ")), TRUE)
for (e in src[keep]) eval(e)
NBOOT <- 200L
cerr <- function(n, sd_) { lam <- rnorm(n, 1, 1); x <- sqrt(W) * lam * rnorm(1) + sqrt(1 - W) * rnorm(n)
  x / sqrt(1 + W) * sd_ }                      # Var(lam F) = E[lam^2] = 2
one_ci <- function(gamma, h, n_rep, cor_S, seed) {
  set.seed(seed)
  sR <- 1; sE <- sqrt(sR^2 * (1 - h) / h * n_rep)
  R <- stats::rnorm(N_GENE, 0, sR); eR <- cerr(N_GENE, sE / sqrt(n_rep)); eA <- cerr(N_GENE, sE / sqrt(n_rep))
  Rhat <- R + eR; e <- gamma * R + eA - Rhat
  k <- sqrt(1 / cor_S^2 - 1); S <- R + cerr(N_GENE, sR * k)
  f <- iv_gain(e, Rhat, S, naive = FALSE, n_boot = NBOOT, seed = seed)
  c(gamma_iv = f$gamma, lo = f$ci[1], hi = f$ci[2])
}
## check the realised correlation on one draw
set.seed(1); M <- sapply(1:200, function(i) cerr(300, 1)); rr <- cor(t(M)); cat(sprintf("target rho %.3f, realised mean pairwise correlation %.4f\n", RHO, mean(rr[upper.tri(rr)])))
t0 <- Sys.time()
RES <- rbindlist(lapply(seq_len(nrow(GRID)), function(i) { g <- GRID[i]
  m <- t(vapply(seq_len(NSIM), function(s) one_ci(g$gamma, g$h, g$n_rep, g$cor_S, seed = i * 1000L + s), numeric(3)))
  data.table(g, rho = RHO, iv_mean = mean(m[, 1]), coverage = mean(m[, 2] <= g$gamma & g$gamma <= m[, 3]),
             below = mean(m[, 3] < g$gamma), above = mean(m[, 2] > g$gamma)) }))
fwrite(RES, file.path(hsd, sprintf("p1_iv_coverage_rho%03d.tsv", round(1000 * RHO))), sep = "\t")
cat(sprintf("=== 83 | rho %.3f | %.0f s ===\npooled coverage %.4f (nominal 0.95); per cell %.3f to %.3f; misses below %.4f above %.4f; mean IV bias %.4f\n",
    RHO, as.numeric(difftime(Sys.time(), t0, units = "secs")), mean(RES$coverage), min(RES$coverage), max(RES$coverage),
    mean(RES$below), mean(RES$above), mean(RES$iv_mean - RES$gamma)))
