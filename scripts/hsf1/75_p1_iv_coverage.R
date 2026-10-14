#!/usr/bin/env Rscript
## 75_p1_iv_coverage.R -- does the instrumented estimator's 95% bootstrap interval cover the true gain?
##
## P1 reports IV intervals; this checks that they mean what they say. No new data: the same 96-cell design grid as 31_p1_design_grid.R, the same
## data-generating function (transplanted by parse, so it cannot drift), the same seeds, and the same estimator from
## the released package, here called with a gene-level percentile bootstrap (iv_gain(n_boot = 200)), which is the
## interval the real-data analysis uses (11_gain_iv.R, 200 resamples over genes).
## Coverage = share of the 200 datasets per cell whose interval contains the true gain. With 200 datasets the Monte
## Carlo standard error of a coverage near 0.95 is about 0.015, so a cell is flagged only below 0.95 - 2 x 0.015.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
hsd <- file.path(root, "data/02_processed/hsf1")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
source(file.path(root, "pkg/sharedstat/R/estimators.R"))
src <- parse(file.path(root, "scripts/hsf1/31_p1_design_grid.R"))
keep <- vapply(src, function(e) { s <- paste(deparse(e), collapse = " ")
  grepl("^(N_GENE|NSIM|one_dataset|GRID) <-", s) }, TRUE)
for (e in src[keep]) eval(e)
stopifnot(N_GENE == 2000L, NSIM == 200L, nrow(GRID) == 96L)
NBOOT <- 200L
one_ci <- function(gamma, h, n_rep, cor_S, seed) {
  ## identical draws to one_dataset(): same seed, same order of rnorm calls
  set.seed(seed)
  sR <- 1; sE <- sqrt(sR^2 * (1 - h) / h * n_rep)
  R <- stats::rnorm(N_GENE, 0, sR); eR <- stats::rnorm(N_GENE, 0, sE / sqrt(n_rep))
  eA <- stats::rnorm(N_GENE, 0, sE / sqrt(n_rep))
  Rhat <- R + eR; e <- gamma * R + eA - Rhat
  k <- sqrt(1 / cor_S^2 - 1); S <- R + stats::rnorm(N_GENE, 0, sR * k)
  f <- iv_gain(e, Rhat, S, naive = FALSE, n_boot = NBOOT, seed = seed)
  c(gamma_iv = f$gamma, lo = f$ci[1], hi = f$ci[2])
}
## fidelity: the point estimates must equal the stored grid's for a sample of cells
LG <- fread(file.path(hsd, "p1_design_grid.tsv"))
t0 <- Sys.time()
RES <- rbindlist(lapply(seq_len(nrow(GRID)), function(i) {
  g <- GRID[i]
  m <- t(vapply(seq_len(NSIM), function(s) one_ci(g$gamma, g$h, g$n_rep, g$cor_S, seed = i * 1000L + s), numeric(3)))
  data.table(g, iv_mean = mean(m[, 1]), coverage = mean(m[, 2] <= g$gamma & g$gamma <= m[, 3]),
             below = mean(m[, 3] < g$gamma), above = mean(m[, 2] > g$gamma), width_median = median(m[, 3] - m[, 2]))
}))
chk <- merge(RES, LG, by = c("gamma", "h", "n_rep", "cor_S"))
stopifnot(nrow(chk) == 96L)
dev <- max(abs(chk$iv_mean - chk$gamma_iv))
fwrite(RES, file.path(hsd, "p1_iv_coverage.tsv"), sep = "\t")
say("=== 75_p1_iv_coverage | %s | %.0f s ===", format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "secs")))
say("96 cells x %d datasets x %d resamples; mean IV estimate reproduces the saved grid to %.2e", NSIM, NBOOT, dev)
say("nominal 95%%: achieved coverage median %.3f, range %.3f to %.3f; cells below 0.92: %d of 96",
    median(RES$coverage), min(RES$coverage), max(RES$coverage), sum(RES$coverage < 0.92))
print(RES[, .(coverage_median = median(coverage), coverage_min = min(coverage), width_median = median(width_median)),
          by = .(h, cor_S)][order(h, cor_S)])
print(RES[, .(coverage_median = median(coverage), coverage_min = min(coverage)), by = gamma])
if (any(RES$coverage < 0.92)) print(RES[coverage < 0.92][order(coverage)])
