#!/usr/bin/env Rscript
## 112_p3_constant_share_null.R -- a second null for the fall in the cis share (Supplementary Note N2, section 17).
##
## Script 108 holds each gene's cis and remainder differences fixed over time. That null cannot tell a change in the
## partition from a constant partition whose overall amplitude changes, and it draws noise per cross, so crosses that
## share a parent and timepoints that share the 0-min sample get independent noise. This null keeps the partition
## constant but lets its amplitude follow the data, and it adds noise to the individual libraries:
##   truth   in each cross, each gene's cis and parental differences are a fixed pattern (the mean over 15, 30 and 60 min
##           of its observed signal divided by that timepoint's noise-corrected parental amplitude) multiplied by the
##           observed parental amplitude at each timepoint. The cis share of the truth is therefore the same at every
##           timepoint, while the amplitude rises or falls as observed.
##   noise   for every library (parent or hybrid allele, timepoint including 0 min) each gene's noise is its own block
##           difference, centred over the four timepoints of that parent or allele and rescaled to one library's variance,
##           multiplied by a random sign drawn per gene, timepoint and block and shared by all libraries of that
##           timepoint and block. The two alleles of a hybrid are one library and the four parents were pooled into one
##           library per sample, so their errors are correlated; a shared sign keeps those correlations (with an
##           independent sign per library the simulated cis null was about twice the observed one). The column
##           sim_over_obs_null_ms_* reports the ratio of simulated to observed null mean squares as a check. The noise of
##           a parental library enters all three crosses that contain that parent, and the noise of a 0-min library
##           enters all three timepoints, as in the data. Each gene keeps its own noise variance, so noise that depends
##           on expression level is kept.
##   variant the same, with the amplitude taken from the mean parental response (amp axis of script 05) instead.
##   variant "gene": each gene keeps its own observed temporal profile. The true parental difference is the observed
##           parental signal at each timepoint and the true cis difference is rho_g times it, with rho_g the gene's
##           least-squares ratio of cis to parental difference over the three timepoints. Every gene's partition is then
##           constant in time, but genes rise and fall with their own kinetics, so the aggregate share can still change
##           if early- and late-responding genes differ in rho_g. This tests whether a change in which genes dominate
##           the response, rather than a change within genes, could produce the fall.
## Statistics and decision rule as in script 108: change in the median share over the six crosses from 15 to 60 min,
## fraction of simulations with a fall at least as large as observed (rule: below 5%), and fraction falling in all six.
suppressPackageStartupMessages({ library(data.table); library(edgeR) })
root <- Sys.getenv("PBL_ROOT"); ext <- file.path(root, "data/03_external/fay2023"); hsd <- file.path(root, "data/02_processed/hsf1")
SPECIES <- c("Sc", "Sp", "Sk", "Su"); HYB <- c("ScxSp", "ScxSk", "ScxSu", "SpxSk", "SpxSu", "SkxSu")
TP <- c(15L, 30L, 60L); TA <- c(0L, TP); PRIOR <- 2; NS <- as.integer(Sys.getenv("NSIM", "1000")); set.seed(20261013L)

## ---- library-level logCPM, built exactly as in script 05 (TMM within species / within hybrid)
dat <- read.csv(file.path(ext, "htseq.ortho.csv")); design <- as.data.table(read.csv(file.path(ext, "design.csv")))
one <- dat$one == 1; cts <- as.matrix(dat[one, -(1:14)]); rownames(cts) <- dat$HOG[one]
design[, col := Name][, block := fifelse(as.integer(sub("^[A-Za-z]+", "", Name)) <= 32L, "block1", "block2")]
P <- design[Strain %in% SPECIES & Allele == Strain & Condition == "heat"][, Time := as.integer(Time)]
H <- design[Strain %in% HYB & Condition == "heat"][, Time := as.integer(Time)]
lcpm <- function(cols) cpm(calcNormFactors(DGEList(counts = cts[, cols, drop = FALSE]), method = "TMM"), log = TRUE, prior.count = PRIOR)
LIB <- list()   # LIB[[unit]][[block]] = genes x 4 timepoints (0, 15, 30, 60); unit = species or hybrid:allele
for (sp in SPECIES) { ps <- P[Strain == sp]; L <- lcpm(ps$col)
  LIB[[sp]] <- lapply(c("block1", "block2"), function(b) sapply(TA, function(t) L[, ps[block == b & Time == t, col]])) }
for (g in HYB) { hs <- H[Strain == g]; L <- lcpm(hs$col)
  for (al in unique(hs$Allele)) LIB[[paste0(g, ":", al)]] <- lapply(c("block1", "block2"), function(b)
    sapply(TA, function(t) L[, hs[Allele == al & block == b & Time == t, col]])) }
POL <- data.table(pair = HYB, s1 = substr(HYB, 1, 2), s2 = substr(HYB, 4, 5))

## ---- gene sets as in script 108 (from the deposited decomposition table)
X <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
SETS <- list(core = X[stratum == "Hsf1_core", unique(hog)], all = X[, unique(hog)])

resp <- function(M) M[, -1, drop = FALSE] - M[, 1]              # response at 15, 30, 60 min against 0 min, same block
axes_from <- function(lib, g) {                                  # per block: cis and par (genes x 3)
  p <- POL[pair == g]
  lapply(1:2, function(b) list(
    cis = resp(lib[[paste0(g, ":", p$s2)]][[b]]) - resp(lib[[paste0(g, ":", p$s1)]][[b]]),
    par = resp(lib[[p$s2]][[b]]) - resp(lib[[p$s1]][[b]]),
    amp = (resp(lib[[p$s2]][[b]]) + resp(lib[[p$s1]][[b]])) / 2))
}
ampl <- function(sig, nul) sqrt(pmax(0, colMeans(sig^2) - colMeans(nul^2)))   # per timepoint
share_cells <- function(sig_c, nul_c, sig_t, nul_t) { a <- ampl(sig_c, nul_c)^2; b <- ampl(sig_t, nul_t)^2; a / (a + b) }

## ---- check: the library-level rebuild reproduces the deposited decomposition exactly
for (g in HYB) { A <- axes_from(lapply(LIB, function(u) lapply(u, function(m) m[SETS$all, , drop = FALSE])), g)
  for (k in seq_along(TP)) { x <- X[pair == g & time == TP[k]][match(SETS$all, hog)]
    stopifnot(max(abs((A[[1]]$cis[, k] + A[[2]]$cis[, k]) / 2 - x$cis)) < 1e-8,
              max(abs((A[[1]]$par[, k] + A[[2]]$par[, k]) / 2 - x$par)) < 1e-8) } }
cat("library-level rebuild reproduces cis and par of the deposited table exactly\n")

## a share is undefined (0/0) when both amplitudes of a cross and timepoint are at the noise floor; such a cross is left
## out of that timepoint's median and the simulation cannot count as falling in all six (floor = 1 records the event)
summ <- function(SH) { m15 <- median(SH[, 1], na.rm = TRUE); m60 <- median(SH[, 3], na.rm = TRUE)
  c(med15 = m15, med60 = m60, drop = m15 - m60, rel = (m15 - m60) / m15,
    all6 = !anyNA(SH[, c(1, 3)]) && all(SH[, 3] < SH[, 1]), floor = anyNA(SH)) }

run_set <- function(set, scale_axis) {
  lib <- lapply(LIB, function(u) lapply(u, function(m) m[SETS[[set]], , drop = FALSE]))
  ## per-gene, per-library noise: the block difference of that library, minus its mean over the four timepoints of the
  ## same unit (a constant block offset cancels in every response t - 0 and so is not noise of any contrast), rescaled
  ## by sqrt(4/3) for the centring, and by 1/sqrt(2) so that one block's draw has the variance of one library
  NOI <- lapply(lib, function(u) { d <- u[[2]] - u[[1]]; (d - rowMeans(d)) * sqrt(4 / 3) / sqrt(2) })
  TR <- list(); OBS <- matrix(NA, 6, 3); MS <- list()
  for (i in seq_along(HYB)) { g <- HYB[i]; A <- axes_from(lib, g)
    sig <- lapply(c("cis", "par", "amp"), function(a) (A[[1]][[a]] + A[[2]][[a]]) / 2)
    nul <- lapply(c("cis", "par", "amp"), function(a) (A[[2]][[a]] - A[[1]][[a]]) / 2); names(sig) <- names(nul) <- c("cis", "par", "amp")
    OBS[i, ] <- share_cells(sig$cis, nul$cis, sig$par - sig$cis, nul$par - nul$cis)
    if (scale_axis == "gene") {
      rho <- rowSums(sig$cis * sig$par) / pmax(rowSums(sig$par^2), 1e-12)
      TR[[g]] <- list(cis = sig$par * rho, par = sig$par)
    } else {
      k <- ampl(sig[[scale_axis]], nul[[scale_axis]])             # observed amplitude at 15, 30, 60 min
      pat_c <- rowMeans(sweep(sig$cis, 2, k, "/")); pat_p <- rowMeans(sweep(sig$par, 2, k, "/"))
      TR[[g]] <- list(cis = outer(pat_c, k), par = outer(pat_p, k)) }
    MS[[g]] <- c(cis = colMeans(nul$cis^2), par = colMeans(nul$par^2)) }
  o <- summ(OBS)
  CAL <- list()
  sims <- t(replicate(NS, {
    SG <- lapply(1:2, function(b) matrix(sample(c(-1, 1), length(NOI[[1]]), replace = TRUE), nrow(NOI[[1]])))
    sim_lib <- lapply(NOI, function(d) lapply(1:2, function(b) d * SG[[b]]))
    SH <- matrix(NA, 6, 3)
    for (i in seq_along(HYB)) { g <- HYB[i]; E <- axes_from(sim_lib, g)
      cb <- lapply(1:2, function(b) TR[[g]]$cis + E[[b]]$cis); pb <- lapply(1:2, function(b) TR[[g]]$par + E[[b]]$par)
      sc <- (cb[[1]] + cb[[2]]) / 2; nc <- (cb[[2]] - cb[[1]]) / 2; sp <- (pb[[1]] + pb[[2]]) / 2; np <- (pb[[2]] - pb[[1]]) / 2
      SH[i, ] <- share_cells(sc, nc, sp - sc, np - nc)   # NaN only if both amplitudes are at the floor
      if (length(CAL) < 6L * 20L) CAL[[length(CAL) + 1L]] <<- data.table(pair = g, time = TP, cis = colMeans(nc^2) / MS[[g]][1:3], par = colMeans(np^2) / MS[[g]][4:6]) }
    summ(SH) }))
  cal <- rbindlist(CAL)[, .(sim_over_obs_null_ms_cis = median(cis), sim_over_obs_null_ms_par = median(par))]
  cbind(data.table(set = set, amplitude_from = scale_axis, obs_med15 = o[["med15"]], obs_med60 = o[["med60"]], obs_drop = o[["drop"]],
                   obs_all6 = as.logical(o[["all6"]]), null_med15 = median(sims[, "med15"]), null_med60 = median(sims[, "med60"]),
                   null_drop_median = median(sims[, "drop"]), null_drop_q95 = quantile(sims[, "drop"], .95, names = FALSE),
                   n_ge_obs = sum(sims[, "drop"] >= o[["drop"]]), p_drop_ge_obs = mean(sims[, "drop"] >= o[["drop"]]),
                   p_all6_fall = mean(sims[, "all6"] == 1), n_sims_with_floor_cell = sum(sims[, "floor"] == 1),
                   ## the same comparison on the proportional scale, fall / share at 15 min, because the simulated share can
                   ## start at a different level from the observed one
                   obs_rel_drop = o[["rel"]], null_rel_drop_median = median(sims[, "rel"]),
                   null_rel_drop_q95 = quantile(sims[, "rel"], .95, names = FALSE),
                   n_rel_ge_obs = sum(sims[, "rel"] >= o[["rel"]])), cal)
}
OUT <- rbindlist(lapply(c("par", "amp", "gene"), function(a) rbindlist(lapply(c("core", "all"), run_set, scale_axis = a))))
fwrite(OUT, file.path(hsd, "p3_constant_share_null.tsv"), sep = "\t"); print(OUT)
