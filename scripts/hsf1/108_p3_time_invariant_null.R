#!/usr/bin/env Rscript
## 108_p3_time_invariant_null.R -- time-invariant null for the 15 -> 60 min fall of the P3 cis share. Each gene's cis
## and trans signal is held at its mean over time, noise is redrawn from the observed null values of each cross x time
## (1000 simulations), and the observed drop in the median share and the fall in all six crosses are compared with the
## null, for the Hsf1 core and all genes. The null was fixed before this analysis was run.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT"); hsd <- file.path(root, "data/02_processed/hsf1"); set.seed(20261004L); NS <- 1000L
X <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
share_of <- function(d) { r <- function(s, n) sqrt(max(0, mean(s^2) - mean(n^2)))
  d[, .(cis = r(cis, cis_nul), trn = r(trans, trans_nul)), by = .(pair, time)][, share := cis^2 / (cis^2 + trn^2)] }
ST <- list(core = X[stratum == "Hsf1_core"], all = unique(X[, !c("stratum"), with = FALSE], by = c("hog", "pair", "time")))
OBS <- lapply(ST, share_of)
## observed floor hits (signal ms below null ms) by cross x time, per component
FL <- rbindlist(lapply(names(ST), function(s) ST[[s]][, .(set = s, cis_floor = mean(cis^2) <= mean(cis_nul^2), trn_floor = mean(trans^2) <= mean(trans_nul^2)), by = .(pair, time)]))
sim1 <- function(d) {
  tr <- d[, .(tcis = mean(cis), ttrn = mean(trans)), by = .(hog, pair)]
  s <- merge(d[, .(hog, pair, time, cis_nul, trans_nul)], tr, by = c("hog", "pair"))
  ## the half-sum noise and the half-difference null are independent with equal variance (script 05, line 39), so each
  ## is drawn independently from the observed null values of that cross x time (keeps its depth/dispersion structure)
  s[, `:=`(ec1 = sample(cis_nul), ec2 = sample(cis_nul), et1 = sample(trans_nul), et2 = sample(trans_nul)), by = .(pair, time)]
  s[, `:=`(cis = tcis + ec1, cis_nul = ec2, trans = ttrn + et1, trans_nul = et2)]
  share_of(s) }
summ <- function(sh) { w <- dcast(sh, pair ~ time, value.var = "share")
  c(med15 = median(w[["15"]]), med60 = median(w[["60"]]), drop = median(w[["15"]]) - median(w[["60"]]), all6 = all(w[["60"]] < w[["15"]])) }
OUT <- rbindlist(lapply(names(ST), function(s) {
  o <- summ(OBS[[s]]); S <- t(replicate(NS, summ(sim1(ST[[s]]))))
  data.table(set = s, obs_med15 = o["med15"], obs_med60 = o["med60"], obs_drop = o["drop"], obs_all6 = as.logical(o["all6"]),
             null_drop_median = median(S[, "drop"]), null_drop_q95 = quantile(S[, "drop"], .95),
             p_drop_ge_obs = mean(S[, "drop"] >= o["drop"]), p_all6_fall = mean(S[, "all6"] == 1)) }))
fwrite(OUT, file.path(hsd, "p3_time_invariant_null.tsv"), sep = "\t"); fwrite(FL, file.path(hsd, "p3_noise_floor_hits.tsv"), sep = "\t")
print(OUT); print(FL[, .(cis_floor = sum(cis_floor), trn_floor = sum(trn_floor), n = .N), by = .(set, time)])
