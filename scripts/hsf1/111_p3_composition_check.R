#!/usr/bin/env Rscript
## 111_p3_composition_check.R -- composition check of the P3 cis-share fall, with the pass criterion
## fixed before this analysis was run.
## Does the 15 -> 60 min fall of the cis share survive removing a shift common to all genes in each contrast?
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT"); hsd <- file.path(root, "data/02_processed/hsf1")
NS <- 1000L
X <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
stopifnot(max(abs(X$par - (X$cis + X$trans))) < 1e-9)

## the share and the time-invariant null, exactly as script 108
share_of <- function(d) { r <- function(s, n) sqrt(max(0, mean(s^2) - mean(n^2)))
  d[, .(cis = r(cis, cis_nul), trn = r(trans, trans_nul)), by = .(pair, time)][, share := cis^2 / (cis^2 + trn^2)] }
sim1 <- function(d) {
  tr <- d[, .(tcis = mean(cis), ttrn = mean(trans)), by = .(hog, pair)]
  s <- merge(d[, .(hog, pair, time, cis_nul, trans_nul)], tr, by = c("hog", "pair"))
  s[, `:=`(ec1 = sample(cis_nul), ec2 = sample(cis_nul), et1 = sample(trans_nul), et2 = sample(trans_nul)), by = .(pair, time)]
  s[, `:=`(cis = tcis + ec1, cis_nul = ec2, trans = ttrn + et1, trans_nul = et2)]
  share_of(s) }
summ <- function(sh) { w <- dcast(sh, pair ~ time, value.var = "share")
  c(med15 = median(w[["15"]]), med60 = median(w[["60"]]), drop = median(w[["15"]]) - median(w[["60"]]), all6 = all(w[["60"]] < w[["15"]])) }

ALL <- unique(X[, !c("stratum"), with = FALSE], by = c("hog", "pair", "time"))
REFG <- unique(X[stratum == "not_induced", .(hog, pair, time)])

## descriptive medians
DESC <- rbind(ALL[, .(genes = "all", med_P = median(par), med_C = median(cis), med_T = median(trans)), by = .(pair, time)],
              merge(ALL, REFG, by = c("hog", "pair", "time"))[, .(genes = "not_induced", med_P = median(par), med_C = median(cis),
                                                                     med_T = median(trans)), by = .(pair, time)])

## centring: offsets per pair x time from all genes (V1) or not-induced genes (V2)
centre <- function(d, ref) {
  off <- ref[, .(oP = median(par), oC = median(cis), oPn = median(cis_nul + trans_nul), oCn = median(cis_nul)), by = .(pair, time)]
  d <- merge(d, off, by = c("pair", "time"))
  d[, `:=`(P2 = par - oP, C2 = cis - oC, Pn2 = (cis_nul + trans_nul) - oPn, Cn2 = cis_nul - oCn)]
  d[, `:=`(cis = C2, trans = P2 - C2, cis_nul = Cn2, trans_nul = Pn2 - Cn2)]
  d[, c("oP", "oC", "oPn", "oCn", "P2", "C2", "Pn2", "Cn2") := NULL]; d }
REF <- list(V0_uncentred = NULL, V1_genome = ALL, V2_not_induced = merge(ALL, REFG, by = c("hog", "pair", "time")))

OUT <- list()
for (v in names(REF)) {
  set.seed(20261004L)
  base_all <- if (is.null(REF[[v]])) ALL else centre(copy(ALL), REF[[v]])
  base_core <- if (is.null(REF[[v]])) X[stratum == "Hsf1_core"] else centre(copy(X[stratum == "Hsf1_core"]), REF[[v]])
  for (s in c("core", "all")) {
    d <- if (s == "core") base_core else base_all
    o <- summ(share_of(d)); S <- t(replicate(NS, summ(sim1(d))))
    OUT[[length(OUT) + 1]] <- data.table(variant = v, set = s, obs_med15 = o["med15"], obs_med60 = o["med60"], obs_drop = o["drop"],
      obs_all6 = as.logical(o["all6"]), null_drop_q95 = quantile(S[, "drop"], .95), p_drop_ge_obs = mean(S[, "drop"] >= o["drop"]))
  }
}
OUT <- rbindlist(OUT)
## the uncentred variant must reproduce the saved script-108 values (p3_time_invariant_null.tsv)
L <- fread(file.path(hsd, "p3_time_invariant_null.tsv"))
chk <- merge(OUT[variant == "V0_uncentred"], L[, .(set, l15 = obs_med15, l60 = obs_med60)], by = "set")
stopifnot(max(abs(chk$obs_med15 - chk$l15)) < 1e-9, max(abs(chk$obs_med60 - chk$l60)) < 1e-9)

ref_drop <- OUT[variant == "V0_uncentred", setNames(obs_drop, set)]
OUT[, pass := obs_drop >= 0.5 * ref_drop[set] & p_drop_ge_obs <= 0.01 & (set != "all" | obs_all6)]
VER <- OUT[variant != "V0_uncentred", .(variant_pass = all(pass)), by = variant]
n_pass <- sum(VER$variant_pass)
verdict <- if (n_pass == 2) "fall survives" else if (n_pass == 1) "partly" else "fall not robust to composition"
fwrite(DESC, file.path(hsd, "p3_composition_medians.tsv"), sep = "\t")
fwrite(OUT, file.path(hsd, "p3_composition_check.tsv"), sep = "\t")
print(DESC[, .(genes, pair, time, med_P = round(med_P, 3), med_C = round(med_C, 3), med_T = round(med_T, 3))])
print(OUT); print(VER); cat("VERDICT:", verdict, "\n")
