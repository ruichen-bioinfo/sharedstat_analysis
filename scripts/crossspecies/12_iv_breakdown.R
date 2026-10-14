#!/usr/bin/env Rscript
## 12_iv_breakdown.R -- how much instrument contamination would change each corrected conclusion in P1 Table 1.
##
## No new model is fitted. This rebuilds exactly the vectors 11_gain_iv.R uses (same inputs, same filters, same
## instrument), asserts that it reproduces the IV slopes saved in gain_iv.tsv, and then does arithmetic on them.
##
## The algebra. The IV slope is b_iv = cov(e, S) / cov(R_hat, S), with e = (gamma - 1) R + v - u and R_hat = R + u.
## If the instrument is contaminated, cov(S, u) = d != 0 (cov(S, v) = 0 kept), write c = cov(S, R) and C = c + d, the
## observed cov(R_hat, S). Then b_iv = ((gamma - 1) c - d) / C, so the true slope is
##     b = (b_iv * C + d) / (C - d) = (b_iv + q) / (1 - q),   q = d / C,
## the share of the instrument's observed covariance with the reference that comes from the reference's error.
## The corrected sign is undone at q* = -b_iv. In correlation units, d = rho sd(S) sd(u) and C = r sd(S) sd(R_hat), so
##     rho* = q* r / sqrt(1 - h),   r = cor(R_hat, S),   h = reliability of R_hat.
## The same closed form reproduces the simulation sweep in Supplementary Figure S2B to three decimals (checked below).
## h: the reliability of the reference arm (h_nofilt in p1_anchor_diagnosis.tsv, see below); q* needs no h.
suppressPackageStartupMessages(library(data.table))
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
outd <- file.path(root, "data/02_processed/crossspecies")
say <- function(...) cat(sprintf(...), "\n", sep = "")
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu"); METHODS <- c("DESeq2","edgeR","limma"); ALPHA <- 0.05
partner <- c(ScxSp = "SkxSu", SkxSu = "ScxSp", ScxSu = "SpxSk", SpxSk = "ScxSu", ScxSk = "SpxSu", SpxSu = "ScxSk")
X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod.tsv.gz"))))[condition == "heat" & group %in% HYB]
Rtab <- X[model == "A" & term == "time", .(gene, group, method, R = eff, padj_t = padj)]
instrument <- function(g, mt) {
  keep <- setdiff(HYB, c(g, partner[[g]]))
  s <- Rtab[group %in% keep & method == mt, .(S = mean(R, na.rm = TRUE), nS = sum(is.finite(R))), by = gene]
  s[nS >= 3]
}
G <- fread(file.path(outd, "gain_iv.tsv"))[model == "A"]
res <- list()
for (g in HYB) for (mt in METHODS) {
  ii <- X[group == g & method == mt & model == "A" & term == "inter", .(gene, ei = eff)]
  D <- merge(merge(ii, Rtab[group == g & method == mt, .(gene, R, padj_t)], by = "gene"), instrument(g, mt), by = "gene")
  D <- D[is.finite(ei) & is.finite(R) & is.finite(S)]   # reading (i), no significance filter: same genes as gain_iv.tsv
  biv <- cov(D$ei, D$S) / cov(D$R, D$S)
  landed <- G[group == g & method == mt, slope_iv]
  stopifnot(length(landed) == 1, abs(biv - landed) < 1e-10)          # same vectors as the saved analysis in gain_iv.tsv
  res[[length(res) + 1]] <- data.table(group = g, method = mt, n = nrow(D), slope_iv = biv, r_RS = cor(D$R, D$S))
}
B <- rbindlist(res)
## Note: h must be the reliability of the reference arm on the genes and timepoint the IV uses (Table 1: model A,
## 15 min, reference-arm padj<0.05). p1_measured_reliability.tsv holds the reliability
## of the allele-paired (DoD) contrast, which is a different quantity.
H <- fread(file.path(root, "data/02_processed/hsf1/p1_anchor_diagnosis.tsv"))[, .(group = hybrid, h_med = h_nofilt)]   # (i) gene set
B <- merge(B, H, by = "group")
B[, q_star := -slope_iv][, rho_star := q_star * r_RS / sqrt(1 - h_med)]
fwrite(B, file.path(outd, "iv_breakdown_by_method.tsv"), sep = "\t")
S <- B[, .(slope_iv = median(slope_iv), r_RS = median(r_RS), h_med = h_med[1], q_star = median(q_star),
           rho_star = median(rho_star)), by = group]
## ranks: a common q applied to all six; the smallest |q| at which the ordering by |b| changes
rk <- function(q) rank(-abs((S$slope_iv + q) / (1 - q)), ties.method = "first")
r0 <- rk(0); qs <- seq(0, 0.5, by = 1e-4)
q_up <- qs[which(vapply(qs, function(q) any(rk(q) != r0), TRUE))[1]]
q_dn <- -qs[which(vapply(qs, function(q) any(rk(-q) != r0), TRUE))[1]]
S[, `:=`(rank_q_pos = q_up, rank_q_neg = q_dn)]
fwrite(S, file.path(outd, "iv_breakdown.tsv"), sep = "\t")
say("=== 12_iv_breakdown | %s ===", format(Sys.time()))
say("saved IV slopes reproduced in all %d hybrid x method cells", nrow(B))
say("%-7s %9s %7s %6s %8s %9s", "hybrid", "slope_iv", "r(R,S)", "h", "q*", "rho*")
for (i in seq_len(nrow(S))) with(S[i], say("%-7s %+9.4f %7.3f %6.3f %+8.4f %+9.4f", group, slope_iv, r_RS, h_med, q_star, rho_star))
say("ordering by |corrected slope| first changes at a common q of %+.4f (positive) or %+.4f (negative)", q_up, q_dn)
## positive control: the closed form against the S2B simulation sweep (gamma 1.3, h 0.75, S0 = 0.8 R + N(0, 0.36))
SW <- fread(file.path(root, "data/02_processed/hsf1/p1_robustness_rho_sweep.tsv"))
vu <- 1 / 3; cS <- 0.8; sS <- sqrt(1)
pred <- vapply(SW$rho, function(rho) { c <- sqrt(1 - rho^2) * cS; d <- rho * sS * sqrt(vu)
  ((0.3 * c - d) / (c + d) - (1.3 * 0.75 - 1)) / (0.3 - (1.3 * 0.75 - 1)) }, 0)
say("closed form vs S2B simulation, surviving share: max |difference| %.4f over %d rho values", max(abs(pred - SW$surviving)), nrow(SW))
stopifnot(max(abs(pred - SW$surviving)) < 0.01)
