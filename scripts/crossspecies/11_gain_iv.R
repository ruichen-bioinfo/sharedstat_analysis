# 11_gain_iv.R — an artefact-free estimate of the allelic response gain, with an evaluation of
# the two analyses that the artefact affects.
#
# Why the OLS slope is biased. Regressing each cross's interaction effect on its own shared
# response (script 10) gives a negative slope in 36 of 36 cells, -0.106 to -0.325, with
# no exception. That cannot be a species-level gain: S. paradoxus is the coefficient allele
# in ScxSp and SpxSk but the reference allele in SpxSu, so a consistent species ordering
# cannot make all six negative. The six slopes are, however, exactly consistent with the
# single ordering Sc > Sk > Sp > Su, which is the reverse of alphabetical order, i.e. with
# "the reference allele always responds more". Alphabetical order carries no biology, so
# this is a reference-level artefact.
#
# The mechanism is in what script 02 stores (eff_col()): for `time` the effect is the
# Time15 main-effect coefficient, which is the reference allele's 15 min response R; for
# `inter` it is the Time15:Allele coefficient at the same timepoint, which is exactly
# C - R. The two share R's measurement error with opposite sign, so with a true gain
# gamma (C_signal = gamma * R_signal) and reliability h = var(signal)/var(total),
#
#     cov(R, C - R) = cov(R,C) - var(R) = (gamma - 1) var(sig_R) - var(noise_R)
#     slope_OLS     = (gamma - 1) h - (1 - h)
#
# At gamma = 1 this is -(1 - h). The observed -0.11 to -0.33 imply h = 0.67 to 0.89, which
# is an unremarkable reliability for a difference-of-differences at n = 2. So the observed
# slopes are fully explained with no gain at all, and the OLS slopes of script 10 section 1
# do not estimate a gain.
#
# Instrumental-variable estimate. The bias comes from one specific thing, R's measurement
# error appearing on both axes. An instrument correlated with R's signal but not its noise
# removes it. The other crosses' shared responses are exactly that: the same genes'
# heat-response programme measured in physically separate libraries. With
# S = mean shared response over the other crosses,
#
#     gamma_hat = 1 + cov(eff_inter, S) / cov(R, S)
#
# is consistent for the true gain, and the gap between slope_OLS and this IV slope is a
# direct measurement of the artefact rather than an assumption about it. Co-pooled partners
# are excluded from S, because those libraries are shared and would reintroduce a common
# noise path (script 09 found no leakage, median excess -0.015 over 36 cells, but the
# instrument should not depend on that result holding).
#
# Bland-Altman gives a second, independent route: regressing C - R on the average response
# (C + R)/2 = eff_time + eff_inter/2 cancels the shared error to first order, because
# cov(C - R, (C + R)/2) = [var(C) - var(R)]/2, which is zero under equal allelic noise.
# It needs no instrument but assumes the two alleles are equally noisy, so it is reported
# alongside IV rather than instead of it. Agreement between two estimators with different
# failure modes is the point.
#
# Analyses evaluated:
#   section 3  the thermal-niche test on the IV gain rather than the OLS slope
#   section 4  the additive four-species model of script 09, which is immune to first
#              order because the artefact enters as a term that is near-constant across
#              the six crosses and the centred design cannot fit a constant (Dc %*% 1 = 0,
#              so a constant goes entirely into the residual and only depresses R2). That
#              is an analytic argument, so it is checked here by injection instead.
#   section 5  the polarity test of script 08, whose responder mask was taken from the
#              same cross's `time` term and therefore selected on a variable sharing noise
#              with the quantity being tested. Re-run with the mask taken from the other
#              crosses, and split by induced versus repressed, since a slope on eff_time
#              and a sign fraction pooled over both directions are not the same statistic.
#
# Usage: Rscript scripts/crossspecies/11_gain_iv.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
outd <- file.path(root, "data/02_processed/crossspecies")
HYB <- c("ScxSk","ScxSp","ScxSu","SkxSu","SpxSk","SpxSu")
METHODS <- c("DESeq2","edgeR","limma")
SPECIES <- c("Sc","Sk","Sp","Su")
ALPHA <- 0.05
partner <- c(ScxSp = "SkxSu", SkxSu = "ScxSp", ScxSu = "SpxSk",
             SpxSk = "ScxSu", ScxSk = "SpxSu", SpxSu = "ScxSk")
CLS <- data.table(sp = SPECIES,
                  thermal = c("thermo","cryo","thermo","cryo"),
                  source  = c("domesticated","wild","wild","domesticated"))
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(20260905)
SELECT <- Sys.getenv("GAIN_IV_SELECT", "none"); stopifnot(SELECT %in% c("none", "padj"))
SUFFIX <- if (SELECT == "padj") "_padj_filter" else ""
say("=== 11_gain_iv | %s | gene selection: %s ===", format(Sys.time()), SELECT)

X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod.tsv.gz"))))
X <- X[condition == "heat" & group %in% HYB]

## R = reference allele response, per cross x method (model A only: `time` is not
## estimable once library is blocked, so model B borrows model A's shared response)
Rtab <- X[model == "A" & term == "time", .(gene, group, method, R = eff, padj_t = padj)]

## instrument: mean shared response over the crosses that are neither the focal cross nor
## its co-pooling partner. Four crosses contribute, so the instrument is itself less noisy
## than any single cross's estimate.
instrument <- function(g, mt) {
  keep <- setdiff(HYB, c(g, partner[[g]]))
  s <- Rtab[group %in% keep & method == mt, .(S = mean(R, na.rm = TRUE),
                                              nS = sum(is.finite(R))), by = gene]
  s[nS >= 3]
}

## ------------------------------------------------------- (1) OLS vs IV per cross ----
say("\n=== (1) the artefact, measured: OLS slope vs IV slope ===")
say("OLS regresses eff_inter on the SAME cross's response, so it carries the shared-error")
say("bias. IV uses the other crosses' response as the instrument and does not.")
say("gain = 1 + slope_IV, the estimand of P1 section 1.")
say("%-8s %-7s %-2s %7s %9s %9s %8s %9s", "cross","method","M","n","slope_OLS","slope_IV",
    "gain_IV","artefact")
resG <- list()
for (g in HYB) for (mt in METHODS) for (mdl in c("A","B")) {
  ii <- X[group == g & method == mt & model == mdl & term == "inter", .(gene, ei = eff)]
  D <- merge(merge(ii, Rtab[group == g & method == mt, .(gene, R, padj_t)], by = "gene"),
             instrument(g, mt), by = "gene")
  ## The primary reading is (i), no significance filter, as fixed before this analysis was run.
  ## Selecting on the regressor (padj of R itself) correlates its
  ## error with the instrument. SELECT=padj gives the padj-filtered reading, written to
  ## gain_iv_padj_filter.tsv for Supplementary Table S5.
  D <- D[is.finite(ei) & is.finite(R) & is.finite(S)]
  if (SELECT == "padj") D <- D[!is.na(padj_t) & padj_t < ALPHA]
  if (nrow(D) < 200) next
  ols <- unname(coef(lm(ei ~ R, D))[2])
  ## IV: cov(y, S) / cov(x, S). Equivalent to two-stage least squares with one instrument.
  iv  <- cov(D$ei, D$S) / cov(D$R, D$S)
  ## bootstrap the IV slope over genes, so the CI reflects gene-level sampling
  bs <- replicate(200, { k <- sample(nrow(D), replace = TRUE)
                         cov(D$ei[k], D$S[k]) / cov(D$R[k], D$S[k]) })
  resG[[length(resG)+1]] <- data.table(group = g, method = mt, model = mdl, n = nrow(D),
    slope_ols = ols, slope_iv = iv, gain_iv = 1 + iv, artefact = ols - iv,   # the estimand is 1 + slope (P1 sec. 1), not 2^slope
    iv_lo = unname(quantile(bs, 0.025)), iv_hi = unname(quantile(bs, 0.975)),
    ba = unname(coef(lm(ei ~ I(R + ei/2), D))[2]))
  say("%-8s %-7s %-2s %7d %+9.4f %+9.4f %8.3f %+9.4f", g, mt, mdl, nrow(D), ols, iv,
      1 + iv, ols - iv)
}
G <- rbindlist(resG)
fwrite(G, file.path(outd, paste0("gain_iv", SUFFIX, ".tsv")), sep = "\t")
say("\nOLS slope negative in %d of %d cells; IV slope negative in %d of %d",
    sum(G$slope_ols < 0), nrow(G), sum(G$slope_iv < 0), nrow(G))
say("artefact (OLS - IV) is negative in %d of %d cells, median %+.4f",
    sum(G$artefact < 0), nrow(G), median(G$artefact))
say("implied reliability h = 1 + artefact ... median %.3f, consistent with n=2 replicates",
    1 + median(G$artefact))

## -------------------------------------------------- (2) the two components ranked ---
say("\n=== (2) coherent gain and gene-specific rewiring rank the crosses differently ===")
NS <- X[model == "A" & term == "inter" & method == "DESeq2",
        .(n_sig = sum(!is.na(padj) & padj < ALPHA)), by = group]
s <- merge(G[model == "A", .(iv = median(slope_iv), ols = median(slope_ols),
                             ba = median(ba)), by = group], NS, by = "group")
s[, `:=`(rank_iv = rank(-abs(iv)), rank_ols = rank(-abs(ols)), rank_sig = rank(-n_sig))]
setorder(s, rank_iv)
say("%-8s %9s %9s %9s %8s %8s %8s", "cross","slope_IV","slope_BA","slope_OLS",
    "rk_gain","n_sig","rk_sig")
for (i in seq_len(nrow(s))) say("%-8s %+9.4f %+9.4f %+9.4f %8.0f %8d %8.0f",
  s$group[i], s$iv[i], s$ba[i], s$ols[i], s$rank_iv[i], s$n_sig[i], s$rank_sig[i])
say("\nSpearman(gain rank, rewiring rank): IV %+.3f | OLS %+.3f",
    cor(s$rank_iv, s$rank_sig, method = "spearman"),
    cor(s$rank_ols, s$rank_sig, method = "spearman"))
say("Spearman(IV rank, OLS rank) = %+.3f -- the artefact reorders the crosses",
    cor(s$rank_iv, s$rank_ols, method = "spearman"))
fwrite(s, file.path(outd, paste0("gain_components", SUFFIX, ".tsv")), sep = "\t")

## ------------------------------------- (3) does the IV gain track thermal niche? ----
say("\n=== (3) thermal niche and isolation source, tested on the IV gain ===")
orient <- function(g, axis) {
  a <- sort(unlist(strsplit(g, "x")))
  k1 <- CLS[sp == a[1]][[axis]]; k2 <- CLS[sp == a[2]][[axis]]
  if (k1 == k2) return(NA_integer_)
  pos <- if (axis == "thermal") "cryo" else "wild"
  if (k2 == pos) 1L else -1L
}
for (axis in c("thermal","source")) {
  say("\n-- axis: %s (positive = %s allele driven harder) --", axis,
      if (axis == "thermal") "cryophilic" else "wild")
  say("%-8s %-30s %9s %9s %-22s", "cross","classes","slope_IV","95% CI","oriented")
  vals <- c()
  for (g in HYB) {
    a <- sort(unlist(strsplit(g, "x")))
    fl <- orient(g, axis)
    r <- G[group == g & model == "A"]
    iv <- median(r$slope_iv)
    say("%-8s %-30s %+9.4f %9s %-22s", g,
        sprintf("%s=%s, %s=%s", a[1], CLS[sp==a[1]][[axis]], a[2], CLS[sp==a[2]][[axis]]),
        iv, sprintf("%+.3f,%+.3f", median(r$iv_lo), median(r$iv_hi)),
        if (is.na(fl)) "matched (null channel)" else sprintf("%+.4f", fl*iv))
    if (!is.na(fl)) vals[g] <- fl*iv
  }
  say("  informative crosses: %d | %d of %d in the predicted direction | median %+.4f",
      length(vals), sum(vals > 0), length(vals), median(vals))
  say("  %s", if (all(vals > 0)) "consistent with this axis"
              else if (all(vals < 0)) "consistent in magnitude but with the opposite sign"
              else "NOT consistent: the sign flips across crosses")
}
say("\nlargest |IV gain| overall: %s", G[model == "A"][which.max(abs(slope_iv)),
    sprintf("%s (slope %+.4f, gain %.3f) -- a %s x %s pair",
            group, slope_iv, gain_iv,
            CLS[sp == sort(unlist(strsplit(group,"x")))[1], thermal],
            CLS[sp == sort(unlist(strsplit(group,"x")))[2], thermal])])

## ------------------------ (4) is the additive species model immune to the artefact? --
## The artefact adds roughly -(1-h) * R to every cross's interaction effect, and R is
## nearly the same quantity in all six crosses, so it enters as a near-constant vector.
## The centred design cannot represent a constant, so the modifiers should be unmoved.
## Tested by injecting an artefact of the observed magnitude and refitting.
say("\n=== (4) additive four-species model: artefact injection test ===")
Dm <- matrix(0, 6, 4, dimnames = list(HYB, SPECIES))
for (g in HYB) { a <- sort(unlist(strsplit(g, "x"))); Dm[g, a[1]] <- -1; Dm[g, a[2]] <- 1 }
Dc <- Dm %*% (diag(4) - matrix(1/4, 4, 4))
Hm <- Dc %*% MASS::ginv(Dc); Gi <- MASS::ginv(Dc)
say("Dc %%*%% 1 = (%s) -- the design cannot fit a constant, so a constant artefact",
    paste(sprintf("%.0f", Dc %*% rep(1, 4)), collapse = ","))
say("goes entirely into the residual and can only depress R2, never bias the modifiers.")
say("%-7s %-10s %7s %8s %8s %8s %8s %8s", "method","data","n_gene","R2_med",
    SPECIES[1], SPECIES[2], SPECIES[3], SPECIES[4])
resI <- list()
for (mt in METHODS) {
  W <- dcast(X[method == mt & model == "A" & term == "inter", .(gene, group, eff)],
             gene ~ group, value.var = "eff")
  Rw <- dcast(Rtab[method == mt, .(gene, group, R)], gene ~ group, value.var = "R")
  M <- merge(W, Rw, by = "gene", suffixes = c("", ".R"))
  Y <- as.matrix(M[, ..HYB]); RR <- as.matrix(M[, paste0(HYB, ".R"), with = FALSE])
  keep <- complete.cases(Y) & complete.cases(RR) & rowSums(Y^2) > 0
  Y <- Y[keep, , drop = FALSE]; RR <- RR[keep, , drop = FALSE]
  lam <- -median(G[method == mt & model == "A", artefact])   # observed artefact size
  for (lab in c("observed","artefact removed","artefact doubled")) {
    Z <- switch(lab, observed = Y, `artefact removed` = Y + lam * RR,
                `artefact doubled` = Y - lam * RR)
    r2 <- 1 - rowSums((Z - Z %*% t(Hm))^2) / rowSums(Z^2)
    Mh <- Z %*% t(Gi); colnames(Mh) <- SPECIES
    ok <- is.finite(r2) & r2 > 0.8
    resI[[length(resI)+1]] <- data.table(method = mt, data = lab, n = sum(ok),
      r2_med = median(r2, na.rm = TRUE), Sc = median(Mh[ok,"Sc"]), Sk = median(Mh[ok,"Sk"]),
      Sp = median(Mh[ok,"Sp"]), Su = median(Mh[ok,"Su"]))
    say("%-7s %-16s %7d %8.3f %8.3f %8.3f %8.3f %8.3f", mt, lab, sum(ok),
        median(r2, na.rm = TRUE), median(Mh[ok,"Sc"]), median(Mh[ok,"Sk"]),
        median(Mh[ok,"Sp"]), median(Mh[ok,"Su"]))
  }
}
II <- rbindlist(resI)
fwrite(II, file.path(outd, paste0("species_modifier_artefact_test", SUFFIX, ".tsv")), sep = "\t")
say("\nmodifier ordering under each data version:")
for (lab in unique(II$data)) {
  s2 <- II[data == lab, lapply(.SD, median), .SDcols = SPECIES]
  o <- SPECIES[order(-unlist(s2))]
  say("  %-16s %s", lab, paste(sprintf("%s(%+.3f)", o, unlist(s2)[o]), collapse = " > "))
}

## --------------- (5) polarity, with the responder mask taken from other crosses ------
## The script 08 mask came from the same cross's `time` padj, which shares noise with the
## interaction effect, so genes entered the test partly because the reference allele's
## response was noise-inflated. Selecting on the other crosses' response removes that
## path. Induced and repressed genes are separated, because a gain multiplies the response
## and therefore flips the sign of the interaction effect between them.
say("\n=== (5) polarity with an independent responder mask, split by direction ===")
say("frac>0 = fraction favouring the cryophilic allele (matched pairs: the 2nd allele).")
say("%-8s %-22s %-7s %-9s %7s %8s %11s", "cross","axis","method","direction","n","frac>0","binom_p")
resP <- list()
for (g in HYB) for (mt in METHODS) {
  a <- sort(unlist(strsplit(g, "x")))
  fl <- orient(g, "thermal"); ax <- if (is.na(fl)) "matched (null channel)" else "cryo vs thermo"
  keep <- setdiff(HYB, c(g, partner[[g]]))
  ## independent mask: significant shared response in at least two other crosses, and its
  ## direction taken from those crosses too, so neither membership nor sign uses this cross
  om <- Rtab[group %in% keep & method == mt,
             .(nsig = sum(!is.na(padj_t) & padj_t < ALPHA), Sm = mean(R, na.rm = TRUE)),
             by = gene][nsig >= 2]
  ii <- X[group == g & method == mt & model == "A" & term == "inter", .(gene, ei = eff)]
  D <- merge(ii, om, by = "gene")[is.finite(ei) & is.finite(Sm) & ei != 0]
  for (dir in c("induced","repressed")) {
    s2 <- D[if (dir == "induced") Sm > 0 else Sm < 0]
    if (nrow(s2) < 50) next
    ## orient to the cryophilic allele; for a gain the sign of eff_inter relative to the
    ## response direction is what carries meaning, so repressed genes are sign-corrected
    e <- s2$ei * (if (is.na(fl)) 1L else fl) * (if (dir == "repressed") -1L else 1L)
    bt <- binom.test(sum(e > 0), length(e), 0.5)
    resP[[length(resP)+1]] <- data.table(group = g, axis = ax, method = mt, direction = dir,
      n = length(e), frac_pos = mean(e > 0), binom_p = bt$p.value, med = median(e))
    say("%-8s %-22s %-7s %-9s %7d %7.1f%% %11.3g", g, ax, mt, dir, length(e),
        100*mean(e > 0), bt$p.value)
  }
}
P <- rbindlist(resP)
fwrite(P, file.path(outd, paste0("polarity_independent_mask", SUFFIX, ".tsv")), sep = "\t")
say("\n--- summary: does the direction split agree? a gain predicts the same sign in both ---")
say("%-8s %-22s %10s %10s %s", "cross","axis","induced","repressed","agree")
for (g in HYB) {
  s2 <- P[group == g]
  if (!nrow(s2)) next
  fi <- median(s2[direction == "induced", frac_pos]); fr <- median(s2[direction == "repressed", frac_pos])
  say("%-8s %-22s %9.1f%% %9.1f%% %s", g, s2$axis[1], 100*fi, 100*fr,
      if (is.na(fi) || is.na(fr)) "-" else if ((fi - 0.5) * (fr - 0.5) > 0) "yes (gain-like)" else "no")
}
for (ax in unique(P$axis)) {
  s2 <- P[axis == ax]
  say("%-22s median frac>0 %.1f%% over %d cells | %d with p<0.05",
      ax, 100*median(s2$frac_pos), nrow(s2), sum(s2$binom_p < 0.05))
}

say("\nwrote gain_iv.tsv, gain_components.tsv, species_modifier_artefact_test.tsv, polarity_independent_mask.tsv")
say("=== done ===")
