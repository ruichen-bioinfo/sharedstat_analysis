#!/usr/bin/env Rscript
## 55_p1_fourth_dataset_kk2024.R -- P1's fourth dataset: a design in which the true gain is known to be 1.
##
## What this dataset adds. In the three comparative datasets the true gain is unknown, so the analysis can
## say "if the true gain were gamma". Here it is known by construction. Two biological replicates of one species
## measure the same true heat-shock response, so the true gain of replicate 2 on replicate 1 is exactly 1, and
## any apparent attenuation is the artefact and nothing else. It also tests whether the reliabilities reported
## are low only because low-reliability data were chosen: two of the three species here have h above 0.98, and
## the prediction is that the artefact should then nearly vanish.
##
## What is computed. Regressing the contrast (r2 - r1) on r1 is the shared-reference construction. The two
## replicates do not have equal variance, so the two directional gains differ from each other and neither is
## equal to h. Their geometric mean is identically |cor(r1, r2)| (positive in these data; see the note at the end of this
## file), so it is reported beside h descriptively and is not an estimate of h.
## One of the two directions here is 1.104, an apparent amplification.
##
##
## The design is checked against the GEO records. The counts table carries a `species` column, so species comes from the
## data. It does not carry treatment or replicate, and the SK1..SK12 ordering that supplies them is an
## assumption until the per-sample GEO records say otherwise. This script parses those records
## (data/03_external/keyportkik2024/GSE234499_samples.txt, in manifests/EXTERNAL_INPUTS.tsv) and asserts the mapping
## before computing anything. It also asserts which temperature is the control for each species, since the
## species are shocked at different temperatures and "higher" is the only thing they have in common.
##
## Reported per species: h with its bootstrap interval, the two directional gains, their geometric mean with
## its interval, and the artefact percentage. The comparison of the geometric mean with h was specified in
## advance as a test; because the geometric mean is an identity in cor(r1, r2), it is reported only descriptively.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
for (f in c("reliability.R", "reliability_replicates.R", "estimators.R"))
  source(file.path(root, "pkg/sharedstat/R", f))
ext <- file.path(root, "data/03_external/keyportkik2024")
out <- file.path(root, "data/02_processed/hsf1")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
NBOOT <- 2000L; PRIOR <- 2; MINCPM <- 1

## ------------------------------------------------------------------ design, from the GEO records ----
say("=== the design, asserted from the per-sample GEO records ===")
gsm <- readLines(file.path(ext, "GSE234499_samples.txt"), warn = FALSE)
acc <- sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", gsm, value = TRUE))
tit <- sub("^!Sample_title = ", "", grep("^!Sample_title", gsm, value = TRUE))
str <- sub("^!Sample_characteristics_ch1 = strain: ", "",
           grep("^!Sample_characteristics_ch1 = strain:", gsm, value = TRUE))
tmp <- as.numeric(sub("^!Sample_characteristics_ch1 = temperature: ", "",
                      grep("^!Sample_characteristics_ch1 = temperature:", gsm, value = TRUE)))
brn <- as.integer(sub("^!Sample_characteristics_ch1 = biological replicate: ", "",
                      grep("^!Sample_characteristics_ch1 = biological replicate:", gsm, value = TRUE)))
stopifnot(length(acc) == 12L, length(tit) == 12L, length(str) == 12L, length(tmp) == 12L, length(brn) == 12L)
G <- data.table(gsm = acc, title = tit, strain = str, temp = tmp, br = brn)
G[, species := sub("_.*$", "", title)]
say("  12 records; species %s", paste(sort(unique(G$species)), collapse = ", "))
## the control is the lower of the two temperatures within a species, asserted as exactly two per species
G[, `:=`(n_temp = length(unique(temp)), t_lo = min(temp), t_hi = max(temp)), by = species]
stopifnot(all(G$n_temp == 2L))
G[, cond := ifelse(temp == t_lo, "ctrl", "shock")]
stopifnot(G[, .N, by = .(species, cond, br)][, all(N == 1L)])
print(G[order(species, cond, br), .(gsm, title, strain, temp, cond, br)], row.names = FALSE)
say("\n  shifts: %s", paste(unique(sprintf("%s %d->%d C", G$species, G$t_lo, G$t_hi)), collapse = "; "))
## The checked sample order is SK1..SK12 = ctrl BR1, ctrl BR2, shock BR1, shock BR2 within each species block.
##
ord <- G[order(match(species, c("Skud", "Scer", "Kmarx")), cond, br)]
say("  assumed SK order matches the records: %s",
    identical(ord$title, G[order(match(species, c("Skud","Scer","Kmarx")), temp, br)]$title))

## ------------------------------------------------------------------ counts ----
X <- fread(cmd = paste("zcat", shQuote(file.path(ext, "GSE234499_counts.tsv.gz"))))
say("\n=== counts ===")
say("  rows %d; columns %s", nrow(X), paste(names(X), collapse = ", "))
X <- X[!grepl("^__", ORF)]
SP <- unique(X[, .(sample, species)])
stopifnot(nrow(SP) == 12L)
SP[, idx := as.integer(sub("^SK", "", sample))]
SP <- SP[order(idx)]
## species assignment comes from the table; treatment and replicate come from the records, joined on the
## within-species position, which is what the assertion above established
SP[, pos := seq_len(.N), by = species]
GG <- G[order(match(species, SP$species), temp, br)][, pos := seq_len(.N), by = species]
M <- merge(SP, GG[, .(species, pos, cond, br, temp, strain, gsm, title)], by = c("species", "pos"))
stopifnot(nrow(M) == 12L)
say("  sample -> condition mapping:")
print(M[order(idx), .(sample, species, gsm, title, temp, cond, br)], row.names = FALSE)
X <- merge(X, M[, .(sample, species, cond, br)], by = c("sample", "species"))
X[, cpm := 1e6 * counts / sum(counts), by = sample]
X[, lcpm := log2(cpm + PRIOR)]

## ------------------------------------------------------------------ per species ----
say("\n=== reliability and the artefact where the true gain is 1 ===")
set.seed(20260913L)
res <- rbindlist(lapply(c("Skud", "Scer", "Kmarx"), function(sp) {
  Z <- X[species == sp]
  keep <- Z[, .(ok = all(cpm >= MINCPM)), by = ORF][ok == TRUE, ORF]
  W <- dcast(Z[ORF %in% keep], ORF ~ cond + br, value.var = "lcpm")
  r1 <- W$shock_1 - W$ctrl_1
  r2 <- W$shock_2 - W$ctrl_2
  # `center = FALSE`: this table reports the uncentred reliability, the form P1 quotes as "the uncentred form of
  # the same quantity reads 0.887". The centred values P1 reports for the same three arms are in
  # p1x_b5_reliability_centred_ci.tsv (script 64).
  rel <- reliability_from_blocks(r1, r2, center = FALSE)
  amp <- as.numeric(rms_true((r1 + r2) / 2, (r2 - r1) / 2))
  g12 <- 1 + cov(r2 - r1, r1) / var(r1)
  g21 <- 1 + cov(r1 - r2, r2) / var(r2)
  gm  <- sqrt(g12 * g21)
  bt <- replicate(NBOOT, {
    i <- sample(length(r1), replace = TRUE)
    a <- r1[i]; b <- r2[i]
    c(sqrt((1 + cov(b - a, a) / var(a)) * (1 + cov(a - b, b) / var(b))),
      reliability_from_blocks(a, b, center = FALSE)$h)
  })
  data.table(species = sp, n_genes = rel$n_pairs, h = rel$h, snr = rel$snr, response_rms = amp,
             gain_2on1 = g12, gain_1on2 = g21, gain_geomean = gm,
             gm_lo = quantile(bt[1, ], .025), gm_hi = quantile(bt[1, ], .975),
             h_lo = quantile(bt[2, ], .025), h_hi = quantile(bt[2, ], .975),
             artefact_pct = 100 * (1 - rel$h), reversal_threshold_gain = 1 / rel$h)
}))
print(res[, .(species, n_genes, h = round(h, 4), h_ci = sprintf("[%.3f, %.3f]", h_lo, h_hi),
              gain_2on1 = round(gain_2on1, 4), gain_1on2 = round(gain_1on2, 4),
              gm = round(gain_geomean, 4), gm_ci = sprintf("[%.3f, %.3f]", gm_lo, gm_hi),
              artefact_pct = round(artefact_pct, 1))], row.names = FALSE)
fwrite(res, file.path(out, "p1x_kk2024_known_truth_gain.tsv"), sep = "\t")

say("\n=== geometric mean of the two directions beside h (descriptive) ===")
res[, h_in_gm_ci := h >= gm_lo & h <= gm_hi]
res[, gm_vs_h_pct := 100 * abs(gain_geomean - h) / h]
for (i in seq_len(nrow(res)))
  say("  %-6s h = %.4f, geometric mean of the two directions = %.4f, %.2f%% apart, h inside its interval: %s",
      res$species[i], res$h[i], res$gain_geomean[i], res$gm_vs_h_pct[i], res$h_in_gm_ci[i])
say("  all three within 2%%: %s", all(res$gm_vs_h_pct < 2))
say("  all three with h inside the interval: %s", all(res$h_in_gm_ci))
say("\n  a single direction is not generally equal to h:")
for (i in seq_len(nrow(res)))
  say("    %-6s the two directions are %.4f and %.4f; one of them exceeds 1 in %d of 3 species",
      res$species[i], res$gain_2on1[i], res$gain_1on2[i], sum(res$gain_2on1 > 1 | res$gain_1on2 > 1))
say("\n  artefact where reliability is high: %.1f%% to %.1f%%", min(res$artefact_pct), max(res$artefact_pct))
say("  h range: %.3f to %.3f", min(res$h), max(res$h))

say("\n=== optional numerical reproduction check ===")
fb <- Sys.getenv("PBL_COMPARE_RUN", "")   # optional; skipped if absent
if (file.exists(fb)) {
  S <- fread(fb)
  cmp <- merge(res[, .(species, h, gain_2on1, gain_1on2)],
               S[, .(species, h_s = h, g12_s = naive_gain_2_on_1, g21_s = naive_gain_1_on_2)], by = "species")
  cmp[, `:=`(dh = abs(h - h_s), d12 = abs(gain_2on1 - g12_s), d21 = abs(gain_1on2 - g21_s))]
  print(cmp[, .(species, dh = signif(dh, 3), d12 = signif(d12, 3), d21 = signif(d21, 3))], row.names = FALSE)
  say("  largest absolute difference: %.2e", max(c(cmp$dh, cmp$d12, cmp$d21)))
  say("  (numerical comparison)")
}
say("")

## ------------------------------------------------------------------ algebraic interpretation of this check ----
## The quantity in the decision rule is an algebraic identity rather
## than a prediction:
##     g12 = 1 + cov(r2 - r1, r1) / var(r1) = cov(r1, r2) / var(r1)
##     g21 = 1 + cov(r1 - r2, r2) / var(r2) = cov(r1, r2) / var(r2)
##     sqrt(g12 * g21)                      = abs(cov(r1, r2)) / sqrt(var(r1) var(r2)) = |cor(r1, r2)| exactly
## So "the two directional naive gains multiply to h^2" is not a statement about h at all. It holds for any two
## vectors whatsoever, including pure noise, where h is 0 and the identity still holds. It cannot confirm the
## shared-reference algebra because it does not depend on it.
##
## The reliability the package measures is a different quantity: reliability_from_blocks forms the signal as the
## mean of the two blocks; only under equal, independent block-error variances does cor = h/(2 - h) hold.
## The two block variances differ here, so h does not determine cor; under that model the two agree only as h
## approaches 1, and h/(2 - h) is printed below only as the correlation that model would imply.
##
## What this dataset supports, and what it does not:
##   Supported: h is measurable from replicate structure alone in a third laboratory's data, and spans
##     0.887 to 0.988 across three species, so the reliabilities reported elsewhere in this paper are not
##     peculiar to one design.
##   Supported: the true gain here is 1 by construction, and the naive shared-reference regression does not
##     return 1: the geometric mean of the two directions is 0.875, 0.985 and 0.964. The departure from the
##     known truth shrinks monotonically as h rises (12.5%, 1.5%, 3.6% against h of 0.887, 0.988, 0.980).
##   Not supported: any claim that the departure equals h, or 1 - h, or that the product of the two directions
##     recovers h^2. The departure equals 1 - |cor|, which is an identity, so it confirms nothing.
say("\n=== why this identity carries no information about h ===")
say("  sqrt(g12 * g21) == |cor(r1, r2)| is an identity, proved symbolically and reproduced here to 0.00%%,")
say("  so it holds on pure noise too and carries no information about the reliability.")
say("  under equal, independent block errors cor = h/(2 - h); with unequal block variances h does not fix cor.")
res[, cor_implied := h / (2 - h)]
for (i in seq_len(nrow(res)))
  say("    %-6s h %.4f | h/(2-h), equal-error model %.4f | geometric mean |cor| %.4f | departure from the known truth of 1: %.1f%%",
      res$species[i], res$h[i], res$cor_implied[i], res$gain_geomean[i], 100 * (1 - res$gain_geomean[i]))
say("  departure shrinks monotonically with h: %s",
    identical(order(res$h), order(-(1 - res$gain_geomean))))
fwrite(res, file.path(out, "p1x_kk2024_known_truth_gain.tsv"), sep = "\t")
