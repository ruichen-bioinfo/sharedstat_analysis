# 31_p1_design_grid.R -- when does a shared reference arm break the naive estimator, and does the
# correction hold across designs rather than only in ours?
#
# The objection this answers. P1's strongest objection is not that the algebra is wrong; it is
# "this is a specific parameterisation problem in one experimental design, and your 19-121% tells me
# nothing about mine." Two of our four artefacts are each demonstrated on one of the two datasets, so
# per artefact the objection stands even though the paper spans two datasets. A design grid answers it in
# the only way that generalises: by mapping the operating regime in which the naive estimator fails.
#
# Generative model, written out so it can be checked.
#   reference arm      R  ~ N(0, sR)            measured as  Rhat = R + eR,  eR ~ N(0, sE/sqrt(n))
#   response arm       A  = gamma * R           measured as  Ahat = A + eA,  eA ~ N(0, sE/sqrt(n))
#   contrast of interest      e = Ahat - Rhat
#   instrument               S  = R + independent noise, tuned to a target cor(S, R)
# The naive slope of e on Rhat is then, with h = var(R) / (var(R) + var(eR)):
#   slope = (gamma - 1) * h - (1 - h)
# which is exactly the decomposition P1 reports: a signal term attenuated by h, and an artefact term
# -(1 - h) that does not vanish with sample size because it comes from the shared error, not from noise.
# The instrumental-variable estimator gamma = 1 + cov(e, S) / cov(Rhat, S) removes it, because S shares
# the reference's signal but not its error.
#
# Negative control. At gamma = 1 there is no divergence at all. The naive estimator
# nevertheless reports slope = -(1 - h): a spurious divergence whose size is set purely by how noisy the
# reference arm is. That is the cleanest statement of the artefact available, and it needs no real data.
#
# This script adds no new estimator and ranks no published method. It maps where one construction fails.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1"); repd <- file.path(root, "doc/05_reports")
figd <- file.path(root, "doc/04_figures/final")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
sink(file.path(repd, "logs/hsf1_p1_design_grid.log"), split = TRUE)
say("=== 31_p1_design_grid | %s ===", format(Sys.time()))

## ---- transplant the estimator from the package so the paper and the package cannot drift ----
PKG <- file.path(root, "pkg/sharedstat/R/estimators.R")
source(PKG)
say("estimator sourced from the released package: %s", sub(paste0(root, "/"), "", PKG))
stopifnot(is.function(iv_gain))

N_GENE <- 2000L
NSIM   <- 200L

one_dataset <- function(gamma, h, n_rep, cor_S, seed) {
  set.seed(seed)
  sR <- 1
  ## h = var(R) / (var(R) + var(eR))  =>  var(eR) = var(R) * (1 - h) / h, and averaging n replicates
  ## divides the error variance by n, so the per-observation error sd is set to reach the target h at
  ## that replicate number. This is what makes n and h separately interpretable.
  sE <- sqrt(sR^2 * (1 - h) / h * n_rep)
  R  <- stats::rnorm(N_GENE, 0, sR)
  eR <- stats::rnorm(N_GENE, 0, sE / sqrt(n_rep))
  eA <- stats::rnorm(N_GENE, 0, sE / sqrt(n_rep))
  Rhat <- R + eR
  Ahat <- gamma * R + eA
  e <- Ahat - Rhat
  ## instrument: shares R's signal, independent of eR and eA
  k <- sqrt(1 / cor_S^2 - 1)
  S <- R + stats::rnorm(N_GENE, 0, sR * k)
  fit <- iv_gain(e, Rhat, S, naive = TRUE)
  list(gamma_iv = fit$gamma, gamma_naive = fit$gamma_naive,
       h_realised = stats::var(R) / stats::var(Rhat))
}

GRID <- CJ(gamma = c(1.0, 1.1, 1.3, 1.6), h = c(0.5, 0.7, 0.9),
           n_rep = c(2L, 3L, 6L, 12L), cor_S = c(0.5, 0.8))
say("grid: %d cells x %d simulations = %d datasets", nrow(GRID), NSIM, nrow(GRID) * NSIM)

RES <- rbindlist(lapply(seq_len(nrow(GRID)), function(i) {
  g <- GRID[i]
  sims <- lapply(seq_len(NSIM), function(s) one_dataset(g$gamma, g$h, g$n_rep, g$cor_S, seed = i * 1000L + s))
  gi <- vapply(sims, function(x) x$gamma_iv, 0)
  gn <- vapply(sims, function(x) x$gamma_naive, 0)
  data.table(g,
    gamma_naive = mean(gn), gamma_iv = mean(gi),
    bias_naive = mean(gn) - g$gamma, bias_iv = mean(gi) - g$gamma,
    rel_bias_naive = (mean(gn) - g$gamma) / max(abs(g$gamma - 1), 1e-9),
    sd_iv = stats::sd(gi),
    ## a sign reversal is when the naive estimate says the response is damped (gamma < 1) while the
    ## truth is that it is amplified (gamma > 1), or vice versa. For gamma = 1 it is the fraction of
    ## simulations reporting any apparent direction at all beyond a negligible band.
    sign_wrong_naive = if (g$gamma > 1) mean(gn < 1) else mean(abs(gn - 1) > 0.05),
    sign_wrong_iv    = if (g$gamma > 1) mean(gi < 1) else mean(abs(gi - 1) > 0.05))
}))
fwrite(RES, file.path(hsd, "p1_design_grid.tsv"), sep = "\t")

## ---------------------------------------------------- the negative control ----
say("\n================ NEGATIVE CONTROL: TRUE gamma = 1, i.e. NO DIVERGENCE AT ALL ================")
NC <- RES[gamma == 1.0][order(h, n_rep, cor_S)]
say("%6s %7s %8s %14s %14s %16s %16s", "h", "n_rep", "cor(S,R)", "naive gamma", "IV gamma", "naive says diff", "IV says diff")
for (i in seq_len(nrow(NC))) say("%6.1f %7d %8.1f %14.4f %14.4f %15.0f%% %15.0f%%",
    NC$h[i], NC$n_rep[i], NC$cor_S[i], NC$gamma_naive[i], NC$gamma_iv[i],
    100 * NC$sign_wrong_naive[i], 100 * NC$sign_wrong_iv[i])
say("\npredicted naive value under no divergence is 1 - (1 - h) = h:")
for (hh in sort(unique(NC$h))) say("  h = %.1f => predicted %.3f, observed %.4f (mean over %d cells)",
    hh, hh, NC[h == hh, mean(gamma_naive)], NC[h == hh, .N])
say("=> the naive estimator reports a divergence that is set ENTIRELY by how noisy the shared reference")
say("   is, on data containing none. The IV estimator returns %.4f on average across every no-divergence",
    NC[, mean(gamma_iv)])
say("   cell, against a truth of 1.")

## -------------------------------------------------------- the operating regime ----
say("\n================ OPERATING REGIME: WHERE DOES THE NAIVE ESTIMATOR FAIL? ================")
say("relative bias of the naive estimator, as a fraction of the true effect (gamma - 1):")
say("%6s %26s %26s %26s", "h", "gamma=1.1", "gamma=1.3", "gamma=1.6")
for (hh in sort(unique(RES$h))) {
  v <- vapply(c(1.1, 1.3, 1.6), function(gg) RES[h == hh & gamma == gg, mean(rel_bias_naive)], 0)
  say("%6.1f %26.2f %26.2f %26.2f", hh, v[1], v[2], v[3])
}
say("\nsame table for the IV estimator:")
say("%6s %26s %26s %26s", "h", "gamma=1.1", "gamma=1.3", "gamma=1.6")
for (hh in sort(unique(RES$h))) {
  v <- vapply(c(1.1, 1.3, 1.6), function(gg) RES[h == hh & gamma == gg, mean((gamma_iv - gg) / (gg - 1))], 0)
  say("%6.1f %26.4f %26.4f %26.4f", hh, v[1], v[2], v[3])
}
say("\nfraction of simulations in which the naive estimator reports the WRONG DIRECTION:")
say("%6s %14s %14s %14s", "h", "gamma=1.1", "gamma=1.3", "gamma=1.6")
for (hh in sort(unique(RES$h))) {
  v <- vapply(c(1.1, 1.3, 1.6), function(gg) RES[h == hh & gamma == gg, mean(sign_wrong_naive)], 0)
  say("%6.1f %13.0f%% %13.0f%% %13.0f%%", hh, 100*v[1], 100*v[2], 100*v[3])
}
say("\nreplicate number does NOT rescue the naive estimator, which is the whole point:")
for (nn in sort(unique(RES$n_rep))) say("  n_rep = %2d => mean relative bias %.2f (IV: %+.4f)",
    nn, RES[n_rep == nn & gamma > 1, mean(rel_bias_naive)], RES[n_rep == nn & gamma > 1, mean(bias_iv)])
say("=> the artefact term -(1 - h) comes from the SHARED ERROR, not from sampling noise, so averaging")
say("   more replicates at fixed h leaves it untouched. This is why it is not a power problem.")

say("\nthe instrument's strength affects the IV estimator's VARIANCE, not its centre:")
for (cc in sort(unique(RES$cor_S))) say("  cor(S,R) = %.1f => IV bias %+.4f, IV sd %.4f",
    cc, RES[cor_S == cc & gamma > 1, mean(bias_iv)], RES[cor_S == cc, mean(sd_iv)])

say("\nWHAT THIS DOES AND DOES NOT SHOW.")
say("It shows the regime in which one construction fails and that one correction removes it across that")
say("regime, on data where the truth is set by us. It does NOT compare against any published method, and")
say("no method is ranked here: the comparison is between two readings of the same data, one of which is")
say("algebraically biased by construction. It also does not establish that any real design has a")
say("particular h; that has to be measured, and the paper reports the measured value where it has one.")
writeLines(capture.output(sessionInfo()), file.path(repd, "logs/sessionInfo_p1_grid.txt"))
say("\n=== done ===")
sink()
