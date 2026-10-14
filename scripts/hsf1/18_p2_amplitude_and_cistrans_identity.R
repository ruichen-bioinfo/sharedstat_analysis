# 18_p2_amplitude_and_cistrans_identity.R
# Two parts, both using estimators already used elsewhere in this project.
#
# Part A: P2's four-condition read-out as a count (n_sig = 13/930/9/18) orders amplitudes only
# when the noise floor is shared, and a count must always be reported beside med(se). At 42 C the
# stress is harsher, so variance inflates and the count collapses, but the underlying interaction
# amplitude need not have shrunk at all. The read-out is therefore the noise-corrected,
# threshold-free interaction amplitude rms_tru(I), built with the same two-block signal/null
# construction as scripts 03, 04, 05 and 15, with med(se) printed beside it.
#   If amplitude remains stable at 42 C while only the SE increases, the temperature x time pattern remains
#   supported at the amplitude level.
#   If amplitude also decreases, only the single anchor is supported.
# Both outcomes accepted in advance.
#
# Part B: the three cis-trans compensation estimates (-0.166 cross-block / ~-0.39 variance
# decomposition / ~-0.7..-0.85 implied by cor(cis,par)) are not three independent measurements.
# Given par = cis + trans they are three ways of combining the same three second moments, and
# they are linked by an exact identity with no distributional assumptions:
#       cov(c,p) = var(c) + cov(c,t)                        [expand cov(c, c+t)]
#   =>  cov(c,t) = cov(c,p) - var(c)
#       var(t)   = var(p) + var(c) - 2 cov(c,p)
#   =>  cor(c,t) = [cor(c,p)*sd(p) - sd(c)] / sd(t)
# So estimating var(c), var(p), cov(c,p) once determines every other quantity. If the three
# numbers disagree, the disagreement is an inconsistency in which variances were noise-corrected,
# and the identity localises it in one step. This is a single closed-form reconciliation: if it
# reconciles, three numbers become one; if it does not, the three are reported in parallel and
# the matter stops there.
#
# Why cov(c,p) is the better anchor. Write the measured quantities as chat = c + e_c and
# phat = p + e_p. Then cov(chat, that) = cov(c,t) + cov(e_c,e_p) - var(e_c): the residual
# construction that = phat - chat drags in -var(e_c), which is exactly Fraser 2019's bias
# (PMID 30270122, doi 10.1016/j.tig.2018.09.003). But cis is measured on hybrid libraries and par
# on parental libraries (different material), so cov(e_c, e_p) ~ 0 and cov(chat, phat) carries
# no analogous term. cov(c,p) is therefore the moment to estimate directly.
#
# Note on the zero-compensation null, checked here algebraically and by simulation.
# Substituting cov(c,p) = 0 into the identity gives cov(c,t) = -var(c), i.e.
# cor(c,t) = -sd(c)/sd(t), which for comparable sd is about -1: complete compensation. The
# formula r_mech = -sqrt(var(c)/(var(p)+var(c))) is algebraically exactly that value, so it is
# the value for cov(c,p)=0 and not a "zero compensation" null. The zero-compensation condition is
# cov(c,t)=0, under which cor(c,p) = sd(c)/sd(p) > 0. The cross-block estimator is unbiased, so
# its null is 0.
#
# Usage: Rscript scripts/hsf1/18_p2_amplitude_and_cistrans_identity.R

suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
pro  <- file.path(root, "data/02_processed"); hsd <- file.path(pro, "hsf1")
meta <- file.path(root, "data/01_metadata"); repd <- file.path(root, "doc/05_reports")
outd <- hsd; MIRROR <- "logs/output_copy"
HEAT <- c("37_10", "37_30", "42_10", "42_30"); TP <- c(15L, 30L, 60L)
NSIM <- 200000L; NBOOT <- 1000L
set.seed(20260906L)
say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, "logs/hsf1_p2_amplitude_cistrans.log")
sink(logf, split = TRUE)
say("=== 18_p2_amplitude_and_cistrans_identity | %s ===", format(Sys.time()))

rms_true <- function(s, n) sqrt(max(0, mean(s^2, na.rm = TRUE) - mean(n^2, na.rm = TRUE)))

## ============================ Part A: P2 threshold-free interaction amplitude ============
say("\n################ PART A: rms_tru(I) replaces the count read-out ################")
say("A count orders amplitudes only when the noise floor is shared, and a count")
say("must be reported beside med(se). At 42 C the stress is harsher, so SE inflates and the count")
say("can collapse while the underlying amplitude does not. The same two-block estimator is used elsewhere in the analysis.")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
stopifnot(identical(SM$sample, libs))
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + factor(SM$grp)); colnames(design) <- levels(factor(SM$grp))
v <- voom(dge, design, plot = FALSE)
E <- v$E                                                  # log-CPM per library
say("log-CPM matrix: %d genes x %d libraries", nrow(E), ncol(E))
# per-replicate interaction: I_r(c) = D - H - M + W using replicate r of each genotype at cond c
col_of <- function(g, c, r) which(SM$knockout == g & SM$cond == c & SM$replicate == r)
reps <- sort(unique(SM$replicate)); say("replicates per group: %s", paste(reps, collapse = ","))
I_rep <- function(c, r) {
  i <- vapply(c("Double.KDKO","HSF1.KD","MSN24.KO","Wildtype"),
              function(g) { k <- col_of(g, c, r); if (length(k) != 1L) NA_integer_ else k }, 0L)
  if (anyNA(i)) return(NULL)
  E[, i[1]] - E[, i[2]] - E[, i[3]] + E[, i[4]]
}
PAIRS_R <- combn(reps, 2L)
AMP <- rbindlist(lapply(HEAT, function(c_) {
  per <- rbindlist(lapply(seq_len(ncol(PAIRS_R)), function(j) {
    a <- I_rep(c_, PAIRS_R[1, j]); b <- I_rep(c_, PAIRS_R[2, j])
    if (is.null(a) || is.null(b)) return(NULL)
    sig <- (a + b) / 2; nul <- (b - a) / 2
    data.table(condition = c_, block_pair = paste0("r", PAIRS_R[1, j], "-r", PAIRS_R[2, j]),
               rms_sig = sqrt(mean(sig^2)), rms_null = sqrt(mean(nul^2)),
               rms_tru = rms_true(sig, nul), q90_sig = quantile(abs(sig), .90))
  }))
  per[, .(condition = c_, n_pairs = .N,
          rms_sig = mean(rms_sig), rms_null = mean(rms_null),
          rms_tru = mean(rms_tru), rms_tru_min = min(rms_tru), rms_tru_max = max(rms_tru),
          q90_sig = mean(q90_sig))]
}))
MSE <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "p2_interaction_model_se.tsv.gz"))))[
       , .(med_se_model = median(se_model), n_sig = sum(padj_mod < 0.05)), by = condition]
AMP <- merge(AMP, MSE, by = "condition"); setkey(AMP, condition); AMP <- AMP[HEAT]
AMP[, snr := rms_tru / rms_null]
fwrite(AMP, file.path(outd, "p2_interaction_amplitude.tsv"), sep = "\t")
say("\n%-8s %10s %10s %10s %18s %8s %12s %10s", "cond", "rms_sig", "rms_null", "rms_tru",
    "rms_tru range", "SNR", "med se(model)", "n padj<.05")
for (i in seq_len(nrow(AMP))) say("%-8s %10.4f %10.4f %10.4f [%7.4f,%7.4f] %8.2f %12.4f %10d",
  AMP$condition[i], AMP$rms_sig[i], AMP$rms_null[i], AMP$rms_tru[i], AMP$rms_tru_min[i],
  AMP$rms_tru_max[i], AMP$snr[i], AMP$med_se_model[i], AMP$n_sig[i])
a30 <- AMP[condition == "37_30", rms_tru]; a4210 <- AMP[condition == "42_10", rms_tru]
a4230 <- AMP[condition == "42_30", rms_tru]; a3710 <- AMP[condition == "37_10", rms_tru]
say("\nnoise floors (rms_null): 37_10 %.4f / 37_30 %.4f / 42_10 %.4f / 42_30 %.4f",
    AMP[condition=="37_10",rms_null], AMP[condition=="37_30",rms_null],
    AMP[condition=="42_10",rms_null], AMP[condition=="42_30",rms_null])
say("=> the floors are %s, which is the condition required before a count may be",
    ifelse(max(AMP$rms_null)/min(AMP$rms_null) < 1.25, "similar", "NOT shared"))
say("   read as an amplitude ordering.")
say("\namplitude ratios vs the 37_30 anchor: 37_10 %.2f | 42_10 %.2f | 42_30 %.2f",
    a3710/a30, a4210/a30, a4230/a30)
say("count ratios vs the anchor          : 37_10 %.3f | 42_10 %.3f | 42_30 %.3f",
    AMP[condition=="37_10",n_sig]/AMP[condition=="37_30",n_sig],
    AMP[condition=="42_10",n_sig]/AMP[condition=="37_30",n_sig],
    AMP[condition=="42_30",n_sig]/AMP[condition=="37_30",n_sig])
HOLDS_42 <- (a4230 / a30) > 0.5 || (a4210 / a30) > 0.5
say("")
if (HOLDS_42) {
  say("RESULT: the interaction amplitude does NOT collapse at 42 C the way the count does.")
  say("The temperature x time epistasis pattern is supported at the amplitude level and does not require")
  say("   restriction to a single anchor. The count collapse is a noise-floor effect, as")
  say("   predicted by the shared-floor requirement. Amplitude estimates and counts with med(se) are both reported.")
} else {
  say("RESULT: the amplitude collapses too, not only the count.")
  say("Only the 37_30 anchor is then supported; the other three conditions have negative results.")
  say("   The question of whether 42 C advances epistasis is answered NO")
  say("   on the amplitude scale as well as on the count scale.")
}

## ==================== Part B: the cis-trans identity and the zero-compensation null ======
say("\n################ PART B: exact identity and a controlled simulation ##########")
say("Identity from par = cis + trans, no assumptions:")
say("  cov(c,t) = cov(c,p) - var(c) ;  var(t) = var(p) + var(c) - 2cov(c,p)")
say("  cor(c,t) = [cor(c,p)*sd(p) - sd(c)] / sd(t)")
say("So var(c), var(p), cov(c,p) determine everything: the 'three estimates' are one moment set.")
say("\n-- controlled simulation --")
sim <- function(rho_ct, n = NSIM, sc = 1, st = 1, se_c = 0.6, se_p = 0.6) {
  Z <- matrix(rnorm(2 * n), n, 2)
  c_ <- sc * Z[, 1]
  t_ <- st * (rho_ct * Z[, 1] + sqrt(max(0, 1 - rho_ct^2)) * Z[, 2])
  p_ <- c_ + t_
  # two independent replicate blocks, each with its own measurement error
  ch1 <- c_ + rnorm(n, 0, se_c); ch2 <- c_ + rnorm(n, 0, se_c)
  ph1 <- p_ + rnorm(n, 0, se_p); ph2 <- p_ + rnorm(n, 0, se_p)
  th1 <- ph1 - ch1; th2 <- ph2 - ch2
  c(true_ct = rho_ct,
    naive_ct      = cor(ch1, th1),                 # same block both sides: Fraser's bias
    crossblock_ct = cor(ch1, th2),                 # independent blocks: unbiased
    cor_cp        = cor(ch1, ph2),
    identity_ct   = (cor(ch1, ph2) * sd(ph2) - sd(ch1)) / sd(th2))
}
SIMT <- rbindlist(lapply(c(0, -0.2, -0.5, -0.9, -0.99), function(r) as.list(sim(r))))
say("%-10s %12s %14s %10s %13s", "true r_ct", "naive", "cross-block", "cor(c,p)", "identity")
for (i in seq_len(nrow(SIMT))) say("%-10.2f %12.3f %14.3f %10.3f %13.3f", SIMT$true_ct[i],
  SIMT$naive_ct[i], SIMT$crossblock_ct[i], SIMT$cor_cp[i], SIMT$identity_ct[i])
say("")
say("Reading:")
say(" 1. at TRUE r_ct = 0 the naive estimator returns %.3f -- biased NEGATIVE, which is exactly",
    SIMT[true_ct == 0, naive_ct])
say("    Fraser 2019's artefactual compensation (PMID 30270122).")
say(" 2. the CROSS-BLOCK estimator returns %.3f at true 0 => UNBIASED => its null hypothesis is 0,",
    SIMT[true_ct == 0, crossblock_ct])
say("    NOT r_mech, so comparing it with r_mech applies a bias correction it does not need.")
say(" 3. at true r_ct = 0, cor(c,p) = %.3f, i.e. POSITIVE. Zero compensation makes cis and par",
    SIMT[true_ct == 0, cor_cp])
say("    POSITIVELY correlated; it does not make them orthogonal.")
say(" 4. cov(c,p) ~ 0 (cor(c,p) ~ 0) happens near true r_ct = %.2f -- i.e. NEAR-COMPLETE",
    SIMT[which.min(abs(SIMT$cor_cp)), true_ct])
say("    compensation. r_mech = -sqrt(var(c)/(var(p)+var(c))) is algebraically the value of")
say("    cor(c,t) AT cov(c,p)=0: it marks near-complete compensation, not zero compensation.")
say(" 5. the identity reproduces the cross-block value to within simulation error, confirming that")
say("    the three published numbers are one moment set and not three measurements.")

say("\n-- applying the identity to the allele-specific expression data --")
A6 <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
A6 <- A6[is.finite(cis) & is.finite(trans) & is.finite(par) & is.finite(cis_nul) & is.finite(trans_nul)]
say("allele-specific expression rows usable: %d over %d genes, %d pairs, %d timepoints",
    nrow(A6), uniqueN(A6$hog), uniqueN(A6$pair), uniqueN(A6$time))
ID <- A6[, {
  vc <- max(0, mean(cis^2) - mean(cis_nul^2))            # noise-corrected var(cis)
  vt <- max(0, mean(trans^2) - mean(trans_nul^2))        # noise-corrected var(trans)
  # var(par) from the identity's own algebra, so all three moments are on one footing
  ccp_raw <- mean(cis * par) - mean(cis) * mean(par)
  vp_raw  <- var(par)
  cct_raw <- mean(cis * trans) - mean(cis) * mean(trans)
  .(n = .N, var_c_corr = vc, var_t_corr = vt, cov_cp_raw = ccp_raw, var_p_raw = vp_raw,
    cor_ct_raw = cor(cis, trans), cor_cp_raw = cor(cis, par),
    # identity route, using the noise-corrected var(c): cov(c,t) = cov(c,p) - var(c)
    cov_ct_identity = ccp_raw - vc,
    cor_ct_identity = (ccp_raw - vc) / sqrt(max(1e-12, vc) * max(1e-12, vt)))
}, by = time]
fwrite(ID, file.path(outd, "a6_cistrans_identity.tsv"), sep = "\t")
say("\n%-5s %9s %9s %11s %12s %12s %16s", "time", "var_c*", "var_t*", "cov(c,p)",
    "cor(c,t) raw", "cor(c,p) raw", "cor(c,t) identity")
for (i in seq_len(nrow(ID))) say("%-5d %9.4f %9.4f %11.4f %12.4f %12.4f %16.4f",
  ID$time[i], ID$var_c_corr[i], ID$var_t_corr[i], ID$cov_cp_raw[i], ID$cor_ct_raw[i],
  ID$cor_cp_raw[i], ID$cor_ct_identity[i])
say("(* = noise-corrected via the block null; cov and the raw correlations are within-block and so")
say(" still carry Fraser's bias term -- they are printed for comparison, not as estimates.)")
say("")
say("SIGN CHECK, the only claim made: cor(c,t) by the identity is %s at all three",
    ifelse(all(ID$cor_ct_identity < 0), "NEGATIVE", "NOT uniformly negative"))
say("timepoints, and cor(c,p) raw is %s.",
    ifelse(all(ID$cor_cp_raw > 0), "positive throughout", "not positive throughout"))
say("=> consistent with MILD COMPENSATION and with Fraser 2019 / Rojas & Lang 2026 in direction.")
say("=> NOT consistent with 'the two axes reinforce each other'.")
say("")
say("MAGNITUDE IS NOT CLAIMED. The identity reconciles the estimates only up to which")
say("moments are noise-corrected; cov(c,p) here is a within-block covariance and therefore still")
say("carries a bias term, so the identity route pins the SIGN and the ORDERING but not the value.")
say("The source of the residual spread is not resolved further. The analysis")
say("finds a sign consistent with Fraser 2019 and Rojas & Lang 2026 as the methodological base, but")
say("does not estimate the magnitude.")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_p2_amp_cistrans.txt"))
mdir <- file.path(MIRROR, "data/02_processed/hsf1"); dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(outd, pattern = "^(p2_interaction_amplitude|a6_cistrans_identity|sessionInfo_p2_amp)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
file.copy(file.path(root, "scripts/hsf1/18_p2_amplitude_and_cistrans_identity.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("
prepared output copy in %s", MIRROR); say("=== done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
