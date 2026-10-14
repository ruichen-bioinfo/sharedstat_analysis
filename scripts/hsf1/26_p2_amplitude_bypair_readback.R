# 26_p2_amplitude_bypair_readback.R
#
# This analysis evaluates the statement that "for two of the three replicate pairings at 42 C the corrected
# amplitude is 0". The saved summary contains the mean, minimum, and maximum of the per-pairing values, and
# therefore establishes only that at least one pairing has a value of zero. The stronger statement requires a trace
# from the per-pairing values; this script saves the
# per-pairing table so the statement can be checked.
#
# It does not re-derive the pairing logic. As in scripts 11/12/13/16/23, it takes I_rep / PAIRS_R /
# rms_true from script 18 by parse()+eval(), then verifies that the
# per-pairing values reproduce the saved mean, minimum, and maximum in
# p2_interaction_amplitude.tsv before output is written. If the check fails, the extracted implementation
# does not match and no output is saved.

suppressPackageStartupMessages({ library(data.table); library(edgeR); library(limma) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
pro  <- file.path(root, "data/02_processed"); meta <- file.path(root, "data/01_metadata")
hsd  <- file.path(pro, "hsf1"); outd <- hsd
repd <- file.path(root, "doc/05_reports")
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
sink(file.path(repd, "logs/hsf1_p2_amplitude_bypair.log"), split = TRUE)
say("=== 26_p2_amplitude_bypair_readback | %s ===", format(Sys.time()))

SRC <- file.path(root, "scripts/hsf1/18_p2_amplitude_and_cistrans_identity.R")
say("\nestimator source: %s", basename(SRC))

## ---- rebuild the same voom matrix script 18 builds, then transplant its estimator -------------
HEAT <- c("37_10", "37_30", "42_10", "42_30")
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
stopifnot(identical(SM$sample, libs))
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(M, group = factor(SM$grp)); keep <- filterByExpr(dge, group = SM$grp)
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
design <- model.matrix(~ 0 + factor(SM$grp)); colnames(design) <- levels(factor(SM$grp))
E <- voom(dge, design, plot = FALSE)$E
say("log-CPM matrix: %d genes x %d libraries", nrow(E), ncol(E))

# use rms_true, col_of, I_rep, and PAIRS_R from script 18
EX <- parse(SRC)
want <- c("rms_true", "col_of", "reps", "I_rep", "PAIRS_R")
got  <- character(0)
for (e in EX) {
  if (is.call(e) && length(e) >= 3L && identical(as.character(e[[1]]), "<-")) {
    nm <- as.character(e[[2]])
    if (length(nm) == 1L && nm %in% want) { eval(e); got <- c(got, nm) }
  }
}
say("functions read from script 18: %s", paste(got, collapse = ", "))
stopifnot(all(want %in% got), is.function(I_rep), is.function(rms_true), ncol(PAIRS_R) == 3L)

## ---- per-pairing values, exactly as script 18 forms them before collapsing --------------------
PER <- rbindlist(lapply(HEAT, function(c_) {
  rbindlist(lapply(seq_len(ncol(PAIRS_R)), function(j) {
    a <- I_rep(c_, PAIRS_R[1, j]); b <- I_rep(c_, PAIRS_R[2, j])
    if (is.null(a) || is.null(b)) return(NULL)
    sig <- (a + b) / 2; nul <- (b - a) / 2
    data.table(condition = c_, block_pair = paste0("r", PAIRS_R[1, j], "-r", PAIRS_R[2, j]),
               rms_sig = sqrt(mean(sig^2)), rms_null = sqrt(mean(nul^2)),
               rms_tru = rms_true(sig, nul))
  }))
}))

## ---- validation against the saved mean, minimum, and maximum ----------
AMP <- fread(file.path(hsd, "p2_interaction_amplitude.tsv"))
CHK <- PER[, .(rms_sig = mean(rms_sig), rms_null = mean(rms_null), rms_tru = mean(rms_tru),
               rms_tru_min = min(rms_tru), rms_tru_max = max(rms_tru)), by = condition]
CMP <- merge(AMP[, .(condition, rms_sig, rms_null, rms_tru, rms_tru_min, rms_tru_max)],
             CHK, by = "condition", suffixes = c("_landed", "_recomputed"))
say("\n---- validation: recomputed vs saved p2_interaction_amplitude.tsv ----")
worst <- 0
for (cc in c("rms_sig", "rms_null", "rms_tru", "rms_tru_min", "rms_tru_max")) {
  d <- max(abs(CMP[[paste0(cc, "_landed")]] - CMP[[paste0(cc, "_recomputed")]]))
  worst <- max(worst, d); say("  %-14s max |saved - recomputed| = %.3e", cc, d)
}
say("worst discrepancy over all five columns: %.3e", worst)
stopifnot(worst < 1e-9)
say("=> the extracted estimator reproduces script 18; the table below uses the same arithmetic.")

## ---- per-condition summary -----------------------------------------------------------
say("\n---- per replicate pairing ----")
say("%-8s %-10s %10s %10s %10s %6s", "cond", "block_pair", "rms_sig", "rms_null", "rms_tru", "==0?")
for (i in seq_len(nrow(PER))) say("%-8s %-10s %10.4f %10.4f %10.4f %6s",
    PER$condition[i], PER$block_pair[i], PER$rms_sig[i], PER$rms_null[i], PER$rms_tru[i],
    PER$rms_tru[i] == 0)
PER[, at_zero := rms_tru == 0]
Z <- PER[, .(n_pairs = .N, n_at_zero = sum(at_zero)), by = condition]
say("\n---- count of pairings whose noise-corrected amplitude is exactly 0 ----")
for (i in seq_len(nrow(Z))) say("  %-8s %d of %d pairings at zero", Z$condition[i], Z$n_at_zero[i], Z$n_pairs[i])
n42 <- Z[condition %in% c("42_10", "42_30"), sum(n_at_zero)]
say("\nat 42 C, summed over its two conditions: %d of %d pairings at zero", n42, Z[condition %in% c("42_10","42_30"), sum(n_pairs)])
say("42_10 alone: %d of 3 | 42_30 alone: %d of 3", Z[condition=="42_10", n_at_zero], Z[condition=="42_30", n_at_zero])
say("\nEvaluated statement: \"for two of the three replicate pairings at 42 C the corrected amplitude is 0\".")
say("The result depends on which of the two 42 C conditions is meant; the")
say("per-condition counts provide the corresponding values.")
say("")

fwrite(PER, file.path(outd, "p2_interaction_amplitude_bypair.tsv"), sep = "\t")
fwrite(Z,   file.path(outd, "p2_interaction_amplitude_zeropairs.tsv"), sep = "\t")
say("\nwrote p2_interaction_amplitude_bypair.tsv (%d rows), p2_interaction_amplitude_zeropairs.tsv", nrow(PER))
writeLines(capture.output(sessionInfo()), file.path(repd, "logs/sessionInfo_p2_bypair.txt"))
say("\n=== done ===")
sink()
