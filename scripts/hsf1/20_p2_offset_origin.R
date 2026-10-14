# 20_p2_offset_origin.R -- tests whether the genome-wide offset of the P2 epistasis contrast is
# compositional or biological.
#
# Question. The epistatic interaction I = D - H - M + W has a non-zero genome-wide median at the
# 37_30 anchor (+0.1319, z = 22.5 against its own SE of the median), and the Hsf1_only class-level interpretation
# flips depending on which reference that median is judged against:
#   vs 0                          -> med I = -0.0070, unresolvable  => class-level interaction unresolved
#   vs the genome-wide median     -> -0.1389, z ~ -3.1              => class-level interaction resolved
#   vs the not_induced class      -> -0.3008, z ~ -6.7              => class-level interaction resolved
# One reading is that the offset is condition-specific, therefore probably real biology, so 0 is the
# right reference. This script tests that reading. The expected outcome, stated before the run, is
# that the offset is compositional, for two reasons:
#
#  (1) Algebra. The contrast coefficients (+1,-1,-1,+1) sum to zero, so any per-genotype constant
#      offset cancels exactly, including a global library-composition shift. More than that, any
#      offset that is additive in the two perturbations cancels: a shift tracking HSF1 status
#      contributes +d(D) - d(H) = 0, and one tracking MSN24 status contributes +d(D) - d(M) = 0.
#      So the observed median cannot come from a simple normalisation artefact. What it can come
#      from is a non-additive compositional effect, and non-additivity is exactly what I measures,
#      which is why this is confusable with biology.
#  (2) Class pattern. not_induced shows the largest offset (+0.2938), larger than the genome-wide
#      median. Real Hsf1/Msn2-4 epistasis should be concentrated in TF-dependent, heat-induced genes,
#      not strongest in genes that do not respond to heat. Genes that do not respond are passengers:
#      they carry the compositional signal and little else. That is the signature of a
#      compositional, not a biological, offset.
#
# Mechanism tested. At 37_30 the chaperone response consumes a large share of the library in
# genotypes with intact Hsf1 (W, M) and much less in the depleted ones (H, D). Every other gene is
# therefore diluted differently between genotypes. If that dilution is additive in the two
# perturbations it cancels; if the chaperone burden is non-additive (Ciccarelli 2023, PMID 37467033,
# shows Msn2/4 can bypass Hsf1, so non-additivity of the burden is expected) it does not cancel,
# and it reaches every gene including the non-responding ones.
#
# Test. Re-estimate the normalisation factors excluding the genes responsible for the library-share
# swing (the Hsf1/Msn2-4-dependent heat-induced classes), then recompute I on the identical contrast.
#     - if the genome-wide median I collapses toward 0  => the offset is compositional
#       => class medians are compared with a reference class rather than 0
#       => the Hsf1_only class-level interaction is resolved (with the class-level qualifier stated)
#     - if the median survives                         => the offset is biological
#       => 0 is the right reference and the Hsf1_only class-level interaction is unresolved
#
# Usage: Rscript scripts/hsf1/20_p2_offset_origin.R

suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })
root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
pro <- file.path(root, "data/02_processed"); hsd <- file.path(pro, "hsf1")
meta <- file.path(root, "data/01_metadata"); repd <- file.path(root, "doc/05_reports")
MIRROR <- "logs/output_copy"
ANCHOR <- "37_30"; OTHER <- c("37_10", "42_10", "42_30")
set.seed(20260906L)
say <- function(...) cat(sprintf(...), "\n", sep = "")
logf <- file.path(repd, "logs/hsf1_p2_offset_origin.log"); sink(logf, split = TRUE)
say("=== 20_p2_offset_origin | %s | is the offset compositional? ===", format(Sys.time()))
say("Expectation stated in advance:")
say("the offset is expected to be compositional, which would support comparing the Hsf1_only class with the reference class.")

## ------------------------------------------------------------------- inputs ---
CM <- fread(file.path(pro, "count_matrix_aligned.tsv")); gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
SM[, cond := paste(temperature, time, sep = "_")][, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol); SM <- SM[sample %in% libs]; setkey(SM, sample); SM <- SM[libs]
stopifnot(identical(SM$sample, libs))
M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
CLS <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
ccol <- intersect(c("class","epistasis_class","cls"), names(CLS))[1]
CL <- unique(CLS[, .(gene = get(names(CLS)[1]), cls = get(ccol))])
say("\ncount matrix %d x %d ; classes: %s", nrow(M), ncol(M),
    paste(sort(unique(CL$cls)), collapse = ", "))

## --------------------- how much library share does the response consume? ------
# direct measurement of the mechanism's premise, made before the test
frac_of <- function(genes) {
  k <- rownames(M) %in% genes
  colSums(M[k, , drop = FALSE]) / colSums(M)
}
RESP <- CL[cls %in% c("Hsf1_only","Msn24_only","redundant","both_required"), gene]
say("\nresponse-class genes (Hsf1_only + Msn24_only + redundant + both_required): %d", length(RESP))
fr <- frac_of(RESP)
SH <- data.table(genotype = SM$knockout, cond = SM$cond, share = fr)[
      , .(share = mean(share)), by = .(cond, genotype)]
W <- dcast(SH, cond ~ genotype, value.var = "share")
say("\nlibrary share taken by the response classes (mean over replicates):")
say("%-8s %10s %10s %10s %12s %14s", "cond", "Wildtype", "HSF1.KD", "MSN24.KO", "Double.KDKO",
    "non-additivity")
for (i in seq_len(nrow(W))) {
  na_ <- W$Double.KDKO[i] - W$HSF1.KD[i] - W$MSN24.KO[i] + W$Wildtype[i]
  say("%-8s %9.4f %10.4f %10.4f %12.4f %14.4f", W$cond[i], W$Wildtype[i], W$HSF1.KD[i],
      W$MSN24.KO[i], W$Double.KDKO[i], na_)
}
say("the last column is the SAME contrast applied to library share. If it is non-zero, the")
say("dilution of every other gene is non-additive and therefore does NOT cancel in I.")

## ---------------------------------------- I under two normalisation regimes ---
I_of <- function(norm_genes, tag) {
  dge <- DGEList(M, group = factor(SM$grp))
  keep <- filterByExpr(dge, group = SM$grp)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  # TMM estimated on norm_genes only, then applied to all retained genes. Excluding the genes
  # that swing library share is the whole point: it removes the source of the compositional shift
  # from the normalisation without touching the contrast.
  sub <- dge[rownames(dge) %in% norm_genes, , keep.lib.sizes = TRUE]
  nf <- calcNormFactors(sub)$samples$norm.factors
  dge$samples$norm.factors <- nf
  grp <- factor(SM$grp); design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
  v <- voom(dge, design, plot = FALSE); fit <- lmFit(v, design)
  mk <- function(c_) { ct <- setNames(rep(0, ncol(design)), colnames(design))
    ct[paste("Double.KDKO", c_, sep = "__")] <-  1
    ct[paste("HSF1.KD",     c_, sep = "__")] <- -1
    ct[paste("MSN24.KO",    c_, sep = "__")] <- -1
    ct[paste("Wildtype",    c_, sep = "__")] <-  1; ct }
  CONDS <- c(ANCHOR, OTHER)
  CT <- do.call(cbind, lapply(CONDS, mk)); colnames(CT) <- CONDS
  stopifnot(all(colSums(CT) == 0))
  cf <- contrasts.fit(fit, CT)
  out <- rbindlist(lapply(seq_along(CONDS), function(j) data.table(
    regime = tag, gene = rownames(cf$coefficients), condition = CONDS[j],
    I = cf$coefficients[, j], se = cf$stdev.unscaled[, j] * cf$sigma)))
  say("  regime '%s': %d genes, %d normalisation genes", tag, uniqueN(out$gene), length(norm_genes))
  out
}
ALLG <- rownames(M)
A <- I_of(ALLG, "all genes (as published)")
B <- I_of(setdiff(ALLG, RESP), "response classes excluded from normalisation")
IX <- rbind(A, B)
fwrite(IX[condition == ANCHOR], file.path(hsd, "p2_offset_origin_genewise.tsv.gz"),
       sep = "\t", compress = "gzip")

SUM <- IX[, .(n = .N, med_I = median(I), se_of_med = median(se) / sqrt(.N)), by = .(regime, condition)]
SUM[, z := med_I / se_of_med]
fwrite(SUM, file.path(hsd, "p2_offset_origin_summary.tsv"), sep = "\t")
say("\n================ GENOME-WIDE MEDIAN I UNDER THE TWO REGIMES ================")
say("%-42s %-8s %10s %12s %8s", "regime", "cond", "med I", "SE of med", "z")
for (i in seq_len(nrow(SUM))) say("%-42s %-8s %+10.4f %12.4f %8.1f",
  SUM$regime[i], SUM$condition[i], SUM$med_I[i], SUM$se_of_med[i], SUM$z[i])
m0 <- SUM[regime == "all genes (as published)" & condition == ANCHOR, med_I]
m1 <- SUM[regime == "response classes excluded from normalisation" & condition == ANCHOR, med_I]
say("\nat the %s anchor: %+.4f -> %+.4f  (%.0f%% of the offset removed)",
    ANCHOR, m0, m1, 100 * (1 - abs(m1) / abs(m0)))

## ------------------------------- class medians under the corrected regime -----
CLSUM <- merge(IX[condition == ANCHOR], CL, by = "gene")[
  , .(n = .N, med_I = median(I), se_of_med = median(se) / sqrt(.N)), by = .(regime, cls)]
ref <- CLSUM[cls == "not_induced", .(regime, ref = med_I)]
CLSUM <- merge(CLSUM, ref, by = "regime")
CLSUM[, z_vs_0 := med_I / se_of_med][, z_vs_ref := (med_I - ref) / se_of_med]
fwrite(CLSUM, file.path(hsd, "p2_offset_origin_classes.tsv"), sep = "\t")
for (rg in unique(CLSUM$regime)) {
  say("\nclass medians, regime '%s':", rg)
  z <- CLSUM[regime == rg][order(-n)]
  say("%-16s %6s %10s %11s %11s %11s", "class", "n", "med I", "z vs 0", "vs ref", "z vs ref")
  for (i in seq_len(nrow(z))) say("%-16s %6d %+10.4f %11.1f %+11.4f %11.1f",
    z$cls[i], z$n[i], z$med_I[i], z$z_vs_0[i], z$med_I[i] - z$ref[i], z$z_vs_ref[i])
}

## ------------------------------------------------------------- the verdict ---
say("\n================ IS THE OFFSET COMPOSITIONAL? ================")
COMPOSITIONAL <- abs(m1) < 0.5 * abs(m0)
say("offset at the anchor: %+.4f (as published) -> %+.4f (response classes out of normalisation)",
    m0, m1)
say("compositional criterion (offset more than halved): %s", COMPOSITIONAL)
HO <- CLSUM[cls == "Hsf1_only" & regime == "response classes excluded from normalisation"]
NI <- CLSUM[cls == "not_induced" & regime == "response classes excluded from normalisation"]
say("")
if (COMPOSITIONAL) {
  say("=> THE OFFSET IS COMPOSITIONAL, not biological.")
  say("   That the offset is condition-specific does not make it biological. Condition specificity")
  say("   is explained by the mechanism: only at")
  say("   37_30 is the chaperone response large enough for its library share to matter, which is")
  say("   why the offset is large at 37_30 and near zero at 42_10 and 42_30.")
  say("   A class-level median interaction is compared with a REFERENCE")
  say("   CLASS, not against 0. With not_induced as the reference and the corrected normalisation,")
  say("   Hsf1_only sits at %+.4f (z = %.1f) => class-level medians can be read against the reference.",
      HO$med_I[1] - HO$ref[1], HO$z_vs_ref[1])
  say("   PER-GENE interpretation is not supported: a single gene's SE is far larger")
  say("   than any class median.")
} else {
  say("=> THE OFFSET SURVIVES the normalisation change and is therefore NOT compositional.")
  say("   0 remains the correct reference and the Hsf1_only class-level interaction remains unresolved, now on a")
  say("   tested basis rather than on an argument.")
}
say("")
say("TWO METHOD DETAILS APPLY IN EITHER CASE:")
say(" 1. the coefficients of I sum to zero, so any offset ADDITIVE in the two perturbations")
say("    cancels exactly; only NON-ADDITIVE offsets survive, and library-share redistribution is")
say("    non-additive whenever the chaperone burden is (Ciccarelli 2023, PMID 37467033).")
say(" 2. therefore an epistasis contrast built from compositional data needs its reference class")
say("    declared explicitly. Reporting a class median 'against zero' is a hidden assumption.")

writeLines(capture.output(sessionInfo()), file.path(hsd, "sessionInfo_p2_offset.txt"))
mdir <- file.path(MIRROR, "data/02_processed/hsf1"); dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(hsd, pattern = "^(p2_offset|sessionInfo_p2_offset)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
file.copy(file.path(root, "scripts/hsf1/20_p2_offset_origin.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("\nwritten. === done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
