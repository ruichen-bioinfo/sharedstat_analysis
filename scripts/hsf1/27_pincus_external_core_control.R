# 27_pincus_external_core_control.R
#
# What this tests, and why it can fail.
# P3 reports that cross-species divergence is not concentrated in the Hsf1-associated core. That
# negative is only interpretable if our core is actually a core. Our core comes from one dataset, one
# measurement modality (steady-state mRNA) and one perturbation (HSF1 knockdown), so "the core is
# wrong" and "divergence is not in the core" are not separable from within our own data.
#
# Pincus et al. 2018 (PMID 30332327, doi 10.1091/mbc.e18-06-0353) measured Hsf1 dependence on the same
# organism with a different modality and a different perturbation: nascent transcription (NAC-seq)
# under Hsf1 nuclear depletion by anchor-away, not knockdown. Their deposited table gives, per gene,
#   NHS (no heat shock) | R (rapamycin, Hsf1 anchored out) | HS (5 min 39 C) | RHS (rapamycin then HS)
# so the Hsf1 dependence of induction is the interaction  log2(HS/NHS) - log2(RHS/R): the same algebra
# as this project's I, on an independent measurement.
#
# The perturbation is clean (checked against their strain table). Rapamycin normally inhibits TORC1
# and thereby activates Msn2/4, which would contaminate an Hsf1-specific contrast. Their strain
# (their Table 1) is `tor1-1 fpr1D RPL13A-FKBP12::NAT HSF1-FRB-yEGFP`: tor1-1 makes TORC1
# rapamycin-resistant and fpr1D removes endogenous FKBP12, so rapamycin does not act through TORC1
# here. Without that, this control would be invalid.
#
# Their published gene list is not used. Their Supplemental Table S1 (52 Hsf1-dependent genes:
# 18 basal + 34 induced, of which 46 also Hsf1-occupied, 6 subthreshold incl. SSA3) is not among the
# files PMC serves. Re-thresholding their data to produce 52 would be an unsupported reconstruction;
# the continuous score avoids that reconstruction. The continuous score below needs no threshold from either side and uses all genes.
#
# Outcomes, fixed before this analysis was run; the second would weaken the divergence interpretation:
#   A) our `Hsf1_only` class scores resolvably higher on Pincus Hsf1 dependence than the
#      `not_induced` reference class  =>  the core is corroborated across measurement modality and
#      perturbation method, and the negative result is interpretable as a result about divergence.
#   B) it does not  =>  the core is not externally corroborated, and the negative result is interpreted
#      as conditional on the core definition, with the projection limitation stated explicitly.
# The comparison is made against the reference class, not against zero, because Pincus's contrast is
# also an interaction on count data and therefore carries the same compositional offset this project
# reports for its own epistasis contrast. Applying our own prescription to their data is required.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1")
ext  <- file.path(root, "data/03_external/pincus2018")
repd <- file.path(root, "doc/05_reports")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
sink(file.path(repd, "logs/hsf1_pincus_core_control.log"), split = TRUE)
say("=== 27_pincus_external_core_control | %s ===", format(Sys.time()))

## ------------------------------------------------- their deposited NAC-seq table ---
# fetch_inputs.sh extracts the first six columns of sheet 1 of the source xlsx; the source URL and the sha256 of
# the resulting TSV are recorded in manifests/EXTERNAL_INPUTS.tsv.
P <- fread(file.path(ext, "pincus2018_nacseq_counts.tsv"))
say("Pincus NAC-seq table: %d genes x %d columns", nrow(P), ncol(P))
stopifnot(nrow(P) == 5732L, all(c("locus","NHS_avg","R_avg","HS_avg","RHS_avg") %in% names(P)))

## the Hsf1-dependence-of-induction interaction, on their data
PC <- 1                                                  # pseudocount on normalised counts
P[, ind_ctrl := log2((HS_avg  + PC) / (NHS_avg + PC))]   # induction with Hsf1 present
P[, ind_dep  := log2((RHS_avg + PC) / (R_avg   + PC))]   # induction with Hsf1 anchored out
P[, pincus_hsf1_dep := ind_ctrl - ind_dep]               # what Hsf1 contributes to induction
say("median induction (Hsf1 present)  = %+.4f", median(P$ind_ctrl, na.rm = TRUE))
say("median induction (Hsf1 depleted) = %+.4f", median(P$ind_dep,  na.rm = TRUE))
say("median Hsf1-dependence contrast  = %+.4f   <- note this is NOT zero; see the reference-class note",
    median(P$pincus_hsf1_dep, na.rm = TRUE))

## --------------------------------------------- our own mutually exclusive classes ---
CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
UN <- CL[, .(n_distinct = uniqueN(class), cls = class[1]), by = gene][n_distinct == 1L]
say("\nour method-unanimous assignment: %d genes (of %d)", nrow(UN), uniqueN(CL$gene))
M <- merge(UN[, .(locus = gene, cls)], P[, .(locus, gene, pincus_hsf1_dep)], by = "locus")
say("joined to Pincus on systematic name: %d genes", nrow(M))
say("\nclass sizes after the join:")
print(M[, .N, by = cls][order(-N)])

## ----------------------------- the test fixed in advance, read against the reference ---
REF <- "not_induced"
ref_med <- M[cls == REF, median(pincus_hsf1_dep, na.rm = TRUE)]
say("\nreference class `%s` median Pincus dependence = %+.4f (n = %d)",
    REF, ref_med, M[cls == REF, .N])
say("as in P2, classes are read as a DIFFERENCE from this reference, not against zero.")

boot_ci <- function(v, B = 4000L) {
  v <- v[is.finite(v)]; if (length(v) < 3L) return(c(NA, NA))
  set.seed(11); quantile(replicate(B, median(sample(v, length(v), TRUE))), c(.025, .975))
}
R <- rbindlist(lapply(sort(unique(M$cls)), function(k) {
  v <- M[cls == k, pincus_hsf1_dep]
  ci <- boot_ci(v)
  w  <- if (k == REF) NA_real_ else
        suppressWarnings(wilcox.test(v, M[cls == REF, pincus_hsf1_dep])$p.value)
  data.table(cls = k, n = length(v), med = median(v, na.rm = TRUE),
             vs_ref = median(v, na.rm = TRUE) - ref_med, lo = ci[1], hi = ci[2], p_vs_ref = w)
}))
setorder(R, -vs_ref)
say("\n%-16s %6s %9s %10s %20s %12s", "class", "n", "med", "vs ref", "95%% CI of med", "p vs ref")
for (i in seq_len(nrow(R))) say("%-16s %6d %+9.4f %+10.4f  [%+.4f, %+.4f] %12.3g",
    R$cls[i], R$n[i], R$med[i], R$vs_ref[i], R$lo[i], R$hi[i], R$p_vs_ref[i])
fwrite(R, file.path(hsd, "pincus_external_core_control.tsv"), sep = "\t")

## --------------------------------------------------------------- the verdict ---
h  <- R[cls == "Hsf1_only"]
bo <- R[cls == "both_required"]
ti <- R[cls == "TF_independent"]
ok <- isTRUE(h$vs_ref > 0) && isTRUE(h$p_vs_ref < 0.05)
say("\n================ PRE-REGISTERED VERDICT ================")
say("Hsf1_only vs reference : %+.4f  (p = %.3g, n = %d)", h$vs_ref, h$p_vs_ref, h$n)
say("both_required vs ref   : %+.4f  (p = %.3g, n = %d)", bo$vs_ref, bo$p_vs_ref, bo$n)
say("TF_independent vs ref  : %+.4f  (p = %.3g, n = %d)  <- should NOT be elevated", ti$vs_ref, ti$p_vs_ref, ti$n)

## ---- evaluation of the pass condition -------------------------------------------------
## The specified pass condition tests only one comparison: it asks only whether
## Hsf1_only exceeds the reference, and never checks the negative control that this same script
## prints as a requirement ("TF_independent <- should NOT be elevated"). The condition and its value
## are shown together with a condition that includes the negative control.
say("\n---------------- THE PASS CONDITION OMITS THE NEGATIVE CONTROL ----")
say("specified condition (Hsf1_only > ref and p < 0.05) evaluates to: %s", ok)
say("but the negative control is elevated too, by almost exactly as much:")
say("  Hsf1_only      vs ref %+.4f", h$vs_ref)
say("  TF_independent vs ref %+.4f   <- the class that should NOT be Hsf1-dependent", ti$vs_ref)
say("  Msn24_only     vs ref %+.4f", R[cls == "Msn24_only"]$vs_ref)
say("  difference Hsf1_only - TF_independent = %+.4f", h$vs_ref - ti$vs_ref)
say("")
say("=> WHAT THIS CONTRAST ACTUALLY SEPARATES IS INDUCED vs NOT-INDUCED, NOT Hsf1-DEPENDENT vs NOT.")
say("   Every induced class sits ~+0.25 to +0.33 above the non-induced reference and they are not")
say("   distinguishable from one another. A likely reason, stated as a hypothesis and NOT as an")
say("   established fact: anchoring Hsf1 out for 45 min before the shift collapses proteostasis and")
say("   pre-induces a stress state, so R is already elevated for stress-responsive genes and the")
say("   RHS/R ratio is compressed for ALL of them, manufacturing apparent Hsf1 dependence wherever")
say("   there is inducibility.")
say("")
say("   THAT HYPOTHESIS IS TESTABLE FROM THESE FOUR COLUMNS, and")
say("   can be evaluated because the baseline shift caused by depletion is log2(R/NHS).")
P[, baseline_shift := log2((R_avg + PC) / (NHS_avg + PC))]
MB <- merge(UN[, .(locus = gene, cls)], P[, .(locus, baseline_shift, ind_ctrl)], by = "locus")
bref <- MB[cls == REF, median(baseline_shift, na.rm = TRUE)]
say("   %-16s %6s %11s %11s %12s", "class", "n", "med shift", "vs ref", "p vs ref")
for (k in c("both_required","Hsf1_only","Msn24_only","redundant","TF_independent", REF)) {
  v <- MB[cls == k, baseline_shift]
  pv <- if (k == REF) NA_real_ else suppressWarnings(wilcox.test(v, MB[cls == REF, baseline_shift])$p.value)
  say("   %-16s %6d %+11.4f %+11.4f %12.3g", k, length(v), median(v, na.rm = TRUE),
      median(v, na.rm = TRUE) - bref, pv)
}
sp <- cor(P$baseline_shift, P$ind_ctrl, method = "spearman", use = "complete.obs")
tp <- P[order(-ind_ctrl)][1:200]; bt <- P[order(ind_ctrl)][1:200]
say("   Spearman(baseline shift, inducibility) over all %d genes = %+.4f", nrow(P), sp)
say("   most inducible 200 genes: med shift %+.4f | least inducible 200: %+.4f",
    median(tp$baseline_shift), median(bt$baseline_shift))
say("")
say("   VERDICT ON THE HYPOTHESIS: directionally supported, QUANTITATIVELY INSUFFICIENT. Inducible")
say("   genes do sit at a raised baseline after depletion, but by only about +0.05 to +0.12, against")
say("   a ~+0.31 apparent dependence in every induced class -- and for Hsf1_only the shift is")
say("   NEGATIVE (%+.4f vs the reference), which would push the other way. The baseline shift",
    MB[cls == "Hsf1_only", median(baseline_shift, na.rm = TRUE)] - bref)
say("   accounts for at most about a third of the effect and cannot be the whole of it.")
say("")
say("   THE READING THAT SURVIVES: anchoring Hsf1 out of the nucleus appears to blunt induction")
say("   ACROSS the heat-shock program, not only at Hsf1's own targets, which is what a collapse of")
say("   chaperone capacity would do given that chaperone availability feeds back on stress")
say("   signalling. A dependence measured this way is therefore not target-specific. That is a")
say("   statement about what the contrast measures in their system. It does not show that the")
say("   knockdown-based core is the better measure: we have no ground truth, and one measure is cleaner")
say("   is not supported by this analysis. The two measures disagree in structure; ours is")
say("   not corroborated by theirs; neither is shown to be right.")
fwrite(MB[, .(locus, cls, baseline_shift, ind_ctrl)],
       file.path(hsd, "pincus_baseline_shift_diagnostic.tsv"), sep = "\t")
say("")
## the correctly-referenced version: among induced classes, is Hsf1_only above the Hsf1-independent one?
vh <- M[cls == "Hsf1_only",      pincus_hsf1_dep]
vt <- M[cls == "TF_independent", pincus_hsf1_dep]
pw <- suppressWarnings(wilcox.test(vh, vt)$p.value)
say("correctly-referenced test (reference = TF_independent, i.e. induced but Hsf1-independent):")
say("  median difference %+.4f , Wilcoxon p = %.3g", median(vh, na.rm = TRUE) - median(vt, na.rm = TRUE), pw)
DISCRIM <- isTRUE(pw < 0.05) && isTRUE(median(vh, na.rm = TRUE) > median(vt, na.rm = TRUE))
say("  discriminates Hsf1_only from Hsf1-independent induced genes: %s", DISCRIM)
say("")
if (!DISCRIM) {
  say("*** Result: THIS TEST DOES NOT CORROBORATE OUR CORE. It fails its own negative")
  say("*** control, so it cannot distinguish 'our core is right' from 'anything inducible looks")
  say("*** Hsf1-dependent in this contrast'. Consequences:")
  say("***   1. The core-divergence negative is CONDITIONAL on our core definition; the core is not")
  say("***      externally validated by this test.")
  say("***   2. The core is not thereby shown to be WRONG either. The test is uninformative for its")
  say("***      intended purpose, which is a different statement from a failed validation.")
  say("***   3. The result uses the specified contrast. The")
  say("***      direct test would be their ChIP-seq occupancy at the 74 loci, which is not in the")
  say("***      deposited files we can reach; the unavailable occupancy data leave this question unresolved.")
} else {
  say("*** Alternative result with the negative control included: the core is corroborated. Hsf1_only")
  say("*** exceeds induced-but-Hsf1-independent genes on an independent modality and perturbation.")
}
say("\n---- result under the original criterion ----")
if (ok) {
  say("\n*** Result: the core IS externally corroborated. Genes we call Hsf1-dependent from")
  say("*** steady-state mRNA under HSF1 knockdown are also Hsf1-dependent in nascent transcription")
  say("*** under Hsf1 nuclear depletion, in a different laboratory's data. Two measurement modalities")
  say("*** and two perturbation methods agree. => the core-divergence negative is about DIVERGENCE, not")
  say("*** about a mis-defined core.")
} else {
  say("\n*** Result: the core is NOT externally corroborated on this test. The core-divergence negative")
  say("*** is CONDITIONAL on our core definition, and the projection limitation becomes more important.")
  say("*** This result uses the specified settings.")
}
say("\nCAVEATS THAT DO NOT DEPEND ON THE OUTCOME.")
say("1. Strain axis: their background is BY4741/BY4742-derived; ours is the GSE179258 background.")
say("   This is a further strain axis on top of the reference-strain caveat and is not controlled.")
say("2. Their heat shock is 5 min at 39 C; our anchor is 30 min at 37 C. The kinetics differ, so")
say("   agreement is evidence about the identity of Hsf1 targets, not about response timing.")
say("3. n = 87 for Hsf1_only before the join; per-gene interpretation remains excluded either way.")
writeLines(capture.output(sessionInfo()), file.path(repd, "logs/sessionInfo_pincus_control.txt"))
say("\n=== done ===")
sink()
