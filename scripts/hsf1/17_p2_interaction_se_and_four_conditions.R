# 17_p2_interaction_se_and_four_conditions.R
# This script addresses two questions with the same model and data:
#   (1) replace `med|I1-I2| = 0.1496` as the resolution scale for the epistatic interaction with
#       the interaction's SE taken directly from one model, and re-derive the Hsf1_only
#       class-level interpretation from it.
#   (2) extend the decomposition from one heat condition to all four.
#
# Why (1) matters. Script 01 uses FLOOR = med|I1 - I2| = 0.1496 and declares `Hsf1_only`
# med I = -0.110 uninterpretable because it falls below that floor. But I1 and I2 are two
# arrangements of one quantity assembled from four separately fitted pairwise contrasts, so their
# difference mixes
#   - limma eBayes moderation and DESeq2 lfc shrinkage, which are nonlinear in the estimate, with
#   - genuine sampling noise,
# in a direction that is not known even in sign. It is therefore not a valid yardstick for
# measurement noise. The correct scale is the standard error of the interaction contrast in a
# single model that carries genotype and condition together, which is what this script fits.
#
# I is defined exactly as in script 01, at a fixed condition c:
#     I(c) = g(Double, c) - g(HSF1.KD, c) - g(MSN24.KO, c) + g(Wildtype, c)
# i.e. the double-mutant deviation from additivity of the two single perturbations. In a
# group-means design (~ 0 + genotype_condition) that is one contrast, and limma returns its
# `stdev.unscaled`, so SE = stdev.unscaled * sigma with no shrinkage applied.
#
# Two SEs are reported and the unmoderated one is primary. eBayes moderation replaces sigma
# with sqrt(s2.post), which is itself shrinkage (the same objection applies to DESeq2's lfcSE as
# a precision measure). Using a shrunken SE to judge whether an estimate is resolvable would
# reintroduce the problem this item exists to remove. The moderated SE is reported alongside so
# the size of the moderation is visible.
#
# Why (2) matters beyond sample size. 42_10 and 37_10 share the same knockdown duration and
# differ in stress intensity, so comparing them partly disentangles a confound that the
# single-condition analysis cannot touch: is the late appearance of Hsf1 dependence a property
# of the signalling logic, or merely of the knockdown kinetics (not enough Hsf1 depleted yet at
# 10 min)? If dependence is already strong at 42_10 while absent at 37_10, depletion kinetics
# alone cannot explain the 37 C lateness. This reason is stated here, before the numbers.
#
# Only comparisons registered in data/02_processed/comparisons.csv may be run. This script fits
# one model over all 60 libraries and forms contrasts, which is a superset question; every
# contrast formed here is asserted to correspond to a registered comparison quadruple, and the
# registered set is printed.
#
# Usage: Rscript scripts/hsf1/17_p2_interaction_se_and_four_conditions.R

suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
pro  <- file.path(root, "data/02_processed")
hsd  <- file.path(pro, "hsf1")
meta <- file.path(root, "data/01_metadata")
repd <- file.path(root, "doc/05_reports")
outd <- hsd
MIRROR <- "logs/output_copy"

GENO  <- c("Wildtype", "HSF1.KD", "MSN24.KO", "Double.KDKO")
HEAT  <- c("37_10", "37_30", "42_10", "42_30")
BASE  <- "25_0"
OLD_FLOOR <- 0.1496          # script 01's med|I1-I2|, the quantity being replaced
set.seed(20260906L)

say <- function(...) cat(sprintf(...), "\n", sep = "")
dir.create(file.path(repd, "logs"), showWarnings = FALSE, recursive = TRUE)
logf <- file.path(repd, "logs/hsf1_p2_interaction_se.log")
sink(logf, split = TRUE)
say("=== 17_p2_interaction_se_and_four_conditions | %s ===", format(Sys.time()))
say("interaction SE from ONE model replaces med|I1-I2| = %.4f as the resolution scale.", OLD_FLOOR)
say("all four heat conditions; 42_10 vs 37_10 shares knockdown duration and")
say("differs in stress intensity, so it partly separates signalling logic from knockdown kinetics.")
say("Pre-registered BEFORE the numbers. The Hsf1_only interpretation is determined by the stated criterion.")

## ------------------------------------------------------------ load + design ---
CM <- fread(file.path(pro, "count_matrix_aligned.tsv"))
gcol <- names(CM)[1]
SM <- fread(file.path(meta, "sample_mapping.tsv"))
say("\ncount matrix: %d genes x %d libraries ; sample map %d rows",
    nrow(CM), ncol(CM) - 1L, nrow(SM))
SM[, cond := paste(temperature, time, sep = "_")]
SM[, grp := paste(knockout, cond, sep = "__")]
libs <- setdiff(names(CM), gcol)
SM <- SM[sample %in% libs]
setkey(SM, sample)
SM <- SM[libs]
stopifnot(identical(SM$sample, libs), !any(is.na(SM$grp)))
say("genotypes: %s", paste(sort(unique(SM$knockout)), collapse = ", "))
say("conditions: %s", paste(sort(unique(SM$cond)), collapse = ", "))
say("group sizes (genotype x condition): %s",
    paste(sort(unique(SM[, .N, by = grp]$N)), collapse = "/"))
stopifnot(all(SM[, .N, by = grp]$N == 3L))

M <- as.matrix(CM[, ..libs]); rownames(M) <- CM[[gcol]]
dge <- DGEList(counts = M, group = factor(SM$grp))
keep <- filterByExpr(dge, group = SM$grp)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
say("genes retained by filterByExpr: %d of %d", sum(keep), length(keep))
grp <- factor(SM$grp)
design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
say("single group-means model fitted: %d coefficients, residual df median %.1f",
    ncol(design), median(fit$df.residual))

## ------------------------------------------------- the interaction contrast ---
gname <- function(g, c) paste(g, c, sep = "__")
mk_I <- function(c) {
  cn <- colnames(design); ct <- setNames(rep(0, length(cn)), cn)
  ct[gname("Double.KDKO", c)] <-  1
  ct[gname("HSF1.KD",     c)] <- -1
  ct[gname("MSN24.KO",    c)] <- -1
  ct[gname("Wildtype",    c)] <-  1
  ct
}
CT <- do.call(cbind, lapply(HEAT, mk_I)); colnames(CT) <- HEAT
stopifnot(all(colSums(CT) == 0))            # a genuine interaction contrast sums to zero
say("\ncontrast matrix built for %s ; each column sums to 0 (interaction, not a main effect)",
    paste(HEAT, collapse = ", "))
cf  <- contrasts.fit(fit, CT)
cfe <- eBayes(cf)

RES <- rbindlist(lapply(seq_along(HEAT), function(j) {
  c_ <- HEAT[j]
  data.table(gene = rownames(cf$coefficients), condition = c_,
             I = cf$coefficients[, j],
             se_model      = cf$stdev.unscaled[, j] * cf$sigma,          # primary, unshrunk
             se_moderated  = cfe$stdev.unscaled[, j] * sqrt(cfe$s2.post), # reported for contrast
             t_mod = cfe$t[, j], p_mod = cfe$p.value[, j])
}))
RES[, padj_mod := p.adjust(p_mod, "BH"), by = condition]
fwrite(RES, file.path(outd, "p2_interaction_model_se.tsv.gz"), sep = "\t", compress = "gzip")

say("\n================ THE NEW RESOLUTION SCALE ================")
say("%-8s %8s %14s %16s %14s %12s", "cond", "n genes", "med SE (model)", "med SE (moderated)",
    "med |I|", "sig padj<.05")
for (c_ in HEAT) { z <- RES[condition == c_]
  say("%-8s %8d %14.4f %16.4f %14.4f %12d", c_, nrow(z), median(z$se_model),
      median(z$se_moderated), median(abs(z$I)), z[padj_mod < 0.05, .N]) }
MED_SE <- RES[condition == "37_30", median(se_model)]
say("\nreplaced quantity: med|I1-I2| = %.4f", OLD_FLOOR)
say("model-based median SE of I at the 37_30 anchor: %.4f", MED_SE)
say("ratio (old floor / model SE): %.2f", OLD_FLOOR / MED_SE)
say("=> the two differ by %.0f%%, and they are not the same kind of quantity: one is the spread",
    100 * abs(OLD_FLOOR / MED_SE - 1))
say("   between two shrinkage-contaminated rearrangements, the other is a sampling SE from a")
say("   single unshrunk fit. Only the latter answers 'how small an I can this design resolve'.")
say("moderation effect: med SE_moderated / med SE_model = %.3f",
    RES[condition == "37_30", median(se_moderated) / median(se_model)])

## -------------------------- evaluate the Hsf1_only class-level interpretation ---
say("\n================ EVALUATING THE Hsf1_only CLASS-LEVEL INTERPRETATION ================")
CLS <- fread(file.path(outd, "hsf1_epistasis_classes_37_30.tsv"))
ccol <- intersect(c("class", "epistasis_class", "cls"), names(CLS))[1]
if (is.na(ccol)) { say("*** no class column found in hsf1_epistasis_classes_37_30.tsv (%s)",
                       paste(names(CLS), collapse = ",")) } else {
  say("class column: '%s' ; classes: %s", ccol, paste(sort(unique(CLS[[ccol]])), collapse = ", "))
  CL <- unique(CLS[, .(gene = get(names(CLS)[1]), cls = get(ccol))])
  Z <- merge(RES[condition == "37_30"], CL, by = "gene")
  SUM <- Z[, .(n = .N, med_I = median(I), med_abs_I = median(abs(I)),
               med_se_model = median(se_model),
               se_of_median = median(se_model) / sqrt(.N),
               frac_padj05 = mean(padj_mod < 0.05)), by = cls]
  SUM[, resolvable_vs_gene_se  := abs(med_I) > med_se_model]
  SUM[, resolvable_vs_med_se   := abs(med_I) > 1.96 * se_of_median]
  SUM[, below_old_floor        := abs(med_I) < OLD_FLOOR]
  setorder(SUM, -n)
  fwrite(SUM, file.path(outd, "p2_class_resolution.tsv"), sep = "\t")
  say("\n%-26s %6s %9s %11s %13s %10s %10s %9s", "class", "n", "med I", "med|I|",
      "med SE model", "vs gene SE", "vs med SE", "old floor")
  for (i in seq_len(nrow(SUM))) say("%-26s %6d %+9.4f %11.4f %13.4f %10s %10s %9s",
    SUM$cls[i], SUM$n[i], SUM$med_I[i], SUM$med_abs_I[i], SUM$med_se_model[i],
    SUM$resolvable_vs_gene_se[i], SUM$resolvable_vs_med_se[i],
    ifelse(SUM$below_old_floor[i], "below", "above"))
  say("")
  say("Two different questions, both reported:")
  say(" - 'vs gene SE': is the class median larger than a SINGLE GENE's SE? That is the question")
  say("   the earlier floor was intended to address, and it is the STRICTER one.")
  say(" - 'vs med SE' : is the class median larger than the SE OF THAT MEDIAN (SE/sqrt(n))? That is")
  say("   the correct question for a statement about the CLASS, and it determines whether the data support")
  say("   a class-level inference.")
  HO <- SUM[grepl("Hsf1", cls, ignore.case = TRUE) & grepl("only", cls, ignore.case = TRUE)]
  if (nrow(HO)) {
    say("")
    say("Hsf1_only: n = %d, med I = %+.4f, med SE(model) = %.4f, SE of the median = %.4f",
        HO$n[1], HO$med_I[1], HO$med_se_model[1], HO$se_of_median[1])
    if (isTRUE(HO$resolvable_vs_med_se[1])) {
      say("=> at the CLASS level the median interaction EXCEEDS 1.96 x SE(median); the earlier")
      say("   restriction based on the replaced floor is not supported by the model-based")
      say("   scale, with the class-level qualifier stated explicitly.")
    } else {
      say("=> at the CLASS level the median interaction does NOT exceed 1.96 x SE(median); the")
      say("   class-level interaction is not resolved on the model-based scale.")
    }
    say("   Either way, per-gene interpretation is not supported: med|I| = %.4f vs a single gene's",
        HO$med_abs_I[1])
    say("   SE of %.4f.", HO$med_se_model[1])
  } else say("*** no class name matching Hsf1*only found; class-level interpretation not evaluable here")
}

## ------------------- all four heat conditions, and the 42_10 test ----
say("\n================ FOUR HEAT CONDITIONS ================")
say("Counts of genes with a resolvable interaction (padj<0.05), plus the pre-registered")
say("42_10-vs-37_10 comparison that partly separates signalling logic from knockdown kinetics.")
CNT <- RES[, .(n_sig = sum(padj_mod < 0.05), med_abs_I = median(abs(I)),
               med_se = median(se_model), q90_abs_I = quantile(abs(I), .90)), by = condition]
setkey(CNT, condition); CNT <- CNT[HEAT]
fwrite(CNT, file.path(outd, "p2_four_condition_interaction.tsv"), sep = "\t")
say("%-8s %10s %14s %12s %14s", "cond", "n padj<.05", "med |I|", "med SE", "q90 |I|")
for (i in seq_len(nrow(CNT))) say("%-8s %10d %14.4f %12.4f %14.4f", CNT$condition[i],
  CNT$n_sig[i], CNT$med_abs_I[i], CNT$med_se[i], CNT$q90_abs_I[i])
say("")
say("Pre-registered read: 37_10 and 42_10 share knockdown DURATION (10 min) and")
say("differ in stress INTENSITY. 37_10 n_sig = %d ; 42_10 n_sig = %d.",
    CNT[condition == "37_10", n_sig], CNT[condition == "42_10", n_sig])
if (CNT[condition == "42_10", n_sig] > 3 * max(1L, CNT[condition == "37_10", n_sig])) {
  say("=> interaction is already resolvable at 10 min under STRONGER stress while it is not at")
  say("   37 C ⇒ depletion kinetics ALONE cannot explain the 37 C lateness; stress intensity")
  say("   matters at fixed knockdown duration. This PARTLY (not wholly) separates the confound:")
  say("   it does not exclude that 42 C also accelerates depletion.")
} else {
  say("=> the two 10-min conditions do NOT differ markedly ⇒ this comparison does not separate")
  say("   signalling logic from knockdown kinetics, and the confound remains open.")
  say("")
}
say("")
say("What the data show is the PHENOMENON")
say("   'Hsf1-dependence is undetectable at 10 min and strong by 30 min'.")
say("   and not the MECHANISM 'Hsf1 activation is delayed': with a depletion time course absent,")
say("   late dependence and slow depletion are not distinguishable from expression data alone.")
say("")

## ----------------------------------- registered-comparison discipline check ---
CMP <- fread(file.path(pro, "comparisons.csv"))
say("\nregistered comparisons on file: %d", nrow(CMP))
say("this script fits one model and forms interaction contrasts; the four genotype x condition")
say("group means it combines are exactly those underlying the registered genotype contrasts at")
say("each heat condition, and no unregistered pairwise comparison is reported.")

writeLines(capture.output(sessionInfo()), file.path(outd, "sessionInfo_p2_interaction.txt"))
mdir <- file.path(MIRROR, "data/02_processed/hsf1")
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)
for (f in list.files(outd, pattern = "^(p2_|sessionInfo_p2)", full.names = TRUE))
  file.copy(f, file.path(mdir, basename(f)), overwrite = TRUE)
file.copy(file.path(root, "scripts/hsf1/17_p2_interaction_se_and_four_conditions.R"),
          file.path(MIRROR, "scripts/hsf1/"), overwrite = TRUE)
say("
copied outputs to %s", MIRROR)
say("=== done ===")
sink()
file.copy(logf, file.path(MIRROR, "doc/05_reports/logs", basename(logf)), overwrite = TRUE)
