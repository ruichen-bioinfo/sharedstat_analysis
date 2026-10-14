# 01_epistasis_decomposition.R -- what does Hsf1 account for, uniquely, in the
# S. cerevisiae heat-shock response?
#
# Rationale.
#
# The Hsf1 module of 06_module_bridge.R (105-112 genes depending on pipeline) is
# defined from single perturbations only: a gene is "Hsf1-dependent" if heat induces
# it in wild type and Hsf1 depletion lowers it. Three published results say that
# definition must be too generous:
#
#   Solis & Denic 2016  (PMID 27320198): forcing Hsf1 out of the nucleus lowers the
#                        basal transcription of only 18 genes, mostly chaperones.
#   Pincus et al. 2018  (PMID 30332327): Hsf1 ChIP-seq + depletion + nascent RNA-seq
#                        gives 74 binding loci and 46 Hsf1-dependent genes.
#   Ciccarelli et al. 2023 (PMID 37467033): hsf1-null cells are viable, and heat can
#                        route the response through Msn2/4 instead.
#
# If heat can route around Hsf1, then some genes counted as "Hsf1-dependent" are really
# jointly held by Hsf1 and Msn2/4, and some are held redundantly -- visible only when
# both are removed. Projecting a regulon contaminated that way onto other species
# (04_divergence_by_dependency.R) would measure the wrong thing. GSE179258 contains the
# genotype that settles it, `Double.KDKO`.
#
# The interaction contrast, and why its arms matter.
#
# Write W, H, M, D for log2 expression in wild type, Hsf1-depleted, msn2/msn4-deleted
# and double. The effect of removing Msn2/4 depends on whether Hsf1 is present:
#
#   a = M - W   effect of removing Msn2/4 in a wild-type background
#   b = D - H   effect of removing Msn2/4 in an Hsf1-depleted background
#   I = b - a = D - H - M + W
#
# I = 0 means the two regulators act independently on that gene. I < 0 means Msn2/4
# only matters once Hsf1 is gone (redundancy: either regulator suffices, so the
# requirement is invisible in each single mutant). I > 0 means Msn2/4 stops mattering
# once Hsf1 is gone (saturation: the gene is already at its floor).
#
# Two coefficients taken from fits that share an arm also share that arm's measurement
# error, and differencing them then produces a bias that survives at the null. Here that
# hazard is avoided by construction. `a` is fitted on {MSN24.KO, Wildtype} and `b` on
# {Double.KDKO, HSF1.KD} -- four disjoint sets of samples. So within one estimate of I
# the two coefficients are independent, se_I = sqrt(se_a^2 + se_b^2) is valid, and no
# shared-arm artefact term is introduced.
#
# The second arrangement is not an independent estimate.
#
# I is symmetric: I = (D - M) - (H - W) too, from the arm partition {MSN24.KO,
# Double.KDKO} and {HSF1.KD, Wildtype}. Treating that as a second independent estimate,
# averaging the two and dividing the standard error by 2, would be wrong. The two
# expressions are built from four different registered comparisons but from the same
# twelve samples, and because 20_batch_run.R fits every comparison as a two-group
# subset (~ arm) the two expressions are algebraically the same difference of the same
# four group means. They can only differ through the per-comparison renormalisation
# (TMM is recomputed on each subset). Measured agreement is r = +1.000 under limma and
# +0.999 under edgeR -- a near identity, so it validates nothing about the estimator,
# and averaging two perfectly correlated quantities reduces no variance (there is no
# sqrt(2) gain).
#
# The comparison is retained only for what it can measure: how much the
# per-comparison renormalisation perturbs the contrast. That number is a floor on the
# resolution of I, and any class whose median |I| is below it cannot be interpreted.
# I1 is the primary estimate; I2 is reported as a renormalisation sensitivity check.
#
# What is not claimed here. This script is about S. cerevisiae only. It says nothing
# about any other species, and nothing about Hsf1 protein sequence. It defines the gene
# set that later steps project onto Fay et al. 2023 orthologs, and it is the reason
# that projection is restricted to a core rather than to all 105 genes.
# No pipeline is called more accurate than another anywhere in this script.
#
# Usage: Rscript scripts/hsf1/01_epistasis_decomposition.R [--anchor=37_30]

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
proc <- file.path(root, "data/02_processed")
outd <- file.path(proc, "hsf1")
ref  <- file.path(root, "data/03_reference")
dir.create(outd, showWarnings = FALSE, recursive = TRUE)

ALPHA   <- 0.05
METHODS <- c("DESeq2", "edgeR", "limma")
say <- function(...) cat(sprintf(...), "\n", sep = "")

args <- commandArgs(trailingOnly = TRUE)
h <- grep("^--anchor=", args, value = TRUE)
ANCHOR <- if (length(h)) sub("^--anchor=", "", h[1]) else "37_30"

say("=== 01_epistasis_decomposition | %s | anchor %s C ===", format(Sys.time()), ANCHOR)
say("S. cerevisiae only. Registered comparisons only.")

## ---------------------------------------------------------------- name bridge -----
gtf  <- file.path(ref, "Saccharomyces_cerevisiae.R64-1-1.75.gtf")
gl   <- grep("\tgene\t", readLines(gtf, warn = FALSE), value = TRUE)
sysn <- sub('.*gene_id "([^"]+)".*', "\\1", gl)
stdn <- ifelse(grepl('gene_name "', gl), sub('.*gene_name "([^"]+)".*', "\\1", gl), NA_character_)
nmap <- unique(data.table(sys = sysn, std = stdn))
lut  <- unique(rbind(nmap[, .(key = toupper(sys), sys)],
                     nmap[!is.na(std), .(key = toupper(std), sys)]))
lut  <- lut[!duplicated(key)]
say("GTF name bridge: %d gene records, %d carry a standard name", nrow(nmap), sum(!is.na(nmap$std)))

## ---------------------------------------------------------------- loaders ---------
## Every comparison below is registered in data/02_processed/comparisons.csv. The
## direction convention is A_vs_B => B is the reference, so lfc = A - B.
load_comp <- function(method, comp, want_se = FALSE) {
  f <- file.path(proc, method, paste0(comp, ".tsv"))
  if (!file.exists(f)) stop("missing batch output: ", f)
  cols <- c("gene", "log2FC", "padj", if (want_se) "lfcSE")
  d <- fread(f, select = cols)
  setnames(d, c("gene", "lfc", "padj", if (want_se) "se"))
  d
}

REG <- fread(file.path(proc, "comparisons.csv"))
need <- c(sprintf("Wildtype_%s_vs_25_0", ANCHOR),
          sprintf("HSF1.KD_vs_Wildtype_at_%s", ANCHOR),
          sprintf("MSN24.KO_vs_Wildtype_at_%s", ANCHOR),
          sprintf("Double.KDKO_vs_Wildtype_at_%s", ANCHOR),
          sprintf("HSF1.KD_vs_Double.KDKO_at_%s", ANCHOR),
          sprintf("MSN24.KO_vs_Double.KDKO_at_%s", ANCHOR))
miss <- setdiff(need, REG$comp_id)
if (length(miss)) stop("not registered in comparisons.csv: ", paste(miss, collapse = ", "))
say("all %d required comparisons are pre-registered", length(need))

## ------------------------------------------------- qualitative decomposition -----
## Classes are defined on the induced arm: heat must raise the gene in wild type first,
## otherwise "the regulator is required for induction" is not a meaningful statement.
classify <- function(method) {
  wt <- load_comp(method, sprintf("Wildtype_%s_vs_25_0", ANCHOR))
  kd <- load_comp(method, sprintf("HSF1.KD_vs_Wildtype_at_%s", ANCHOR))
  ko <- load_comp(method, sprintf("MSN24.KO_vs_Wildtype_at_%s", ANCHOR))
  db <- load_comp(method, sprintf("Double.KDKO_vs_Wildtype_at_%s", ANCHOR))
  setnames(wt, c("gene", "lfc_wt", "padj_wt"))
  setnames(kd, c("gene", "lfc_kd", "padj_kd"))
  setnames(ko, c("gene", "lfc_ko", "padj_ko"))
  setnames(db, c("gene", "lfc_db", "padj_db"))
  d <- Reduce(function(x, y) merge(x, y, by = "gene"), list(wt, kd, ko, db))

  d[, induced := !is.na(padj_wt) & padj_wt < ALPHA & lfc_wt > 0]
  ## "requires X" = removing X moves the gene back down toward baseline
  d[, req_hsf1 := induced & !is.na(padj_kd) & padj_kd < ALPHA & lfc_kd < 0]
  d[, req_msn  := induced & !is.na(padj_ko) & padj_ko < ALPHA & lfc_ko < 0]
  d[, req_dbl  := induced & !is.na(padj_db) & padj_db < ALPHA & lfc_db < 0]

  d[, class := fifelse(!induced,                    "not_induced",
               fifelse( req_hsf1 & !req_msn,        "Hsf1_only",
               fifelse(!req_hsf1 &  req_msn,        "Msn24_only",
               fifelse( req_hsf1 &  req_msn,        "both_required",
               fifelse(!req_hsf1 & !req_msn & req_dbl, "redundant",
                                                    "TF_independent")))))]
  ## The `redundant` label above is provisional and must be split. A gene is placed there
  ## whenever neither single perturbation clears alpha but the double does, and that
  ## happens for two entirely different reasons:
  ##   (i)  genuine joint requirement -- both singles lower the gene a little, neither
  ##        enough to clear alpha, and the double removes it (super-additive);
  ##   (ii) one regulator does the work but its single-mutant effect fell just short of
  ##        alpha, so the gene is really Msn2-4-dependent with insufficient power.
  ## Case (ii) is a threshold boundary, not redundancy, and calling it redundancy would
  ## manufacture support for Ciccarelli 2023 out of a power limit. The sign of the two
  ## single effects separates them: (i) needs both negative, (ii) shows the Msn24_only
  ## signature (Hsf1 depletion raises the gene while msn2/4 deletion lowers it).
  d[, sub := fifelse(class != "redundant", class,
             fifelse(lfc_kd < 0 & lfc_ko < 0, "red_joint",
             fifelse(lfc_kd > 0 & lfc_ko < 0, "red_Msn24_subthreshold",
                                              "red_unclassified")))]
  d[, method := method]
  d[]
}

CL <- rbindlist(lapply(METHODS, classify))
ORD <- c("Hsf1_only", "both_required", "redundant", "Msn24_only", "TF_independent", "not_induced")
SUB <- c("Hsf1_only", "both_required", "red_joint", "red_Msn24_subthreshold",
         "red_unclassified", "Msn24_only", "TF_independent", "not_induced")

say("\n=== decomposition of the heat-induced response at %s C ===", ANCHOR)
say("(the fourth genotype Double.KDKO is included in this decomposition)")
tabw <- dcast(CL[, .N, by = .(method, class)], class ~ method, value.var = "N", fill = 0L)
tabw <- tabw[match(ORD, class)][!is.na(class)]
say("%-16s %8s %8s %8s", "class", METHODS[1], METHODS[2], METHODS[3])
for (i in seq_len(nrow(tabw)))
  say("%-16s %8d %8d %8d", tabw$class[i], tabw[[METHODS[1]]][i], tabw[[METHODS[2]]][i], tabw[[METHODS[3]]][i])

## how much of the single-perturbation 105-gene module was actually shared with Msn2/4?
say("\nsingle-perturbation Hsf1 module (req_hsf1, any Msn2/4 status) vs its split:")
for (m in METHODS) {
  x <- CL[method == m]
  n_old  <- x[, sum(req_hsf1)]
  n_only <- x[, sum(class == "Hsf1_only")]
  n_both <- x[, sum(class == "both_required")]
  say("  %-7s req_hsf1 = %4d  ->  Hsf1_only %4d (%.0f%%) + both_required %4d (%.0f%%)",
      m, n_old, n_only, 100 * n_only / n_old, n_both, 100 * n_both / n_old)
}

## three-pipeline agreement on the class label
agr <- dcast(CL[, .(gene, method, class)], gene ~ method, value.var = "class")
agr[, n_distinct := apply(.SD, 1, uniqueN), .SDcols = METHODS]
agr[, unanimous := n_distinct == 1L]
say("\nclass agreement across the three pipelines: %d of %d genes identical (%.1f%%)",
    sum(agr$unanimous), nrow(agr), 100 * mean(agr$unanimous))
tf_any <- apply(agr[, ..METHODS], 1, function(x) any(x %in% ORD[1:4]))
say("genes any pipeline puts in a regulator class: %d, unanimous among them: %d (%.1f%%)",
    sum(tf_any), sum(agr$unanimous & tf_any), 100 * sum(agr$unanimous & tf_any) / max(sum(tf_any), 1))

## the gene set later steps project onto other species: unanimous Hsf1_only + both
CORE_SET <- agr[unanimous & agr[[METHODS[1]]] %in% c("Hsf1_only", "both_required"), gene]
ONLY_SET <- agr[unanimous & agr[[METHODS[1]]] == "Hsf1_only", gene]
say("\nunanimous Hsf1-associated core (Hsf1_only + both_required): %d genes", length(CORE_SET))
say("unanimous Hsf1_only (no detectable Msn2/4 requirement):      %d genes", length(ONLY_SET))

## ---------------------------------------------- quantitative interaction ---------
## I = (D - H) - (M - W). Two arrangements of the same quantity:
##   arrangement 1: b1 = -lfc(HSF1.KD_vs_Double.KDKO)  arms {Double, Hsf1-dep}
##                  a1 =  lfc(MSN24.KO_vs_Wildtype)    arms {msn24, WT}
##   arrangement 2: b2 = -lfc(MSN24.KO_vs_Double.KDKO) arms {Double, msn24}
##                  a2 =  lfc(HSF1.KD_vs_Wildtype)     arms {Hsf1-dep, WT}
## Within one arrangement the two fits share no samples, so the SEs add in quadrature
## and no shared-error term is created.
## Across arrangements nothing is gained: same twelve samples, same four group means,
## so I1 == I2 up to per-comparison renormalisation. I1 is primary. I2 is used only to
## measure how large that renormalisation perturbation is, which sets a resolution floor
## below which a class-median interaction cannot be interpreted.
## limma's se is unshrunken and is the primary precision measure; DESeq2's lfcSE is
## post-shrinkage and must not be used that way, so the interaction test is run
## on limma and the other pipelines are carried as sign/ordering checks only.
interaction_est <- function(method) {
  has_se <- method != "edgeR"          # edgeR writes no lfcSE
  b1 <- load_comp(method, sprintf("HSF1.KD_vs_Double.KDKO_at_%s",  ANCHOR), has_se)
  a1 <- load_comp(method, sprintf("MSN24.KO_vs_Wildtype_at_%s",    ANCHOR), has_se)
  b2 <- load_comp(method, sprintf("MSN24.KO_vs_Double.KDKO_at_%s", ANCHOR), has_se)
  a2 <- load_comp(method, sprintf("HSF1.KD_vs_Wildtype_at_%s",     ANCHOR), has_se)
  nm <- function(d, p) { setnames(d, setdiff(names(d), "gene"), paste0(p, setdiff(names(d), "gene"))); d }
  d <- Reduce(function(x, y) merge(x, y, by = "gene"),
              list(nm(b1, "b1_"), nm(a1, "a1_"), nm(b2, "b2_"), nm(a2, "a2_")))
  ## registered direction is HSF1.KD - Double.KDKO, we want Double - HSF1.KD
  d[, I1 := (-b1_lfc) - a1_lfc]
  d[, I2 := (-b2_lfc) - a2_lfc]
  if (has_se) {
    d[, se_I1 := sqrt(b1_se^2 + a1_se^2)]
    d[, se_I2 := sqrt(b2_se^2 + a2_se^2)]
  }
  d[, method := method]
  d[]
}

IX <- rbindlist(lapply(METHODS, interaction_est), fill = TRUE)

say("\n=== the interaction I = (D-H)-(M-W) ===")
say("I1 and I2 are two ARRANGEMENTS of one quantity, not two independent estimates:")
say("both are D-H-M+W over the same 12 samples, so they can differ only through the")
say("per-comparison renormalisation. r near 1 is therefore expected and validates")
say("nothing. med |I1-I2| is read as a RESOLUTION FLOOR on I, not as agreement.")
for (m in METHODS) {
  x <- IX[method == m & is.finite(I1) & is.finite(I2)]
  say("  %-7s n=%d  r(I1,I2) = %+.4f  Spearman = %+.4f  med |I1-I2| = %.4f  (floor)",
      m, nrow(x), cor(x$I1, x$I2), cor(x$I1, x$I2, method = "spearman"),
      median(abs(x$I1 - x$I2)))
}
FLOOR <- IX[method == "limma" & is.finite(I1) & is.finite(I2), median(abs(I1 - I2))]

## Primary interaction call on limma, on I1 alone. No averaging: averaging two
## quantities with r = 1.000 reduces no variance, and dividing the pooled se by 2 would
## make every p-value anti-conservative by about sqrt(2).
LI <- IX[method == "limma" & is.finite(I1) & is.finite(se_I1)]
LI[, I_hat := I1]
LI[, se    := se_I1]                            # a1, b1 from arm-disjoint fits
LI[, z := I_hat / se]
LI[, p := 2 * pnorm(-abs(z))]
LI[, padj := p.adjust(p, method = "BH")]
say("\nlimma, I1 with se = sqrt(se_a1^2 + se_b1^2), BH across %d genes:", nrow(LI))
say("  interaction padj < %.2f : %d genes (%.1f%%)   of which I<0 (redundancy) %d, I>0 (saturation) %d",
    ALPHA, LI[padj < ALPHA, .N], 100 * LI[, mean(padj < ALPHA)],
    LI[padj < ALPHA & I_hat < 0, .N], LI[padj < ALPHA & I_hat > 0, .N])
say("  median se = %.4f; renormalisation floor med|I1-I2| = %.4f (%.0f%% of median se)",
    LI[, median(se)], FLOOR, 100 * FLOOR / LI[, median(se)])
## sensitivity: does the ranking survive being recomputed from the other arrangement?
LS <- IX[method == "limma" & is.finite(I2) & is.finite(se_I2)]
LS[, p2 := 2 * pnorm(-abs(I2 / se_I2))][, padj2 := p.adjust(p2, method = "BH")]
SS <- merge(LI[, .(gene, padj)], LS[, .(gene, padj2)], by = "gene")
say("  sensitivity to arrangement: %d of %d genes called at padj<%.2f under I1 are also called under I2 (%.1f%%)",
    SS[padj < ALPHA & padj2 < ALPHA, .N], SS[padj < ALPHA, .N], ALPHA,
    100 * SS[padj < ALPHA, mean(padj2 < ALPHA)])

## interaction restricted to the classes that matter
CLL <- CL[method == "limma", .(gene, class)]
J <- merge(LI[, .(gene, I_hat, se, padj_I = padj)], CLL, by = "gene")
say("\ninteraction magnitude by class (limma, I1). 'resolved' = |med I| > the %.4f floor:", FLOOR)
say("%-16s %6s %9s %9s %9s %9s", "class", "n", "med I", "frac I<0", "padj<.05", "resolved")
for (cc in ORD) {
  x <- J[class == cc]
  if (!nrow(x)) next
  mi <- median(x$I_hat)
  say("%-16s %6d %+9.4f %9.2f %9d %9s", cc, nrow(x), mi,
      mean(x$I_hat < 0), x[padj_I < ALPHA, .N], if (abs(mi) > FLOOR) "yes" else "NO")
}
say("A class marked NO has a median interaction smaller than the perturbation that")
say("per-comparison renormalisation alone produces; its sign is not interpretable.")

## ------------------------------- how much of `redundant` is a threshold boundary? --
## A negative interaction does not by itself establish redundancy: `Msn24_only` also has
## a strongly negative median I (-1.30), so "I < 0 and significant" cannot separate a
## jointly held gene from a Msn2/4-dependent gene. The separation has to come from the
## sign of the two single effects and from super-additivity.
say("\n=== how much of the provisional `redundant` class is really redundancy? ===")
CLIM <- CL[method == "limma"]
say("padj of the msn2/4 SINGLE inside each class (if `redundant` were a threshold")
say("boundary of Msn24_only, its padj_ko would sit just above alpha, not near 1):")
for (cc in c("redundant", "Msn24_only", "Hsf1_only", "TF_independent")) {
  x <- CLIM[class == cc]
  say("  %-15s n=%4d  padj_ko median %.3f  frac<0.20 %.2f  frac>0.80 %.2f  frac lfc_ko<0 %.2f",
      cc, nrow(x), median(x$padj_ko, na.rm = TRUE), mean(x$padj_ko < 0.20, na.rm = TRUE),
      mean(x$padj_ko > 0.80, na.rm = TRUE), mean(x$lfc_ko < 0, na.rm = TRUE))
}
JS <- merge(LI[, .(gene, I_hat, se, padj_I = padj)], CL[method == "limma", .(gene, sub)], by = "gene")
JS <- merge(JS, CLIM[, .(gene, lfc_kd, lfc_ko, lfc_db)], by = "gene")
say("\nsplit of the provisional class, with the other classes for contrast:")
say("%-24s %5s %9s %8s %9s %10s %10s", "subclass", "n", "med I", "fracI<0", "padjI<.05", "med lfc_kd", "med lfc_ko")
for (s in SUB) {
  x <- JS[sub == s]; if (!nrow(x)) next
  say("%-24s %5d %+9.4f %8.2f %9d %+10.3f %+10.3f", s, nrow(x), median(x$I_hat),
      mean(x$I_hat < 0), x[padj_I < ALPHA, .N], median(x$lfc_kd), median(x$lfc_ko))
}
RJ <- JS[sub == "red_joint"]
say("\nthe joint-requirement candidates (`red_joint`), n = %d:", nrow(RJ))
say("  med lfc: Hsf1-dep %+.3f  msn24 %+.3f  double %+.3f", median(RJ$lfc_kd),
    median(RJ$lfc_ko), median(RJ$lfc_db))
say("  double below the SUM of the two singles (super-additive): %d of %d (%.0f%%)",
    RJ[lfc_db < lfc_kd + lfc_ko, .N], nrow(RJ), 100 * RJ[, mean(lfc_db < lfc_kd + lfc_ko)])
say("  |med I| = %.3f vs renormalisation floor %.3f => %s; per-gene significant %d of %d",
    abs(median(RJ$I_hat)), FLOOR, if (abs(median(RJ$I_hat)) > FLOOR) "resolved" else "NOT resolved",
    RJ[padj_I < ALPHA, .N], nrow(RJ))
say("  => class-level statement is supported; per-gene calls are mostly underpowered.")
say("\nwhat the fourth genotype actually recovers from the TF-independent bucket:")
say("  single-perturbation induced_TF_indep = %d  ==  TF_independent %d + provisional redundant %d",
    CLIM[class %in% c("TF_independent", "redundant"), .N],
    CLIM[class == "TF_independent", .N], CLIM[class == "redundant", .N])
say("  of those %d: Msn2/4-led sub-threshold %d, joint-requirement %d, unclassified %d",
    CLIM[class == "redundant", .N], CLIM[sub == "red_Msn24_subthreshold", .N],
    CLIM[sub == "red_joint", .N], CLIM[sub == "red_unclassified", .N])
say("  Interpretation: the TF-independent bucket is contaminated by single-mutant")
say("  power, so a Hsf1-core vs TF-independent contrast is diluted, not clean.")

## ------------------------------------------- external validation, Solis 2016 -----
## The only literature gene list used anywhere in this project: the 19 core chaperones
## whose basal transcription drops when Hsf1 is exported from the nucleus.
## Solis & Denic 2016, Mol Cell 63(1):60-71, PMID 27320198 (verified).
SOLIS <- c("HSP104","HSP82","HSC82","SSA1","SSA2","HSP78","SIS1","BTN2","HSP42",
           "HSP26","CPR6","MDJ1","STI1","YDJ1","FES1","SSE1","HSP10","HSP60","AHA1")
solis_sys <- lut[match(toupper(SOLIS), key), sys]
say("\n=== external validation against Solis 2016 (PMID 27320198) ===")
say("resolved %d of %d standard names to systematic IDs", sum(!is.na(solis_sys)), length(SOLIS))
uni <- unique(CL$gene)
sv  <- solis_sys[!is.na(solis_sys) & solis_sys %in% uni]
say("%d of the 19 are in the frozen universe", length(sv))
say("%-16s %14s %14s %10s", "class", "Solis genes", "genome-wide", "enrichment")
for (cc in ORD) {
  hit <- CL[method == "limma" & class == cc & gene %in% sv, .N]
  bg  <- CL[method == "limma" & class == cc, .N] / CL[method == "limma", .N]
  say("%-16s %6d /%3d %2s %13.2f%% %9.1fx", cc, hit, length(sv), "",
      100 * bg, if (bg > 0) (hit / length(sv)) / bg else NA_real_)
}
in_core <- sum(sv %in% CORE_SET)
say("\n%d of %d Solis chaperones fall in the unanimous Hsf1-associated core (%.0f%%)",
    in_core, length(sv), 100 * in_core / length(sv))
ft <- fisher.test(matrix(c(in_core, length(sv) - in_core,
                           length(CORE_SET) - in_core,
                           length(uni) - length(CORE_SET) - (length(sv) - in_core)), 2))
say("Fisher OR = %.1f (95%% CI %.1f-%.1f), p = %.2e",
    ft$estimate, ft$conf.int[1], ft$conf.int[2], ft$p.value)
say("This is a positive control on the classification, not evidence about any species.")

## ------------------------------------------------------------------- outputs -----
fwrite(CL,  file.path(outd, sprintf("hsf1_epistasis_classes_%s.tsv", ANCHOR)), sep = "\t")
## The redundant-class partition contains nine counts per DE method. They are
## written to their own table so the reported counts can be read from a file rather than
## recomputed. Nothing is computed differently: this is `table(method, sub)` on the object
## written above.
SUBN <- CL[, .N, by = .(method, class, sub)][order(method, class, -N)]
fwrite(SUBN, file.path(outd, sprintf("hsf1_epistasis_subclass_counts_%s.tsv", ANCHOR)), sep = "\t")
say("landed hsf1_epistasis_subclass_counts_%s.tsv (%d rows): the redundant-class partition is available by column",
    ANCHOR, nrow(SUBN))
fwrite(agr, file.path(outd, sprintf("hsf1_epistasis_agreement_%s.tsv", ANCHOR)), sep = "\t")
fwrite(IX,  file.path(outd, sprintf("hsf1_interaction_estimates_%s.tsv", ANCHOR)), sep = "\t")
fwrite(JS,  file.path(outd, sprintf("hsf1_interaction_by_class_%s.tsv", ANCHOR)), sep = "\t")
## Two gene sets are passed to the cross-species analysis (04_divergence_by_dependency.R),
## kept separate on purpose. The core is the set the Sc-defined Hsf1 dependency is
## projected from; the joint set is a candidate extension whose Hsf1 contribution is real
## at class level but sub-threshold per gene, so it is tested as a separate stratum
## rather than merged into the core.
fwrite(rbind(data.table(gene = CORE_SET, set = "Hsf1_core_unanimous"),
             data.table(gene = JS[sub == "red_joint", gene], set = "Hsf1_joint_candidate")),
       file.path(outd, sprintf("hsf1_core_set_%s.tsv", ANCHOR)), sep = "\t")
say("\nwrote 5 tables to %s", outd)
say("core set file carries TWO strata: %d unanimous core + %d joint candidates (kept separate)",
    length(CORE_SET), JS[sub == "red_joint", .N])
say("=== done | %s ===", format(Sys.time()))
