# 06_module_bridge.R — does cross-species divergence of the heat-shock response sit
# in the Hsf1 module or the Msn2/4 module?
#
# This is the script that joins the two halves of the project.
#
# Fay et al. 2023 can see that the heat-shock response has diverged in cis between
# Saccharomyces species, and their GO enrichment points at de novo protein folding
# (chaperones, n=8, P=7.0e-5) and cytoplasmic translation. But with only cross-species
# data they must attribute genes to pathways through annotation. Their stated
# conclusion goes the other way from that enrichment: because the more heat-sensitive
# species mount the stronger response, they read the heat-shock response as "a signal
# of heat stress rather than a driver of thermotolerance". Chaperone induction being
# the most diverged module and the response being a passive stress readout are in
# tension, and annotation-based attribution cannot resolve it.
#
# GSE179258 (Muhlhofer et al. 2024 FEBS Lett, PMID 38366111)
# is S. cerevisiae with Hsf1 depleted, with msn2/msn4 deleted, and with both, across
# 25/37/42 C. That lets each gene's regulator be assigned from a perturbation rather
# than from an annotation:
#
#   induced       : up in wild type after the shift
#   Hsf1-depend.  : induced, and significantly lower when Hsf1 is depleted
#   Msn2/4-depend.: induced, and significantly lower when msn2/msn4 are deleted
#
# and then asked whether the genes whose response has diverged between species are
# drawn preferentially from one module.
#
# Anchor condition. Fay shift cells 25 to 37 C and sample at 15/30/60 min, so the
# comparable yeast contrast is 37 C, not 42 C. Within 37 C the 30-min timepoint is
# used: at 10 min the Hsf1-depleted strain differs from wild type in only 65-103 genes
# depending on method, versus 2792-2863 at 30 min, i.e. the Hsf1 effect at mild heat is
# late-onset. Defining the module at 10 min would define it out of existence. 42 C/30
# min is carried as a sensitivity analysis.
#
# Both halves are DE calls, so both inherit the method dependence measured in 02-04.
# The module assignment is therefore made under each of DESeq2, edgeR and limma
# separately, and the cross-species divergence set likewise, giving a 3x3 grid. If the
# answer to "which module" is the same in all nine cells, that evolutionary conclusion
# is robust to method choice; if it flips, the flip is the result. Both outcomes are
# scientifically interpretable.
#
# Agreement and disagreement among methods are
# summarized without treating any method as ground truth.
#
# Usage: Rscript scripts/crossspecies/06_module_bridge.R [--anchor=37_30]

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
proc <- file.path(root, "data/02_processed")
outd <- file.path(proc, "crossspecies")
ref  <- file.path(root, "data/03_reference")
ALPHA <- 0.05
METHODS <- c("DESeq2","edgeR","limma")
say <- function(...) cat(sprintf(...), "\n", sep = "")
read_maybe_gz <- function(p) if (grepl("\\.gz$", p)) fread(cmd = paste("zcat", shQuote(p))) else fread(p)

args <- commandArgs(trailingOnly = TRUE)
h <- grep("^--anchor=", args, value = TRUE)
ANCHOR <- if (length(h)) sub("^--anchor=", "", h[1]) else "37_30"
say("=== 06_module_bridge | %s | anchor %s C ===", format(Sys.time()), ANCHOR)

## ------------------------------------------------ systematic and standard names ---
## Fay's ortholog table names genes by the S. cerevisiae standard name where one exists
## (PDR3, IRA1) and falls back to the systematic name (YLR407W). Our count matrix is
## keyed on systematic names throughout, so the Ensembl GTF supplies the bridge.
gtf <- file.path(ref, "Saccharomyces_cerevisiae.R64-1-1.75.gtf")
gl  <- grep("\tgene\t", readLines(gtf, warn = FALSE), fixed = FALSE, value = TRUE)
sysn <- sub('.*gene_id "([^"]+)".*', "\\1", gl)
stdn <- ifelse(grepl('gene_name "', gl), sub('.*gene_name "([^"]+)".*', "\\1", gl), NA_character_)
nmap <- unique(data.table(sys = sysn, std = stdn))
say("GTF: %d gene records, %d with a standard name", nrow(nmap), sum(!is.na(nmap$std)))

## key lookup accepts either spelling, case-insensitively
lut <- rbind(nmap[, .(key = toupper(sys), sys)],
             nmap[!is.na(std), .(key = toupper(std), sys)])
lut <- unique(lut)[!duplicated(key)]

## ------------------------------------------------------- yeast module definition ---
load_comp <- function(method, comp) {
  f <- file.path(proc, method, paste0(comp, ".tsv"))
  if (!file.exists(f)) stop("missing batch output: ", f)
  d <- fread(f, select = c("gene","log2FC","padj"))
  setnames(d, c("gene","lfc","padj"))
  d
}

## classify every gene under one method; `sig_dn` means the perturbation lowered it
classify <- function(method, anchor) {
  wt  <- load_comp(method, sprintf("Wildtype_%s_vs_25_0", anchor))
  kd  <- load_comp(method, sprintf("HSF1.KD_vs_Wildtype_at_%s", anchor))
  ko  <- load_comp(method, sprintf("MSN24.KO_vs_Wildtype_at_%s", anchor))
  setnames(wt, c("gene","lfc_wt","padj_wt"))
  setnames(kd, c("gene","lfc_kd","padj_kd"))
  setnames(ko, c("gene","lfc_ko","padj_ko"))
  d <- merge(merge(wt, kd, by = "gene"), ko, by = "gene")
  d[, induced  := !is.na(padj_wt) & padj_wt < ALPHA & lfc_wt > 0]
  d[, repressed := !is.na(padj_wt) & padj_wt < ALPHA & lfc_wt < 0]
  ## dependence = the perturbation moves the gene back toward baseline
  d[, hsf1_dep := induced & !is.na(padj_kd) & padj_kd < ALPHA & lfc_kd < 0]
  d[, msn_dep  := induced & !is.na(padj_ko) & padj_ko < ALPHA & lfc_ko < 0]
  d[, hsf1_dep_rep := repressed & !is.na(padj_kd) & padj_kd < ALPHA & lfc_kd > 0]
  d[, msn_dep_rep  := repressed & !is.na(padj_ko) & padj_ko < ALPHA & lfc_ko > 0]
  d[, module := fifelse(hsf1_dep & msn_dep, "both",
                fifelse(hsf1_dep, "Hsf1", fifelse(msn_dep, "Msn24",
                fifelse(induced, "induced_TF_indep", "not_induced"))))]
  d[, method := method]
  d[]
}

MOD <- rbindlist(lapply(METHODS, classify, anchor = ANCHOR))
say("\n=== yeast module sizes at %s C (per DE method) ===", ANCHOR)
tab <- dcast(MOD[, .N, by = .(method, module)], module ~ method, value.var = "N", fill = 0)
say("%-18s %8s %8s %8s", "module", METHODS[1], METHODS[2], METHODS[3])
for (i in seq_len(nrow(tab))) say("%-18s %8d %8d %8d", tab$module[i],
  tab[[METHODS[1]]][i], tab[[METHODS[2]]][i], tab[[METHODS[3]]][i])
say("\nrepression arm (RPGs live here): Hsf1-dep %s | Msn2/4-dep %s",
    paste(MOD[, sum(hsf1_dep_rep), by = method]$V1, collapse = "/"),
    paste(MOD[, sum(msn_dep_rep),  by = method]$V1, collapse = "/"))

## robustness of the assignment itself: how often do the three methods agree?
agr <- dcast(MOD[, .(gene, method, module)], gene ~ method, value.var = "module")
agr[, n_distinct := apply(.SD, 1, function(x) uniqueN(x)), .SDcols = METHODS]
agr[, gold := n_distinct == 1]
say("\nmodule assignment agreement across the three methods: %d of %d genes identical (%.1f%%)",
    sum(agr$gold), nrow(agr), 100*mean(agr$gold))
any_tf  <- apply(agr[, ..METHODS], 1, function(x) any(x %in% c("Hsf1","Msn24","both")))
gold_tf <- agr$gold & agr[[METHODS[1]]] %in% c("Hsf1","Msn24","both")
say("among genes any method assigns to a TF module: %d, of which %d get the same assignment from all three (%.1f%%)",
    sum(any_tf), sum(gold_tf), 100*sum(gold_tf)/max(sum(any_tf), 1))
fwrite(MOD, file.path(outd, sprintf("yeast_modules_%s.tsv", ANCHOR)), sep = "\t")
fwrite(agr, file.path(outd, sprintf("yeast_module_agreement_%s.tsv", ANCHOR)), sep = "\t")

## ------------------------------------------------------------------- the bridge ---
X <- read_maybe_gz(file.path(outd, "fay_multimethod.tsv.gz"))
hyb <- c("ScxSp","ScxSk","ScxSu","SpxSk","SpxSu","SkxSu")
hogs <- unique(X[, .(gene, Name)])
hogs[, sys := lut$sys[match(toupper(Name), lut$key)]]
say("\n=== bridge: Fay HOG -> S. cerevisiae systematic name ===")
say("HOGs: %d | mapped: %d (%.1f%%) | unmapped example: %s",
    nrow(hogs), sum(!is.na(hogs$sys)), 100*mean(!is.na(hogs$sys)),
    paste(head(hogs[is.na(sys)]$Name, 6), collapse = ", "))
## and how many of those land on a gene our yeast run actually tested
tested <- unique(MOD$gene)
hogs[, in_yeast := sys %in% tested]
say("of the mapped HOGs, %d are also in the yeast shared universe (%.1f%% of all HOGs)",
    sum(hogs$in_yeast, na.rm = TRUE), 100*mean(hogs$in_yeast, na.rm = TRUE))
fwrite(hogs, file.path(outd, "hog_to_scer.tsv"), sep = "\t")

## -------------------------------------------------------- the hypothesis test ----
## Background is deliberately narrow: HOGs that map to a yeast gene we tested, that
## respond to heat in this hybrid (otherwise there is no response to diverge in), and
## that are heat-induced in wild-type yeast (otherwise the Hsf1/Msn2-4 module
## assignment does not apply). Foreground is the genes called as cis-diverged in
## response. Any enrichment is therefore relative to other heat-induced responders,
## not to the genome, which is the comparison the question is actually about.
fisher_row <- function(a, b, c_, d) {
  ft <- fisher.test(matrix(c(a, b, c_, d), 2, 2))
  list(or = unname(ft$estimate), lo = ft$conf.int[1], hi = ft$conf.int[2], p = ft$p.value)
}

res <- list()
for (ym in METHODS) {
  md <- MOD[method == ym, .(sys = gene, module)]
  key <- merge(hogs[!is.na(sys) & in_yeast == TRUE, .(gene, sys)], md, by = "sys")
  for (xm in METHODS) for (mdl in c("A","B")) for (g in hyb) {
    S <- X[condition == "heat" & group == g & method == xm & model == mdl]
    if (nrow(S) == 0) next
    tt <- S[term == "time",  .(gene, padj)]; setnames(tt, "padj", "padj_time")
    ii <- S[term == "inter", .(gene, padj)]; setnames(ii, "padj", "padj_inter")
    D <- merge(merge(key, ii, by = "gene"), tt, by = "gene", all.x = TRUE)
    ## model B has no estimable `time` term; fall back to model A's responder set
    if (all(is.na(D$padj_time))) {
      ttA <- X[condition == "heat" & group == g & method == xm & model == "A" & term == "time",
               .(gene, padj_time = padj)]
      D[, padj_time := NULL]; D <- merge(D, ttA, by = "gene", all.x = TRUE)
    }
    bg <- D[!is.na(padj_time) & padj_time < ALPHA & module %in% c("Hsf1","Msn24","both","induced_TF_indep")]
    if (nrow(bg) < 50) next
    bg[, div := !is.na(padj_inter) & padj_inter < ALPHA]
    ## (a) each module against all other heat-induced responders
    for (m0 in c("Hsf1","Msn24","both")) {
      a <- bg[module == m0 & div == TRUE, .N]; b <- bg[module == m0 & div == FALSE, .N]
      c_ <- bg[module != m0 & div == TRUE, .N]; d <- bg[module != m0 & div == FALSE, .N]
      if (a + b < 10) next
      f <- fisher_row(a, b, c_, d)
      res[[length(res)+1]] <- data.table(test = "vs_rest", yeast_method = ym,
        xs_method = xm, model = mdl, group = g, module = m0,
        n_mod = a + b, n_div = a, rate = a/(a+b),
        rate_rest = c_/max(c_+d,1), or = f$or, lo = f$lo, hi = f$hi, p = f$p)
    }
    ## (b) head to head: Hsf1-only versus Msn2/4-only, which is the actual question
    a <- bg[module == "Hsf1" & div == TRUE, .N]; b <- bg[module == "Hsf1" & div == FALSE, .N]
    c_ <- bg[module == "Msn24" & div == TRUE, .N]; d <- bg[module == "Msn24" & div == FALSE, .N]
    if (a + b >= 10 && c_ + d >= 10) {
      f <- fisher_row(a, b, c_, d)
      res[[length(res)+1]] <- data.table(test = "Hsf1_vs_Msn24", yeast_method = ym,
        xs_method = xm, model = mdl, group = g, module = "Hsf1_over_Msn24",
        n_mod = a + b, n_div = a, rate = a/(a+b), rate_rest = c_/max(c_+d,1),
        or = f$or, lo = f$lo, hi = f$hi, p = f$p)
    }
  }
}
R <- rbindlist(res)
fwrite(R, file.path(outd, sprintf("module_divergence_%s.tsv", ANCHOR)), sep = "\t")

## ----------------------------------------------------------------- reporting -----
say("\n=== (a) divergence enrichment of each module, vs other heat-induced responders ===")
say("(model A, pooled over hybrids: median odds ratio and how many hybrids exceed 1)")
say("%-14s %-7s %-7s %6s %8s %8s %10s", "module","yeast","xspec","n_hyb","OR_med","OR_rng","p<0.05")
for (m0 in c("Hsf1","Msn24","both")) for (ym in METHODS) for (xm in METHODS) {
  s <- R[test == "vs_rest" & module == m0 & yeast_method == ym & xs_method == xm & model == "A"]
  if (nrow(s) == 0) next
  say("%-14s %-7s %-7s %6d %8.2f %8s %10d", m0, ym, xm, nrow(s), median(s$or),
      sprintf("%.1f-%.1f", min(s$or), max(s$or)), sum(s$p < 0.05))
}

say("\n=== (b) THE TEST: Hsf1 module versus Msn2/4 module, head to head ===")
say("OR > 1 means the Hsf1-dependent genes are the more diverged; OR < 1 means Msn2/4.")
say("%-7s %-7s %-2s %-8s %7s %7s %9s %9s %8s %11s", "yeast","xspec","M","hybrid",
    "n_Hsf1","n_Msn", "rate_Hsf1","rate_Msn","OR","p")
HH <- R[test == "Hsf1_vs_Msn24"]
setorder(HH, model, yeast_method, xs_method, group)
for (i in seq_len(nrow(HH))) say("%-7s %-7s %-2s %-8s %7d %7s %8.1f%% %8.1f%% %8.2f %11.3g",
  HH$yeast_method[i], HH$xs_method[i], HH$model[i], HH$group[i], HH$n_mod[i], "-",
  100*HH$rate[i], 100*HH$rate_rest[i], HH$or[i], HH$p[i])

say("\n--- verdict grid: sign of the Hsf1-vs-Msn2/4 comparison in each of the 9 method pairs ---")
say("%-7s %-7s %-2s %8s %8s %8s %s", "yeast","xspec","M","n_hyb","OR>1","OR<1","reading")
for (mdl in c("A","B")) for (ym in METHODS) for (xm in METHODS) {
  s <- HH[yeast_method == ym & xs_method == xm & model == mdl]
  if (nrow(s) == 0) next
  up <- sum(s$or > 1); dn <- sum(s$or < 1)
  say("%-7s %-7s %-2s %8d %8d %8d %s", ym, xm, mdl, nrow(s), up, dn,
      if (up > dn) "Hsf1 module more diverged" else if (dn > up) "Msn2/4 module more diverged" else "tied")
}
sg <- HH[model == "A"]
if (nrow(sg)) say("\nmodel A overall: %d of %d hybrid x method-pair cells favour Hsf1 (median OR %.2f); %d cells reach p<0.05",
                  sum(sg$or > 1), nrow(sg), median(sg$or), sum(sg$p < 0.05))

say("\nwrote yeast_modules_%s.tsv, yeast_module_agreement_%s.tsv, hog_to_scer.tsv, module_divergence_%s.tsv",
    ANCHOR, ANCHOR, ANCHOR)
say("=== done ===")
