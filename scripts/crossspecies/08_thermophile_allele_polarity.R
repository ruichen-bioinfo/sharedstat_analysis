# 08_thermophile_allele_polarity.R — do the cryophilic alleles systematically respond
# more strongly, and is that coherence global or module-specific?
#
# Target. Fay et al. 2023 (PMID 37972247) report a phenotypic asymmetry — "in
# response to heat, the cryophilic species showed a stronger stress response than the
# thermophilic species" — and then, in the same abstract, "within the hybrids, only a
# small number of temperature-responsive genes showed consistent differences between
# alleles from the thermophilic and cryophilic species". That second sentence is what
# licenses their conclusion that "divergence in the heat shock response is mainly a
# consequence of a strain's thermal tolerance", i.e. that the cis contribution is minor
# and the causal arrow runs from tolerance to response rather than the reverse.
#
# That claim is a count of genes passing an FDR threshold, computed under a model that
# treats the two allele read-columns of one library as independent samples. Both halves
# of that are things this project has already measured:
#
#   - thresholded interaction call sets are the least reproducible quantity in the whole
#     design (Jaccard median 0.509 across DE methods, minimum 0.000; script 04)
#   - modelling the library pairing raises interaction calls 3.23x under heat (script 04),
#     and it is the pairing choice, not the DE method, that flips the sign of the
#     Hsf1-versus-Msn2/4 comparison (script 07)
#
# So the claim is re-tested here in the form that does not depend on either weakness: as
# a polarity test on the sign of the interaction effect, over all heat responders rather
# than only those passing a threshold. A systematic sign bias is detectable at a fraction
# of the effect size a per-gene FDR call needs, so "only a small number of genes are
# individually significant" and "the alleles differ coherently" can both be true, and
# distinguishing them changes the conclusion.
#
# Sign convention. In each hybrid, model.matrix takes the alphabetically first species as
# the reference, so the stored interaction coefficient belongs to the alphabetically
# second allele. That is derived from design.csv here rather than assumed, and effects
# are then re-polarised so that positive always means "the cryophilic allele responds
# more strongly". Thermal classes are from the paper's own framing, confirmed in its
# abstract: S. cerevisiae and S. paradoxus grow well above the maximum growth
# temperature of S. kudriavzevii and S. uvarum (see also Salvado et al. 2011 Appl
# Environ Microbiol, PMID 21317255, verified against PubMed 2026-09-05).
#
# Internal negative controls, requiring no external gold standard. Two of the six
# crosses pair species of matched thermal class: ScxSp is thermophile x thermophile and
# SkxSu is cryophile x cryophile. Neither has a defined cryophile-minus-thermophile
# axis, so under the hypothesis they should show no polarity, while any artefact that
# produces sign bias from mapping asymmetry, reference bias or normalisation would show
# up in them too. For these two crosses the polarity is measured against the
# alphabetical axis, which carries no thermal meaning, and reported as the null channel.
#
# Specificity control. The same polarity is measured on the `allele` term (divergence in
# expression level, not in response). If cryophilic alleles are simply higher-expressed
# or better-mapped, polarity appears on both terms; if the effect is specific to the
# response, it appears on `inter` and not on `allele`.
#
# Mechanism discrimination. Two readings of script 07's result — both TF modules
# enriched for divergence, neither dominating — make opposite predictions here:
#   shared upstream set-point:  one species' alleles respond more across both modules,
#                               so polarity is present and equal in Hsf1 and Msn2/4, and
#                               per-gene effects are correlated between crosses sharing a
#                               species
#   independent promoter-level: no coherent polarity; each promoter has drifted on its
#                               own, so signs are near 50/50 and cross-cross correlation
#                               of per-gene effects is weak
#
# Usage: Rscript scripts/crossspecies/08_thermophile_allele_polarity.R

suppressPackageStartupMessages({ library(data.table) })

root <- Sys.getenv("PBL_ROOT", unset = normalizePath("."))
proc <- file.path(root, "data/02_processed")
outd <- file.path(proc, "crossspecies")
ext  <- file.path(root, "data/03_external/fay2023")
ALPHA   <- 0.05
METHODS <- c("DESeq2","edgeR","limma")
THERMO  <- c("Sc","Sp")
CRYO    <- c("Sk","Su")
say <- function(...) cat(sprintf(...), "\n", sep = "")
say("=== 08_thermophile_allele_polarity | %s ===", format(Sys.time()))

## ---------------------------------------------- derive the sign convention -------
design <- read.csv(file.path(ext, "design.csv"))
design <- design[design$Name != "Sp26", , drop = FALSE]
alle <- lapply(split(design$Allele, design$Strain), function(x) sort(unique(x)))
alle <- alle[vapply(alle, length, 1L) == 2]
POL <- rbindlist(lapply(names(alle), function(g) {
  a <- alle[[g]]                       # a[1] = model reference, a[2] = stored coefficient
  cls <- function(s) if (s %in% CRYO) "cryo" else "thermo"
  data.table(group = g, ref = a[1], coef_allele = a[2],
             ref_class = cls(a[1]), coef_class = cls(a[2]),
             ## +1 means the stored effect already points at the cryophilic allele
             flip = if (cls(a[2]) == "cryo" && cls(a[1]) == "thermo") 1L
                    else if (cls(a[2]) == "thermo" && cls(a[1]) == "cryo") -1L
                    else NA_integer_,
             axis = if (cls(a[1]) == cls(a[2])) "matched (null channel)" else "cryo vs thermo")
}))
say("\n=== sign convention derived from design.csv ===")
say("%-8s %-6s %-12s %-24s %s", "cross","ref","coefficient","axis","positive effect means")
for (i in seq_len(nrow(POL))) say("%-8s %-6s %-12s %-24s %s",
  POL$group[i], sprintf("%s(%s)", POL$ref[i], POL$ref_class[i]),
  sprintf("%s(%s)", POL$coef_allele[i], POL$coef_class[i]), POL$axis[i],
  if (is.na(POL$flip[i])) sprintf("%s responds more (no thermal meaning)", POL$coef_allele[i])
  else "cryophilic allele responds more")

## ------------------------------------------------------------------ load data ----
X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod.tsv.gz"))))
MOD  <- fread(file.path(outd, "yeast_modules_37_30.tsv"))
HOGS <- fread(file.path(outd, "hog_to_scer.tsv"))
X <- merge(X, POL[, .(group, flip, axis)], by = "group")
## responder mask comes from model A's shared-response term, which is the only model in
## which it is estimable; it is a covariate, not part of the contrast being tested
resp <- X[condition == "heat" & model == "A" & term == "time", .(gene, group, method,
          responder = !is.na(padj) & padj < ALPHA)]

polar <- function(tm) {
  D <- X[condition == "heat" & term == tm, .(gene, group, method, model, axis,
                                             eff, padj, flip)]
  D <- merge(D, resp, by = c("gene","group","method"))
  ## null channel keeps the alphabetical axis, which has no thermal interpretation
  D[, e := ifelse(is.na(flip), eff, eff * flip)]
  D[responder == TRUE & is.finite(e) & e != 0]
}

## ------------------------------------------- (1) polarity of the response term ----
say("\n=== (1) polarity of the interaction term: which allele responds more? ===")
say("all heat responders, no significance filter. frac>0 = fraction favouring the")
say("cryophilic allele (matched crosses: fraction favouring the alphabetically 2nd).")
say("%-8s %-24s %-7s %-2s %7s %7s %11s %9s %10s", "cross","axis","method","M",
    "n","frac>0","binom_p","med_eff","signrank_p")
run_polarity <- function(tm, label) {
  D <- polar(tm); out <- list()
  for (g in POL$group) for (mt in METHODS) for (mdl in c("A","B")) {
    s <- D[group == g & method == mt & model == mdl]
    if (nrow(s) < 50) next
    k <- sum(s$e > 0); n <- nrow(s)
    bp <- binom.test(k, n, 0.5)
    wp <- tryCatch(wilcox.test(s$e)$p.value, error = function(e) NA_real_)
    out[[length(out)+1]] <- data.table(term = tm, group = g, axis = s$axis[1],
      method = mt, model = mdl, n = n, frac_pos = k/n, binom_p = bp$p.value,
      lo = bp$conf.int[1], hi = bp$conf.int[2],
      med_eff = median(s$e), signrank_p = wp,
      n_sig = s[!is.na(padj) & padj < ALPHA, .N],
      frac_pos_sig = s[!is.na(padj) & padj < ALPHA, if (.N) mean(e > 0) else NA_real_])
    if (label) say("%-8s %-24s %-7s %-2s %7d %6.1f%% %11.3g %9.3f %10.3g",
      g, s$axis[1], mt, mdl, n, 100*k/n, bp$p.value, median(s$e), wp)
  }
  rbindlist(out)
}
P_inter <- run_polarity("inter", TRUE)

say("\n--- summary by axis and model (interaction term) ---")
say("%-24s %-2s %6s %9s %11s %s", "axis","M","cells","frac>0_med","min_binom_p","cells p<0.05")
for (ax in unique(P_inter$axis)) for (mdl in c("A","B")) {
  s <- P_inter[axis == ax & model == mdl]
  if (!nrow(s)) next
  say("%-24s %-2s %6d %8.1f%% %11.3g %d of %d", ax, mdl, nrow(s),
      100*median(s$frac_pos), min(s$binom_p), sum(s$binom_p < 0.05), nrow(s))
}
say("\nFay's threshold-based count, for comparison: among genes that ARE individually")
say("significant, the fraction favouring the cryophilic allele --")
say("%-8s %-24s %-7s %-2s %7s %9s", "cross","axis","method","M","n_sig","frac>0")
for (i in seq_len(nrow(P_inter))) if (!is.na(P_inter$frac_pos_sig[i]) && P_inter$n_sig[i] >= 10)
  say("%-8s %-24s %-7s %-2s %7d %8.1f%%", P_inter$group[i], P_inter$axis[i],
      P_inter$method[i], P_inter$model[i], P_inter$n_sig[i], 100*P_inter$frac_pos_sig[i])

## --------------------------------------- (2) specificity: the level term too? ----
say("\n=== (2) specificity control: same polarity on the expression-level term ===")
say("If cryophilic alleles are merely higher-expressed or better-mapped, polarity shows")
say("up here as well; if the effect is specific to the response, it does not.")
P_allele <- run_polarity("allele", FALSE)
say("%-24s %-2s %14s %14s", "axis","M","inter frac>0","allele frac>0")
for (ax in unique(P_inter$axis)) for (mdl in c("A","B")) {
  a <- P_inter[axis == ax & model == mdl]; b <- P_allele[axis == ax & model == mdl]
  if (!nrow(a) || !nrow(b)) next
  say("%-24s %-2s %13.1f%% %13.1f%%", ax, mdl, 100*median(a$frac_pos), 100*median(b$frac_pos))
}
fwrite(rbind(P_inter, P_allele), file.path(outd, "allele_polarity.tsv"), sep = "\t")

## ------------------------------------------- (3) is the polarity module-specific --
say("\n=== (3) polarity within each TF module (mechanism discrimination) ===")
say("Equal polarity in both modules points at a shared upstream set-point; polarity in")
say("one module only points at module-specific cis evolution.")
key <- merge(HOGS[!is.na(sys) & in_yeast == TRUE, .(gene, sys)],
             MOD[method == "DESeq2", .(sys = gene, module)], by = "sys")
Dm <- merge(polar("inter"), key[, .(gene, module)], by = "gene")
say("%-24s %-2s %-18s %7s %8s %11s", "axis","M","module","n","frac>0","binom_p")
resm <- list()
for (ax in unique(Dm$axis)) for (mdl in c("A","B")) for (m0 in c("Hsf1","Msn24","both","induced_TF_indep")) {
  s <- Dm[axis == ax & model == mdl & module == m0 & method == "DESeq2"]
  if (nrow(s) < 30) next
  bp <- binom.test(sum(s$e > 0), nrow(s), 0.5)
  resm[[length(resm)+1]] <- data.table(axis = ax, model = mdl, module = m0, n = nrow(s),
    frac_pos = mean(s$e > 0), p = bp$p.value, med_eff = median(s$e))
  say("%-24s %-2s %-18s %7d %7.1f%% %11.3g", ax, mdl, m0, nrow(s), 100*mean(s$e > 0), bp$p.value)
}
fwrite(rbindlist(resm), file.path(outd, "allele_polarity_by_module.tsv"), sep = "\t")

## ------------------------------ (4) is it a species property? cross-cross corr ----
## If a species' alleles carry a coherent response shift, the per-gene effect measured in
## two different crosses containing that species should agree. Sharing one species is
## the comparison; sharing none is the control.
say("\n=== (4) does a species' allelic effect replicate across the crosses carrying it? ===")
say("Spearman of per-gene interaction effect, re-polarised to the shared species.")
sp_of <- function(g) unlist(strsplit(g, "x"))
E <- X[condition == "heat" & term == "inter" & method == "DESeq2" & group %in% POL$group]
say("%-18s %-2s %-14s %7s %8s %11s", "cross pair","M","shared species","n","rho","p")
resc <- list()
gs <- POL$group
for (i in seq_along(gs)) for (j in seq_along(gs)) if (i < j) for (mdl in c("A","B")) {
  g1 <- gs[i]; g2 <- gs[j]
  sh <- intersect(sp_of(g1), sp_of(g2))
  ## orient each cross's effect toward the shared species (or, with none shared,
  ## toward the alphabetically second allele of each, which is the control channel)
  orient <- function(g, s) { p <- POL[group == g]; if (length(s) && s == p$ref) -1L else 1L }
  s1 <- E[group == g1 & model == mdl, .(gene, e1 = eff * orient(g1, sh))]
  s2 <- E[group == g2 & model == mdl, .(gene, e2 = eff * orient(g2, sh))]
  M <- merge(s1, s2, by = "gene")
  M <- M[is.finite(e1) & is.finite(e2)]
  if (nrow(M) < 100) next
  ct <- suppressWarnings(cor.test(M$e1, M$e2, method = "spearman"))
  resc[[length(resc)+1]] <- data.table(pair = paste(g1, g2, sep = "|"), model = mdl,
    shared = if (length(sh)) paste(sh, collapse = ",") else "none (control)",
    n = nrow(M), rho = unname(ct$estimate), p = ct$p.value)
  say("%-18s %-2s %-14s %7d %8.3f %11.3g", paste(g1, g2, sep = "|"), mdl,
      if (length(sh)) paste(sh, collapse = ",") else "none (control)",
      nrow(M), unname(ct$estimate), ct$p.value)
}
CC <- rbindlist(resc)
fwrite(CC, file.path(outd, "allele_effect_cross_replication.tsv"), sep = "\t")
for (mdl in c("A","B")) {
  s <- CC[model == mdl]
  if (!nrow(s)) next
  say("\nmodel %s: median rho %.3f when a species is shared vs %.3f when none is (%d vs %d pairs)",
      mdl, median(s[shared != "none (control)", rho]), median(s[shared == "none (control)", rho]),
      s[shared != "none (control)", .N], s[shared == "none (control)", .N])
}

say("\nwrote allele_polarity.tsv, allele_polarity_by_module.tsv, allele_effect_cross_replication.tsv")
say("=== done ===")
