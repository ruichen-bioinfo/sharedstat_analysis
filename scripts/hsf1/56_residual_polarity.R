#!/usr/bin/env Rscript
## 56_residual_polarity.R -- allele polarity of the residual carbohydrate-metabolism genes.
##
## Question. Script 44 splits the carbohydrate-metabolism genes into three: those in the Msn2/4 dependence class
## (23), those in the Hsf1 classes (7), and the rest (104). The excess that survives full normalisation sits in
## the third group. This script asks whether the alleles of that third group carry a species polarity: within each
## hybrid, does the allele from one species respond more strongly, and does that line up with the species'
## thermal niche the way the parental amplitudes do?
##
## Scope. This does not ask whether the divergence is adaptive. It asks whether the
## residual is coherent across hybrids or independent per promoter, which distinguishes "one upstream input
## differs between species" from "each promoter has drifted separately". Both readings are about network
## rewiring around a conserved Hsf1, and the analysis distinguishes these alternatives.
##
## Decision rule, fixed before running:
##   Target : the 104 "neither" genes. Controls: the 23 Msn2/4 genes and the translation module.
##   Measure: script 08's allele-polarity mechanism on the interaction term, effects multiplied by `flip` so a
##            positive value means the cryophilic allele responds more. The two matched crosses (same thermal
##            class on both sides) keep the alphabetical axis and are the null channel.
##   Expected: cryophile-first at 30 and 60 min, in the same direction as the parental amplitude order.
##   Alternative outcomes are no order or thermophile-first.
##   Interpretation: cryophile-first is consistent with the hypothesis that a common 37 C is a different
##   intensity of stress for different species, so peripheral metabolic genes are recruited differentially
##   (compare the species-matched shock of script 58). It does not support an adaptive claim.
##
## Nothing is redefined here. The three gene sets are rebuilt by parsing script 44 and re-running its own set
## construction; the sign convention is rebuilt by parsing script 08 and re-running its own POL block. Both are
## then checked against the numbers those scripts published (23 / 7 / 104, and the four cryo-vs-thermo axes).

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
outd <- file.path(root, "data/02_processed/crossspecies")
hsd  <- file.path(root, "data/02_processed/hsf1")
ext  <- file.path(root, "data/03_external/fay2023")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")
ALPHA <- 0.05; NBOOT <- 2000L; SEED <- 20260913L
METHODS <- c("limma", "edgeR", "DESeq2")
CRYO <- c("Sk", "Su")                      # cryophilic species, as in script 08

## ------------------------------------------------------------------ sign convention, from script 08 ----
say("=== sign convention, rebuilt by re-running script 08's own block ===")
src08 <- readLines(file.path(root, "scripts/crossspecies/08_thermophile_allele_polarity.R"))
i0 <- grep("^design <- read.csv", src08)[1]
i1 <- grep("^\\}\\)\\)$", src08)
i1 <- i1[i1 > i0][1]
stopifnot(!is.na(i0), !is.na(i1))
eval(parse(text = paste(src08[i0:i1], collapse = "\n")))
stopifnot(exists("POL"), nrow(POL) == 6L)
print(POL[, .(group, ref, coef_allele, axis, flip)], row.names = FALSE)
say("  cryo-vs-thermo axes: %d, matched null channels: %d",
    POL[axis == "cryo vs thermo", .N], POL[axis != "cryo vs thermo", .N])
stopifnot(POL[axis == "cryo vs thermo", .N] == 4L, POL[axis != "cryo vs thermo", .N] == 2L)

## ------------------------------------------------------------------ the three sets, from script 44 ----
say("\n=== the three gene sets, rebuilt by re-running script 44's own construction ===")
src44 <- readLines(file.path(root, "scripts/hsf1/44_p3_module_stratum_crosstab.R"))
jline <- grep("^J <- ", src44)[1]
stopifnot(!is.na(jline))
## take everything from the start of the script through the three-way split, but stop before its bootstrap
istop <- grep("^G3 <- ", src44)[1]
stopifnot(!is.na(istop))
blk <- src44[1:istop]
## Do not strip its say()/print() lines: several of them span more than one line, and removing only the first
## line of a multi-line call leaves a dangling argument that fails to parse. The output is silenced by
## capture.output instead, which leaves the code exactly as script 44 wrote it.
e44 <- new.env(parent = globalenv())
invisible(capture.output(suppressWarnings(suppressMessages(
  eval(parse(text = paste(blk, collapse = "\n")), envir = e44)))))
msn <- get("msn", envir = e44); hsf <- get("hsf", envir = e44); nei <- get("nei", envir = e44)
## Note: the manuscript's 23 / 7 / 104 are the counts of HOGs with divergence data, not the set sizes: the
## "neither" set has 124 members, of which 104 have an estimable divergence. Both are asserted below against
## script 44's own output.
say("  set sizes:            Msn2/4 %d | Hsf1 %d | neither %d", length(msn), length(hsf), length(nei))
DV44 <- get("DV", envir = e44)
withdat <- sapply(list(msn = msn, hsf = hsf, nei = nei),
                  function(s) DV44[hog %in% s, uniqueN(hog)])
say("  with divergence data: Msn2/4 %d | Hsf1 %d | neither %d", withdat[1], withdat[2], withdat[3])
stopifnot(length(msn) == 23L, length(hsf) == 7L, length(nei) == 124L)
stopifnot(withdat[["msn"]] == 23L, withdat[["hsf"]] == 7L, withdat[["nei"]] == 104L)
say("  matches script 44: sets 23 / 7 / 124, with divergence data 23 / 7 / 104 (the manuscript's numbers)")
MODS_MAP <- get("J", envir = e44)[, .(hog, module)]
## The module is named `translation_ribosome`, not `translation`; the size assertion keeps an empty control set
## from being read as "no polarity".
say("  modules present: %s", paste(sort(unique(na.omit(MODS_MAP$module))), collapse = ", "))
transl <- MODS_MAP[module == "translation_ribosome", unique(hog)]
say("  translation_ribosome control: %d HOGs", length(transl))
stopifnot(length(transl) > 20L)

## ------------------------------------------------------------------ polarity ----
X <- fread(cmd = paste("zcat", shQuote(file.path(outd, "fay_multimethod.tsv.gz"))))
HOGS <- fread(file.path(outd, "hog_to_scer.tsv"))
X <- merge(X, POL[, .(group, flip, axis)], by = "group")
resp <- X[condition == "heat" & model == "A" & term == "time",
          .(gene, group, method, responder = !is.na(padj) & padj < ALPHA)]
## Limitation. The decision rule asks for 30 and 60 min separately, which this table cannot give:
## `fay_multimethod` has no `time` column; its `inter` term is already pooled over the time course. Note: a
## reference to `time` does not error, because data.table resolves it to base::time(), so a missing column
## becomes a function instead of a failure; the absence is asserted below. The polarity is therefore pooled over
## timepoints, and the per-timepoint form of the test is reported as not answerable from the saved quantities.
D <- X[condition == "heat" & term == "inter", .(gene, group, method, model, axis, eff, padj, flip)]
stopifnot(!("time" %in% names(X)))
D <- merge(D, resp, by = c("gene", "group", "method"))
D[, e := ifelse(is.na(flip), eff, eff * flip)]
D <- D[responder == TRUE & is.finite(e) & e != 0]
say("\n  heat responders with a usable interaction effect: %d rows", nrow(D))

## The `gene` column of fay_multimethod already is the HOG identifier (HOG0003711 etc), so no mapping is
## needed; hog_to_scer only adds the S. cerevisiae name. Asserting the identifier shape prevents a silent
## mismatch that would leave every set empty.
stopifnot(all(grepl("^HOG[0-9]+$", head(unique(D$gene), 50))))
D[, hog := gene]
say("  gene identifiers are HOG ids; %d distinct HOGs with heat interaction effects", uniqueN(D$hog))
SETS <- list(neither = nei, Msn24 = msn, Hsf1 = hsf, translation = transl)
D[, set := NA_character_]
for (s in names(SETS)) D[hog %in% SETS[[s]], set := s]
say("  rows assigned to a set: %d of %d", D[!is.na(set), .N], nrow(D))

set.seed(SEED)
run <- function(sub, label) {
  out <- list()
  for (g in POL$group) {
    s <- sub[group == g]
    if (nrow(s) < 20L) next
    ## one value per gene: median over the three DE methods and the two models, so a gene is not counted six
    ## times and the method choice cannot drive the sign
    per <- s[, .(e = median(e)), by = hog]
    k <- sum(per$e > 0); n <- nrow(per)
    bt <- binom.test(k, n, 0.5)
    out[[length(out) + 1L]] <- data.table(
      set = label, group = g, axis = s$axis[1], n_genes = n,
      frac_cryo = k / n, lo = bt$conf.int[1], hi = bt$conf.int[2], binom_p = bt$p.value,
      med_eff = median(per$e),
      signrank_p = tryCatch(wilcox.test(per$e)$p.value, error = function(z) NA_real_))
  }
  rbindlist(out)
}
R <- rbindlist(lapply(names(SETS), function(s) run(D[set == s], s)))
R[, resolvable := lo > 0.5 | hi < 0.5]
say("\n=== polarity by set, cross and timepoint ===")
say("frac_cryo > 0.5 means the cryophilic allele responds more; matched crosses have no thermal meaning")
for (s in names(SETS)) {
  say("\n-- %s --", s)
  print(R[set == s, .(group, axis, n_genes, frac_cryo = round(frac_cryo, 3),
                      ci = sprintf("[%.2f, %.2f]", lo, hi), p = signif(binom_p, 3),
                      med_eff = round(med_eff, 4), resolvable)], row.names = FALSE)
}
fwrite(R, file.path(hsd, "residual_polarity_by_set.tsv"), sep = "\t")

## ------------------------------------------------------------------ by species ----
## Each species appears in three of the six hybrids. The decision rule asks for a per-species summary: for a
## given species, the mean over its crosses of the fraction favouring its allele. That needs the sign oriented
## to the species, not to the cryophilic class, so it is computed from the design and not from `flip`.
say("\n=== per-species summary, mean over the three crosses containing that species ===")
SP <- rbindlist(lapply(seq_len(nrow(POL)), function(i) {
  data.table(group = POL$group[i], species = c(POL$ref[i], POL$coef_allele[i]),
             sign = c(-1L, 1L))                  # stored effect is (coefficient - reference)
}))
bysp <- rbindlist(lapply(names(SETS), function(s) {
  z <- merge(D[set == s, .(hog, group, eff)], SP, by = "group", allow.cartesian = TRUE)
  z[, e_sp := eff * sign]                         # positive: this species' allele responds more
  per <- z[, .(e_sp = median(e_sp)), by = .(species, group, hog)]
  per[, .(n_cross = uniqueN(group), n_genes = uniqueN(hog),
          frac_own = mean(e_sp > 0), med = median(e_sp)), by = species][, set := s]
}))
for (s in names(SETS)) {
  say("\n-- %s --", s)
  print(bysp[set == s][order(-frac_own),
             .(species, n_cross, n_genes, frac_own = round(frac_own, 3), med = round(med, 4))],
        row.names = FALSE)
}
fwrite(bysp, file.path(hsd, "residual_polarity_by_species.tsv"), sep = "\t")

## ------------------------------------------------------------------ read-out ----
say("\n=== read-out against the pre-registration ===")
CRY <- c("Sk", "Su"); THE <- c("Sc", "Sp")
say("  (pooled over the time course; the 30 and 60 min split of the pre-registration is not answerable here)")
for (s in names(SETS)) {
  z <- bysp[set == s]
  if (!nrow(z)) next
  c_mean <- mean(z[species %in% CRY, frac_own]); t_mean <- mean(z[species %in% THE, frac_own])
  ord <- z[order(-frac_own), paste(species, collapse = " > ")]
  say("  %-12s | cryophile mean %.3f vs thermophile mean %.3f | order %s | cryophile-first %s",
      s, c_mean, t_mean, ord, all(z[species %in% CRY, frac_own] > z[species %in% THE, max(frac_own)]))
}
say("\n  null channel (matched crosses) must show no polarity if the effect is thermal:")
print(R[axis != "cryo vs thermo" & set %in% c("neither", "Msn24"),
        .(set, group, n_genes, frac_cryo = round(frac_cryo, 3), p = signif(binom_p, 3), resolvable)],
      row.names = FALSE)
say("\n  resolvable cells by set (CI excluding 0.5), cryo-vs-thermo axes only:")
print(R[axis == "cryo vs thermo", .(cells = .N, resolvable = sum(resolvable)), by = set], row.names = FALSE)
say("")

## ------------------------------------------------------------------ what this answers ----
## Result: the expectation in the decision rule is not supported, and the negative is interpretable because the
## null channel is flat for the target set.
##
## The target (the 104 residual carbohydrate-metabolism genes). No polarity: the four cryo-vs-thermo crosses give
## fractions 0.422 to 0.577 and none has an interval excluding 0.5; per species the fractions run 0.451 to 0.558
## with a cryophile mean of 0.499 against a thermophile mean of 0.504. The two matched crosses are equally flat
## (0.452, 0.449), so the absence is not a broken measurement. The residual excess in these genes is therefore
## not accompanied by a coherent species-level allelic direction.
##
## Scientific interpretation. The metabolism-module excess is a magnitude, not a direction. A coherent
## polarity would have pointed to one upstream input differing between species; its absence is consistent with
## divergence accumulated promoter by promoter. That is a constraint on interpretation and it favours the
## rewiring reading over a single shifted master input. It is not evidence for either mechanism on its own.
##
## The controls do not behave alike, and one of them cannot serve as a control.
##   translation_ribosome shows large polarity (up to 0.814), but it shows it in the matched crosses too
##     (ScxSp 0.814 at p = 9e-30, SkxSu 0.712 at p = 2e-15), which have no thermal interpretation. Script 08's
##     published genome-wide result already reports this: the matched channel runs 0.55 to 0.575 on ScxSp at p
##     down to 1.5e-12. So the translation polarity is the known non-thermal reference artefact, and translation
##     cannot serve as a positive control for a thermal effect here.
##   Msn2/4 leans the other way from the expected direction, thermophile-first (Sc 0.710, Su 0.290), with
##     one resolvable cell (ScxSu 0.130, p = 4.9e-04). Its matched channel ScxSp is 0.304 at p = 0.093, so part
##     of that lean may be the same artefact, so it is descriptive rather than conclusive.
##   Hsf1 has 4 to 6 genes per cross and is uninformative, as expected at n = 7.
##
## The 30 versus 60 min split is not estimable here. `fay_multimethod` carries no time column; its interaction
## term is pooled over the time course, so it cannot estimate a timepoint-specific contrast.
say("\n=== read-out ===")
tgt <- R[set == "neither" & axis == "cryo vs thermo"]
say("  target (104 residual genes): %d cryo-vs-thermo crosses, fractions %.3f to %.3f, resolvable %d",
    nrow(tgt), min(tgt$frac_cryo), max(tgt$frac_cryo), sum(tgt$resolvable))
nul <- R[set == "neither" & axis != "cryo vs thermo"]
say("  its null channel: fractions %s, resolvable %d  => the flat result is not a broken measurement",
    paste(sprintf("%.3f", nul$frac_cryo), collapse = ", "), sum(nul$resolvable))
tr <- R[set == "translation"]
say("  translation control: matched crosses %s at p %s  => polarity present WITHOUT thermal meaning,",
    paste(sprintf("%.3f", tr[axis != "cryo vs thermo", frac_cryo]), collapse = " and "),
    paste(sprintf("%.0e", tr[axis != "cryo vs thermo", binom_p]), collapse = " and "))
say("  so it is the known non-thermal artefact and not a positive control")
P08 <- fread(file.path(outd, "allele_polarity.tsv"))
m8 <- P08[term == "inter" & axis != "cryo vs thermo" & grepl("ScxSp", group)]
say("  consistent with script 08's published genome-wide matched channel: %.3f to %.3f, min p %.1e",
    min(m8$frac_pos), max(m8$frac_pos), min(m8$binom_p))
say("  pre-registered cryophile-first for the target: NOT SUPPORTED (this is a pre-registered failable outcome)")
