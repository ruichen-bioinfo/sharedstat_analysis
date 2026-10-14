#!/usr/bin/env Rscript
## 46_p2_sequence_orthogonal_validation.R -- check P2's expression-derived classes against evidence that
## never saw the expression data.
##
## Background. P2's central reinterpretation is that most of the "redundantly controlled" class is sub-threshold
## Msn2/4 dependence rather than joint control. That claim is currently supported only by the same counts it
## reinterprets. Two independent evidence sources are already on disk and neither is derived from these
## libraries: YEASTRACT's literature-curated regulator-target records, and the promoter sequences of the four
## reference genomes. If the reinterpretation is right, the genes it moves should look like Msn2/4 targets by
## sequence and by curation, not like jointly controlled genes.
##
## Decision rule, written before any number was computed:
##   P1. The `redundant` class, and specifically its sub-threshold-Msn2/4 part, should carry STRE (CCCCT)
##       counts and YEASTRACT Msn2/4 target fractions closer to `Msn24_only` than to `both_required`.
##   P2. `Hsf1_only` should be the most HSE-enriched class.
##   Failure condition: if `redundant` sits with `both_required` on both sequence measures, the
##       reinterpretation is not supported by independent evidence. If the
##       classes are indistinguishable on all three measures, the test is uninformative rather than
##       confirmatory, and that is what gets reported.
## No measure will be added, and no class regrouped, after seeing these numbers.
##
## Independence of the evidence. YEASTRACT records include microarray-based "Indirect" evidence from deletion
## experiments, which is the same kind of experiment P2 analyses even though it is not the same data. Direct
## (binding-based) and Indirect records are therefore reported separately, and the Direct column is the one
## the decision rule is judged on. STRE and HSE counts are pure sequence and share nothing with expression.
##
## Fidelity check. Only the per-stratum aggregate of the HSE table is saved, not the per-gene table, so the
## scanner is transplanted from script 10 and must reproduce the saved positive-control table exactly before anything
## new is computed.

suppressPackageStartupMessages({ library(data.table) })
root <- Sys.getenv("PBL_ROOT", "<PBL_ROOT>")
hsd  <- file.path(root, "data/02_processed/hsf1"); csd <- file.path(root, "data/02_processed/crossspecies")
refd <- file.path(root, "data/03_reference")
say  <- function(fmt, ...) cat(sprintf(fmt, ...), "\n", sep = "")

STRE <- "CCCCT"          # Msn2/4 stress response element, motif fixed in the decision rule
NPERM <- 2000L
SEED <- 20260912L

## ------------------------------------------------------------------ transplant scanner from script 10 ----
SRC <- file.path(root, "scripts/hsf1/10_promoter_hse_architecture_v2.R")
EX  <- parse(SRC)
need <- c("count_hse", "PRc")
have <- function() all(vapply(need, exists, TRUE, envir = globalenv()))
n_eval <- 0L
for (e in EX) {
  s <- paste(deparse(e), collapse = " ")
  if (grepl("^sink\\(|^fwrite\\(|quit\\(|^writeLines\\(", s)) next
  if (grepl("^say\\(", s)) next
  try(eval(e, envir = globalenv()), silent = TRUE)
  n_eval <- n_eval + 1L
  if (have()) break
}
say("transplant: evaluated %d of %d expressions from %s", n_eval, length(EX), basename(SRC))
stopifnot(have(), is.function(count_hse), is.data.table(PRc), all(c("sys", "species", "prom") %in% names(PRc)))
say("scanned promoter universe: %d rows, %d genes, %d species",
    nrow(PRc), uniqueN(PRc$sys), uniqueN(PRc$species))

HSg <- cbind(PRc[, .(sys, species)],
             rbindlist(lapply(PRc$prom, function(s) as.data.table(count_hse(s)))))
stopifnot(nrow(HSg) == nrow(PRc), all(c("n_classic", "n_gapped", "max_run") %in% names(HSg)))

## ------------------------------------------------------------------ fidelity check ----
say("\n=== fidelity check: reproduce the saved positive-control table per-stratum carriage ===")
## script 10 reads the per-gene cis/trans table here, not the stratum-aggregated amplitude table; using
## the latter loses the hog column and with it the gene-to-stratum map the check needs.
A6 <- fread(cmd = paste("zcat", shQuote(file.path(hsd, "a6_cis_trans_37_30.tsv.gz"))))
HB <- fread(file.path(csd, "hog_to_scer.tsv")); sysmap <- unique(HB[, .(hog = gene, sys)])
ST <- unique(A6[, .(hog, stratum)])
ST[is.na(stratum) | stratum == "", stratum := "unassigned"]
ST <- merge(ST, sysmap, by = "hog")
HSm <- merge(HSg, ST, by = "sys", allow.cartesian = TRUE)
PC <- HSm[, .(n = .N, any_classic = mean(n_classic > 0), any_gapped = mean(n_gapped > 0),
              any_either = mean(n_classic > 0 | n_gapped > 0), mean_maxrun = mean(max_run)), by = stratum]
REF <- fread(file.path(hsd, "b3v2_positive_control.tsv"))
CMP <- merge(REF, PC, by = "stratum", suffixes = c("_pub", "_our"))
stopifnot(nrow(CMP) == nrow(REF))
w <- max(abs(c(CMP$n_pub - CMP$n_our, CMP$any_classic_pub - CMP$any_classic_our,
               CMP$any_gapped_pub - CMP$any_gapped_our, CMP$any_either_pub - CMP$any_either_our)))
say("  strata compared: %d; worst absolute difference: %.3g", nrow(CMP), w)
stopifnot(w < 1e-12)
say("  reproduced. The scanner and the promoter universe are script 10's, not a reimplementation.")

## ------------------------------------------------------------------ STRE, same promoters ----
## Counted on the same extracted windows the HSE scan used, so the two sequence measures share a universe and
## a window definition. Both strands, because the element is not palindromic and either orientation binds.
revcomp <- function(s) chartr("ACGTacgt", "TGCAtgca", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
count_fixed <- function(s, motif) {
  n1 <- length(gregexpr(motif, s, fixed = TRUE)[[1]])
  if (regexpr(motif, s, fixed = TRUE) == -1L) n1 <- 0L
  rc <- revcomp(motif)
  n2 <- length(gregexpr(rc, s, fixed = TRUE)[[1]])
  if (regexpr(rc, s, fixed = TRUE) == -1L) n2 <- 0L
  n1 + n2
}
PRc[, n_stre := vapply(prom, count_fixed, 0L, motif = STRE)]
PRc[, plen_bp := nchar(prom)]
say("\nSTRE (%s, both strands) counted on the same %d promoter windows; median %d per promoter, max %d",
    STRE, nrow(PRc), median(PRc$n_stre), max(PRc$n_stre))

## ------------------------------------------------------------------ YEASTRACT ----
Y <- fread(file.path(refd, "yeastract2019-flat-file.tsv"), header = FALSE, quote = "", sep = "\t",
           select = 1:9, col.names = c("tf_sys", "tf_name", "tgt_sys", "tgt_name", "pmid", "date",
                                       "condition", "effect", "evidence"))
Y[, tf_name := toupper(tf_name)]
tgt <- function(tf, ev = NULL) {
  z <- Y[tf_name == tf]
  if (!is.null(ev)) z <- z[evidence == ev]
  unique(z$tgt_sys)
}
MSN_D <- union(tgt("MSN2", "Direct"), tgt("MSN4", "Direct"))
MSN_A <- union(tgt("MSN2"), tgt("MSN4"))
HSF_D <- tgt("HSF1", "Direct"); HSF_A <- tgt("HSF1")
say("YEASTRACT: Msn2/4 direct %d targets, any %d; Hsf1 direct %d, any %d",
    length(MSN_D), length(MSN_A), length(HSF_D), length(HSF_A))

## ------------------------------------------------------------------ P2 classes ----
CL <- fread(file.path(hsd, "hsf1_epistasis_classes_37_30.tsv"))
UN <- CL[, .(nd = uniqueN(class), cls = class[1]), by = gene][nd == 1L, .(sys = gene, cls)]
say("\nP2 method-unanimous classes: %d genes over %d classes", nrow(UN), uniqueN(UN$cls))

## the sub-threshold split inside the redundant class, from P2's own `sub` column
subcol <- intersect(c("sub", "subclass"), names(CL))[1]
if (!is.na(subcol)) {
  SB <- unique(CL[class == "redundant", .(sys = gene, sub = get(subcol))])
  SB <- SB[, .(sub = if (uniqueN(sub) == 1L) sub[1] else NA_character_), by = sys]
  say("redundant-class sub-labels available: %s", paste(sort(unique(na.omit(SB$sub))), collapse = ", "))
} else { SB <- data.table(sys = character(), sub = character()); say("no sub-label column found") }

SC <- PRc[species == "Sc", .(sys, n_stre, plen_bp)]
SC <- merge(SC, HSg[species == "Sc", .(sys, n_classic, n_gapped, max_run)], by = "sys")
D <- merge(UN, SC, by = "sys")
D[, `:=`(hse_any = n_classic > 0 | n_gapped > 0,
         stre_per_kb = 1000 * n_stre / plen_bp,
         msn_direct = sys %in% MSN_D, msn_any = sys %in% MSN_A,
         hsf_direct = sys %in% HSF_D, hsf_any = sys %in% HSF_A)]
say("P2 classes with a scanned Sc promoter: %d of %d genes", nrow(D), nrow(UN))

ORD <- c("Hsf1_only", "both_required", "redundant", "Msn24_only", "TF_independent", "unassigned",
         "not_induced")
D[, cls := factor(cls, levels = intersect(ORD, unique(cls)))]
TAB <- D[, .(n = .N,
             hse_pct = 100 * mean(hse_any),
             stre_med = as.numeric(median(n_stre)), stre_per_kb = mean(stre_per_kb),
             msn_direct_pct = 100 * mean(msn_direct), msn_any_pct = 100 * mean(msn_any),
             hsf_direct_pct = 100 * mean(hsf_direct), hsf_any_pct = 100 * mean(hsf_any)), by = cls]
setorder(TAB, cls)
say("\n=== the three orthogonal measures, by P2 class ===")
print(TAB[, .(cls, n, hse_pct = round(hse_pct, 1), stre_med, stre_per_kb = round(stre_per_kb, 2),
              msn_direct_pct = round(msn_direct_pct, 1), msn_any_pct = round(msn_any_pct, 1),
              hsf_direct_pct = round(hsf_direct_pct, 1))], row.names = FALSE)
fwrite(TAB, file.path(hsd, "p2x_class_sequence_curation.tsv"), sep = "\t")

## ------------------------------------------------------------------ the comparisons fixed in the rule ----
say("\n=== P1: is `redundant` closer to `Msn24_only` or to `both_required`? ===")
gap <- function(col) {
  r <- TAB[cls == "redundant"][[col]]; m <- TAB[cls == "Msn24_only"][[col]]; b <- TAB[cls == "both_required"][[col]]
  data.table(measure = col, redundant = r, Msn24_only = m, both_required = b,
             closer_to = if (abs(r - m) < abs(r - b)) "Msn24_only" else "both_required")
}
G <- rbindlist(lapply(c("stre_per_kb", "msn_direct_pct", "msn_any_pct", "hse_pct"), gap))
print(G[, .(measure, redundant = round(redundant, 2), Msn24_only = round(Msn24_only, 2),
            both_required = round(both_required, 2), closer_to)], row.names = FALSE)

set.seed(SEED)
say("\n  permutation test, redundant against both_required, %d draws:", NPERM)
for (v in c("stre_per_kb", "msn_direct")) {
  a <- D[cls == "redundant"][[v]]; b <- D[cls == "both_required"][[v]]
  obs <- mean(a) - mean(b); pool <- c(a, b); na <- length(a)
  nul <- replicate(NPERM, { i <- sample(length(pool), na); mean(pool[i]) - mean(pool[-i]) })
  say("    %-14s redundant %.3f vs both_required %.3f, diff %+.3f, p_emp %.3f",
      v, mean(a), mean(b), obs, (1 + sum(abs(nul) >= abs(obs))) / (1 + NPERM))
}

if (nrow(SB) && any(!is.na(SB$sub))) {
  say("\n=== the same measures inside the redundant class, split by P2's own sub-label ===")
  DS <- merge(D[cls == "redundant"], SB, by = "sys")
  T2 <- DS[!is.na(sub), .(n = .N, hse_pct = 100 * mean(hse_any), stre_per_kb = mean(stre_per_kb),
                          msn_direct_pct = 100 * mean(msn_direct)), by = sub]
  print(T2[order(-n), .(sub, n, hse_pct = round(hse_pct, 1), stre_per_kb = round(stre_per_kb, 2),
                        msn_direct_pct = round(msn_direct_pct, 1))], row.names = FALSE)
  fwrite(T2, file.path(hsd, "p2x_redundant_sublabel_sequence.tsv"), sep = "\t")
}

say("\n=== P2: is `Hsf1_only` the most HSE-enriched class? ===")
o <- TAB[order(-hse_pct)]
say("  HSE carriage ranking: %s", paste(sprintf("%s %.1f%%", o$cls, o$hse_pct), collapse = " > "))
say("  Hsf1_only is top: %s", identical(as.character(o$cls[1]), "Hsf1_only"))
say("")
